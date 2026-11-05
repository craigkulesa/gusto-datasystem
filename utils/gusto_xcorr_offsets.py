#!/usr/bin/env python
"""
GUSTO mixer offset pipeline (skeleton).

Measure on-sky offsets between GUSTO mixers by cross-correlating
per-mixer sub-images against the FIDUCIAL mixer of the band over
small sky regions and velocity ranges, then convert the sky offsets
to focal-plane offsets using the gondola attitude (AZ/EL/roll) of the
Level-1 OTF scans recorded in each cube's scan-list extension.
All offsets (and focal-plane angles) are defined with respect to the
fiducial mixer, whose offset is zero by definition.

Dependencies: astropy, scipy, numpy, matplotlib, bayesicfitting, PyYAML.
No spectral_cube is used: cubes are handled as numpy arrays + astropy WCS.
"""

import argparse
import logging
import os
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import yaml
from scipy import ndimage

import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.io import fits
from astropy.table import QTable
from astropy.time import Time
from astropy.wcs import WCS

from BayesicFitting import (Gauss, Kernel2dModel, LevenbergMarquardtFitter,
                            PolySurfaceModel)

log = logging.getLogger("gusto_xcorr")


def load_config(path):
    """Parse the general configuration file (gusto_xcorr.yaml).

    Structure:
        data:
          level1_root: ...     # Level-1 root; per-target products:
                              #   <level1_root>/<target>/BAND_{seq}_{scanid}_L10.fits
          level2_root: ...     # Level-2 root; per-mixer cubes:
                              #   <level2_root>/<target>/...
        bands:
          <BAND>:
            fiducial: <mixer>          # offset reference, zero by definition
            mixers: [<mixer>, ...]     # usable mixers for this band
            cubes:                     # optional explicit filenames
              <mixer>: <file>          # default: {band}_L2_{mixer}.fits
        options:
          min_snr / collapse / patch / beam / kernel /
          roll / pix_band / cal_file / selfcheck /
          launch_lat / launch_lon

    Returns (data_cfg, bands, options).
    """
    global _SYS_FLOOR_DEG, _SYS_FLOOR_BY_BAND
    with open(path) as f:
        cfg = yaml.safe_load(f)
    data_cfg = cfg["data"]
    # read options BEFORE the band loop: per-band values may fall
    # back to the global option (pix_band)
    options = cfg.get("options", {})
    # systematic floor on measurement sigmas (deg), in quadrature
    # everywhere inputs are combined; default 0.005 deg = 18"
    _SYS_FLOOR_DEG = float(options.get("sys_floor_deg", 0.005))
    _SYS_FLOOR_BY_BAND = {}
    # top-level analysis-time exclusion list (see apply_exclusions):
    # entries match band/target/region/mixer/vel subsets, each with
    # a 'reason'; the ECSV measurement record is NEVER edited
    options["exclude"] = cfg.get("exclude", []) or []
    bands = {}
    for band, b in cfg["bands"].items():
        fid = b["fiducial"]
        # mixers is OPTIONAL: when absent, the mixer set is
        # auto-discovered from the cube naming convention at run time
        mixers = list(b.get("mixers", []))
        if mixers and fid not in mixers:
            mixers.append(fid)
        # pix_band: which PIXEL ARRAY (B1/B2) this band's mixers
        # belong to, for the calibration file naming. PER-BAND,
        # not global -- e.g. NII is the B1 array while CII is B2;
        # a global setting mislabels one band's rows (NII mixers
        # written as B2M* when they are really B1M*).
        pix_band = b.get("pix_band", options.get("pix_band", "B2"))
        # per-band systematic floor override (deg): e.g. NII
        # typically needs a larger floor than the CII default
        if b.get("sys_floor_deg") is not None:
            _SYS_FLOOR_BY_BAND[band] = float(b["sys_floor_deg"])
        # optional explicit cube filenames (see main(); default is
        # the naming convention <target>_<band>_<mixer>_at_<Beam>_<kernel>.fits,
        # which needs the target/beam/kernel resolved at run time)
        cubes = dict(b.get("cubes", {}))
        bands[band] = {"fiducial": fid, "mixers": mixers,
                       "pix_band": pix_band,
                       "cubes": cubes}
    log.info("Config: bands %s, level1_root %s, level2_root %s",
             list(bands), data_cfg.get("level1_root"),
             data_cfg.get("level2_root"))
    return data_cfg, bands, options


# ----------------------------------------------------------------------
# 1. Inputs and region definitions
# ----------------------------------------------------------------------

@dataclass
class Region:
    name: str
    frame: str                  # 'galactic' or 'icrs'
    lon: u.Quantity             # deg ( Galactic l or RA )
    lat: u.Quantity             # deg ( Galactic b or Dec )
    radius: u.Quantity          # arcsec; ANGULAR half-size of the sky
                                # box (never pixels -- generic regions)
    vel_ranges: list = field(default_factory=list)   # [(vlo, vhi), ...] km/s
    min_snr: float = 5.0


def load_regions(path):
    """Parse the YAML region file -> list[Region]."""
    with open(path) as f:
        cfg = yaml.safe_load(f)
    regions = []
    for r in cfg["regions"]:
        regions.append(Region(
            name=r["name"],
            frame=r.get("frame", "galactic"),
            lon=u.Quantity(r["lon"]),
            lat=u.Quantity(r["lat"]),
            radius=u.Quantity(r.get("radius", 300 * u.arcsec)),
            vel_ranges=[tuple(v) for v in r.get("vel_ranges", [])],
            min_snr=float(r.get("min_snr", 5.0)),
        ))
    log.info("Loaded %d regions from %s", len(regions), path)
    return regions


def load_cube(path):
    """Open a GUSTO cube: (data, header, wcs, seq_map).

    data: 3D array, spatial axes + spectral axis as given by the WCS.
    seq_map: 2D SEQ_ID extension (HDU 3) -- same spatial WCS as the
    data; at each spatial pixel it records the SEQ_ID of the Level-1
    sequence (calibratable unit: OTF scan + REF/HOT) that contributed
    there. The SEQ_ID is the FIRST number in the L1 filename, so the
    filename lookup is exact.
    """
    with fits.open(path) as hdul:
        hdu = hdul[0] if hdul[0].data is not None else hdul[1]
        data = np.asarray(hdu.data, dtype=float)
        header = hdu.header.copy()
        # SEQ_ID is HDU 3 (verify the name as a sanity check)
        try:
            seq_hdu = hdul[3]
            seq_map = seq_hdu.data
            if seq_hdu.name and seq_hdu.name.upper() != "SEQ_ID":
                log.warning("HDU 3 named %r, expected SEQ_ID "
                            "(using it anyway)", seq_hdu.name)
        except IndexError:
            # fallback: search by extension name
            seq_map = None
            for h in hdul[1:]:
                if h.name and h.name.upper() == "SEQ_ID":
                    seq_map = h.data
                    break
        if seq_map is None:
            log.warning("No SEQ_ID extension found in %s", path)
    wcs = WCS(header)
    log.info("Loaded %s: %s, SEQ_ID ext: %s", path, data.shape,
             "n/a" if seq_map is None else str(seq_map.shape))
    return data, header, wcs, seq_map


def scanids_in_box(scan_map, pixbox):
    """Unique SEQ_IDs contributing inside the pixel box (x0, x1, y0, y1).
    scan_map is indexed [y, x] to match the spatial grid."""
    x0, x1, y0, y1 = pixbox
    ids = np.unique(scan_map[y0:y1, x0:x1])
    return ids[ids >= 0]        # drop padding / no-data values


# ----------------------------------------------------------------------
# 2. Sub-cube extraction (small sky boxes)
# ----------------------------------------------------------------------

def skybox_to_pixbox(wcs, region, shape2d):
    """Convert the region center + ANGULAR radius to a rectangular
    pixel box (x0, x1, y0, y1), snapped inside the cube. The angular
    size is divided by the pixel scale of the cube's grid, so the
    same region definition works for any grid resolution."""
    # spatial SkyCoord in the cube's frame
    center = SkyCoord(region.lon, region.lat, frame=region.frame)
    # assumes cube world axes: lon, lat, velocity -> adapt if order differs
    xp, yp = wcs.celestial.world_to_pixel(center)
    nx, ny = shape2d
    # convert each pixel scale to a plain float FIRST, then abs:
    # np.abs() on the Quantity array raises
    # "only dimensionless scalar quantities can be converted"
    cd = wcs.celestial.proj_plane_pixel_scales()
    cd0 = abs(float(cd[0].to_value(u.deg)))
    cd1 = abs(float(cd[1].to_value(u.deg)))
    sx = int(2 * region.radius.to_value(u.deg) / cd0)
    sy = int(2 * region.radius.to_value(u.deg) / cd1)
    x0 = max(0, int(round(xp - sx / 2))); x1 = min(nx, x0 + sx)
    y0 = max(0, int(round(yp - sy / 2))); y1 = min(ny, y0 + sy)
    if (x1 - x0) < sx or (y1 - y0) < sy:
        log.warning("Region %s clipped at cube edge", region.name)
    return (x0, x1, y0, y1)


def vel_to_slice(wcs, vel_range):
    """km/s range -> (i0, i1) spectral pixel indices (i1 exclusive).
    Handles VRAD/VOPT/FREQ CTYPE3 and the sign of CDELT3. The given
    range is NORMALIZED to (min, max): a reversed range like (18, 10)
    is treated as 10..18 km/s (with a warning), not as an empty slice.
    Ranges falling outside the cube's velocity span are clipped and
    warned about; a fully outside range raises instead of silently
    producing an empty (all-NaN/zero) image."""
    orig = (float(vel_range[0]), float(vel_range[1]))
    vlo, vhi = sorted(orig)
    if (vlo, vhi) != orig:
        log.warning("velocity range %s reversed, using %.1f..%.1f km/s",
                    orig, vlo, vhi)
    # TODO: handle FREQ via equivalencies (u.spectral)
    with u.set_enabled_equivalencies(u.spectral()):
        v = u.Quantity([vlo, vhi], u.km / u.s)
        ip = wcs.spectral.world_to_pixel(v)
    lo, hi = int(np.floor(min(ip))), int(np.ceil(max(ip)))
    lo = max(lo, 0)
    hi = min(hi, wcs.spectral.array_shape[0])
    if hi <= lo:
        raise ValueError(
            f"velocity range {vlo}..{vhi} km/s has no overlap with the "
            f"cube's spectral axis (pixel indices {min(ip):.1f}, "
            f"{max(ip):.1f})")
    if (lo, hi) != (int(np.floor(min(ip))), int(np.ceil(max(ip))) + 1):
        log.warning("velocity range %.1f..%.1f km/s partially outside "
                    "the cube; clipped to pixels %d:%d", vlo, vhi, lo, hi)
    return (lo, hi)


def extract_subcube(data, wcs, pixbox, vslice):
    """Slice the 3D cube and return (subcube, subwcs).

    subwcs keeps valid world coordinates: copy the WCS and shift
    CRPIX by the box origin. Assumes spatial axes are (x, y) = last
    two array axes, spectral axis first; adapt to actual layout."""
    x0, x1, y0, y1 = pixbox
    i0, i1 = vslice
    sub = data[i0:i1, y0:y1, x0:x1]           # <- adjust axis order to data layout
    subwcs = wcs.deepcopy()
    subwcs.wcs.crpix[0] -= x0
    subwcs.wcs.crpix[1] -= y0
    return sub, subwcs


def collapse(subcube, method="moment0"):
    """Collapse the velocity axis of the sub-cube to a 2D sub-image."""
    if method == "moment0":
        # moment-0: sum * |dv| (skip NaNs)
        img = np.nansum(subcube, axis=0)
    elif method == "mean":
        img = np.nanmean(subcube, axis=0)
    else:
        raise ValueError(f"unknown collapse method: {method}")
    return img


def preprocess_image(img, opts):
    """Clean a collapsed 2D image BEFORE cross-correlation.

    Three optional steps, applied IDENTICALLY to every mixer image
    (same recipe, same parameters), so they suppress noise without
    biasing the measured offsets:

    1. OUTLIER REJECTION (outlier_clip, robust sigma; 0 = off):
       MAD-based sigma clip around the image median; outlier pixels
       (spikes, bad detector samples) are replaced by the local
       3x3 median.
    2. DESTRIPING (destripe_kernel px, destripe_axis 'x' or 'y';
       0 = off): subtract a background estimated by a MEDIAN FILTER
       ALONG destripe_axis with width destripe_kernel:
       out = img - median_filter(img, size). The filter follows
       structure that is EXTENDED along its direction, so stripes
       running PARALLEL to destripe_axis are estimated and
       subtracted away; compact sources (narrower than the kernel
       along that direction) do not enter the median and SURVIVE
       the subtraction. Rule: set destripe_axis to the direction
       the stripes RUN along (scan direction), NOT across them.
       CRITICAL WHY: the stripes otherwise sit in the correlation
       surface as spurious bumps that can beat a weak real peak
       (random lock-on offsets, not zero-shift -- the stripes are
       per-mixer, not common-mode).
    3. SMOOTHING (smooth_fwhm_px; 0 = off): small Gaussian; raises
       the correlation SNR, broadens (never shifts) the peak.

    Config (options:): outlier_clip: 5.0, destripe_kernel: 25,
    destripe_axis: "x", smooth_fwhm_px: 0 -- tune per band; NII
    typically benefits from all three. NaNs are filled before
    filtering and re-masked after, so coverage holes stay holes.
    """
    from scipy import ndimage
    if not np.any(np.isfinite(img)):
        return img
    out = np.array(img, dtype=float)
    mask = ~np.isfinite(out)
    fill = float(np.nanmedian(out))
    out[mask] = fill

    clip = float(opts.get("outlier_clip", 0.0))
    if clip > 0:
        med = np.median(out)
        mad = np.median(np.abs(out - med)) * 1.4826
        if mad > 0:
            bad = np.abs(out - med) > clip * mad
            if np.any(bad):
                local = ndimage.median_filter(out, size=3)
                out[bad] = local[bad]

    k = int(opts.get("destripe_kernel", 0))
    if k > 2:
        axis = str(opts.get("destripe_axis", "x")).lower()
        size = (1, k) if axis == "x" else (k, 1)
        out = out - ndimage.median_filter(out, size=size,
                                          mode="nearest")

    fwhm = float(opts.get("smooth_fwhm_px", 0.0))
    if fwhm > 0:
        out = ndimage.gaussian_filter(out, fwhm / 2.3548,
                                      mode="nearest")

    out[mask] = np.nan
    log.info("preprocess: outlier_clip=%s destripe=%spx/%s "
             "smooth=%spx  (0 = step off)",
             clip if clip else "off",
             k if k else "off",
             str(opts.get("destripe_axis", "x")).lower(),
             fwhm if fwhm else "off")
    return out


def baseline_correct(subcube, vel_ranges, wcs, pix0, min_pts=6):
    """Per-pixel SPECTRAL baseline correction -- the FIRST image step.

    The subcube must span the UNION of the region's velocity ranges
    EXTENDED by `baseline_pad` spectral elements on each side (the
    caller cuts it that way; see main). LINE channels -- inside ANY
    of the region's vel_ranges, overlapping/adjacent ranges combine
    into one exclusion window -- are EXCLUDED from the fit. A LINEAR
    (deg-1) baseline is fit per SPATIAL PIXEL on the remaining
    continuum channels and subtracted over ALL channels
    (extrapolated across the line window).

    Pathological cubes (NII: data dropouts limiting the valid
    spectral coverage per pixel): fall back to an OFFSET (deg 0)
    when fewer than `min_pts` valid continuum channels remain;
    pixels with < 2 valid points are left UNCORRECTED (counted and
    logged). vel_ranges are km/s; pix0 is the absolute spectral
    index of the subcube's first channel, so the line window is
    located in world velocity, robust against per-cube grid
    offsets."""
    nch, ny, nx = subcube.shape
    with u.set_enabled_equivalencies(u.spectral()):
        v = wcs.spectral.pixel_to_world(
            np.arange(pix0, pix0 + nch)).to_value(u.km / u.s)
    line = np.zeros(nch, dtype=bool)
    for lo, hi in vel_ranges:
        line |= (v >= min(lo, hi)) & (v <= max(lo, hi))
    base = np.where(~line)[0]
    if base.size < 2:
        log.warning("baseline: no continuum channels outside the "
                    "line window (pad too small?) -- NOT correcting")
        return subcube
    out = np.array(subcube, dtype=float)
    n_lin = n_off = n_skip = 0
    for iy in range(ny):
        for ix in range(nx):
            spec = out[:, iy, ix]
            good = base[np.isfinite(spec[base])]
            if good.size >= min_pts:
                x = v[good]
                A = np.vstack([x, np.ones_like(x)]).T
                sol, *_ = np.linalg.lstsq(A, spec[good], rcond=None)
                out[:, iy, ix] -= sol[0] * v + sol[1]
                n_lin += 1
            elif good.size >= 2:
                out[:, iy, ix] -= float(np.mean(spec[good]))
                n_off += 1
            else:
                n_skip += 1
    log.info("baseline: %d px linear, %d px offset-only, %d px "
             "skipped (%d line / %d continuum channels)",
             n_lin, n_off, n_skip, int(line.sum()), base.size)
    return out


def snr_gate(img, min_snr):
    """Rough SNR gate on the sub-image. Returns (ok, snr).

    Detects degenerate images explicitly: all-NaN (coverage holes or
    an empty slice) or a constant image give a clear warning instead
    of a silent SNR 0.0 skip."""
    if not np.any(np.isfinite(img)):
        log.warning("SNR gate: image is entirely non-finite "
                    "(empty slice or no coverage)")
        return False, 0.0
    med = np.nanmedian(img)
    mad = 1.4826 * np.nanmedian(np.abs(img - med))
    peak = np.nanmax(img)
    if mad <= 0:
        log.warning("SNR gate: image is (near-)constant, MAD=0 -- "
                    "likely an empty or bad slice")
        return False, 0.0
    snr = (peak - med) / mad
    return (snr >= min_snr), snr


# ----------------------------------------------------------------------
# 3. Cross-correlation and 2D peak fit (bayesicfitting)
# ----------------------------------------------------------------------

def fft_xcorr(img_a, img_b):
    """FFT cross-correlation of two equal-shaped images.
    Returns (corr, shift_axes) where corr is the correlation surface
    and integer peak = argmax. Zero-shift is at the center index.
    NaNs are filled with the image median before correlating."""
    def prep(img):
        img = np.nan_to_num(img, nan=np.nanmedian(img))
        return img - np.nanmean(img)
    a, b = prep(img_a), prep(img_b)
    ny, nx = a.shape
    corr = np.fft.irfft2(np.fft.rfft2(a) * np.conj(np.fft.rfft2(b)),
                         s=a.shape)
    # circular shift so zero lag is at (ny//2, nx//2)
    corr = np.roll(corr, (ny // 2, nx // 2), axis=(0, 1))
    return corr


def integer_peak(corr):
    """Coarse (integer) peak, as pixel shift relative to zero lag."""
    ny, nx = corr.shape
    j, i = np.unravel_index(np.argmax(corr), corr.shape)
    dy = j - ny // 2
    dx = i - nx // 2
    return (dy, dx), corr[j, i]


def fit_peak(corr, int_peak, patch=7):
    """Subpixel 2D elliptical-Gaussian + constant fit around the
    integer peak, with bayesicfitting. Returns (dx, dy, sig_dx,
    sig_dy, rho, fit_result) or None if the fit fails.

    Model: Kernel2dModel(kernel=Gauss(), shape='E') -- an elliptical
    Gaussian kernel -- summed with PolySurfaceModel(0), a CONSTANT
    offset level for the correlation surface. (A linear surface
    (degree 1) proved too flexible: its tilt dragged the fitted peak
    away from the correlation maximum.)
    """
    ny, nx = corr.shape
    dy, dx = int_peak
    jc, ic = ny // 2 + dy, nx // 2 + dx
    x0w, x1w = max(0, ic - patch // 2), min(nx, ic - patch // 2 + patch)
    y0w, y1w = max(0, jc - patch // 2), min(ny, jc - patch // 2 + patch)
    win = corr[y0w:y1w, x0w:x1w]
    if not np.all(np.isfinite(win)):
        win = np.nan_to_num(win, nan=0.0)

    # fit domain as 1D flattened (x, y) pairs, shifts relative to
    # zero lag at (nx//2, ny//2)
    xx, yy = np.meshgrid(np.arange(x0w, x1w) - nx // 2,
                         np.arange(y0w, y1w) - ny // 2)
    x = np.stack([xx.ravel(), yy.ravel()], axis=1)
    yd = win.ravel()

    # elliptical Gaussian kernel model + constant offset
    # NOTE: kernel takes a Gauss() INSTANCE, not the Gauss class
    kernel = Kernel2dModel(kernel=Gauss(), shape='E')
    model = kernel + PolySurfaceModel(0)
    # initial guesses: Kernel2dModel parameters (amplitude, xcenter,
    # ycenter, width1, width2 -- 5 params) + PolySurfaceModel(0)
    # (offset -- 1 param) = 6 total; offset = window median
    model.parameters = [np.nanmax(win), dx, dy, 2.0, 2.0,
                        float(np.nanmedian(win))]

    fitter = LevenbergMarquardtFitter(x, model)
    try:
        fitter.fit(yd)
    except Exception as exc:
        log.warning("peak fit failed: %s", exc)
        return None
    p = fitter.parameters
    log.info("peak fit params [amp, x0, y0, w1, w2, off] = %s",
             np.array2string(np.asarray(p), precision=3))
    # parameter covariance: fitter.covariance (confirmed present).
    # peak position = parameters [1] (xcenter) and [2] (ycenter) --
    # verify the Kernel2dModel + PolySurfaceModel parameter order in
    # the installed version.
    cov_all = np.asarray(fitter.covariance)
    cov = cov_all[1:3, 1:3]          # 2x2 covariance of (x0, y0)
    sig_dx = np.sqrt(cov[0, 0])
    sig_dy = np.sqrt(cov[1, 1])
    rho = (cov[0, 1] / (sig_dx * sig_dy)
           if (sig_dx > 0 and sig_dy > 0) else 0.0)
    # sanity: fitted peak must stay near the integer peak, else the
    # fit slid off-peak (bad initial widths, clipped window, or the
    # surface absorbing the Gaussian) -- warn loudly
    if abs(p[1] - dx) > patch / 2 or abs(p[2] - dy) > patch / 2:
        log.warning("fitted peak (%.2f, %.2f) far from integer peak "
                    "(%d, %d) -- fit may have slid off-peak",
                    p[1], p[2], dx, dy)
    # sanity: sigma of the peak position much larger than the window
    # means the covariance is inflated (bad fit or unweighted LM);
    # the error ellipse would dwarf the correlation peak
    if max(sig_dx, sig_dy) > patch / 2:
        log.warning("peak-position sigma (%.2f, %.2f) larger than "
                    "patch/2 (%.1f) -- covariance looks inflated",
                    sig_dx, sig_dy, patch / 2)
    return (float(p[1]), float(p[2]), float(sig_dx), float(sig_dy),
            float(rho), fitter)


def peak_snr(corr, int_peak):
    """Peak significance: peak / rms of the correlation away from the
    peak (mask +- 3 pixels around it)."""
    ny, nx = corr.shape
    j, i = ny // 2 + int_peak[0], nx // 2 + int_peak[1]
    mask = np.ones_like(corr, dtype=bool)
    mask[max(0, j - 3):j + 4, max(0, i - 3):i + 4] = False
    rms = np.std(corr[mask])
    return corr[j, i] / rms if rms > 0 else 0.0


# ----------------------------------------------------------------------
# 4. Sky -> AltAz -> focal plane conversion (astropy)
# ----------------------------------------------------------------------

def calculate_az_el_offsets(observer_time, observer_lat, observer_lon,
                            observer_alt, gal_lon, gal_lat,
                            ref_lon=351.0, ref_lat=0.5):
    """Tangential AZ/EL offsets (deg) of a Galactic position relative
    to a reference (boresight) position, at a given epoch.

    Uses AltAz + spherical_offsets_to: handles the cos(elevation)
    factor in the azimuth difference correctly (a raw AZ difference
    shrinks on the sky by cos(el); spherical_offsets_to does this
    properly and stays well-behaved away from the zenith).
    """
    location = EarthLocation(lat=observer_lat * u.deg,
                             lon=observer_lon * u.deg,
                             height=observer_alt * u.m)
    frame = AltAz(location=location, obstime=Time(observer_time))

    target = SkyCoord(l=gal_lon * u.deg, b=gal_lat * u.deg,
                      frame='galactic').transform_to(frame)
    reference = SkyCoord(l=ref_lon * u.deg, b=ref_lat * u.deg,
                         frame='galactic').transform_to(frame)

    daz, dalt = reference.spherical_offsets_to(target)
    return daz.to(u.deg).value, dalt.to(u.deg).value


def galactic_to_altaz_matrix(observer_time, location, ref_lon, ref_lat,
                             step=1 * u.arcsec):
    """2x2 rotation matrix: Galactic tangential (dlon, dlat) -> AltAz
    tangential (daz, del) at the reference position, at this epoch.

    Two probe offsets (+step along each Galactic tangential axis) form
    the numerical Jacobian; this replaces hand-derived parallactic
    -angle algebra and absorbs all frame conventions. Use the boresight
    or region center as the reference; the rotation varies by ~0.01 deg
    over a 0.5 deg box, so either is fine.
    """
    frame = AltAz(location=location, obstime=Time(observer_time))
    ref = SkyCoord(l=ref_lon * u.deg, b=ref_lat * u.deg, frame='galactic')
    ref_altaz = ref.transform_to(frame)

    step_as = step.to_value(u.arcsec)   # q.arcsec shorthand REMOVED in
    M = np.zeros((2, 2))                # astropy >= 7: use to_value()
    for j, probe in enumerate([(step, 0 * u.arcsec),
                               (0 * u.arcsec, step)]):
        p = SkyCoord(ref.spherical_offsets_by(*probe))
        daz, dalt = ref_altaz.spherical_offsets_to(
            p.transform_to(frame))
        M[:, j] = [daz.to_value(u.arcsec), dalt.to_value(u.arcsec)]
    M /= step_as
    return M


def skyoffset_to_focal(dlon, dlat, M_frame, roll):
    """Convert tangential sky offsets to focal-plane offsets.

    M_frame: Galactic->AltAz rotation matrix from
    galactic_to_altaz_matrix at the region/boresight position and
    scan epoch. roll: gondola roll about the boresight (radian).

    Returns (f_dx, f_dy, R) with R the total 2x2 rotation matrix
    (frame rotation * roll), for covariance propagation.
    """
    c, s = np.cos(roll), np.sin(roll)
    R_roll = np.array([[c, s], [-s, c]])
    R = R_roll @ M_frame
    f = R @ np.array([dlon.to(u.arcsec).value, dlat.to(u.arcsec).value])
    return f[0] * u.arcsec, f[1] * u.arcsec, R


def cov_rotate(C, R):
    """Propagate the full 2x2 covariance through rotation R."""
    return R @ C @ R.T


def seqid_at_center(seq_map, pixbox, wcs_celestial, region):
    """SEQ_ID at the region center.

    The region center sky position identifies the pixel that contains
    it; the SEQ_ID extension (HDU 3) of the mixer cube records, at
    each pixel, the Level-1 sequence that contributed there. The
    center pixel's SEQ_ID identifies the Level-1 product whose primary
    FITS header holds the gondola attitude for the region.
    """
    center = SkyCoord(region.lon, region.lat, frame=region.frame)
    xp, yp = wcs_celestial.world_to_pixel(center)
    x0, x1, y0, y1 = pixbox
    ix = int(np.clip(round(float(xp)) - x0, 0, x1 - x0 - 1))
    iy = int(np.clip(round(float(yp)) - y0, 0, y1 - y0 - 1))
    return int(seq_map[iy, ix])


def find_l1_product(l1_dir, band, seq_id):
    """Locate the Level-1 product for a SEQ_ID.

    Filename convention: BAND_{seq}_{scanid}_L10.fits, e.g.
    NII_10097_24232_L10.fits -- first number is the unique sequence
    id (a calibratable unit: OTF scan + REF/HOT), second is the
    SCANID. The gridder records SEQ_ID, which matches the FIRST
    filename field, so the lookup is EXACT.

    IMPORTANT: the SEQ_ID in the filename is zero-padded to 5
    digits (e.g. seq 42 -> '00042'). The SEQ_ID extension in the
    cube stores a plain integer, so it must be formatted with
    leading zeros before globbing.
    """
    seq_str = f"{int(seq_id):05d}"
    matches = sorted(Path(l1_dir).glob(f"{band}_{seq_str}_*_L10.fits"))
    if len(matches) == 1:
        return matches[0]
    if not matches:
        # fall back to the unpadded form, just in case the
        # convention is mixed across products
        matches = sorted(Path(l1_dir).glob(f"{band}_{int(seq_id)}_*_L10.fits"))
    if len(matches) == 1:
        return matches[0]
    if not matches:
        raise FileNotFoundError(f"no L1 product for {band} sequence "
                                f"{seq_id} (tried {seq_str}) in {l1_dir}")
    raise RuntimeError(f"ambiguous L1 products for {band} sequence "
                       f"{seq_id}: {[m.name for m in matches]}")


def get_scan_attitude(l1_dir, band, seq_id):
    """Gondola POSITION from the Level-1 product -- the key to the
    focal-plane conversion for a FLOATING telescope.

    The gondola flew a circumpolar (near-)orbit around the South
    Pole: over the mission its position swept through 360 deg of
    longitude. The AltAz frame for the Galactic->focal rotation is
    defined by WHERE and WHEN the telescope observed, so the
    location must come from this scan's L1 product:

      GON_LAT  : gondola geographic latitude (deg)
      GON_LON  : gondola geographic longitude (deg)
      GON_ALT  : gondola altitude (m)
      epoch    : median of the UNIXTIME column of the scan

    Returns (lat, lon, alt, obstime, az, el) -- az/el (boresight
    Az/El) are OPTIONAL and only used as a pointing sanity check;
    the offset math needs location + epoch + target position only.
    """
    l1_path = find_l1_product(l1_dir, band, seq_id)
    hdr = fits.getheader(l1_path, 0)
    lat = hdr["GON_LAT"] * u.deg
    lon = hdr["GON_LON"] * u.deg
    alt = hdr["GON_ALT"] * u.m
    # optional boresight az/el (pointing sanity check only)
    az = hdr.get("GON_AZ") or hdr.get("AZ")
    el = hdr.get("GON_EL") or hdr.get("EL")
    if az is not None and el is not None:
        az = az * u.deg
        el = el * u.deg
    else:
        log.warning("no boresight AZ/EL keywords in %s -- continuing "
                     "(not needed for the offset math)", l1_path.name)
    # unixtime is a COLUMN in the level-1 data (table extension):
    # median of the scan is sufficient (legs < 3 min)
    with fits.open(l1_path) as hdul:
        tdata = None
        for h in hdul[1:]:
            if h.data is not None and "UNIXTIME" in h.data.names:
                tdata = h.data["UNIXTIME"]
                break
        if tdata is None:
            raise KeyError(f"no UNIXTIME column in {l1_path}")
        obstime = Time(np.median(tdata), format='unix')
    return lat, lon, alt, obstime, az, el


# NOTES_ON_ROLL
# -------------
# True model: d_sky(t) = M_frame(t) . R(theta0) . f, where
# M_frame(t) is the Galactic->AltAz rotation (from astropy, given
# GON_LON/GON_LAT + obstime) and theta0 is the constant focal-plane
# roll about the boresight (not reported by the gondola).
# - With theta0 = 0 the derived offsets live in a "sky-defined focal
#   frame": internally consistent, correct magnitudes, but rotated
#   by theta0 w.r.t. the physical instrument axes. Fine if only
#   relative offsets within the array matter.
# - theta0 is SOLVABLE from the data: M_frame(t) varies between
#   epochs, so a wrong theta0 causes epoch-to-epoch scatter of the
#   focal offsets, while the true theta0 makes all epochs agree.
#   Solve by minimizing epoch scatter over theta0 (1-parameter scan
#   over the weighted-mean chi2 of all mixer offsets).


# ----------------------------------------------------------------------
# 5. Combination: weighted mean + chi2 consistency
# ----------------------------------------------------------------------

def weighted_mean(xs, sigs):
    """Inverse-variance weighted mean. Returns (mean, sig_int, chi2, dof)."""
    w = 1.0 / np.asarray(sigs) ** 2
    m = np.sum(w * np.asarray(xs)) / np.sum(w)
    sig_int = 1.0 / np.sqrt(np.sum(w))
    chi2 = np.sum(w * (np.asarray(xs) - m) ** 2)
    dof = len(xs) - 1
    return m, sig_int, chi2, dof


def scaled_uncertainty(sig_int, chi2, dof):
    """PDG-style: never smaller than the internal estimate."""
    if dof <= 0:
        return sig_int
    return sig_int * np.sqrt(max(1.0, chi2 / dof))


_SYS_FLOOR_DEG = 0.005     # systematic floor on per-measurement
                           # sigmas (deg), added in quadrature in
                           # combine_measurements; overruled by the
                           # config option sys_floor_deg in
                           # load_config. Rationale: the formal LM
                           # sigmas are ~1 order optimistic (patch
                           # pixels are correlated AND the Gaussian
                           # model-shape error is systematic); without
                           # a floor, targets with ONE measurement
                           # (chi2 undefined, no PDG inflation) and
                           # targets whose slices "agree" within
                           # optimistic errors keep absurdly small
                           # sigmas and dominate every weighted mean.
_SYS_FLOOR_BY_BAND = {}    # per-band override, set in load_config
                           # from bands.<BAND>.sys_floor_deg -- e.g.
                           # NII typically needs a larger floor than
                           # CII (noisier data, weaker lines)


def _floor_of(row):
    """systematic floor applying to one measurement row (deg):
    its band's override when configured, else the global floor."""
    band = str(row.get("band", "")) if "band" in row else ""
    return _SYS_FLOOR_BY_BAND.get(band, _SYS_FLOOR_DEG)


def combine_measurements(rows, floor_cap=False):
    """Combine all measurement rows of one mixer (both axes
    independently). Returns summary dict.

    Per mixer also the inverse-variance weighted mean sky offset
    from the fiducial in DEGREES (dlon_deg/dlat_deg + sigma), for
    calibration use. Each input sigma is inflated in quadrature by
    the systematic floor (_SYS_FLOOR_DEG) BEFORE the weighted mean,
    so chi2 and the PDG rescaling operate on honest errors: with
    floor >> formal sigma, a single-slice target's sigma IS the
    floor (no more fake chi2=0 dominance), and multi-slice chi2
    measures real disagreement.

    floor_cap=True for combines WITHIN one target/epoch (per-run
    summaries, per-region rows, per-target summaries): overlapping
    or continuous velocity ranges make those slices CORRELATED, so
    averaging cannot beat the systematic floor -- the mean's sigma
    is capped at the floor (it otherwise shrinks as floor/sqrt(n),
    which is exactly what drives chi2/dof << 1 for correlated
    near-duplicate slices). Cross-target combines (floor_cap=False,
    default) keep the 1/sqrt(N) shrinkage: different fields/epochs
    are genuinely independent."""
    out = {}
    # focal-plane offsets (fx, fy) in DEGREES: the calibration
    # quantities. With roll=0, fx ~ azimuthal angular offset,
    # fy ~ elevation angular offset from the fiducial.
    for axis in ("fx", "fy"):
        xs = [r[axis].to_value(u.deg) for r in rows
              if r["sig_" + axis].value > 0]
        ss = [np.hypot(r["sig_" + axis].to_value(u.deg),
                       _floor_of(r))
              for r in rows if r["sig_" + axis].value > 0]
        if not xs:
            continue
        m, sig_int, chi2, dof = weighted_mean(xs, ss)
        if dof == 0:
            # single-measurement combine: chi2 undefined; zero it
            # EXACTLY (weighted_mean returns ~1e-32 float noise)
            chi2 = 0.0
        sig = scaled_uncertainty(sig_int, chi2, dof)
        if floor_cap:
            # correlated slices share the systematic: the mean is
            # no better known than one measurement
            sig = max(sig, _floor_of(rows[0]))
        out[axis] = m * u.deg
        out["sig_" + axis] = sig * u.deg
        out["chi2_" + axis] = chi2
        out["dof_" + axis] = dof
    out["n_meas"] = len(rows)
    return out


# ----------------------------------------------------------------------
# 6. QA plots
# ----------------------------------------------------------------------

def plot_region_montage(images, labels, extent, region,
                        vel_range, outfile, snr=None):
    """QC montage of the region images BEFORE cross-correlation.

    Panels: one per image (fiducial first) on a shared color scale +
    a difference panel (mixer - fiducial). Colorbars, region center
    marker, and SNR in the title. Written before the x-corr runs so
    problems (noisy region, coverage holes, wrong velocity range,
    structure mismatch between mixers) are visible immediately.

    extent: world-coordinate box from the ORIGINAL cube WCS, in
    the cube's NATIVE frame, plus the frame name (montage_extent
    returns a 5-tuple). The region center marker is drawn at the
    region's world coordinates TRANSFORMED INTO THAT FRAME -- e.g.
    an ICRS-defined region against a Galactic cube is transformed
    before plotting, so marker and image always agree.
    """
    import matplotlib.pyplot as plt
    n = len(images)
    fig, axes = plt.subplots(1, n + 1, figsize=(4 * (n + 1), 4))
    axes = np.atleast_1d(axes).ravel()
    # imshow wants a 4-tuple; the 5th element (frame name) is for
    # markers and axis labels
    ext4 = tuple(extent[:4])
    vmax = np.nanpercentile(np.nanmax(images, axis=0), 99)
    # region center in the MONTAGE frame (cube-native frame)
    frame = extent[4] if len(extent) == 5 else region.frame
    rc = SkyCoord(region.lon, region.lat, frame=region.frame)
    rc = rc.galactic if frame == "galactic" else rc.icrs
    cx = float(rc.l.to_value(u.deg)) if frame == "galactic" \
        else float(rc.ra.to_value(u.deg))
    cy = float(rc.b.to_value(u.deg)) if frame == "galactic" \
        else float(rc.dec.to_value(u.deg))
    # WRAP: the extent may be display-normalized (a box straddling
    # l=0 / RA=0 shows as negative longitudes, montage_extent), so
    # shift the marker into the SAME convention -- else it lands
    # 360 deg away from the panel (marker "missing" near the wrap)
    while cx < ext4[0] - 180.0:
        cx += 360.0
    while cx > ext4[1] + 180.0:
        cx -= 360.0
    lonlab = "Gal. lon (deg)" if frame == "galactic" else "RA (deg)"
    latlab = "Gal. lat (deg)" if frame == "galactic" else "Dec (deg)"
    from matplotlib.ticker import MaxNLocator
    for ax, im, lab in zip(axes, images, labels):
        pcm = ax.imshow(im, origin="lower", extent=ext4, cmap="inferno",
                        vmin=0, vmax=vmax)
        ax.plot(cx, cy, "w+", ms=12, mew=1.5)
        ax.set_title(f"Mixer {lab}", fontsize=10)
        ax.set_xlabel(lonlab)
        ax.set_ylabel(latlab)
        # few, well-spaced ticks: lon values are long (e.g. 348.63)
        # and overlap at default tick density
        ax.xaxis.set_major_locator(MaxNLocator(3))
        ax.yaxis.set_major_locator(MaxNLocator(4))
        fig.colorbar(pcm, ax=ax, shrink=0.8)
    # difference panel: measured - fiducial (first image)
    diff = images[-1] - images[0]
    ax = axes[n]
    dmax = np.nanpercentile(np.abs(diff), 99)
    pcm = ax.imshow(diff, origin="lower", extent=ext4, cmap="RdBu_r",
                    vmin=-dmax, vmax=dmax)
    ax.plot(cx, cy, "k+", ms=12, mew=1.5)
    ax.set_title(f"Difference (Mixer {labels[-1]} - "
                 f"Mixer {labels[0]})", fontsize=10)
    ax.set_xlabel(lonlab)
    ax.xaxis.set_major_locator(MaxNLocator(3))
    fig.colorbar(pcm, ax=ax, shrink=0.8)
    title = f"{region.name}  v = {vel_range} km/s"
    if snr is not None:
        title += f"   SNR = {snr:.1f}"
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(outfile, dpi=120)
    plt.close(fig)
    log.info("wrote %s", outfile)


def plot_xcorr_panel(corr, int_peak, fit, subwcs, label, outfile):
    """Correlation surface + integer peak + fitted (subpixel) peak.

    Shows ONLY the sigma fit contours: the 1-, 2- and 3-sigma
    confidence contours of the fitted peak position (from the
    parameter covariance), centered on the fitted peak. A zoom inset
    around the peak makes the (usually small) contours visible."""
    import matplotlib.pyplot as plt
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))
    ny, nx = corr.shape
    extent = (-(nx // 2), nx - nx // 2, -(ny // 2), ny - ny // 2)
    im = ax1.imshow(corr, origin="lower", extent=extent, cmap="viridis")
    ax1.plot(int_peak[1], int_peak[0], "rx", ms=10, label="integer")
    if fit is not None:
        ax1.plot(fit[0], fit[1], "w+", ms=12, label="fit")
        # sigma fit contours: 1-, 2-, 3-sigma confidence contours of
        # the fitted peak position, from (sig_dx, sig_dy, rho)
        sigx, sigy, rho = fit[2], fit[3], fit[4]
        cov = np.array([[sigx ** 2, rho * sigx * sigy],
                        [rho * sigx * sigy, sigy ** 2]])
        vals, vecs = np.linalg.eigh(cov)
        ang = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
        from matplotlib.patches import Ellipse
        for nsig, lw in ((1, 1.4), (2, 0.8), (3, 0.5)):
            ax1.add_patch(Ellipse((fit[0], fit[1]),
                                  2 * nsig * np.sqrt(vals[0]),
                                  2 * nsig * np.sqrt(vals[1]),
                                  angle=ang, fill=False, color="w",
                                  lw=lw))
        # zoom inset around the fitted peak, same sigma contours
        pad = max(3.0, 4 * max(np.sqrt(vals)))
        axi = ax1.inset_axes([0.62, 0.08, 0.35, 0.35])
        axi.imshow(corr, origin="lower", extent=extent, cmap="viridis")
        axi.plot(int_peak[1], int_peak[0], "rx", ms=7)
        axi.plot(fit[0], fit[1], "w+", ms=9)
        for nsig, lw in ((1, 1.4), (2, 0.8), (3, 0.5)):
            axi.add_patch(Ellipse((fit[0], fit[1]),
                                  2 * nsig * np.sqrt(vals[0]),
                                  2 * nsig * np.sqrt(vals[1]),
                                  angle=ang, fill=False, color="w",
                                  lw=lw))
        axi.set_xlim(fit[0] - pad, fit[0] + pad)
        axi.set_ylim(fit[1] - pad, fit[1] + pad)
        axi.set_xticks([])
        axi.set_yticks([])
        ax1.indicate_inset_zoom(axi, edgecolor="w")
        # annotate the actual numbers so a huge/tiny ellipse is
        # diagnosable from the plot alone
        ax1.text(0.02, 0.02,
                 f"$x_0={fit[0]:.2f}\\pm{fit[2]:.2f}$  "
                 f"$y_0={fit[1]:.2f}\\pm{fit[3]:.2f}$  "
                 f"$\\rho={fit[4]:.2f}$",
                 transform=ax1.transAxes, color="w", fontsize=8,
                 va="bottom")
    ax1.legend()
    ax1.set_title(f"x-corr: {label}")
    fig.colorbar(im, ax=ax1)
    ax2.hist(corr.ravel(), bins=50)
    ax2.set_title("correlation values")
    fig.tight_layout()
    fig.savefig(outfile, dpi=120)
    plt.close(fig)
    log.info("wrote %s", outfile)


def read_calibration(path):
    """Read the [offsets] calibration file -> ordered list of rows.

    Format:
        [offsets] version = 0.5
        PIX	AZ		EL		COMMENT
        B2M8	0.000000	0.000000	FIDUCIAL
        B2M5	0.031436	-0.032198	THEORY
        B2M5	0.040013	-0.033674	AS_MEASURED

    NOTE: the SAME PIX can appear on MULTIPLE rows (e.g. a THEORY
    row and an AS_MEASURED row) -- so the result is an ordered LIST
    of (pix, {'az', 'el', 'comment'}), NOT a dict. The downstream
    pipeline reads all rows but only uses the AS_MEASURED ones.
    """
    rows = []
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("[") or line.startswith("PIX"):
            continue
        parts = line.split()
        if len(parts) < 3:
            log.warning("skipping unparseable calibration line: %r",
                        line)
            continue
        rows.append((parts[0],
                     {"az": float(parts[1]), "el": float(parts[2]),
                      "comment": parts[3] if len(parts) > 3 else ""}))
    return rows


def update_calibration(cal_rows, summaries, fid_label, pix_band,
                       drop=()):
    """Fold the xcorr-measured RESIDUAL offsets into the calibration.

    The L2 cubes were built from L1 data with the calibration-file
    offsets already APPLIED, so each xcorr measurement is a RESIDUAL
    (measured-minus-assumed). The updated AS_MEASURED offset is:

        az_new = az_applied + fx        (both decimal degrees)
        el_new = el_applied + fy

    FILE POLICY (per the pipeline's convention):
    - THEORY rows are KEPT UNCHANGED. The downstream pipeline reads
      all rows but only ever USES the AS_MEASURED ones; theory
      values stay in the file as reference/documentation.
    - For each mixer with a current measurement: any OLD
      AS_MEASURED row is dropped and a new AS_MEASURED row is
      appended (az/el = applied value + measured residual).
    - FIDUCIAL and other rows pass through unchanged.
    - `drop`: PIX names whose AS_MEASURED rows are REMOVED even
      without a current measurement -- the escape hatch for stale
      rows that crept into an earlier version (e.g. a fiducial
      measured against itself and written out as AS_MEASURED);
      config key 'cal_drop:'. THEORY/FIDUCIAL rows of that PIX
      pass through: the pixel reverts to its theory value.

    cal_rows: ordered list of (pix, {az, el, comment}) from
    read_calibration_rows; returns the new ordered list.
    """
    # measured residuals keyed by pix name
    meas = {}
    for label, s in summaries.items():
        if label == fid_label or "fx" not in s or "fy" not in s:
            continue
        meas[f"{pix_band}M{label}"] = (s["fx"].to_value(u.deg),
                                       s["fy"].to_value(u.deg))
    # applied values (the offsets the L1 was made with) come from
    # the THEORY rows; fall back to the last AS_MEASURED row
    applied = {}
    for pix, c in cal_rows:
        if c["comment"] == "THEORY" or pix not in applied:
            applied[pix] = c
    out = []
    for pix, c in cal_rows:
        if pix in meas and c["comment"] == "AS_MEASURED":
            continue          # superseded: new row appended below
        if c["comment"] == "AS_MEASURED" and pix in drop:
            log.info("cal_drop: removed stale AS_MEASURED row "
                     "for %s (reverts to its THEORY/FIDUCIAL row)",
                     pix)
            continue
        out.append((pix, dict(c)))
    # appended: fresh AS_MEASURED rows for currently measured mixers
    for pix, (dfx, dfy) in sorted(meas.items()):
        a = applied.get(pix, {"az": 0.0, "el": 0.0})
        out.append((pix, {"az": a["az"] + dfx, "el": a["el"] + dfy,
                          "comment": "AS_MEASURED"}))
    return out


def read_cal_version(path):
    """Calibration data version from the [offsets] header line
    ('[offsets] version = 0.5' -> 0.5). This version identifies the
    calibration DATA and is carried into the L1 products; it is
    bumped on every combined run of record (see aggregate_results).
    Defaults to 0.5 if the header carries no version."""
    for line in Path(path).read_text().splitlines():
        if line.strip().lower().startswith("[offsets]"):
            parts = line.replace("=", " ").split()
            try:
                return float(parts[-1])
            except ValueError:
                break
    return 0.5


def write_calibration_rows(rows, outfile, version="0.5"):
    """Write [(pix, {az, el, comment}), ...] in the [offsets] format.
    Used for the updated calibration file. `version` is the
    calibration DATA version (see read_cal_version)."""
    lines = [f"[offsets] version = {version}",
             "PIX\tAZ\t\tEL\t\tCOMMENT"]
    for pix, c in rows:
        lines.append(f"{pix}\t{c['az']:>10.6f}\t{c['el']:>10.6f}"
                     f"\t{c['comment']}")
    outfile.write_text("\n".join(lines) + "\n")
    log.info("wrote %s", outfile)
def aggregate_results(root, band=None, fid_by_band=None,
                      pix_band_by_band=None, roll=None, opts=None):
    """Sweep a tree of per-run output directories and combine.

    Each run directory (created by main()) holds an
    offsets_measurements.ecsv with the FULL per-region, per-slice
    inventory -- including obstime and gondola position. This:

      1. reads every */offsets_measurements.ecsv under <root>
         (optionally filtered to one band via the run-id prefix),
      2. stacks them into one master table,
         <root>/all_measurements.ecsv,
      3. combines per mixer across ALL runs (regions, velocity
         slices, epochs) with the same weighted-mean machinery,
         -> <root>/all_summary.ecsv (+ the per-mixer chi2 flags
         show whether fields agree),
      4. writes the master calibration file,
         <root>/offsets_<band>_deltas.txt (DELTA/residual format) and,
         if the env var GUSTO_CAL_FILE points at the applied
         calibration file, the UPDATED file.

    The combined table keeps region + obstime per row, so epoch
    structure survives for solving the roll theta0 later.
    """
    root = Path(root)
    rows = []
    n_runs = 0
    for mfile in sorted(root.glob(f"*/offsets_measurements.ecsv")):
        rdir = mfile.parent.name
        t = QTable.read(mfile)
        # target per run: the table's own column when present
        # (newer runs), else parsed from the run-dir name
        # <BAND>_<TARGET>_<TIMESTAMP>
        rows_ = []
        for r in t:
            d = dict(zip(t.colnames, [r[c] for c in t.colnames]))
            if "target" not in d or not str(d["target"]):
                parts = rdir.split("_")
                d["target"] = parts[1] if len(parts) > 1 else rdir
            rows_.append(d)
        # (the band filter below selects one band in -b mode; with
        # no -b, rows keep EVERY band -- the cross-band run of
        # record -- and all combines are per (band, mixer))
        if band and not rdir.startswith(str(band) + "_"):
            continue
        rows.extend(rows_)
        n_runs += 1
        log.info("aggregate: read %s (%d rows)", mfile, len(t))
    if not rows:
        raise SystemExit(f"no offsets_measurements.ecsv found under {root}")
    log.info("aggregate: %d runs, %d measurement rows", n_runs, len(rows))
    # analysis-time exclusions (config 'exclude:'): ECSVs stay the
    # complete record; excludes only what enters the combines
    rows_all = rows
    excl = (opts or {}).get("exclude", [])
    suffix = f"_{band}" if band else ""
    rows = apply_exclusions(rows, excl,
                            root / f"exclusions{suffix}.txt")
    # analysis-time mixer filter ('include_mixers: [5]' in config):
    # slice WHO enters the combine WITHOUT re-running measurements
    # -- restricting at measurement time OVERWRITES the per-run
    # record (mixers 2/3 were lost that way once)
    only = [str(m) for m in (opts or {}).get("include_mixers", [])]
    if only:
        n0 = len(rows)
        rows = [r for r in rows if str(r["mixer"]) in only]
        log.info("include_mixers %s: kept %d/%d rows",
                 only, len(rows), n0)
    # EXACT-duplicate safety net: the same run present TWICE under
    # the root (same band/target/region/mixer/vel/obstime AND
    # identical values to 1e-12 deg). Independent measurements are
    # never bit-identical; double-counting makes COMBINED sigmas
    # over-confident by sqrt(2), so the copy is DROPPED with a
    # loud warning pointing at the tree, not silently merged.
    seen = {}
    dedup = []
    n_dup = 0
    for r in rows:
        key = (str(r.get("band", "")), str(r["target"]),
               str(r["region"]), str(r["mixer"]),
               float(r["vel_lo"]), float(r["vel_hi"]),
               str(r.get("obstime", "")))
        sig = (round(float(r["fx"].to_value(u.deg)), 12),
               round(float(r["fy"].to_value(u.deg)), 12))
        if seen.get(key) == sig:
            n_dup += 1
            continue
        seen[key] = sig
        dedup.append(r)
    if n_dup:
        log.warning("dropped %d EXACT duplicate rows: the same run "
                    "is present twice under the aggregate root "
                    "-- check the 'aggregate: read ...' log lines "
                    "and clean the tree", n_dup)
        rows = dedup
    master = QTable(rows=rows_all)
    master.write(root / f"all_measurements{suffix}.ecsv",
                 overwrite=True, format="ascii.ecsv",
                 delimiter=",")
    opts = opts or {}
    cal_file = opts.get("cal_file") or \
        os.environ.get("GUSTO_CAL_FILE")
    # stale AS_MEASURED rows to REMOVE from the carried-over
    # calibration (pixels that crept in historically -- e.g. a
    # pixel later adopted as fiducial -- and have no current
    # measurement to supersede them). Their THEORY/FIDUCIAL rows
    # pass through: the pixel reverts to that value.
    cal_drop = set(opts.get("cal_drop", []) or [])

    # ---- PER-BAND processing --------------------------------------
    # The same mixer number in different bands is a DIFFERENT
    # physical pixel (B1M5 vs B2M5), so summaries, tables, plots,
    # roll solutions and calibration deltas are built PER BAND, in
    # banded (-b) AND cross-band (no -b) mode. All root-level
    # tables are band-suffixed, so a second band's -b aggregate
    # can no longer overwrite the first band's outputs.
    bands_present = sorted({str(r["band"]) for r in rows})
    fid_by_band = fid_by_band or {}
    pix_band_by_band = pix_band_by_band or {}
    band_summaries = {}         # band -> {mixer: summary}
    band_roll_sols = {}         # band -> roll solution
    per_target_all = {}         # tgt -> {(band, mixer): summary}
    for b in bands_present:
        brows = [r for r in rows if str(r["band"]) == b
                 and str(r["region"]) != "-"]
        if not brows:
            continue
        fid = str(fid_by_band.get(b, "")) or None
        pixb = pix_band_by_band.get(b, "B2")
        # combined summaries for THIS band's mixers. The FIDUCIAL
        # is skipped: it is the reference (zero by definition), and
        # any self-measurement rows (fiducial configured in a run's
        # usable list, or self-check rows) must NOT become
        # AS_MEASURED calibration rows or pull the roll solve.
        summ = {}
        for label in sorted({r["mixer"] for r in brows}):
            if str(label) == str(fid):
                continue
            mrows = [r for r in brows if r["mixer"] == label]
            if mrows:
                summ[label] = combine_measurements(mrows)
        band_summaries[b] = summ
        # summary WITH full input breakdown (ECSV + eyeball text)
        srows = summary_rows(brows, summ, fid)
        st = QTable(rows=srows)
        st.write(root / f"all_summary_{b}.ecsv", overwrite=True,
                 format="ascii.ecsv", delimiter=",")
        write_text_table(srows, root / f"all_summary_{b}.txt")
        # per-target summaries (within one target: slices are
        # correlated -- the mean's sigma is capped at the floor)
        per_target = {}
        for tgt in sorted({str(r["target"]) for r in brows}):
            trows = [r for r in brows
                     if str(r["target"]) == tgt]
            if not trows:
                continue
            sm = {}
            for label in sorted({r["mixer"] for r in trows}):
                if str(label) == str(fid):
                    # fiducial: reference, not a measurement --
                    # keep it out of the roll solve and the plot
                    continue
                mrows = [r for r in trows if r["mixer"] == label]
                if mrows:
                    sm[label] = combine_measurements(mrows,
                                                     floor_cap=True)
            per_target[tgt] = sm
            per_target_all[tgt] = {
                **per_target_all.get(tgt, {}),
                **{(b, label): v for label, v in sm.items()}}
        if not fid:
            log.warning("band %s: no fiducial in config "
                        "(bands.%s.fiducial) -- roll solve, plot "
                        "and calibration deltas skipped", b, b)
            continue
        band_roll_sol = solve_roll(
            per_target, root / f"roll_solution_{b}.txt")
        band_roll_sols[b] = band_roll_sol
        # the CONFIG roll, when set, is AUTHORITATIVE for the
        # aggregate outputs: the solved theta0 stays in
        # roll_solution*.txt as a diagnostic but is NOT applied.
        # Measurement and aggregate stay in the SAME convention
        # (e.g. roll: 0.0 everywhere).
        if roll is not None and band_roll_sol is not None:
            band_roll_sols[b] = {**band_roll_sol,
                                 "solved_theta0_deg":
                                     band_roll_sol["theta0_deg"],
                                 "theta0_deg": float(roll) * u.deg,
                                 "config_roll": True}
        plot_aggregate_summary(
            per_target, fid, b,
            root / f"summary_offsets_{b}.png")
        # combined RESIDUALS per band (measured-minus-applied): the
        # provenance record for what offsets_cal.txt folds in
        if pixb:
            write_calibration(summ, fid, pixb,
                              root / f"offsets_{b}_deltas.txt")
    # per-target combined summaries: built PER BAND in the loop
    # above (band_summaries / per_target_all)

    # -- band-INDEPENDENT roll solve ----------------------------------
    # CII and NII are in the SAME focal plane, so the solved roll
    # is one shared parameter across bands: each (band, mixer)
    # pixel keeps its OWN true offset, all share the ONE theta0.
    # Mixers are keyed (band, mixer) because the same number in
    # different bands is a different physical pixel. per_target_all
    # is built in the per-band loop above from the EXCLUSION-
    # FILTERED rows (the same population as every other output).
    cross_sol = None
    if len(bands_present) > 1:
        log.info("cross-band roll solve: bands %s, %d targets",
                 bands_present, len(per_target_all))
        cross_sol = solve_roll(per_target_all,
                               root / "roll_solution.txt")
        # config roll overrides the cross-band solution as well
        if roll is not None and cross_sol is not None:
            cross_sol = {**cross_sol,
                         "solved_theta0_deg": cross_sol["theta0_deg"],
                         "theta0_deg": float(roll) * u.deg,
                         "config_roll": True}

    # -- per-band roll-corrected files (only for NONZERO rotation) --
    # The roll de-rotation is LINEAR, so rotating the COMBINED
    # offsets by -theta0 at the end is mathematically identical to
    # re-running every per-target measurement with roll set in the
    # config. Sigmas rotate through the covariance (C' = R C R^T).
    # Prefer the cross-band solution (same focal plane); fall back
    # to the per-band one. With roll = 0 (the standard setup) the
    # _rot files would be BYTE-IDENTICAL duplicates -- skipped.
    for b, summ in band_summaries.items():
        fid = str(fid_by_band.get(b, "")) or None
        pixb = pix_band_by_band.get(b, "B2")
        theta_sol = cross_sol if cross_sol is not None \
            else band_roll_sols.get(b)
        if not (theta_sol and fid and pixb
                and "theta0_deg" in theta_sol):
            continue
        th = np.radians(
            float(theta_sol["theta0_deg"].to_value(u.deg)))
        if abs(th) <= 1e-6:
            continue
        c, s_ = np.cos(th), np.sin(th)   # de-rotation R(-theta0)
        Rm = np.array([[c, s_], [-s_, c]])
        summ_corr = {}
        for label, s in summ.items():
            if "fx" not in s or "fy" not in s:
                continue
            v = np.array([s["fx"].to_value(u.deg),
                          s["fy"].to_value(u.deg)])
            sg = np.array([s["sig_fx"].to_value(u.deg),
                           s["sig_fy"].to_value(u.deg)])
            u0 = Rm @ v
            sp = np.sqrt(np.diag(Rm @ np.diag(sg ** 2) @ Rm.T))
            # the roll solve's PDG scale carries the cross-target
            # scatter into these corrected sigmas too
            pdg = theta_sol.get("pdg_scale", 1.0)
            if pdg > 1.0:
                sp = sp * pdg
            summ_corr[label] = {"fx": u0[0] * u.deg,
                                "fy": u0[1] * u.deg,
                                "sig_fx": sp[0] * u.deg,
                                "sig_fy": sp[1] * u.deg}
        write_calibration(summ_corr, fid, pixb,
                          root / f"offsets_{b}_deltas_rot.txt")
        if cal_file:
            cal_rows = read_calibration(cal_file)
            upd = update_calibration(cal_rows, summ_corr,
                                     fid, pixb, drop=cal_drop)
            write_calibration_rows(
                upd, root / f"offsets_{b}_cal_updated_rot.txt")
        log.info("wrote roll-corrected calibration files "
                 "(theta0 = %+.3f deg)",
                 float(theta_sol["theta0_deg"].to_value(u.deg)))

    # -- THE combined calibration file --------------------------------
    # offsets_cal.txt is the run-of-record product of the CROSS-
    # BAND aggregate (no -b): ONE file for the whole array, each
    # band's residuals folded in SEQUENTIALLY (applied + delta ->
    # AS_MEASURED). A -b aggregate writes only that band's deltas
    # file and never touches offsets_cal.txt -- two -b runs would
    # otherwise clobber each other's bands.
    if band:
        log.info("offsets_cal.txt is written only by the "
                 "cross-band aggregate (--aggregate without -b); "
                 "this -b run wrote offsets_%s_deltas.txt", band)
    elif cal_file and band_summaries:
        cal_rows = read_calibration(cal_file)
        upd = cal_rows
        for b in bands_present:
            if b in band_summaries and fid_by_band.get(b):
                upd = update_calibration(
                    upd, band_summaries[b], str(fid_by_band[b]),
                    pix_band_by_band.get(b, "B2"),
                    drop=cal_drop)
        # VERSION BUMP: source + 0.1 per run of record, or an
        # explicit 'cal_version:' (e.g. 1.0 for the flight-ready
        # file); the DATA version lives on offsets_cal.txt (git)
        src_v = read_cal_version(cal_file)
        new_v = float(opts.get("cal_version", 0.0)) or \
            round(src_v + 0.1, 1)
        write_calibration_rows(upd, root / "offsets_cal.txt",
                               version=new_v)
        log.info("calibration version: %.1f -> %.1f "
                 "(offsets_cal.txt -> L1 products)", src_v, new_v)
    else:
        log.warning("no cal_file (config 'cal_file:' or "
                    "GUSTO_CAL_FILE env) -- offsets_cal.txt not "
                    "written; combined residuals are in "
                    "offsets_<band>_deltas.txt")
    if band:
        return band_summaries.get(str(band), {})
    return band_summaries


def write_calibration(summaries, fid_label, pix_band, outfile):
    """Write the measured mixer offsets in the calibration file
    format used by the GUSTO pointing solution, e.g.:

        [offsets] version = 0.5
        PIX	AZ		EL		COMMENT
        B2M8	0.000000	0.000000	FIDUCIAL
        B2M5	0.040013	-0.033674	AS_MEASURED

    AZ/EL are the focal-plane offsets fx/fy in DEGREES (roll=0:
    fx along +Az, fy along +El). PIX names combine the band tag
    (B1/B2) with the mixer number; the fiducial pixel is zero by
    definition.

    IMPORTANT: these are DELTA offsets. The L1 (hence the L2 cubes)
    were made WITH the calibration-file offsets already applied, so
    each value here is the residual (measured-minus-assumed), NOT
    the total offset. Rows are flagged AS_DELTA so the downstream
    pipeline (which USES only AS_MEASURED rows) cannot mistake
    deltas for totals. The TOTALS are produced by
    update_calibration(): applied + delta -> AS_MEASURED.
    """
    def fmt(v):
        return f"{v:>{8}.6f}\t"

    lines = ["[offsets] version = 0.5",
             "PIX\tAZ\t\tEL\t\tCOMMENT"]
    fid_pix = f"{pix_band}M{fid_label}"
    lines.append(f"{fid_pix}\t0.000000\t0.000000\tFIDUCIAL")
    for label in sorted(summaries):
        s = summaries[label]
        if label == fid_label or "fx" not in s or "fy" not in s:
            continue
        lines.append(f"{pix_band}M{label}\t"
                     f"{fmt(s['fx'].to_value(u.deg))}"
                     f"{fmt(s['fy'].to_value(u.deg))}"
                     f"AS_DELTA")
    outfile.write_text("\n".join(lines) + "\n")
    log.info("wrote %s (%d measured pixels)",
             outfile, len(lines) - 2)


def plot_summary(summary_by_mixer, outfile):
    """Focal-plane offset vectors with error ellipses per mixer."""
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 7))
    for label, s in summary_by_mixer.items():
        if "fx" not in s or "fy" not in s:
            continue
        # strip units: fx/fy are now Quantities in DEGREES, but
        # matplotlib needs plain floats
        x, y = s["fx"].to_value(u.deg), s["fy"].to_value(u.deg)
        ax.annotate("", xy=(x, y), xytext=(0, 0),
                    arrowprops=dict(arrowstyle="->", lw=1.5))
        ax.plot(x, y, "o")
        sx = float(s.get("sig_fx", 1 * u.deg).to_value(u.deg))
        sy = float(s.get("sig_fy", 1 * u.deg).to_value(u.deg))
        from matplotlib.patches import Ellipse
        ax.add_patch(Ellipse((x, y), 2 * sx, 2 * sy, color="k", lw=1,
                             fill=False))
        ax.text(x, y + 2 * sy, label, fontsize=9)
    ax.axhline(0, color="gray", lw=0.5); ax.axvline(0, color="gray", lw=0.5)
    ax.set_xlabel("focal-plane x (deg)")
    ax.set_ylabel("focal-plane y (deg)")
    ax.set_aspect("equal")
    ax.set_title("Mixer offsets in focal plane")
    fig.tight_layout()
    fig.savefig(outfile, dpi=120)
    plt.close(fig)
    log.info("wrote %s", outfile)


def plot_aggregate_summary(per_target, fid_label, band, outfile):
    """Aggregate offsets, one PANEL PER MIXER, all targets overlaid.

    per_target: {target: {mixer: summary}}. One panel per measured
    mixer (relative to the fiducial, number only in the title); all
    targets in LIGHT BLACK (alpha 0.5 -- vectors + 1/2-sigma
    ellipses, floored sigmas) so the cluster reads at a glance, the
    COMBINED weighted mean in FULL-ALPHA crimson on top (no text
    label needed -- the color carries it), which MUST sit inside
    the cluster. Axes in ARCSEC (degrees are unreadable at these
    scales). The COMBINED ellipse is clipped at the axes, not
    accommodated, so the per-target data sets the scale."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse
    mixers = sorted({m for tgt in per_target for m in per_target[tgt]})
    # drop the fiducial from the panels (zero by definition)
    mixers = [m for m in mixers if str(m) != str(fid_label)]
    if not mixers:
        log.info("aggregate plot: no measured mixers, skipping")
        return
    ncol = min(3, len(mixers))
    nrow = int(np.ceil(len(mixers) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(5 * ncol, 5 * nrow),
                             squeeze=False)
    for ax in axes.ravel()[len(mixers):]:
        ax.set_visible(False)          # unused panels
    for ax, label in zip(axes.ravel(), mixers):
        pts = []                       # (target, x, y, sx, sy) in ARCSEC
        for tgt, sm in per_target.items():
            s = sm.get(label)
            if not s or "fx" not in s:
                continue
            pts.append((tgt,
                        s["fx"].to_value(u.deg) * 3600.0,
                        s["fy"].to_value(u.deg) * 3600.0,
                        float(s.get("sig_fx", 1 * u.deg).to_value(u.deg))
                        * 3600.0,
                        float(s.get("sig_fy", 1 * u.deg).to_value(u.deg))
                        * 3600.0))
        # combined over targets (crimson): weighted mean of the
        # per-target summaries (select sm[label], NOT the dict)
        tgt_sm = [sm[label] for sm in per_target.values()
                  if label in sm and "fx" in sm[label]]
        combined = combine_measurements(tgt_sm) if len(tgt_sm) > 1 \
            else None
        # scale from the per-target points only; tight margin so a
        # near-zero cluster doesn't drown in empty space
        xmax = ymax = 0.0
        for tgt, x, y, sx, sy in pts:
            xmax = max(xmax, abs(x) + 2.5 * sx)
            ymax = max(ymax, abs(y) + 2.5 * sy)
        lim = max(xmax, ymax) * 1.08 if (xmax or ymax) else 20.0

        def inlim(v, frac=0.94):
            """clamp a label position inside the axes"""
            return float(np.clip(v, -lim * frac, lim * frac))

        for tgt, x, y, sx, sy in pts:
            # individuals: LIGHT so the cluster is readable and
            # the crimson COMBINED stands out on top
            ax.annotate("", xy=(x, y), xytext=(0, 0),
                        arrowprops=dict(arrowstyle="->", lw=1.2,
                                        alpha=0.5))
            ax.plot(x, y, "o", ms=4, alpha=0.5)
            for k, a in ((1, 0.5), (2, 0.25)):
                ax.add_patch(Ellipse((x, y), 2 * k * sx, 2 * k * sy,
                                     color="k", lw=1, fill=False,
                                     alpha=a))
            ax.text(inlim(x), inlim(y + 2 * sy), tgt, fontsize=7,
                    alpha=0.6)
        if combined and "fx" in combined:
            xc = combined["fx"].to_value(u.deg) * 3600.0
            yc = combined["fy"].to_value(u.deg) * 3600.0
            scx = float(combined.get("sig_fx",
                                     1 * u.deg).to_value(u.deg)) * 3600.0
            scy = float(combined.get("sig_fy",
                                     1 * u.deg).to_value(u.deg)) * 3600.0
            # COMBINED: full alpha, no text label (color says it)
            for k, a in ((1, 1.0), (2, 0.6)):
                ax.add_patch(Ellipse((xc, yc), 2 * k * scx, 2 * k * scy,
                                     color="crimson", lw=2, fill=False,
                                     alpha=a))
            ax.annotate("", xy=(xc, yc), xytext=(0, 0),
                        arrowprops=dict(arrowstyle="->", lw=2.5,
                                        color="crimson"))
            ax.plot(xc, yc, "o", ms=6, color="crimson")
        ax.axhline(0, color="gray", lw=0.5)
        ax.axvline(0, color="gray", lw=0.5)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
        ax.set_xlabel("focal-plane x (arcsec)")
        ax.set_ylabel("focal-plane y (arcsec)")
        ax.set_title(f"Mixer {label} vs {fid_label}", fontsize=11)
    fig.suptitle(f"{band} offsets by target "
                 f"(fiducial = mixer {fid_label})")
    fig.tight_layout()
    fig.savefig(outfile, dpi=120)
    plt.close(fig)
    log.info("wrote %s", outfile)


def solve_roll(per_target, outfile=None):
    """Solve the constant gondola roll theta0 from per-target offsets.

    MODEL: the measured focal offset of mixer m in target/epoch i is
    the TRUE focal offset v0_m rotated by the constant roll theta0:
        v_im = R(theta0) . v0_m
    All targets and mixers share the SAME theta0 (theta0 is a
    property of the gondola; the epoch-dependent sky geometry is
    already handled by M_frame in the measurement pipeline -- what
    leaks through when roll=0 is ASSUMED is this single common
    rotation). Per-target offsets that disagree in direction but
    agree once rotated by a common theta0 are exactly this
    signature (e.g. one target's fy flip vs the others).

    FIT: 1 shared theta0 + 2 parameters (v0_m) per mixer. For fixed
    theta0 each mixer's v0_m is the weighted mean of the
    de-rotated offsets (weights 1/sigma^2, sigmas rotated with the
    offsets); chi2(theta0) is profiled on a +/-45 deg grid with
    parabolic refinement. sigma_theta0 from the delta-chi2=1
    crossings; chi2/dof flags whether ONE common roll explains the
    target-to-target scatter (chi2/dof ~ 1: yes -- and then the
    per-mixer v0_m are the roll-corrected offsets; >> 1: something
    else varies between targets -- field systematics, not roll).

    Needs >=2 targets with measured mixers; returns None otherwise.
    """
    data = {}                     # mixer -> [(target, vec, sig), ...]
    for tgt, sm in per_target.items():
        for label, s in sm.items():
            if "fx" not in s or "fy" not in s:
                continue
            data.setdefault(label, []).append((
                tgt,
                np.array([s["fx"].to_value(u.deg),
                          s["fy"].to_value(u.deg)]),
                np.array([float(s["sig_fx"].to_value(u.deg)),
                          float(s["sig_fy"].to_value(u.deg))])))
    if len(per_target) < 2 or not data:
        log.info("roll solve: need >=2 targets with measured "
                 "mixers -- skipping")
        return None

    def chi2_of(theta):
        """total chi2 of the one-roll model at theta (radians)."""
        # de-rotate into the true focal frame: u = R(-theta) v
        c, s_ = np.cos(theta), np.sin(theta)
        Rm = np.array([[c, s_], [-s_, c]])
        tot, nd = 0.0, 0
        for m, items in data.items():
            us, sv = [], []
            for tgt, v, sg in items:
                us.append(Rm @ v)
                # sigmas rotate with the offsets: per-axis variance
                # of the rotated diagonal covariance
                Cr = Rm @ np.diag(sg ** 2) @ Rm.T
                sv.append(np.sqrt(np.diag(Cr)))
            us, sv = np.array(us), np.array(sv)
            w = 1.0 / sv ** 2
            u0 = np.sum(us * w, axis=0) / np.sum(w, axis=0)
            tot += float(np.sum((us - u0) ** 2 / sv ** 2))
            nd += 2 * len(us)
        return tot, nd

    # profile scan
    thetas = np.radians(np.arange(-45.0, 45.001, 0.25))
    chi2 = np.array([chi2_of(t)[0] for t in thetas])
    k = int(np.argmin(chi2))
    t0 = thetas[k]
    if 0 < k < len(thetas) - 1:            # parabolic refinement
        d = np.gradient(chi2, thetas)
        t0 = thetas[k] - d[k] / (2.0 * (chi2[k + 1] - 2 * chi2[k]
                                        + chi2[k - 1]) /
                                 (thetas[k + 1] - thetas[k]) ** 2)
    c2_min, nd = chi2_of(t0)
    npar = 1 + 2 * len(data)
    dof = max(nd - npar, 1)

    # sigma from delta-chi2 = 1 crossings
    def crossing(sign):
        for i in range(k, len(thetas) - 1) if sign > 0 \
                else range(k, 0, -1):
            j = i if sign > 0 else i - 1
            if sign < 0:
                i, j = j, j + 1
            if chi2[j] - c2_min >= 1.0:
                # interpolate between thetas[i], thetas[j]
                d0 = chi2[i] - c2_min
                d1 = chi2[j] - c2_min
                return abs(thetas[j] - thetas[i]) * (1.0 - d0) / \
                    (d1 - d0) if d1 > d0 else abs(thetas[j] - thetas[i])
        return np.radians(45.0)
    sig_t = max(crossing(+1), crossing(-1))

    # PDG scaling: with the (known-optimistic) formal per-target
    # sigmas, chi2/dof >> 1 does NOT mean the roll is wrong -- the
    # FIT MINIMUM (theta0 value) is insensitive to the sigma scale,
    # but its reported UNCERTAINTY is not. Inflate sig_theta0 and
    # the per-mixer offset sigmas by sqrt(chi2/dof), the standard
    # honest-error rescale. chi2/dof ~ 1 -> no change.
    pdg = np.sqrt(c2_min / dof)
    sig_t = sig_t * pdg if pdg > 1.0 else sig_t

    # per-mixer roll-corrected offsets at the solution
    c, s_ = np.cos(t0), np.sin(t0)
    Rm = np.array([[c, s_], [-s_, c]])
    mixers = {}
    detgt_rows = []               # per-target de-rotated values + residuals
    for m, items in data.items():
        us, sv, tgs = [], [], []
        for tgt, v, sg in items:
            us.append(Rm @ v)
            Cr = Rm @ np.diag(sg ** 2) @ Rm.T
            sv.append(np.sqrt(np.diag(Cr)))
            tgs.append(tgt)
        us, sv = np.array(us), np.array(sv)
        w = 1.0 / sv ** 2
        u0 = np.sum(us * w, axis=0) / np.sum(w, axis=0)
        su = np.sqrt(1.0 / np.sum(w, axis=0))
        if pdg > 1.0:
            su = su * pdg
        mixers[m] = (u0, su)
        for tgt, ui, sgi in zip(tgs, us, sv):
            detgt_rows.append(
                (m, tgt, ui[0], ui[1],
                 float(np.linalg.norm(ui - u0)) * 3600.0))
    sol = {"theta0_deg": np.degrees(t0) * u.deg,
           "sig_theta0_deg": np.degrees(sig_t) * u.deg,
           "chi2": c2_min, "dof": dof,
           "chi2_per_dof": c2_min / dof, "pdg_scale": pdg,
           "n_targets": len(per_target), "mixers": mixers}
    log.info("roll solution: theta0 = %.3f +/- %.3f deg "
             "(chi2/dof = %.1f -> sigmas PDG-inflated by %.1f, "
             "%d targets)",
             sol["theta0_deg"].value, sol["sig_theta0_deg"].value,
             sol["chi2_per_dof"], pdg, len(per_target))
    if outfile:
        lines = [f"roll solution from {len(per_target)} targets",
                 f"theta0 = {sol['theta0_deg'].value:+.4f} "
                 f"+/- {sol['sig_theta0_deg'].value:.4f} deg",
                 f"chi2 = {c2_min:.1f} / {dof} dof",
                 "",
                 "roll-corrected focal offsets (true frame):",
                 "mixer  fx(deg)      sig     fy(deg)      sig"]
        for m, (u0, su) in sorted(mixers.items()):
            # m may be a plain mixer number (per-band solve) or a
            # (band, mixer) tuple (cross-band solve): str() so the
            # alignment spec works for both
            lines.append(f"{str(m):>5}  {u0[0]:+.6f}  {su[0]:.6f}  "
                         f"{u0[1]:+.6f}  {su[1]:.6f}")
        lines += ["",
                  "interpretation: chi2/dof ~ 1 -> one common roll",
                  "explains the target scatter, and the offsets above",
                  "are the roll-corrected values. If instead chi2/dof",
                  ">> 1, targets disagree beyond a single rotation --",
                  "check per-field systematics before trusting these.",
                  "",
                  "to APPLY the roll: set  roll: <theta0>  in the",
                  "config options and re-run the measurement pipeline",
                  "(per-target runs), so the measured offsets come out",
                  "in the true focal frame directly."]
        lines += ["",
                  "per-target DE-ROTATED offsets (true focal frame):",
                  "the individual deltas corrected by the solved roll;",
                  "|u - u0| is the residual scatter (arcsec) around",
                  "the mixer's combined offset -- this is where to",
                  "check whether ONE roll really aligns the targets:",
                  "residuals ~ half a pixel (10-20 arcsec) = roll does",
                  "the job (chi2 excess = sigma bookkeeping); a",
                  "beam-fraction outlier flags a per-field systematic.",
                  "",
                  "mixer  target               fx(deg)     fy(deg)   |u-u0|(\")"]
        for m, tgt, ux, uy, res in detgt_rows:
            lines.append(f"{str(m):>5}  {str(tgt):<18}  "
                         f"{ux:+.6f}  {uy:+.6f}  {res:8.1f}")
        # chi2(theta) profile: the shape that matters for judging
        # constraint. A sharp basin ~1 deg wide = constrained; a
        # flat curve (sigma hits the scan-boundary fallback, e.g.
        # "theta0 = -42 +/- 45") = UNCONSTRAINED -- the offsets are
        # too short relative to their sigmas for any rotation to
        # matter (tip displacement |v|*theta0 is below the noise);
        # theta0 is then consistent with 0 and should NOT be applied.
        lines += ["",
                  "chi2(theta) profile (dof = %d):" % dof,
                  "theta(deg)   chi2    chi2/dof"]
        for tt in np.arange(-45.0, 45.01, 5.0):
            c2t, _ = chi2_of(np.radians(tt))
            lines.append(f"{tt:+7.1f}  {c2t:8.1f}   {c2t / dof:7.2f}")
        Path(outfile).write_text("\n".join(lines) + "\n")
        log.info("wrote %s", outfile)
    return sol


def montage_extent(wcs_celestial, pixbox):
    """World-coordinate extent of the pixel box, computed from the
    ORIGINAL cube WCS, in the cube's NATIVE frame (no frame
    conversion!). Returns (lon0, lon1, lat0, lat1, frame) in degrees
    for imshow(extent=), where frame is 'galactic' or 'icrs' --
    callers must draw markers and axis labels in THAT frame too
    (transform region coordinates into it; e.g. an ICRS-defined
    LMC region plotted against a Galactic cube axis previously
    produced nonsense extents).

    Implementation via pixel_to_skycoord: SkyCoord knows which
    coordinate is longitude/latitude regardless of the CTYPE axis
    order, eliminating the corner/axis bookkeeping that previously
    mixed lat values into the lon range (axis showed ~-1..348).
    Corners normalized to the NATIVE frame. lon range: corners
    differing in x only; lat: in y."""
    from astropy.wcs.utils import pixel_to_skycoord
    x0, x1, y0, y1 = pixbox
    ctype = wcs_celestial.wcs.ctype[0].upper()
    frame = "galactic" if ctype.startswith("GLON") else "icrs"
    c00 = pixel_to_skycoord(x0, y0, wcs_celestial)
    c10 = pixel_to_skycoord(x1, y0, wcs_celestial)
    c01 = pixel_to_skycoord(x0, y1, wcs_celestial)
    if frame == "galactic":
        lon_a, lon_b = c00.galactic.l.deg, c10.galactic.l.deg
        lat_a, lat_b = c00.galactic.b.deg, c01.galactic.b.deg
    else:
        lon_a, lon_b = c00.icrs.ra.deg, c10.icrs.ra.deg
        lat_a, lat_b = c00.icrs.dec.deg, c01.icrs.dec.deg
    # ascending ranges (lon may decrease along +x; for RA the
    # extent handed to imshow then runs descending, which is
    # CORRECT: RA increases to the left)
    # WRAP: near l=0 / RA=0 the two lon corners sit on either side
    # of the wrap (e.g. 359.7 and 0.3); naive min/max would give
    # an extent covering ~all 360 deg. Unwrap to the SHORT way
    # around, then DISPLAY-NORMALIZE: a box straddling the wrap
    # is shown as NEGATIVE longitudes (359.7..360.3 -> -0.3..+0.3)
    # instead of an axis that runs across 360 deg -- and callers
    # must normalize markers with the same convention.
    if lon_b - lon_a > 180.0:
        lon_a += 360.0
    elif lon_a - lon_b > 180.0:
        lon_b += 360.0
    lo, hi = min(lon_a, lon_b), max(lon_a, lon_b)
    if lo > 180.0:
        lo -= 360.0
        hi -= 360.0
    log.info("montage_extent: corners lon %.4f,%.4f lat %.4f,%.4f "
             "-> extent (%.4f, %.4f, %.4f, %.4f) %s",
             lon_a, lon_b, lat_a, lat_b,
             lo, hi, min(lat_a, lat_b), max(lat_a, lat_b), frame)
    return (lo, hi,
            min(lat_a, lat_b), max(lat_a, lat_b), frame)


# ----------------------------------------------------------------------
# 7. Output
# ----------------------------------------------------------------------

def summary_rows(measurements, summaries, fid_label):
    """Per-mixer summary rows WITH full input breakdown.

    For each mixer (EXCLUDING the fiducial, whose offset is zero by
    definition -- including its self-check row -- so the summary only
    lists real measurements):
      - one row per (region, velocity range): the individual slice,
        input labeled '<region> [vel_lo,vel_hi]'
      - one row per region: weighted mean over that region's slices,
        input labeled '<region>'
      - one COMBINED row: weighted mean over everything.
    The COMBINED row is what write_calibration / update_calibration
    consume. chi2 flags on slice rows are meaningless (n=1, dof=0);
    on region rows they show slice agreement; on COMBINED they show
    region agreement.
    """
    rows = []
    for label in sorted({r["mixer"] for r in measurements}):
        if str(label) == str(fid_label):
            continue          # fiducial: reference, not a measurement
        m_all = [r for r in measurements
                 if r["mixer"] == label and str(r["region"]) != "-"]
        if not m_all:
            continue
        # per (region, velocity range): individual slice rows
        for r in m_all:
            rows.append({"mixer": label,
                         "input": f"{r['region']} "
                                  f"[{r['vel_lo']},{r['vel_hi']}]",
                         **combine_measurements([r])})
        # per region: weighted mean over that region's slices
        # (correlated when the vel ranges overlap/are continuous --
        # sigma capped at the floor)
        for rg in sorted({str(r["region"]) for r in m_all}):
            mrows = [r for r in m_all if str(r["region"]) == rg]
            rows.append({"mixer": label, "input": rg,
                         **combine_measurements(mrows,
                                                floor_cap=True)})
        if label in summaries:
            rows.append({"mixer": label, "input": "COMBINED",
                         **summaries[label]})
    return rows


def write_text_table(rows, outfile):
    """Human-readable FIXED-WIDTH text version of summary rows.

    The ECSV files remain the FULL-PRECISION machine record
    (reverting earlier rounding: a data file should carry what the
    fit produced, and re-analysis must not depend on display
    choices). This is the EYEBALL copy: degrees at 5 decimals
    (0.036 arcsec), floats at 6 significant digits, aligned
    columns -- readable directly, no other program needed."""
    if not rows:
        return
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)

    def fmt(v):
        if v is None:
            return "-"
        if hasattr(v, "unit"):
            x = float(v.value)
            if v.unit == u.deg:
                return f"{x:+.5f}"
            return f"{x:.6g}"
        if isinstance(v, (float, np.floating)):
            return f"{v:.6g}"
        return str(v)

    headers = []
    for k in keys:
        sample = next((r[k] for r in rows
                       if k in r and r[k] is not None), None)
        unit = getattr(sample, "unit", None)
        headers.append(f"{k} [{unit}]" if unit is not None else k)
    cols = []
    for k, h in zip(keys, headers):
        vals = [fmt(r.get(k)) for r in rows]
        w = max([len(h)] + [len(v) for v in vals])
        cols.append((h.ljust(w), [v.ljust(w) for v in vals]))
    lines = ["  ".join(h for h, _ in cols),
             "  ".join("-" * len(h) for h, _ in cols)]
    for i in range(len(rows)):
        lines.append("  ".join(v[i] for _, v in cols))
    Path(outfile).write_text("\n".join(lines) + "\n")
    log.info("wrote %s", outfile)


def _row_matches(r, e):
    """Does measurement row r match one config 'exclude:' entry?
    An entry matches on EVERY field it specifies: band, target,
    region, mixer, vel ([lo,hi], matched on overlap)."""
    for key in ("band", "target", "region", "mixer"):
        if key in e and str(r[key]) != str(e[key]):
            return False
    if "vel" in e:
        vlo, vhi = sorted(float(x) for x in e["vel"])
        if not (float(r["vel_lo"]) <= vhi and float(r["vel_hi"]) >= vlo):
            return False
    return True


def apply_exclusions(rows, excl, outfile):
    """Apply the analysis-time exclusion list; write an AUDIT file.

    Outliers identified in the figures are excluded by CONFIG
    (top-level 'exclude:' list), never by editing the ECSV record.
    Matching rows are kept out of ALL combines (summaries, roll
    solve, calibration); the measurement ECSVs stay complete, and
    every exclusion is documented with its reason in
    <outfile> -- the record the manual line-deletion never left."""
    if not excl:
        return rows
    kept, out = [], []
    for r in rows:
        hit = next((e for e in excl if _row_matches(r, e)), None)
        if hit is not None:
            out.append((r, hit))      # audit pairs (row, entry)
        else:
            kept.append(r)            # BARE rows -- not pairs!
    if out:
        lines = ["# exclusions applied at aggregate time",
                 "# (config 'exclude:'; ECSV record unchanged)",
                 "# target  region  mixer  vel_range(deg)",
                 "# fx(deg)   fy(deg)   reason"]
        for r, e in out:
            lines.append(f"{r['target']}  {r['region']}  "
                         f"{r['mixer']}  "
                         f"[{float(r['vel_lo']):.1f},"
                         f"{float(r['vel_hi']):.1f}]  "
                         f"{float(r['fx'].to_value(u.deg)):+.5f}  "
                         f"{float(r['fy'].to_value(u.deg)):+.5f}  "
                         f"{e.get('reason', '')}")
        Path(outfile).write_text("\n".join(lines) + "\n")
        log.info("exclusions: %d rows removed (audit: %s)",
                 len(out), outfile)
    else:
        log.info("exclusions: config list matched no rows")
    return kept


def save_results(measurements, summaries, fid_label, outdir):
    """Save the per-measurement table and the summary table (ECSV).

    The summary table shows ALL inputs -- one row per (mixer,
    region, velocity range), then per region, then the COMBINED
    weighted average per mixer (summary_rows). The measurements
    table remains the full inventory with epoch/attitude per row."""
    mt = QTable(rows=measurements)
    # ECSV: comma delimiter (astropy default is space, which makes
    # the data block awkward to eyeball / ingest in other tools)
    mt.write(outdir / "offsets_measurements.ecsv", overwrite=True,
             format="ascii.ecsv", delimiter=",")
    # eyeball copy: fixed-width text table, readable precision
    srows = summary_rows(measurements, summaries, fid_label)
    st = QTable(rows=srows)
    st.write(outdir / "offsets_summary.ecsv", overwrite=True,
             format="ascii.ecsv", delimiter=",")
    write_text_table(srows, outdir / "offsets_summary.txt")
    log.info("wrote results tables to %s", outdir)


# ----------------------------------------------------------------------
# Main pipeline
# ----------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("-c", "--config", required=True,
                   help="general config YAML: root directories, "
                        "per-band fiducial + usable mixers, options")
    p.add_argument("-b", "--band", default=None,
                   help="band label, e.g. CII; optional with "
                        "--aggregate (no -b = cross-band run of "
                        "record over all bands)")
    p.add_argument("-t", "--target", default=None,
                   help="target name; data live under "
                        "<level1_root>/<target>/ and "
                        "<level2_root>/<target>/")
    p.add_argument("-o", "--output", default=None,
                   help="output directory; defaults to the current "
                        "directory when not given")
    p.add_argument("--ref-cube", default=None,
                   help="combined (mixer 0) reference cube; only sets "
                        "the common WCS grid; default resolved from "
                        "the naming convention")
    p.add_argument("-r", "--regions", default=None,
                   help="YAML region file (short form: -r)")
    p.add_argument("--beam", default=None,
                   help="beam size used in the regridding (part of "
                        "the cube filename); default from config")
    p.add_argument("--kernel", default=None, choices=["B", "G", "N"],
                   help="gridding kernel: B (gaussbessel), G (gauss), "
                        "N (nearest); default from config")
    p.add_argument("--collapse", default=None,
                   choices=["moment0", "mean"],
                   help="override config options.collapse")
    p.add_argument("--patch", type=int, default=None,
                   help="override config options.patch")
    p.add_argument("--min-snr", type=float, default=None,
                  help="override min SNR of the collapsed region image")
    p.add_argument("--selfcheck", action="store_true",
                   help="also correlate the fiducial with itself "
                        "(sanity check: must give zero offset)")
    p.add_argument("--aggregate", default=None, metavar="ROOT",
                   help="AGGREGATION MODE: sweep all per-run output "
                        "directories under ROOT, combine all "
                        "measurements per mixer, and write the master "
                        "summary + calibration files. Other run "
                        "arguments are then optional.")
    p.add_argument("--fid", default=None,
                   help="aggregate mode: fiducial mixer for the "
                        "calibration output (default from config)")
    p.add_argument("--pix-band", default=None,
                   help="aggregate mode: B1/B2 pixel band tag "
                        "(default from config or B2)")
    args = p.parse_args(argv)
    if not args.band and not args.aggregate:
        p.error("-b/--band is required for a measurement run; "
                "--aggregate may be run with or without -b "
                "(without = cross-band run of record)")
    return args


def main(argv=None):
    args = parse_args(argv)
    # -o absent -> current directory
    # per-run output dir: <output_root>/<BAND>_<TARGET>_<timestamp>.
    # -o sets the ROOT; each run gets its own timestamped directory,
    # so results are never overwritten and the aggregate mode can
    # sweep the whole tree of runs.
    run_id = f"{args.band}_{args.target}_" \
             f"{Time.now().strftime('%Y%m%dT%H%M%S')}"
    outdir = (Path(args.output) / run_id) if args.output \
        else Path.cwd() / "offsets_runs" / run_id
    outdir.mkdir(parents=True, exist_ok=True)
    log.info("run output directory: %s", outdir)

    # -- AGGREGATION MODE ----------------------------------------------
    # combine the results of many per-field runs; no new correlation
    # is computed. Example:
    #   ./gusto_xcorr_offsets.py -c cfg.yaml --aggregate offsets_runs/ \
    #       --fid 8 --pix-band B2
    # (GUSTO_CAL_FILE env var: path of the APPLIED calibration file;
    # when set, an updated [offsets] file is written too)
    if args.aggregate:
        # Aggregate combines per-run results. The CROSS-BAND run of
        # record (no -b) needs EVERY band's fiducial + pix_band tag
        # -- all from the config; CLI --fid/--pix-band override only
        # the single band of a -b run.
        data_cfg, bands, opts = load_config(args.config) \
            if args.config else ({}, {}, {})
        fid_by_band = {b: str(bc.get("fiducial"))
                       for b, bc in bands.items()
                       if bc.get("fiducial") is not None}
        pix_band_by_band = {
            b: bc.get("pix_band", opts.get("pix_band", "B2"))
            for b, bc in bands.items()}
        if args.band:
            if args.fid is not None:
                fid_by_band[args.band] = str(args.fid)
            if args.pix_band is not None:
                pix_band_by_band[args.band] = args.pix_band
        if not fid_by_band:
            raise SystemExit("no band fiducials in config "
                             "(bands.<band>.fiducial) -- required "
                             "for --aggregate")
        aggregate_results(args.aggregate, band=args.band,
                          fid_by_band=fid_by_band,
                          pix_band_by_band=pix_band_by_band,
                          roll=opts.get("roll", None),
                          opts=opts)
        return 0

    # measurement mode: target + regions ARE required (argparse can't
    # express "required unless --aggregate", so validate here)
    if args.target is None:
        raise SystemExit("measurement mode requires -t/--target")
    if args.regions is None:
        raise SystemExit("measurement mode requires --regions")
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.StreamHandler(),
                  logging.FileHandler(outdir / "run.log", mode="w")])

    # General config: root directories, bands, options. The target is
    # a CLI input; all paths resolve as <root>/<target>/...
    data_cfg, bands_cfg, opts = load_config(args.config)
    band_cfg = bands_cfg[args.band]
    fid_label = band_cfg["fiducial"]
    target = args.target
    l1_dir = Path(data_cfg["level1_root"]) / target
    l2_dir = Path(data_cfg["level2_root"]) / target
    if not l1_dir.is_dir():
        raise SystemExit(f"level-1 target dir not found: {l1_dir}")
    if not l2_dir.is_dir():
        raise SystemExit(f"level-2 target dir not found: {l2_dir}")

    # CLI overrides of config options
    # NOTE: named collapse_method, NOT collapse -- a local variable
    # named `collapse` would shadow the collapse() FUNCTION and
    # cause "TypeError: 'str' object is not callable"
    collapse_method = args.collapse or opts.get("collapse", "moment0")
    patch = args.patch or int(opts.get("patch", 7))
    beam = args.beam or opts.get("beam")
    kernel = args.kernel or opts.get("kernel", "B")
    if beam is None:
        raise SystemExit("beam size not set: use --beam or "
                        "options.beam in the config")

    # Cube naming convention:
    #   <target>_<band>_<mixer>_at_<Beam>_<kernel>.fits
    # mixer: 0 = all mixers combined, else the mixer number.
    # kernel: B (gaussbessel), G (gauss), N (nearest neighbor).
    # Beam: arcsec, plain decimal; trailing zeros are NOT stripped
    # (e.g. 0.72 -> "at_0.72_", 1.0 -> "at_1.0_").
    def beam_str(b):
        s = str(float(b))
        return s + ".0" if "." not in s else s

    def cube_name(mixer_num):
        return (f"{target}_{args.band}_{mixer_num}"
                f"_at_{beam_str(beam)}_{kernel}.fits")

    # Launch-site reference (lat/lon); the EarthLocation height is
    # taken per scan from GON_ALT once attitude lookup is active.
    launch_site = EarthLocation(
        lat=float(opts.get("launch_lat", -77.85)) * u.deg,
        lon=float(opts.get("launch_lon", 166.67)) * u.deg,
        height=36 * u.km)

    def make_location(lat, lon, alt):
        """Per-scan EarthLocation from the gondola's OWN GPS position
        and altitude (circumpolar flight: position sweeps 360 deg in
        longitude over the mission -- a fixed launch site is wrong).
        The launch-site config is only a FALLBACK if the L1 header
        lacks the position keywords."""
        return EarthLocation(lat=lat, lon=lon, height=alt)

    regions = load_regions(args.regions)
    # combined reference cube (mixer 0) only sets the common WCS grid;
    # default resolved from the naming convention unless --ref-cube
    ref_wcs = None
    ref_cube_path = (Path(args.ref_cube) if args.ref_cube
                     else l2_dir / cube_name(0))
    if ref_cube_path.exists():
        _, _, ref_wcs, _ = load_cube(ref_cube_path)
    else:
        log.warning("reference (mixer 0) cube not found: %s -- "
                    "falling back to the fiducial mixer's WCS",
                    ref_cube_path)

    # Per-band mixer cubes: resolved from the naming convention
    # <target>_<band>_<mixer>_at_<Beam>_<kernel>.fits in l2_dir.
    # All cube names share this pattern, so the mixer set can be
    # auto-discovered by globbing the target's level-2 directory
    # (excluding mixer 0, the combined cube). An explicit mixers:
    # list in the config restricts/overrides the discovery; explicit
    # cubes: entries (absolute paths used as-is) override names.
    mixers = band_cfg.get("mixers")
    if not mixers:
        # auto-discover: glob <target>_<band>_*_at_<beam>_<kernel>.fits
        pat = (f"{target}_{args.band}_*_at_{beam_str(beam)}_{kernel}.fits")
        discovered = set()
        for f in l2_dir.glob(pat):
            mstr = f.name.split(f"{args.band}_", 1)[1].split("_at_")[0]
            if mstr.isdigit() and int(mstr) != 0:   # 0 = combined cube
                discovered.add(int(mstr))
        mixers = sorted(discovered)
        log.info("auto-discovered mixers for %s %s: %s",
                 target, args.band, mixers)
        if not mixers:
            raise SystemExit(f"no per-mixer cubes matching {l2_dir}/{pat}")
        if band_cfg["fiducial"] not in mixers:
            raise SystemExit(f"fiducial {band_cfg['fiducial']} not "
                             f"among discovered mixers {mixers}")
    mixer_paths = {}
    for m in mixers:
        fname = band_cfg["cubes"].get(m, cube_name(m))
        fp = Path(fname)
        mixer_paths[m] = fp if fp.is_absolute() else l2_dir / fname
    for m, pth in mixer_paths.items():
        if not pth.exists():
            raise SystemExit(f"missing cube for mixer {m}: {pth}")
    if fid_label not in mixer_paths:
        raise SystemExit(f"fiducial '{fid_label}' not configured "
                         f"for band {args.band}")
    fid_path = mixer_paths[fid_label]
    # the fiducial cube's WCS is the shared grid (fallback: ref cube)
    fid_data, fid_hdr, fid_wcs, fid_seqmap = load_cube(fid_path)
    if ref_wcs is None:
        ref_wcs = fid_wcs

    # -- grid sampling: error inflation --------------------------------
    # PIX_BEAM: pixels per beam (e.g. 3); BEAMFWHM: beam FWHM in
    # arcmin. The pixels are NOT independent at this sampling: the
    # correlation surface has ~N_pix / PIX_BEAM^2 effective
    # independent samples, so the LM fit's formal covariance
    # (which assumes N_pix independent points) is optimistic by
    # ~PIX_BEAM. Inflate the covariance by PIX_BEAM^2 (sigmas by
    # PIX_BEAM) to make the reported uncertainties honest.
    pix_beam = int(fid_hdr.get("PIX_BEAM", 3))
    beam_fwhm = float(fid_hdr.get("BEAMFWHM", np.nan)) * u.arcmin
    log.info("grid: %d pixels/beam, beam FWHM %s -- sigmas inflated "
             "by %d", pix_beam, beam_fwhm, pix_beam)

    # self-correlation of the fiducial (mixer vs ITSELF) is a
    # sanity check: peak must land at zero lag, dx = dy = 0, and the
    # error ellipse shows the method's floor sigma. Enable via
    # config `selfcheck: true` or --selfcheck; off by default.
    selfcheck = args.selfcheck or bool(opts.get("selfcheck", False))

    measurements = []

    for label, mcube_path in mixer_paths.items():
        if label == fid_label and not selfcheck:
            continue          # fiducial is the reference, not a target
        m_data, m_hdr, m_wcs, m_seqmap = load_cube(mcube_path)

        for region in regions:
            shape2d = (fid_data.shape[2], fid_data.shape[1])  # (nx, ny); adapt
            pixbox = skybox_to_pixbox(ref_wcs.celestial, region, shape2d)
            # -- coordinate-chain diagnostic -----------------------
            # region center -> cube pixel -> montage extent; all from
            # the ORIGINAL cube WCS, so a wrong frame/order/axis
            # mapping is visible in the log immediately
            (rxp,), (ryp,) = ref_wcs.celestial.wcs_world2pix(
                [region.lon.to_value(u.deg)],
                [region.lat.to_value(u.deg)], 0)
            log.info("region %s: center (%.4f, %.4f) deg -> cube pixel "
                      "(%.1f, %.1f), pixbox %s, cube shape %s",
                      region.name, region.lon.to_value(u.deg),
                      region.lat.to_value(u.deg), rxp, ryp, pixbox,
                      shape2d)
            if not (0 <= rxp < shape2d[0] and 0 <= ryp < shape2d[1]):
                raise SystemExit(
                    f"region {region.name} center falls OUTSIDE the "
                    f"cube footprint (pixel {rxp:.1f}, {ryp:.1f} vs "
                    f"shape {shape2d}) -- check frame/coordinates")

            # -- optional SPECTRAL BASELINE correction (FIRST) -----
            # cut the UNION of the region's velocity ranges,
            # extended by baseline_pad spectral elements on BOTH
            # sides, correct per pixel EXCLUDING all line channels
            # (overlapping/adjacent ranges combine), then slice
            # each range from the corrected cube below.
            bpad = int(opts.get("baseline_pad", 0))
            bmin = int(opts.get("baseline_min_pts", 6))
            fid_base = mix_base = None
            pad_lo = 0
            if bpad > 0 and region.vel_ranges:
                vsl = [vel_to_slice(ref_wcs, (lo, hi))
                       for (lo, hi) in region.vel_ranges]
                u0 = min(s[0] for s in vsl)
                u1 = max(s[1] for s in vsl)
                pad_lo = max(u0 - bpad, 0)
                fid_base, _ = extract_subcube(
                    fid_data, fid_wcs, pixbox,
                    (pad_lo, min(u1 + bpad, fid_data.shape[0])))
                fid_base = baseline_correct(
                    fid_base, region.vel_ranges, fid_wcs,
                    pad_lo, min_pts=bmin)
                mix_base, _ = extract_subcube(
                    m_data, m_wcs, pixbox,
                    (pad_lo, min(u1 + bpad, m_data.shape[0])))
                mix_base = baseline_correct(
                    mix_base, region.vel_ranges, m_wcs,
                    pad_lo, min_pts=bmin)

            for (vlo, vhi) in region.vel_ranges:
                vslice = vel_to_slice(ref_wcs, (vlo, vhi))

                # -- attitude: SEQ_ID at the region center ---------
                # center sky position -> pixel -> SEQ_ID (HDU 3) of
                # the MEASURED mixer's cube (not the combined or the
                # fiducial's -- e.g. mixer 5 vs 8 uses mixer 5's map,
                # mixer 2 vs 3 uses mixer 2's) -> exact L1 filename
                # match via glob (SEQ_ID is the first field in
                # BAND_{seq}_{scanid}_L10.fits)
                seq_id = seqid_at_center(m_seqmap, pixbox,
                                         ref_wcs.celestial, region)
                # az, el, alt, obstime = get_scan_attitude(
                #     l1_dir, args.band, seq_id)
                # location = make_location(alt)
                # M_frame = galactic_to_altaz_matrix(
                #     obstime, location, region.lon.value,
                #     region.lat.value)
                # theta0 handled as a constant -- see NOTES_ON_ROLL
                # (set to 0, or solved from epoch scatter)

                # -- extract the same sub-cube from fiducial & mixer --
                if fid_base is not None:
                    # slice the range out of the BASELINE-CORRECTED
                    # padded cube: spatial pass-through, relative
                    # spectral window (extract_subcube keeps the
                    # subwcs contract for downstream QA panels)
                    ref_sub, subwcs = extract_subcube(
                        fid_base, fid_wcs,
                        (0, fid_base.shape[2], 0, fid_base.shape[1]),
                        (vslice[0] - pad_lo, vslice[1] - pad_lo))
                    mix_sub, _ = extract_subcube(
                        mix_base, m_wcs,
                        (0, mix_base.shape[2], 0, mix_base.shape[1]),
                        (vslice[0] - pad_lo, vslice[1] - pad_lo))
                else:
                    ref_sub, subwcs = extract_subcube(fid_data,
                                                      fid_wcs,
                                                      pixbox, vslice)
                    mix_sub, _ = extract_subcube(m_data, m_wcs,
                                                pixbox, vslice)

                ref_img = preprocess_image(
                    collapse(ref_sub, collapse_method), opts)
                mix_img = preprocess_image(
                    collapse(mix_sub, collapse_method), opts)

                # -- QA montage: region images BEFORE the SNR gate ---
                # ALWAYS written, also for skipped measurements, so
                # a zero-SNR skip can be diagnosed from the images.
                # Extent from the ORIGINAL cube WCS + pixbox (never a
                # CRPIX-shifted sub-WCS); center marker relative to
                # the box origin.
                reg_extent = montage_extent(ref_wcs.celestial, pixbox)
                ok, snr = snr_gate(mix_img,
                                   args.min_snr or region.min_snr)
                plot_region_montage(
                    [ref_img, mix_img], [fid_label, label],
                    reg_extent, region, (vlo, vhi),
                    outdir / f"img_{label}_{region.name}_{vlo}_{vhi}.png",
                    snr=snr)
                if not ok:
                    log.warning("%s %s %s: SNR %.1f below threshold, skip",
                                label, region.name, (vlo, vhi), snr)
                    continue

                # -- cross-correlation vs the fiducial mixer --------
                corr = fft_xcorr(mix_img, ref_img)
                int_peak, peak_val = integer_peak(corr)
                snr_corr = peak_snr(corr, int_peak)

                # -- subpixel fit (bayesicfitting) ------------------
                fit = fit_peak(corr, int_peak, patch=patch)

                # -- QA x-corr panel --------------------------------
                plot_xcorr_panel(
                    corr, int_peak, fit, subwcs,
                    f"{label} {region.name} {vlo}-{vhi}",
                    outdir / f"xcorr_{label}_{region.name}_{vlo}_{vhi}.png")

                if fit is None:
                    continue
                dx, dy, sdx, sdy, rho, _ = fit
                log.info("%s %s %s: dx=%.3f+-%.3f dy=%.3f+-%.3f "
                         "rho=%.2f (pixels)",
                         label, region.name, (vlo, vhi),
                         dx, sdx, dy, sdy, rho)

                # -- sky -> focal plane -----------------------------
                # the FLOATING telescope's location is the key: the
                # gondola's own position (GON_LAT, GON_LON, GON_ALT)
                # + epoch from this region's SEQ_ID define the AltAz
                # frame; boresight az/el are optional sanity checks
                gon_lat, gon_lon, alt, obstime, az_chk, el_chk = \
                    get_scan_attitude(l1_dir, args.band, seq_id)
                location = make_location(gon_lat, gon_lon, alt)
                log.info("%s %s: gondola at (%.3f, %.3f) deg, "
                         "%.0f m, %s", region.name, label,
                         gon_lat.value, gon_lon.value, alt.value,
                         obstime.isot)
                M_frame = galactic_to_altaz_matrix(
                    obstime, location,
                    region.lon.to_value(u.deg),
                    region.lat.to_value(u.deg))
                # gondola roll about the boresight: fixed relative to
                # the focal plane -> one unknown constant theta0
                # (config `roll` in deg; 0 = focal x along +Az)
                roll = np.radians(float(opts.get("roll", 0.0)))
                # tangential sky offsets (east, north), cos(lat)-projected
                dlon, dlat = pixshift_to_skyoffset(
                    dx, dy, ref_wcs, pixbox,
                    region.lat.to_value(u.deg))
                fx, fy, R = skyoffset_to_focal(dlon, dlat, M_frame,
                                               roll)
                C_sky = np.array([[sdx ** 2, rho * sdx * sdy],
                                  [rho * sdx * sdy, sdy ** 2]])
                C_f = cov_rotate(C_sky, R)
                # PIX_BEAM oversampling: pixels are not independent
                # (see load of fid_hdr) -> inflate the formal
                # covariance so the reported sigmas are honest
                C_f = C_f * float(pix_beam ** 2)
                # fx, fy and C_f are computed above from the attitude
                # chain (placeholder removed)

                measurements.append(dict(
                    band=args.band, target=args.target,
                    mixer=label, region=region.name,
                    vel_lo=vlo, vel_hi=vhi,
                    # epoch + gondola position: essential bookkeeping
                    # when COMBINING runs across fields/epochs (and
                    # for solving the roll theta0 later)
                    obstime=obstime.isot,
                    gon_lat=gon_lat, gon_lon=gon_lon, gon_alt=alt,
                    dx_pix=dx, dy_pix=dy,
                    dlon=dlon, dlat=dlat,
                    # offsets from the fiducial in DEGREES, for
                    # calibration use
                    dlon_deg=dlon.to(u.deg), dlat_deg=dlat.to(u.deg),
                    # focal-plane offsets from the fiducial in DEGREES
                    fx=fx.to(u.deg), fy=fy.to(u.deg),
                    snr_corr=snr_corr,
                    sig_fx=(np.sqrt(C_f[0, 0]) * u.arcsec).to(u.deg),
                    sig_fy=(np.sqrt(C_f[1, 1]) * u.arcsec).to(u.deg),
                    rho_f=C_f[0, 1] / np.sqrt(C_f[0, 0] * C_f[1, 1]),
                ))

    # -- fiducial mixer: zero offset by definition -------------------
    # (angles in the focal plane trivially zero as well)
    # -- fiducial mixer: zero offset by definition -------------------
    # obstime/position columns are placeholders ("-"/NaN) here:
    # the fiducial is the REFERENCE, so no AzEl conversion is done
    # for it -- the epoch + gondola position apply to the MIXER
    # measurement rows only. (Kept as columns so the table is
    # rectangular; treat "-" as "not applicable", not "unknown".)
    measurements.append(dict(
        band=args.band, target=args.target, mixer=fid_label,
        region="-", vel_lo=0, vel_hi=0,
        obstime="-", gon_lat=np.nan * u.deg, gon_lon=np.nan * u.deg,
        gon_alt=np.nan * u.m,
        dx_pix=0.0, dy_pix=0.0, dlon=0 * u.arcsec, dlat=0 * u.arcsec,
        dlon_deg=0 * u.deg, dlat_deg=0 * u.deg,
        fx=0 * u.deg, fy=0 * u.deg,
        snr_corr=np.nan,
        sig_fx=0.0 * u.deg, sig_fy=0.0 * u.deg, rho_f=0.0,
    ))

    # -- combine per mixer ------------------------------------------
    summaries = {}
    mixers = sorted({m["mixer"] for m in measurements})
    for label in mixers:
        rows = [m for m in measurements if m["mixer"] == label]
        # within one run/epoch: slices are correlated (overlapping
        # vel ranges) -- cap the mean's sigma at the floor
        summaries[label] = combine_measurements(rows,
                                                floor_cap=True)
        log.info("%s: %s", label, summaries[label])

    save_results(measurements, summaries, fid_label, outdir)
    plot_summary(summaries, outdir / "summary_offsets.png")
    # measured offsets in the calibration-file format (PIX/AZ/EL).
    # pix_band: which array this run's band belongs to (B1/B2) --
    # from config `pix_band` (e.g. "B2" when the fiducial is B2M8)
    # measured DELTA offsets in the calibration-file format
    # (PIX/AZ/EL). These are residual corrections to the APPLIED
    # calibration (rows flagged AS_DELTA); the pipeline consumes
    # only AS_MEASURED rows, which appear in the updated file.
    # pix_band comes from the BAND config (NII->B1, CII->B2, ...)
    pix_band = band_cfg.get("pix_band", opts.get("pix_band", "B2"))
    write_calibration(summaries, fid_label, pix_band,
                      outdir / f"offsets_{args.band}_deltas.txt")
    # -- update the master calibration file with the residuals ------
    # cal_file: path to the [offsets] file whose offsets were
    # APPLIED when the L1 (hence these L2 cubes) were made. The
    # xcorr offsets are residuals; the updated file is written to
    # the output dir (the original is never modified in place).
    cal_file = opts.get("cal_file")
    # NOTE: no per-run updated-cal file is written: the combined
    # calibration (offsets_cal.txt, applied + combined residuals)
    # is produced ONCE by --aggregate over all runs; per-region
    # cal files were not needed.
    cal_file = opts.get("cal_file")
    if cal_file:
        log.info("cal_file set: run --aggregate to fold this "
                 "run's residuals into offsets_cal.txt")
    return 0


def pixshift_to_skyoffset(dx, dy, wcs, pixbox, lat_deg):
    """Convert a (subpixel) pixel shift on the shared grid to
    TANGENTIAL sky offsets (east, north; arcsec) at the region center.

    dlon is the tangential EAST offset: the raw longitude
    difference is projected by cos(lat), so (dlon, dlat) are true
    angular offsets on the sky -- matching what
    galactic_to_altaz_matrix (spherical-offset probes) and
    skyoffset_to_focal expect. TODO: use the full WCS (CD matrix,
    projection) instead of CDELT for CD/rotated grids.
    """
    cd = wcs.celestial.proj_plane_pixel_scales()
    # convert to plain floats before use (np.abs on Quantity arrays
    # raises -- see skybox_to_pixbox)
    dlon = (dx * abs(cd[0].to_value(u.arcsec))
            * np.cos(np.radians(lat_deg))) * u.arcsec
    dlat = dy * abs(cd[1].to_value(u.arcsec)) * u.arcsec
    return dlon, dlat


if __name__ == "__main__":
    sys.exit(main())