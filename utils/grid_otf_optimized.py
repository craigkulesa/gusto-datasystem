# Copyright (C) 2015 Associated Universities, Inc. Washington DC, USA.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
# Correspondence concerning GBT software should be addressed as follows:
#       GBT Operations
#       National Radio Astronomy Observatory
#       P. O. Box 2
#       Green Bank, WV 24944-0002 USA

import math
import numpy as np
import numpy.ma as ma
import sys
import scipy.special


def grid_otf(data, xsky, ysky, wcsObj, nchan, xsize, ysize, pix_scale, beam_fwhm,
             weight=None, kern="gaussbessel", gauss_fwhm=None, verbose=4):
    """
    Grid individual spectra onto a specified regular grid following the
    recommendations of Mangum et al. (2007).  Written to be of general use
    in gridding OTF data.

    Adapted from IDL code provided by Adam Leroy (aleroy@nrao.edu).

    Inputs:
    data        - nspec × nchan 2-D array of spectra
    xsky        - nspec vector of X sky positions (deg)
    ysky        - nspec vector of Y sky positions (deg)
    wcsObj      - WCS object suitable for gridding xsky/ysky
    weight      - (optional) nspec weight vector; equal weights assumed if None
    beam_fwhm   - telescope beam FWHM in decimal degrees
    gauss_fwhm  - Gaussian kernel FWHM (deg); used only when kern="gauss";
                  defaults to beam_fwhm/3
    kern        - gridding kernel: "gaussbessel" | "gauss" | "nearest"
    verbose     - verbosity level

    Returns: (cube, weight_cube, beam_fwhm) or (None, None, None) on failure.

    Optimisation notes vs. the original:
    - The O(nx*ny*nspec) Python double loop is replaced by vectorised NumPy
      operations.  For each spectrum we accumulate its contribution to all
      pixels within r_support in one shot using boolean masks and array
      arithmetic.
    - The per-pixel interp1d call is replaced by a direct integer-index look-up
      into the pre-computed convolution table.
    - NaN masking is hoisted outside the loop; a NaN-weight copy of `data` is
      prepared once so the inner loop never pays the cost of np.isnan.
    - For the "nearest" kernel the inner body is a single scatter-add with
      np.round, avoiding all distance calculations.
    """

    result = (None, None, None)

    # ------------------------------------------------------------------ #
    # Argument checking                                                    #
    # ------------------------------------------------------------------ #
    if len(data.shape) != 2 or len(xsky.shape) != 1 or len(ysky.shape) != 1:
        if verbose > 1:
            print("data, sky coordinates have unexpected shapes")
            print("data : ", data.shape)
            print("xsky : ", xsky.shape)
            print("ysky : ", ysky.shape)
        return result

    nspec, nchan_data = data.shape

    if nspec == 0 or nchan_data == 0:
        if verbose > 1:
            print("no data given")
        return result

    if nspec != len(xsky) or nspec != len(ysky):
        if verbose > 1:
            print("Number of sky position values does not match number of spectra in data")
        return result

    if kern not in ("gaussbessel", "gauss", "nearest"):
        if verbose > 1:
            print("kern must be one of gaussbessel, gauss, or nearest")
        return result

    cubeShape = (nchan, ysize, xsize)

    if cubeShape[0] != nchan_data:
        if verbose > 1:
            print("Frequency axis in target header and spectra length do not match")
        return result

    print('Masked arrays', ma.is_masked(data))
    print('Total spectra to grid:', nspec)

    if weight is None:
        weight = np.ones(nspec, dtype=np.float32)
    else:
        weight = np.asarray(weight, dtype=np.float32)

    # ------------------------------------------------------------------ #
    # Kernel parameters                                                    #
    # ------------------------------------------------------------------ #
    if kern == "gauss":
        if gauss_fwhm is None:
            gauss_fwhm = beam_fwhm / 3.0
        r_support    = 5.0 * pix_scale
        max_conv_fn  = 1.0
        cutoff_conv_fn = 0.0
        scale_fwhm   = math.sqrt(beam_fwhm**2 + gauss_fwhm**2) / beam_fwhm

    elif kern == "gaussbessel":
        a = 1.55 * beam_fwhm / 3.0
        b = 2.52 * beam_fwhm / 3.0
        r_support    = 3.0 * pix_scale
        max_conv_fn  = 0.5
        cutoff_conv_fn = 0.0
        scale_fwhm   = 1.09

    else:  # nearest
        r_support      = 1.0          # not used in distance sense
        scale_fwhm     = 1.0
        cutoff_conv_fn = 0.0
        max_conv_fn    = 1.0

    # ------------------------------------------------------------------ #
    # Pre-compute convolution look-up table                                #
    # ------------------------------------------------------------------ #
    # We index into this table with  idx = int(dist_sqrd / pre_delta_sqrd)
    # avoiding per-call interp1d overhead.
    if kern != "nearest":
        n_pre         = 10000
        r_support_sqrd   = r_support ** 2
        pre_delta_sqrd   = r_support_sqrd / n_pre
        pre_dist_sqrd    = np.arange(int(n_pre * 1.01 + 1)) * pre_delta_sqrd
        pre_dist         = np.sqrt(pre_dist_sqrd)

        if kern == "gauss":
            pre_conv_fn = max_conv_fn * np.exp(
                -0.5 * (pre_dist / (gauss_fwhm / 2.354)) ** 2
            )
        else:  # gaussbessel
            x = np.pi * pre_dist / a
            x[x == 0.0] = 1e-5          # avoid divide-by-zero
            pre_conv_fn = (scipy.special.j1(x) / x) * np.exp(-(pre_dist / b) ** 2)
            pre_conv_fn[0] = max_conv_fn

        # Convert to float32 for memory efficiency
        pre_conv_fn = pre_conv_fn.astype(np.float32)
        pre_conv_fn_len = len(pre_conv_fn)

        r_support_pix     = r_support / pix_scale
        r_support_pix_sqrd = r_support_pix ** 2
        # Index of the "cap" distance (first non-zero pre_dist_sqrd entry)
        cap_idx = 1

    # ------------------------------------------------------------------ #
    # Pixel coordinates for every spectrum                                 #
    # ------------------------------------------------------------------ #
    zeros = np.zeros(len(xsky))
    x_pix, y_pix, _ = wcsObj.wcs_world2pix(xsky, ysky, zeros, 0)
    x_pix = x_pix.astype(np.float32)
    y_pix = y_pix.astype(np.float32)

    # ------------------------------------------------------------------ #
    # Output arrays                                                        #
    # ------------------------------------------------------------------ #
    data_cube   = np.zeros(cubeShape, dtype=np.float32)
    weight_cube = np.zeros(cubeShape, dtype=np.float32)

    nx = cubeShape[2]
    ny = cubeShape[1]

    # ------------------------------------------------------------------ #
    # Pre-process data: replace NaN/masked values with 0                  #
    # combined_mask[s, c] = True if channel c of spectrum s is bad        #
    # ------------------------------------------------------------------ #
    data_work = np.array(data, dtype=np.float32)          # strips ma mask, gets fill values
    nan_mask  = np.isnan(data_work)                        # (nspec, nchan)

    # Honor masked array mask if present (ma.getmaskarray always returns a
    # full boolean array, never the scalar False that ma.getmask can return)
    if ma.is_masked(data):
        combined_mask = nan_mask | np.asarray(ma.getmaskarray(data), dtype=bool)
    else:
        combined_mask = nan_mask

    data_work[combined_mask] = 0.0

    # per-spectrum-channel weight: starts at weight[s], zeroed where bad
    # shape: (nspec, nchan)
    spec_weight = np.broadcast_to(
        weight[:, np.newaxis], (nspec, nchan_data)
    ).copy()                                               # writeable copy
    spec_weight[combined_mask] = 0.0

    # ------------------------------------------------------------------ #
    # "nearest" kernel — fast scatter path                                 #
    # ------------------------------------------------------------------ #
    if kern == "nearest":
        # Round pixel coords to integer pixel centres
        xi = np.round(x_pix).astype(np.int32)
        yi = np.round(y_pix).astype(np.int32)

        # Valid spectra: pixel inside cube footprint
        valid = (xi >= 0) & (xi < nx) & (yi >= 0) & (yi < ny)

        for s in np.where(valid)[0]:
            i, j = xi[s], yi[s]
            w_vec = spec_weight[s]                # (nchan,)
            wsum  = w_vec.sum()
            if wsum <= 0.0:
                continue
            data_cube[:, j, i]   += data_work[s] * w_vec
            weight_cube[:, j, i] += w_vec

        # Normalise accumulated numerator
        nonzero = weight_cube > 0.0
        data_cube[nonzero] /= weight_cube[nonzero]

        beam_fwhm = scale_fwhm * beam_fwhm
        return (data_cube, weight_cube, beam_fwhm)

    # ------------------------------------------------------------------ #
    # Convolution kernels: iterate over spectra (outer loop = nspec)      #
    #                                                                      #
    # For each spectrum we compute the distance to every pixel within its  #
    # support radius and splat its weighted contribution onto those pixels. #
    # This is O(nspec * n_pix_in_support) instead of O(nx*ny*nspec).      #
    # ------------------------------------------------------------------ #
    if verbose > 3:
        print("Gridding %d spectra ..." % nspec)

    for s in range(nspec):
        if verbose > 3 and nspec > 1:
            sys.stdout.write("\rSpectrum %d / %d" % (s + 1, nspec))
            sys.stdout.flush()

        xs = float(x_pix[s])
        ys = float(y_pix[s])
        ws = weight[s]
        if ws <= 0.0:
            continue

        # Bounding box of support in pixel space
        i_lo = max(0,  int(math.floor(xs - r_support_pix)))
        i_hi = min(nx, int(math.ceil( xs + r_support_pix)) + 1)
        j_lo = max(0,  int(math.floor(ys - r_support_pix)))
        j_hi = min(ny, int(math.ceil( ys + r_support_pix)) + 1)

        if i_lo >= i_hi or j_lo >= j_hi:
            continue

        # Pixel coordinate arrays for this bounding box
        ii = np.arange(i_lo, i_hi, dtype=np.float32)   # (ni,)
        jj = np.arange(j_lo, j_hi, dtype=np.float32)   # (nj,)

        # Squared pixel distances:  shape (nj, ni)
        dist2 = (ii[np.newaxis, :] - xs) ** 2 + (jj[:, np.newaxis] - ys) ** 2

        # Mask pixels inside support circle
        inside = dist2 <= r_support_pix_sqrd
        if not inside.any():
            continue

        # Compute convolution weight for each pixel inside support
        # index into pre_conv_fn: idx = dist2 / (pre_delta_sqrd / pix_scale^2)
        pix_delta_sqrd = pre_delta_sqrd / (pix_scale ** 2)
        raw_idx = (dist2[inside] / pix_delta_sqrd).astype(np.int32)
        raw_idx = np.clip(raw_idx, 0, pre_conv_fn_len - 1)

        conv_vals = pre_conv_fn[raw_idx]                 # (n_inside,)
        # Cap very close pixels
        conv_vals[raw_idx < cap_idx] = max_conv_fn

        # Combined weight: conv * per-spectrum weight scalar
        cw = conv_vals * ws                               # (n_inside,)

        # Flat indices of pixels inside the bounding box that are inside support
        jj_idx, ii_idx = np.where(inside)
        j_abs = jj_idx + j_lo
        i_abs = ii_idx + i_lo

        # For each inside pixel accumulate contribution across all channels
        # data_work[s]: (nchan,), spec_weight[s]: (nchan,)
        # cw[k] is the spatial weight for pixel k
        # We add cw[k] * data_work[s] * spec_weight[s] to data_cube[:,j,i]
        # and    cw[k] * spec_weight[s]              to weight_cube[:,j,i]
        #
        # Vectorised: reshape cw to (n_inside, 1), broadcast with (1, nchan)
        cw2d    = cw[:, np.newaxis]                       # (n_inside, 1)
        sw_row  = spec_weight[s][np.newaxis, :]           # (1, nchan)
        d_row   = data_work[s][np.newaxis, :]             # (1, nchan)

        wt_contrib   = cw2d * sw_row                      # (n_inside, nchan)
        data_contrib = cw2d * sw_row * d_row              # (n_inside, nchan)

        # Scatter-add into output arrays
        np.add.at(weight_cube, (slice(None), j_abs, i_abs),
                  wt_contrib.T)
        np.add.at(data_cube,   (slice(None), j_abs, i_abs),
                  data_contrib.T)

    if verbose > 3:
        print('')

    # ------------------------------------------------------------------ #
    # Normalise: data / weight                                             #
    # ------------------------------------------------------------------ #
    nonzero = weight_cube > 0.0
    data_cube[nonzero] /= weight_cube[nonzero]

    beam_fwhm = scale_fwhm * beam_fwhm
    return (data_cube, weight_cube, beam_fwhm)
