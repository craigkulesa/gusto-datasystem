import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
from scipy.optimize import minimize
from scipy.ndimage import map_coordinates

def load_fits(file_path):
    """Load a FITS file and return data and WCS."""
    with fits.open(file_path) as hdul:
        data = hdul[0].data
        wcs = WCS(hdul[0].header)
    return data, wcs

def click_rectangle(ax, title):
    """Let the user select two corners of a rectangle."""
    plt.title(title)
    coords = plt.ginput(2, timeout=60)  # Wait for user to click two points
    plt.close()
    return coords

def crop_image(data, wcs, coords):
    """Crop the image and WCS to the selected rectangular region."""
    x1, y1 = coords[0]
    x2, y2 = coords[1]
    x_min, x_max = int(min(x1, x2)), int(max(x1, x2))
    y_min, y_max = int(min(y1, y2)), int(max(y1, y2))
    
    cropped_data = data[y_min:y_max, x_min:x_max]
    cropped_wcs = wcs.deepcopy()
    cropped_wcs.wcs.crpix[0] -= x_min
    cropped_wcs.wcs.crpix[1] -= y_min
    
    return cropped_data, cropped_wcs

from scipy.ndimage import map_coordinates

def align_images(data1, wcs1, data2, wcs2, initial_shift=(0, 0)):
    """Optimize the alignment of two images by minimizing their subtraction."""

    def misalignment_error(shift):
        """Compute error as the sum of absolute differences after shifting."""
        lon_shift, lat_shift = shift
        wcs2_shifted = wcs2.deepcopy()
        wcs2_shifted.wcs.crval[0] += lon_shift
        wcs2_shifted.wcs.crval[1] += lat_shift

        # Generate pixel coordinates for the second image
        y_indices, x_indices = np.indices(data2.shape)
        world_coords = wcs2_shifted.wcs_pix2world(x_indices, y_indices, 0)

        # Map these world coordinates to pixel coordinates in the first image
        pix_coords = wcs1.wcs_world2pix(world_coords[0], world_coords[1], 0)

        # Interpolate data1 at the mapped pixel coordinates
        aligned_data = map_coordinates(data1, [pix_coords[1], pix_coords[0]], order=1, mode='constant', cval=0.0)

        # Compute the sum of absolute differences as the alignment error
        error = np.sum(np.abs(data2 - aligned_data))
        return error


    result = minimize(misalignment_error, initial_shift, method='Nelder-Mead')
    return result.x

def calculate_az_el_offsets(observer_time, observer_lat, observer_lon, observer_alt, gal_lon, gal_lat):
    """Calculate azimuth and elevation offsets."""
    observer_location = EarthLocation(lat=observer_lat*u.deg, lon=observer_lon*u.deg, height=observer_alt*u.m)
    observer_time = Time(observer_time)

    lon = (351+gal_lon)*u.deg
    print(lon)
    lat = (0.5+gal_lat)*u.deg
    print(lat)
    gal_coords1 = SkyCoord(l=lon, b=lat, frame='galactic')
    altaz_frame1 = AltAz(location=observer_location, obstime=observer_time)
    altaz_coords1 = gal_coords1.transform_to(altaz_frame1)

    lon = (351)*u.deg
    print(lon)
    lat = (0.5)*u.deg
    print(lat)
    gal_coords2 = SkyCoord(l=lon, b=lat, frame='galactic')
    altaz_frame2 = AltAz(location=observer_location, obstime=observer_time)
    altaz_coords2 = gal_coords2.transform_to(altaz_frame2)

    return (altaz_coords1.az.deg-altaz_coords2.az.deg), (altaz_coords1.alt.deg-altaz_coords2.alt.deg)

def main():
    # Load FITS files
    file1 = input("Enter the path to the first FITS file: ") or "../level2/src/align_b2m8.fits"
    file2 = input("Enter the path to the second FITS file: ") or "../level2/src/align_b2m5.fits"
    data1, wcs1 = load_fits(file1)
    data2, wcs2 = load_fits(file2)

    # Plot the first image for rectangle selection
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.imshow(data1, origin='lower', cmap='gray')
    coords = click_rectangle(ax, "Click two corners to define a rectangle on the first image")
    print("Selected rectangle coordinates:", coords)

    # Crop both images
    cropped_data1, cropped_wcs1 = crop_image(data1, wcs1, coords)
    cropped_data2, cropped_wcs2 = crop_image(data2, wcs2, coords)

    # Perform alignment
    shift = align_images(cropped_data1, cropped_wcs1, cropped_data2, cropped_wcs2)
    print("Optimal Galactic longitude/latitude shifts (degrees):", shift)

    # Get observer info
    observer_time = input("Enter observer time (YYYY-MM-DD HH:MM:SS): ") or "2024-02-04 13:02:32"
    print(observer_time)

    try:
        observer_lat = float(input("Enter observer latitude (degrees): "))
    except ValueError:
        observer_lat = -74.735938
    print(observer_lat)

    try:
        observer_lon = float(input("Enter observer longitude (degrees): "))
    except ValueError:
        observer_lon = 61.630017
    print(observer_lon)

    try:
        observer_alt = float(input("Enter observer altitude (meters): "))
    except ValueError:
        observer_alt = 35172
    print(observer_alt)

    # Calculate azimuth and elevation offsets
    az, el = calculate_az_el_offsets(observer_time, observer_lat, observer_lon, observer_alt, shift[0], shift[1])
    print(f"Azimuth offset: {az:.4f} degrees")
    print(f"Elevation offset: {el:.4f} degrees")

if __name__ == "__main__":
    main()

