#!/usr/bin/env python3

from astropy.io import fits
from astropy import units as u
import numpy as np
from astropy.coordinates import Angle


# Open FITS file
hdu_list = fits.open('NGC7538_CII_subcube.fits')
data = hdu_list[0].data
header = hdu_list[0].header

# Check for NaN values in the data and handle them
intensity = np.nan_to_num(data)  # Convert NaNs to zero

# Extract frequency axis information
frequency = header['CRVAL3'] + np.arange(header['NAXIS3']) * header['CDELT3']  # Frequency axis
pixel_area_sr = (Angle(header['CDELT1'], unit='deg').rad * Angle(header['CDELT2'], unit='deg').rad)  # Convert to steradians

# Integrate intensity over the frequency axis to get flux density
flux_density = np.nansum(intensity, axis=0)  # Summing over the frequency axis

# Calculate the total flux (flux density * pixel area in steradians)
total_flux = np.nansum(flux_density * pixel_area_sr)

# Debugging output for flux
print(f"Total Flux: {total_flux}")

# Convert total flux to luminosity
distance = 2.65 * 10**3 * u.pc  # Distance to NGC 7538 in parsecs
distance_cm = distance.to(u.cm).value  # Convert parsecs to cm

# Ensure total_flux is positive
if total_flux <= 0:
	print("Total flux is zero or negative. Please check the intensity values.")
	hdu_list.close()
	exit()
	
luminosity = total_flux * 4 * np.pi * distance_cm**2  # in erg/s

# Convert luminosity from erg/s to watts
luminosity_watts = luminosity * 1e-7

# Convert watts to solar luminosities (1 L_sun = 3.828 x 10^26 W)
solar_luminosity = luminosity_watts / (3.828 * 10**26)

# Print the result
print(f'[C II] Luminosity: {solar_luminosity:.2e} L_sun')

# Close the FITS file
hdu_list.close()