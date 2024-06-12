#!/usr/bin/env python3

from astropy.io import fits
from astropy import units as u
import numpy as np

# Open FITS file
hdu_list = fits.open('NGC7538_CII_moment0_small.fits')
data = hdu_list[0].data
header = hdu_list[0].header

# Extract intensity (Data is in K km/s)
intensity = data  # Replace with correct data extraction if necessary

# Convert intensity to flux density in erg/s/cm²/sr using the conversion factor
conversion_factor = 7.0e-6  # K km/s to erg/s/cm²/sr from Goicoechea et al. 2015
flux_density_erg = intensity * conversion_factor  # Now in erg/s/cm²/sr

# Convert flux density to W/m²/sr (to apply 1e-3)
flux_density_wm2_sr = (flux_density_erg * u.erg / (u.s * u.cm**2 * u.sr)).to(u.W / (u.m**2 * u.sr)) 

# Check for NaN or infinite values
flux_density_wm2_sr = np.nan_to_num(flux_density_wm2_sr)

# Beam area in deg²
beam_area_deg2 = (header['BMAJ'] * header['BMIN'] * (np.pi / (4 * np.log(2)))) * u.deg**2
# Convert beam area to steradians
beam_area_sr = beam_area_deg2.to(u.sr)

# Sum flux density over all pixels to get total flux density in W/m²/sr
total_flux_density_sr = np.sum(flux_density_wm2_sr)

# Total flux in W/m²
total_flux_wm2 = total_flux_density_sr * beam_area_sr

# Calculate [C II] Luminosity in Watts
distance = 2.65 * 10**3 * u.parsec  # Distance to NGC 7538
luminosity_watts = total_flux_wm2 * 4 * np.pi * (distance.to(u.m))**2  # Luminosity in Watts
#print('[C II] Luminosity in Watt:',luminosity_watts)

# Convert to solar luminosities
solar_luminosity = 3.828 * 10**26 * u.W  # Solar luminosity in Watts
luminosity_solar = luminosity_watts / solar_luminosity

# Print the [C II] Luminosity in solar luminosities
print(f'[C II] Luminosity: {luminosity_solar.to(u.dimensionless_unscaled)} L_sun')
