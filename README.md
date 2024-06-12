
# [C II] Luminosity Calculation for NGC 7538

This repository contains a Python script to calculate the [C II] luminosity of the NGC 7538 region using data from a FITS file. The script performs unit conversions, handles the data, and calculates the luminosity in both Watts and solar luminosities.

## Requirements

- Python 3.x
- Astropy
- Numpy

You can install the necessary Python packages using:

```sh
pip install astropy numpy
```

## Script Overview

The script `Cplus_luminosity.py` performs the following steps:

1. **Open the FITS File**: Reads the FITS file and extracts the data and header.
2. **Extract Intensity Data**: Intensity data is assumed to be in units of K km/s.
3. **Convert Intensity to Flux Density**: Converts intensity to flux density in erg/s/cm²/sr using a conversion factor.
4. **Convert Flux Density to W/m²/sr**: Uses Astropy's unit conversions.
5. **Handle NaN or Infinite Values**: Replaces NaNs or infinite values with zeros.
6. **Calculate Beam Area**: Calculates the beam area in steradians.
7. **Sum Flux Density Over All Pixels**: Sums the flux density to get the total flux density in W/m²/sr.
8. **Calculate Total Flux in W/m²**: Multiplies the total flux density by the beam area.
9. **Calculate [C II] Luminosity**: Converts the total flux to luminosity in Watts using the distance to NGC 7538.
10. **Convert to Solar Luminosities**: Converts the luminosity to solar luminosities.

## Detailed Steps

### 1. Open the FITS File

```python
from astropy.io import fits

# Open the FITS file
hdu_list = fits.open('NGC7538_CII_moment0_small.fits')
data = hdu_list[0].data
header = hdu_list[0].header
```

### 2. Extract Intensity Data

The intensity data in the FITS file is assumed to be in units of K km/s.

```python
intensity = data  # The data is in K km/s
```

### 3. Convert Intensity to Flux Density

```python
conversion_factor = 7.0e-6  # K km/s to erg/s/cm²/sr
flux_density_erg = intensity * conversion_factor  # Now in erg/s/cm²/sr
```

### 4. Convert Flux Density to W/m²/sr

```python
from astropy import units as u

flux_density_wm2_sr = (flux_density_erg * u.erg / (u.s * u.cm**2 * u.sr)).to(u.W / (u.m**2 * u.sr))
```

### 5. Handle NaN or Infinite Values

```python
flux_density_wm2_sr = np.nan_to_num(flux_density_wm2_sr.value) * (u.W / (u.m**2 * u.sr))
```

### 6. Calculate Beam Area

The beam area in steradians needs to be calculated from the beam major (`BMAJ`) and minor (`BMIN`) axes provided in the header.

```python
beam_area_deg2 = (header['BMAJ'] * header['BMIN'] * (np.pi / (4 * np.log(2)))) * u.deg**2
beam_area_sr = beam_area_deg2.to(u.sr)  # Convert to steradians
```

### 7. Sum Flux Density Over All Pixels

```python
total_flux_density_sr = np.sum(flux_density_wm2_sr)
```

### 8. Calculate Total Flux in W/m²

```python
total_flux_wm2 = total_flux_density_sr * beam_area_sr
```

### 9. Calculate [C II] Luminosity

Using the distance to NGC 7538, convert the total flux to luminosity in Watts.

```python
distance = 2.65 * 10**3 * u.parsec  # Distance to NGC 7538
luminosity_watts = total_flux_wm2 * 4 * np.pi * (distance.to(u.m))**2  # Luminosity in Watts
print('[C II] Luminosity in Watt:', luminosity_watts)
```

### 10. Convert to Solar Luminosities

```python
solar_luminosity = 3.828 * 10**26 * u.W  # Solar luminosity in Watts
luminosity_solar = luminosity_watts / solar_luminosity

# Print the [C II] Luminosity in solar luminosities
print(f'[C II] Luminosity: {luminosity_solar.to(u.dimensionless_unscaled)} L_sun')
```

## Running the Script

Ensure that the FITS file `NGC7538_CII_moment0_small.fits` is in the same directory as the script. Run the script using:

```sh
python Cplus_luminosity.py
```

## References

- [Goicoechea et al. 2015](https://iopscience.iop.org/article/10.1088/0004-637X/812/1/75/pdf)
- Astropy Documentation: [https://docs.astropy.org/]

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
