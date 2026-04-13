# ALMA Molecular Spectra Analysis — AGB & post-AGB Stars

Analysis of molecular spectral lines from AGB and post-AGB stars observed with ALMA, with the goal of deriving physical parameters such as mass-loss rates and characterizing the optical depth regime of the circumstellar envelope.

---

## Overview

This project processes and analyzes calibrated ALMA spectra for a sample of evolved stars. The workflow covers the full data analysis pipeline: from raw spectral data ingestion and unit conversion, through line profile fitting, to the derivation of physical parameters.

The main target used for development and testing is **IRAS 06531-0216**, with the analysis designed to be extensible to the full sample.

---

## Scientific Goals

- Identify and characterize molecular emission lines (CO, ¹³CO, SO, CN, HC₃N) in ALMA Band 3 (~3mm) and Band 6 (~1.3mm) spectra
- Convert spectra from flux density (Jy/beam) to main beam brightness temperature (T_mb)
- Fit spectral line profiles to extract peak temperature, systemic velocity and expansion velocity
- Calculate line intensity ratios to assess the optical depth regime of each source
- Derive mass-loss rates for the stellar sample using both optically thin and optically thick prescriptions

---

## Repository Structure

```
.
├── all_spectra/                  # Raw ALMA spectra (.dat files, one per source)
├── Spectra_ALMA.ipynb            # Main analysis notebook
├── Mass_loss_IRAS_06531.ipynb    # Mass-loss rate derivation for IRAS 06531
└── README.md
```

---

## Methods

### Spectral Line Identification
Spectral windows are defined in GHz based on known molecular rest frequencies. The data covers two ALMA bands separated by a gap (~116–211 GHz):

| Band | Frequency Range | Wavelength | Lines |
|------|----------------|------------|-------|
| Band 3 | 85 – 116 GHz | ~3 mm | CO(1-0), SO, CN, HC₃N |
| Band 6 | 211 – 274 GHz | ~1.3 mm | CO(2-1), ¹³CO(2-1) |

### T_mb Conversion
Flux density is converted to main beam brightness temperature using:

$$T_{mb} = \frac{1.222 \times 10^{3} \cdot S_\nu}{\nu^2 \cdot \theta_{maj} \cdot \theta_{min}}$$

where $S_\nu$ is in Jy/beam, $\nu$ in GHz, and $\theta$ in arcseconds. Each spectral line uses its own beam parameters from the FITS header.

### Line Profile Fitting
A generalized flat-top parabolic profile is fitted to each line:

$$T(v) = T_{peak} \cdot \left(1 - \left(\frac{v - V_c}{V_e}\right)^2\right)^\beta$$

The exponent $\beta$ is a free parameter that controls the shape of the profile top: $\beta = 1$ recovers the classical parabola (optically thin shell), while $\beta < 1$ produces a flat-top profile consistent with optically thick emission or resolved shells.

### Intensity Ratios
Line intensities are computed by numerical integration of T_mb over velocity:

$$I = \int T_{mb} \, dv \quad [\text{K·km/s}]$$

Two diagnostic ratios are computed:
- **R = I(CO 2-1) / I(CO 1-0)**: traces the excitation and optical depth of ¹²CO
- **R = I(¹²CO 2-1) / I(¹³CO 2-1)**: traces the isotopic line opacity

### Mass-Loss Rate
Two prescriptions are implemented depending on the optical depth regime:

**Optically thin:**
$$\dot{M} = 4.55 \times 10^{-19} \left(\frac{T_{mb}}{\log(W/0.04) \cdot s(J)}\right)^{5/6} f_{CO}^{-1} \, V_{exp}^{11/6} \, (D\theta)^{5/3}$$

**Optically thick:**
$$\dot{M} = 1.4 \, \frac{T_{mb} \, V_e^2 \, D^2 \, \theta^2}{2 \times 10^{19} \, f_{CO}^{0.85} \, s(J)}$$

---

## Dependencies

```
python >= 3.10
numpy
matplotlib
scipy
astropy
pandas
jupyter
```

Install with:
```bash
pip install numpy matplotlib scipy astropy pandas jupyter
```

---

## Stellar Sample

| Source | Distance (pc) | V_lsr (km/s) |
|--------|--------------|--------------|
| IRAS 06531-0216 | 1158.35 | 39 |
| IRAS 07145-1428 | — | 88 |
| IRAS 07180-1314 | — | 58 |
| IRAS 08500-3254 | — | 58 |
| R Hor | — | 38 |
| RV Aqr | — | 2 |
| RZ Sgr | — | -29.5 |
| UU For | — | -7 |
| V510 Pup | — | 45 |
| V743 Mon | — | — |

---

## Status

- [x] Spectral line visualization (Band 3 & Band 6)
- [x] T_mb conversion with per-line beam correction
- [x] Flat-top profile fitting with β parameter
- [x] Intensity ratio calculation (¹²CO and ¹³CO)
- [x] W factor (dilution factor) calculation
- [x] Mass-loss rate formulas (thin & thick)
- [ ] Full sample analysis
- [ ] Final parameter table

---

## Author

Undergraduate research project in observational astrophysics.
Data: ALMA calibrated spectra.
