# Biomass Internship Work

Research log and codebase for my MSc internship minithesis: investigating whether HH-VV
copolar phase difference (CPD) from ESA's Biomass P-band SAR can serve as a proxy for sea
surface salinity, with secondary workstreams in wave/sea-state characterisation and
cryospheric change detection. Maintained as a working record for the report and for my ESA
supervisor's reference.

### Contents
- [`notebooks/`](notebooks) — case study processing (Amazon, Ganges-Meghna, Guiana Shelf, El Niño Pacific)
- [`figures/`](figures) — output plots and figures by case study
- [`src/`](src) — reusable processing scripts (ocean masking, calibration, MAAP access, orientation fixes)
- [`docs/`](docs) — UCD internship expectations and supervisor information
- [`report/`](report) — internship report and poster

### Data
The Biomass satellite, launched April 2025, carries a P-band SAR (435 MHz) capable of deep
penetration through vegetation, dry sand and ice, making it uniquely suited to secondary
science investigations across diverse Earth surface types. All data is sourced from the ESA
MAAP Explorer platform (explorer.maap.eo.esa.int). Products are mostly Biomass Level 1a
S1_SCS__1S (Single Look Complex Stripmap), full polarimetry (HH, HV, VH, VV), P-band 435 MHz,
processed by the Biomass Central Processing Facility. Raw data files (TIFFs, binary products)
are not tracked in this repository due to size; acquisition guidance is in
[`src/MAAP_product_acquisition.md`](src/MAAP_product_acquisition.md).

### Code sources
All code in this repository is based on and taken directly from:
- Björn's L1a processing script (received 22 May)
- ESA Pol-InSAR course content and notebooks (DLR/ESA, 13 April – 17 July 2026)
- https://github.com/satim-co/PolSARpro.git

### Author
Orlaith Doyle, UCD MSc Space Science & Technology student.
Intern, ESA Climate Action, Sustainability and Science Department (1 May – 30 October 2026).
ESA Supervisor: Björn Rommen.

### References
- ESA Biomass Mission: https://www.esa.int/Applications/Observing_the_Earth/FutureEO/Biomass
- PolSARpro Python: https://github.com/satim-co/PolSARpro
- ESA MAAP Explorer: https://explorer.maap.eo.esa.int
