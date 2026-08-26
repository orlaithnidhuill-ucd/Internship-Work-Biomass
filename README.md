# Biomass Internship Work

Research log and codebase for my MSc internship minithesis: investigating whether HH-VV
copolar phase difference (CPD) from ESA's Biomass P-band SAR can serve as a proxy for sea
surface salinity, with secondary workstreams in wave/sea-state characterisation and
cryospheric change detection. Maintained as a working record for the report and for my ESA
supervisor's reference.

### Contents
- [`2. Notebooks and Study Case Cells`](2.%20Worked%20Notebooks) Case study processing (Amazon, Ganges-Meghna, Guiana Shelf, El Niño Pacific)
- [`3. Output Plots & Figures`](3.%20Output%20Plots%20&%20Figures) Output plots and figures by case study
- [`4. Reusable Processing Scripts`](4.%20Reusable%20Processing%20Scripts) Handy scripts for future work; Ocean masking, calibration, MAAP access, orientation fixes etc
- [`5. Full Report & Poster`](5.%20Full%20Report%20&%20Poster) Final internship report and poster
- [`6. UCD Internship Expectation Information`](6.%20UCD%20Internship%20Expectation%20Information) This is a dedicated module/placement detail folder for my supervisor's access

### Data
The Biomass satellite, launched April 2025, carries a P-band SAR (435 MHz) capable of deep penetration through vegetation, dry sand and ice, making it uniquely suited to secondary
science investigations across diverse Earth surface types. All data is sourced from the ESA MAAP Explorer platform (explorer.maap.eo.esa.int). Products are mostly Biomass Level 1a
S1_SCS__1S (Single Look Complex Stripmap), full polarimetry (HH, HV, VH, VV), P-band 435 MHz, processed by the Biomass Central Processing Facility. Raw data files (TIFFs, binary products)
are not tracked in this repository due to size, though acquisition guidance is in [`4. Reusable Processing Scripts`](4.%20Reusable%20Processing%20Scripts) under "MAAP Product Acquisition".

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
