# neonPLFA

**Continental-scale synthesis of NEON phospholipid fatty acid (PLFA) data**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21463376.svg)](https://doi.org/10.5281/zenodo.21463376)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This repository provides an open-source R workflow and interactive Shiny application for accessing, processing, and analyzing soil microbial community data from the National Ecological Observatory Network (NEON). The workflow processes phospholipid fatty acid (PLFA) measurements from 11,399 soil samples collected across 47 NEON sites between 2017 and 2024, spanning diverse biomes from tropical forests to Arctic tundra.

## Features

- **Automated data acquisition**: Download and harmonize NEON PLFA data (DP1.10104.001) and soil chemistry data (DP1.10086.001)
- **Quality control**: Apply standardized QC procedures and correct for known contaminants (C18:0)
- **Microbial community metrics**: Calculate total biomass, fungal:bacterial ratios, Gram-positive:Gram-negative ratios, cyclopropyl:precursor stress indices, and more
- **Interactive visualization**: Explore temporal trends, spatial patterns, and environmental drivers through a Shiny app
- **Reproducible workflows**: Fully documented code with transparent processing steps

## Data Access

The processed dataset (11,399 samples × 129 variables) is available on Zenodo:

**Hagen, C.J. & SanClements, M.D. (2026).** Continental-scale synthesis of NEON phospholipid fatty acid data (2017-2024). *Zenodo*. [https://doi.org/10.5281/zenodo.21463376](https://doi.org/10.5281/zenodo.21463376)

## Installation

### Requirements

- R version ≥ 4.0.0
- Required R packages:

```r
# Install required packages
install.packages(c(
  "neonUtilities",  # NEON data access
  "dplyr",          # Data manipulation
  "tidyr",          # Data tidying
  "ggplot2",        # Visualization
  "shiny",          # Interactive app
  "DT",             # Interactive tables
  "leaflet",        # Interactive maps
  "broom",          # Statistical summaries
  "stringr",        # String manipulation
  "purrr",          # Functional programming
  "lubridate",      # Date handling
  "jsonlite",       # JSON parsing
  "tibble"          # Modern data frames
))
```

### Download

```bash
# Clone this repository
git clone https://github.com/CedricHagen/neonPLFA.git
cd neonPLFA
```

Or download as ZIP from the [GitHub repository](https://github.com/CedricHagen/neonPLFA).

## Quick Start

### Option 1: Interactive Shiny Application

Launch the NEON PLFA Explorer app to interactively download, process, and visualize data:

```r
# Set working directory to the repository
setwd("path/to/neonPLFA")

# Launch the Shiny app
shiny::runApp("R")
```

The app allows you to:
- Select NEON sites and date ranges
- Download and process data automatically
- Visualize temporal trends and spatial patterns
- Explore relationships between microbial metrics and environmental drivers
- Export processed datasets and plots

### Option 2: Process Data with R Scripts

Use the processing functions directly:

```r
# Load functions
source("R/neon_plfa_functions.R")
source("R/neon_plfa.R")

# Download and process PLFA data for specific sites
# See R/run_neon_plfa_functions.R for examples
```

### Option 3: Use Pre-processed Data

Download the processed dataset directly from Zenodo and analyze in your preferred software:

```r
# Read the processed dataset
plfa_data <- read.csv("neon_plfa_synthesis_v1.2.csv")

# Explore the data
head(plfa_data)
summary(plfa_data)
```

## Data Products Used

This workflow integrates the following NEON data products:

- **DP1.10104.001**: Soil microbe marker genes, metagenomes and PLFA
- **DP1.10086.001**: Soil physical and chemical properties, periodic
- **DP1.00096.001**: Soil physical and chemical properties, Megapit (for soil classification)

## Calculated Metrics

The workflow calculates six key microbial community metrics:

1. **Total microbial biomass** (nmol/g dry soil): Total scaled lipid concentration, C18:0-corrected
2. **Fungal PLFA** (nmol/g): Fungal biomarker (18:2ω6,9)
3. **Bacterial PLFA** (nmol/g): Sum of 12 bacterial biomarkers
4. **Fungal:bacterial (F:B) ratio**: Relative dominance of fungal vs. bacterial decomposers
5. **Gram-positive:Gram-negative ratio**: Relative abundance of bacterial groups
6. **Cyclopropyl:precursor stress index**: Ratio of cyclopropyl fatty acids to their monoenoic precursors (primary: cy19:0/18:1ω7c), a physiological-stress indicator

See the [data dictionary](tables/TableS2_data_dictionary.csv) for complete variable descriptions.

## Repository Structure

The repository root is the public code + dataset release (Shiny app, processed
data, release tables/figures). Each manuscript built on this dataset lives under
`projects/`, and every project keeps its submission-ready deliverables in a
`manuscript/` subfolder.

```
neonPLFA/
├── README.md                          # This file
├── LICENSE                            # MIT License
├── CITATION.cff                       # Citation metadata
├── zenodo_metadata_dataset.json       # Dataset release metadata
├── zenodo_metadata_code.json          # Code release metadata
├── R/
│   ├── app.R                          # Shiny application
│   ├── neon_plfa.R                    # Core functions for data download
│   ├── neon_plfa_functions.R          # PLFA metric calculations
│   └── run_neon_plfa_functions.R      # Example usage
├── data/
│   └── neon_plfa_synthesis_v1.2.csv   # Processed dataset (129 variables)
├── tables/                            # Release tables (incl. TableS2_data_dictionary.csv)
├── figures/                           # Release figures (Figure1–6)
├── docs/
│   └── USER_GUIDE.md                  # Shiny app user guide
├── projects/                          # Manuscripts built on this dataset
│   ├── descriptor/                    # Scientific Data descriptor (this dataset)
│   │   └── manuscript/                #   docx + supplement + tables/ + figures/ + data/
│   │                                  #   + generate_manuscript_outputs_v1.0.R
│   ├── mane/                          # MANE continental test (New Phytologist)
│   │   ├── manuscript/                #   PLFA_MANE_MS.docx + supplement + figures/
│   │   └── 0X_*.R, output/            #   analysis pipeline
│   ├── climate/                       # Climate-regime × microbial (GCB)
│   │   ├── manuscript/                #   PLFA_climate_MS.docx + supplement + figures/
│   │   └── 0X_*.R, output/            #   analysis pipeline
│   └── flux/                          # Microbial → soil respiration (Project D; in progress)
│       ├── manuscript/                #   (placeholder — no draft yet)
│       └── 0X_*.R, output/            #   analysis pipeline
└── archive/                           # Superseded data & older drafts
```

## Documentation

- **Shiny App User Guide**: See `docs/USER_GUIDE.md` for detailed instructions
- **Data Dictionary**: `tables/TableS2_data_dictionary.csv` describes all 129 variables
- **Manuscript**: Full methods and validation are described in the associated Scientific Data publication

## Citation

If you use this dataset or code, please cite:

**Hagen, C.J. & SanClements, M.D. (2026).** A continental-scale synthesis of phospholipid fatty acid data from the National Ecological Observatory Network. *Scientific Data* (in review).

**Data:**  
Hagen, C.J. & SanClements, M.D. (2026). Continental-scale synthesis of NEON phospholipid fatty acid data (2017-2024). *Zenodo*. https://doi.org/10.5281/zenodo.20299460

**Code:**  
Hagen, C.J. & SanClements, M.D. (2026). neonPLFA: R workflow for processing NEON PLFA data (v1.0.0). *Zenodo*. https://doi.org/10.5281/zenodo.YYYYYYY

BibTeX:
```bibtex
@article{hagen2026neon,
  author = {Hagen, Cedric J. and SanClements, Michael D.},
  title = {A continental-scale synthesis of phospholipid fatty acid data from the National Ecological Observatory Network},
  journal = {Scientific Data},
  year = {2026},
  note = {in review}
}

@dataset{hagen2026data,
  author = {Hagen, Cedric J. and SanClements, Michael D.},
  title = {Continental-scale synthesis of NEON phospholipid fatty acid data (2017-2024)},
  year = {2026},
  publisher = {Zenodo},
  doi = {10.5281/zenodo.21463376}
}
```

## Contributing

We welcome contributions! Please:
- Report bugs or request features via [GitHub Issues](https://github.com/CedricHagen/neonPLFA/issues)
- Submit improvements via pull requests
- Follow existing code style and documentation practices

## Support

For questions or issues:
- Open an [issue](https://github.com/CedricHagen/neonPLFA/issues) on GitHub
- Contact: Cedric J. Hagen (cedric.hagen@colorado.edu)

## Acknowledgments

This work was supported by the National Ecological Observatory Network (NEON), a program sponsored by the U.S. National Science Foundation and operated under cooperative agreement by Battelle. This work was supported by the NSF-funded CO-WY ASCEND Engine. We thank H. Cross and S. Weintraub-Leff for discussions and feedback.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

The processed dataset is available under a CC-BY-4.0 license.

## Version History

- **v1.0.0** (2026-06): Initial release
  - 11,399 samples from 47 NEON sites (2017-2024)
  - Complete workflow and Shiny application
  - Manuscript submission to Scientific Data
