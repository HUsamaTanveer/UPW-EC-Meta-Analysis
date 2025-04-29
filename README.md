# UPW-EC-Meta-Analysis

This repository contains the data and code for the paper "Evaluating the Impact of Emerging Contaminants on Membrane Performance in Ultrapure-Water Production for Semiconductor Manufacturing: A PRISMA-Directed Meta-Analysis"

## Repository Structure

UPW-EC-Meta-Analysis/
│
├── README.md                 # Repository overview
├── LICENSE                   # License information
├── requirements.txt          # Python dependencies
│
├── data/                     # Datasets
│   ├── raw/                  # Raw extracted data
│   │   ├── study_metadata.csv
│   │   ├── rejection_data.csv
│   │   ├── flux_decline_data.csv
│   │   └── operating_conditions.csv
│   │
│   └── excel/                # Excel datasets
│       ├── Meta-Analysis_Complete_Dataset.xlsx
│       ├── Quality_Assessment_Results.xlsx
│       └── Figure_Data.xlsx
│
└── code/                     # Analysis code
    ├── 1_data_preprocessing/ # Scripts for data cleaning and preparation
    │   ├── 01_extract_data.py
    │   ├── 02_clean_data.py
    │   └── 03_calculate_effect_sizes.py
    │
    ├── 2_meta_analysis/      # Meta-analysis scripts
    │   ├── 01_pooled_estimates.py
    │   ├── 02_heterogeneity.py
    │   ├── 03_meta_regression.py
    │   └── 04_publication_bias.py
    │
    └── 5_visualization/      # Visualization scripts
        ├── 01_forest_plots.py
        ├── 02_box_plots.py
        ├── 03_bubble_plots.py
        ├── 04_violin_plots.py
        └── 05_impact_ranking.py

## Setup

1. Install required packages: pip install -r requirements.txt

2. The data analysis workflow follows these steps:
   - Data preprocessing (code/1_data_preprocessing/)
   - Meta-analysis and meta-regression (code/2_meta_analysis/)
   - Visualization (code/5_visualization/)


## Note

This repository contains the code and data associated with a manuscript that is currently under review. Additional documentation will be provided upon publication.
