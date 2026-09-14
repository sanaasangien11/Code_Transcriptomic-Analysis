Transcriptomic Analysis Pipeline

This repository contains code for transcriptomic data analysis, covering data processing, exploratory data analysis, sequence extraction, and downstream analysis.

Overview

The project provides a partially reproducible workflow for processing and exploring transcriptomic data and extracting biological sequences to identify the candidate genes. 

Main Components

* Data Processing – quality control, trimming, alignment and organization of transcriptomic datasets
* Exploratory Data Analysis – statistical summaries, data exploration, and visualization
* Sequence Extraction – extraction and processing of relevant nucleotide/transcript sequences
* Parsing Annotations - Parsing the gff/gtf annotations for relevant transcript information
* Transcriptomic Analysis – analysis of transcript-level data and biological patterns 
* Data Visualization – plots and visual summaries of the processed data

Repository Structure

.
├── Scripts/               # Analysis and processing scripts
├── Sequence_Extraction_Code/   # Sequence extraction code
├── Parsing_Annotation_Code/   # Code for extracting relevant transcript-level information
└── README.md              # Project documentation

Workflow

The general workflow is:

1. Data Input
2. Data Processing & Quality Control
3. Exploratory Data Analysis
4. Differential Gene Expression Analysis
5. Visualization and Results
6. Sequence Extraction
7. Correlation Analysis

Requirements

The analysis is implemented using R and python.

Main packages and tools used in this project may include:

* Python
* Pandas
* NumPy
* Biopython
* R / Bioconductor 

Reproducibility

The code is organized to make the transcriptomic analysis workflow reproducible. Some of the parts of the code are reproducible. 

Author

Sanaa Sangien


This project is intended for research purposes.
