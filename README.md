Differential Expression Analysis Using DESeq2

This script is designed to perform differential expression analysis using the R package DESeq2. It processes count data from RNA sequencing experiments, compares treated and control conditions, and outputs the results as a CSV file.

Prerequisites

Before running this script, ensure that the following are installed and available in your R environment:
1. R version 4.0 or later
2. Required R packages:
`data.table`
`dplyr`
`DESeq2`

Script Workflow

1.	Setup
	•	The script begins by clearing the R environment (`rm(list=ls())`) and loading the required libraries (`data.table`, `dplyr`, and `DESeq2`).
	•	The working directory is set to the folder containing the count matrices.
2.	Data Loading
	•	Two sets of count matrices are loaded:
	•	Control condition: `counts` (processed) and `counts2` (unprocessed).
	•	Treated condition: `counts3` (processed) and `counts4` (unprocessed).
3.	Data Integration
	•	The processed and unprocessed counts for both conditions are combined into two metadata tables:
	•	`countsMeta`: Control samples.
	•	`countsMeta2`: Treated samples.
	•	These two tables are merged into a single table (`countsMeta3`) using a full join on the “SuperPos” column, ensuring all positions are included.
	•	Missing values are replaced with zeros, and the “SuperPos” column is removed after setting it as row names.
4.	Sample Metadata
	•	A sample metadata table (`sampleTable`) is created, specifying:
	•	Sample names.
	•	Experimental conditions (`ctr` for control, `trt` for treated).
5.	Differential Expression Analysis
	•	A DESeq2 dataset (`dds`) is created using the combined count matrix (`countsMeta3`) and sample metadata (`sampleTable`).
	•	DESeq2 is run to perform differential expression analysis with a design formula based on the experimental condition (`~ condition`).
6.	Results Export
	•	The results of the differential expression analysis are saved as a CSV file at the specified output path (`outputfile`). The file includes:
	•	Gene/feature IDs as row names.
	•	Log fold changes, p-values, adjusted p-values, and other statistics.

How to Run
1.	Place all required input files in the specified working directory.
2.	Update the working directory path in the script if necessary:
	setwd("~/path/to/your/data")
3.	Update the output file path if needed:
	outputfile <- "~/path/to/output/results.csv"
4.	Run the script in an R environment.

