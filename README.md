# 3D'omics | Micro-scale spatial metagenomics (MSSM)
This webbook provides all the code used in the analyses for the study "Micro-scale spatial metagenomics: revealing high-resolution spatial biogeography of gut microbiomes". The project 3D'omics (http://www.3domics.eu) was Funded by European Union’s Horizon 2020 Research and Innovation programme. 
<br> 

## Input data
The data used for the analyses are can be downloaded from zenodo at: 
https://zenodo.org/records/17091770,
https://zenodo.org/records/17091731, and
https://zenodo.org/records/17091749.
The repository also contains an R object containing the necessary data to run Step 2-Step 4.
<br>


## Analysis

#### 1 - Data Preparation
The code in 01-data_preparation.Rmd loads and tidies the data, generates a filtered and normalized genome count table for downstream quantitative analyses, and saves all necessary data frames along with the customized color palettes in an R object (working_data_object.Rdata).

#### 2 – Main manuscript
The R code in 02-main_manuscript.Rmd reproduces all analyses described in the manuscript. This includes the technical validation of microscale microbiome profiling across intestinal regions; comparison of de novo and reference-based genome catalogues for bacterial quantification in caecal MSSM datasets; assessment of environmental and host contamination; evaluation of inter-individual variation and technical reproducibility; validation of MSSM against FISH and conventional metagenomics; analysis of within-strain SNP-level microdiversity across individuals; and assessment of non-random spatial structuring of microbial communities.
#### 3 – Supporting information
The R code in 03-supporting information.Rmd executes all analyses described in the Supporting Information Text. This includes the validation of lysis condition (S1: Lysis Conditions) and resource optimisation (S2: Resource optimisation).
## Analysis output

## Analysis output

The **bookdown-rendered webbook** containing all the above code and its output is available at:

https://3d-omics.github.io/MSSM

While the webbook provides a user-friendly overview of the procedures, analyses can be directly reproduced using the **.Rmd** files stored in the root directory of this repository. Note that some code chunks that require heavy computation might have been turned off using 'eval=FALSE' or cached using 'cache=TRUE'. To re-render the webbook, you can use the following code:

```r
library(bookdown)
library(htmlwidgets)
library(webshot)

render_book(input = ".", output_format = "bookdown::gitbook", output_dir = "docs")
```


