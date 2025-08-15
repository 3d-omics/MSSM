# 3D'omics | Micro-scale spatial metagenomics (MSSM)
This webbook contains all the code used for the analyses of in "Micro-scale spatial metagenomics: revealing high-resolution spatial biogeography of gut microbiomes".
<br> 

## Input data
The main datasets used for the microbiome analyses, and included in the data folder, are: 
<br>


## Analysis steps

#### Step 1 - Data Preparation
The code in 01-data_preparation.Rmd loads the input data files, produces a filtered and normalized genome count table for further quantitative analyses and groups all the required dataframes and the customised colour palettes in an R object (xxx.Rdata) for downstream analyses.

#### Step 2 – Manuscript Method Development
The R code in 02-method_development.Rmd reproduces all analyses described in the Method Development section of the manuscript. This includes the comparison of macroscale and microscale MAG catalogues, as well as analyses focused on optimizing lysis conditions, evaluating library performance across intestinal sections and positive control reactions, and exploring resource optimization and throughput strategies to improve efficiency and sample processing capacity.

#### Step 3 – Manuscript Method Validation
The R code in 03-method_validation.Rmd executes all analyses for the Method Validation section. It applies the developed MSSM method to two animals, assessing the discriminative power and replicability of MSSM data compared to FISH and macroscale metagenomics.

#### Step 4 – Manuscript Method Implementation
The R code in 04-method_implementation.Rmd runs all analyses described in the Method Implementation section of the manuscript. This includes applying the validated MSSM method to caecum and colon cryosections, enabling the study of spatial patterns in taxonomic and functional variation. It also investigates the spatial and host distribution of Lawsonibacter strains resolved by MSSM and demonstrates the method’s ability to recover SNP-level microdiversity within strains.

## Analysis output

The **bookdown-rendered webbook** containing all the above code and its output is available at:

[.....](...)

While the webbook provides a user-friendly overview of the procedures, analyses can be directly reproduced using the **.Rmd** files stored in the root directory of this repository. Note that some code chunks that require heavy computation might have been turned off using 'eval=FALSE' or cached using 'cache=TRUE'. To re-render the webbook, you can use the following code:

```r
library(bookdown)
library(htmlwidgets)
library(webshot)

render_book(input = ".", output_format = "bookdown::gitbook", output_dir = "docs")
```


