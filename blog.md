---
author: David Watkins
date: YYYY-MM-DD
slug: loadflex-batch
draft: True
title: Loadflex Batch Mode
type: post
categories: Data Science
 
 
 
 
 
 
 

tags: 
  - R
 
 
description: 
keywords:
  - R
 
 
 
 
---
### Loadflex Batch Script

This script automates running various load models for different sites and consitutents, using the `loadflex` package along with some features of `rloadest`. Output is collated across all sites for easy analysis of input data, predicted loads, and model metrics. This post will go through basic setup and use of the script as it exists now. In the future it may be hardened into a part of the actual `loadflex` package.

Installation and setup
----------------------

First, go to the Github repository [USGS-R/loadflexBatch](https://github.com/USGS-R/loadflexBatch) and download the zip file (using the big green button) to your preferred directory for R projects, and unzip it. Open RStudio, start a new project (File -&gt; New Project), select the "Existing directory" option, select the `loadflexBatch-master` folder, and click "Create Project". You will now be inside the `loadflexBatch-master` folder, and have access to the batch script.

``` r
list.files()
```

    ##  [1] "batch.R"                "batchHelperFunctions.R"
    ##  [3] "blog.html"              "blog.md"               
    ##  [5] "blog.pdf"               "blog.Rmd"              
    ##  [7] "loadflexBatch.Rproj"    "NO3.csv"               
    ##  [9] "output"                 "PT.csv"                
    ## [11] "README.html"            "README.md"             
    ## [13] "README.Rmd"             "three_ANA_sites"

Next, we need to install the packages that the script depends. In your console, run

``` r
install.packages('dplyr', 'rloadest', 'devtools', repos = c('https://owi.usgs.gov/R', 'https://cloud.r-project.org'))
```

We will also install the main `loadflex` package, directly from Github to ensure we have the very latest version:

``` r
devtools::install_github("USGS-R/loadflex")
```

The most up-to-date installation instructions can be found at <https://github.com/USGS-R/loadflex#installation>.

Now we are ready to look at the user inputs, file structure, and run the script.

Input parameters and directory setup
------------------------------------

Open the main script, `batch.R`. There are some basic instructions at the top. Below that are the user inputs, set up for the included example data. The user supplies information about input/output folder names and locations, constituents, load units, and load rate units. The constituent names need to match the names of the input subfolders that contain the input data (paired water quality and discharge measurements). The discharge folder, containing the discharge measurements used to make the load predictions, works the same way. `siteInfo` is a .csv file (inside `inputFolder`) that contains metadata for water quality and discharge sites. Look at the example file (`three_ANA_sites/siteInfo.csv`) for reference:

    ##   matching.site   site.id site.name lat lon basin.area constituent   units
    ## 1     RONC02800 RONC02800      Ronc   1   1        100         NO3 mg L^-1
    ## 2     RONC02800 RONC02800      Ronc   1   1        100          PT mg L^-1
    ## 3     RONC02800 RONC02800      Ronc   1   1        100           Q     cms
    ## 4     MOGU02900 MOGU02900      Mogu   2   2        200         NO3 mg L^-1
    ## 5     MOGU02900 MOGU02900      Mogu   2   2        200          PT mg L^-1
    ## 6     MOGU02900 MOGU02900      Mogu   2   2        250           Q     cms
    ## 7     ORIZ02900 ORIZ02900      Oriz   3   3        300         NO3 mg L^-1
    ## 8     ORIZ02900 ORIZ02900      Oriz   3   3        300          PT mg L^-1
    ## 9     ORIZ02900 ORIZ02900      Oriz   3   3        300           Q     cms

Here is input consituent data formatting:

``` r
head(read.csv('three_ANA_sites/NO3/MOGU02900.csv'), 5)
```

    ##         date         Q  NO3 CODIGO_ESTACAO status
    ## 1 2001-02-06 285.00034 0.27      MOGU02900      1
    ## 2 2001-04-03 305.36021 0.28      MOGU02900      1
    ## 3 2001-06-20 115.44747 0.47      MOGU02900      2
    ## 4 2001-08-07  82.08674 0.58      MOGU02900      1
    ## 5 2001-10-02 104.30206 0.58      MOGU02900      2

Once the input parameters are set correctly, you can source the script. It is configured to run with the included example data to start.

``` r
source('batch.R')
```

Behind the scenes
-----------------

The script reads, processes, and writes output for each site/consituent combination individually. Currently, it runs three models for each iteration — the `rloadest` 5-parameter regression model, an interpolation model from `loadflex`, and a composite interpolation-regression model also from `loadflex`. Additional models can be added with some modifications to the script.

All the site metadata is stored in `loadflex::metadata` objects, where it is referenced throughout the script.

If the drainage basin areas for paired water quality and discharge sites are different (in the site metadata file, e.g `siteInfo.csv`), discharge is scaled by the appropriate ratio. The resulting estimates descibe fluxes at the water quality monitoring site.

Output
------

Like the input files, the output files are written to a separate folder for each constituent. Output contains 5 files, which cover all the sites for that consituent:

-   four .csvs:
    -   input data summary
    -   annual predicted loads for each model
    -   multiyear predicted loads for each model
    -   model metrics
-   PDF of plots of input data, predicted loads, and model diagnostics

Additionally, there are four folders corresponding to the four csvs that contain individual .csv files for each site.

Conventions to remember
-----------------------

-   Column names and metadata slots for site information are assumed to refer to the water quality site, unless the name specifies they refer to discharge.

-   Constituent and discharge symbols need to be consistent throughout, including directory names and in the site metadata csv.
