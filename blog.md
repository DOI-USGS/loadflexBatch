Loadflex batch script
---------------------

The script at <https://github.com/USGS-R/loadflexBatch/blob/master/batch.R> automates running multiple load models for multiple sites and consitutents, using the `loadflex` package along with some features of `rloadest`. Output is collated across all sites for easy analysis of input data, predicted loads, and model metrics. This post will go through basic setup and use of the script as it exists now. In the future it may be hardened into a part of the actual `loadflex` package.

Installation and setup
----------------------

First, go to the Github repository [USGS-R/loadflexBatch](https://github.com/USGS-R/loadflexBatch) and download the zip file (using the big green button) to your preferred directory for R projects, and unzip it. Open RStudio, start a new project (File -&gt; New Project), select the "Existing directory" option. select the `loadflexBatch-master` folder, and click "Create Project". Note that on Windows there will be two `loadflexBatch-master` folders — select the lower-level one. You will now be inside the `loadflexBatch-master` folder, and have access to the batch script.

``` r
list.files()
```

    ##  [1] "batch.R"                "batchHelperFunctions.R"
    ##  [3] "blog.md"                "blog.pdf"              
    ##  [5] "blog.Rmd"               "Hirsch_sites"          
    ##  [7] "Hirsch_sites.yml"       "loadflexBatch.Rproj"   
    ##  [9] "README.md"              "three_ANA_sites"       
    ## [11] "three_ANA_sites.yml"

Next, we need to install the packages that the script depends on. In your console, run

``` r
install.packages(
  c('dplyr', 'rloadest', 'devtools', 'yaml'), 
  repos = c('https://owi.usgs.gov/R', 'https://cloud.r-project.org'))
```

We will also install the main `loadflex` package directly from Github to ensure we have the very latest version:

``` r
devtools::install_github("USGS-R/loadflex")
```

The most up-to-date installation instructions can always be found at <https://github.com/USGS-R/loadflex#installation>.

Now we are ready to look at the user inputs, run the script, and inspect the output.

Input parameters and files
--------------------------

Open the main script, `batch.R`, by clicking on it in the 'Files' pane in the lower right of your RStudio window. There is a basic description of the file at the top. Below that is a command to read in the user inputs file for one of the included example datasets. It looks like this:

``` r
inputs <- yaml::yaml.load_file('three_ANA_sites.yml')
```

The default file is `three_ANA_sites.yml`. To run another example (e.g., `Hirsch_sites.yml`), edit the file name in the above line of the `batch.R` script. The user inputs file is in the YAML language (<http://www.yaml.org/>). In YAML, each line contains a `key: value` pair, except blank lines and those lines beginning with `#`, which are comments. In your user input YAML file you should supply information about input/output folder names and locations, constituents, load units, and load rate units. `three_ANA_sites.yml` looks like this:

    # This YAML file contains high-level user options for running batch.R on a
    # specific dataset. YAML files are made of key: value pairs like the ones below.
    # Edit the values (the text to the right of each colon) to describe your data.

    # input file information
    inputFolder: "three_ANA_sites/input" # folder containing all input files and subfolders
    constituents: ["NO3", "PT"] # names of folders inside inputFolder containing constituent data, and the column names for constituents within those data files
    discharge: "Q" # name of the folder inside inputFolder containing daily discharge data, and the column name for discharge within those data files
    date: "date" # column name for dates within the constituent and discharge data files
    siteInfo: "siteInfo.csv" # name of the csv file inside inputFolder containing site and constituent metadata

    # analysis specifications
    minDaysPerYear: 345 # number of days required for a year to be included in the multi-year average

    # desired units for output
    loadUnits: "kg"
    loadRateUnits: "kg/yr"

    # output folder where results files and subfolders will be written
    outputFolder: "three_ANA_sites/output"

The `constituents` named in the user inputs file need to match the names of the input subfolders that contain the input data (paired water quality and discharge measurements). Similarly, `dischargeFolder` should name a subfolder containing the discharge measurements used to make the load predictions.

The file named by `siteInfo` should be a comma-separated (.csv) file inside `inputFolder`. This file contains a table of metadata for individual water quality and discharge sites. Look at the example file (`three_ANA_sites/input/siteInfo.csv`) for reference on the required column names and format:

    ##   matching.site   site.id site.name    lat    lon basin.area constituent
    ## 1     RONC02800 RONC02800      Ronc -20.25 -48.43        100         NO3
    ## 2     RONC02800 RONC02800      Ronc -20.25 -48.43        100          PT
    ## 3     RONC02800 RONC02800      Ronc -20.25 -48.43        100           Q
    ## 4     MOGU02900 MOGU02900      Mogu -20.15 -48.63        200         NO3
    ## 5     MOGU02900 MOGU02900      Mogu -20.15 -48.63        200          PT
    ## 6     MOGU02900 MOGU02900      Mogu -20.15 -48.63        200           Q
    ## 7     ORIZ02900 ORIZ02900      Oriz -21.41 -47.87        300         NO3
    ## 8     ORIZ02900 ORIZ02900      Oriz -21.41 -47.87        300          PT
    ## 9     ORIZ02900 ORIZ02900      Oriz -21.41 -47.87        350           Q
    ##     units
    ## 1 mg L^-1
    ## 2 mg L^-1
    ## 3     cms
    ## 4 mg L^-1
    ## 5 mg L^-1
    ## 6     cms
    ## 7 mg L^-1
    ## 8 mg L^-1
    ## 9     cms

For each row of the `siteInfo` file, the value in the `constituent` column should match both (1) a folder name given as a `constituent` or `dischargeFolder` in the user inputs file, and (2) an actual folder within the `inputFolder`. And then within each of those folders, there should be a separate data file (.csv format) for each site, with a file name equal to the names given in the `site.id` column of the `siteInfo` file.

For this example, the input folder structure is therefore:

    - three_ANA_sites
      - input
        - siteInfo.csv
        - NO3
          - RONC02800.csv
          - MOGU02900.csv
          - ORIZ02800.csv
        - PT
          - RONC02800.csv
          - MOGU02900.csv
          - ORIZ02800.csv
        - Q
          - RONC02800.csv
          - MOGU02900.csv
          - ORIZ02800.csv

The batch script loops over the constituents listed in the user inputs file, finds the corresponding rows in `siteInfo.csv`, and identifies the one water quality file and one discharge file for each site-constituent combination. Both files bear the name of a `site.id` and have the `.csv` suffix, and they appear in a constituent folder (e.g., `NO3`) and the discharge folder (`Q`), respectively.

Water quality and discharge are usually measured at the exact same location, but sometimes water quality is sampled somewhat downstream or upstream of the flow gage. The script handles this possible discrepancy by expecting two columns for site identifiers in the `siteInfo` file. The column called `site.id` describes precise site locations. If water quality is measured downstream or upstream of the flow gage, the concentration and discharge rows in the `siteInfo` file should have different `site.id` values. The column called `matching.site` is used to match up pairs of water quality and discharge sites that you want to combine into a load model. The `matching.site` should always be the same for the concentration row and discharge row to be combined.

Here is an example of how the input consituent data should be formatted, with columns for date of the observation (`date`), discharge (`Q`), the concentration of the constituent (in this case `NO3`), additional columns that are ignored by this script (e.g., `CODIGO_ESTACAO`), and the censoring and data quality code (`status`), which follows the Brazillian ANA's convention of 0=bad value, 1=normal value, 2=value known to be less than or equal to the number given in the constituent (`NO3`) column.

``` r
head(read.csv('three_ANA_sites/input/NO3/MOGU02900.csv'), 5)
```

    ##         date         Q  NO3 CODIGO_ESTACAO status
    ## 1 2001-02-06 285.00034 0.27      MOGU02900      1
    ## 2 2001-04-03 305.36021 0.28      MOGU02900      1
    ## 3 2001-06-20 115.44747 0.47      MOGU02900      1
    ## 4 2001-08-07  82.08674 0.58      MOGU02900      1
    ## 5 2001-10-02 104.30206 0.58      MOGU02900      1

And here is how the input discharge data should be formatted, with columns for date of the observation (`date`) and mean daily discharge (`Q`), plus optional additional columns that will be ignored by `batch.R`. Whereas the constituent data file only has rows for those dates on which concentration was measured, the discharge data file has rows for every date on which flux is to be estimated.

``` r
head(read.csv('three_ANA_sites/input/Q/MOGU02900.csv'), 5)
```

    ##         date        Q CODIGO_ESTACAO
    ## 1 2001-01-01 503.0083      MOGU02900
    ## 2 2001-01-02 478.1621      MOGU02900
    ## 3 2001-01-03 441.5538      MOGU02900
    ## 4 2001-01-04 392.3504      MOGU02900
    ## 5 2001-01-05 344.9383      MOGU02900

Running the script
------------------

Once the input parameters and files are set correctly, you can source the script. Recall that it is configured to run with the `three_ANA_sites` example by default, but you can edit the .yml filename to pick a different sample.

``` r
source('batch.R')
```

Behind the scenes
-----------------

The script reads, processes, and writes output for each site/consituent combination individually. Currently, it runs three models for each iteration — the `rloadest` 5-parameter regression model, an interpolation model from `loadflex`, and a composite interpolation-regression model also from `loadflex`. Additional models can be added with some modifications to the script.

If the drainage basin areas for paired water quality and discharge sites are different (in the site metadata file, e.g `siteInfo.csv`), discharge is scaled by the appropriate ratio. The resulting estimates describe fluxes at the water quality monitoring site.

If a constituent dataset includes censored data (i.e., `status` flags of 2), those data are passed to the `loadReg2` model in the censored data format required by `rloadest`. The composite and interpolation models are then excluded, because these models are not equipped to handle censored data.

All the site metadata is stored in a `loadflex::metadata` R object, which is referenced throughout the script.

The script uses several built-in `loadflex` functions to create the outputs in R, which the script then writes to files. Those functions include `summarizeInputs`, `summarizeModel`, `predictSolute`, and `aggregateSolute`. Descriptions of these functions can be found by typing `?` followed by the function name at the R prompt.

Output files
------------

Like the input files, the output files are written to a separate folder for each constituent. Most important are the five summary files, each of which covers all the sites for that consituent. These are prefixed with the constituent name, so for NO3, these are:

    - output
      - NO3
        - NO3_inputs.csv       # input data summary
        - NO3_annual.csv       # annual predicted loads from each model
        - NO3_multiYear.csv    # multiyear predicted loads from each model
        - NO3_modelMetrics.csv # diagnostics and descriptive metrics for each model
        - NO3_plots.pdf        # plots of input data, predicted loads, and model diagnostics

Additionally, there are four folders containing individual .csv files for each site. The four .csv files listed above combine the contents of those folders.

Conventions to remember
-----------------------

-   Column names for site information (e.g., in the inputs summary) are assumed to refer to the water quality site, unless the name explicitly indicates that they refer to the flow site.

-   IDs for constituents, discharge, and sites need to be consistent across the user inputs file, the `siteInfo` file, column names within the files, and the folder and file names in the `siteInputs` folder.

What's next?
------------

Create your own input files following the conventions above. Edit line 9 of `batch.R` to point to your own user inputs YAML file, e.g.,

``` r
inputs <- yaml::yaml.load_file('my_own_sites.yml')
```

Then run `source('batch.R')` to produce outputs for your sites.

Inspect the output files to learn about the size and quality of the input data, how well each model performed, and how the model predictions varied across models and sites. We'll create a separate document with suggestions for how to use these outputs to assess the quality of the data, models, and predictions.
