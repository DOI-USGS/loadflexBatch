#### Overview of loadflexBatch ####

# loadflexBatch is an R script that runs several loadflex models for many sites 
# and solutes in one big batch. There are three steps to using it.

# 1. Prepare the input data files. Marcelo has done this for many ANA sites.
# Examples are in three_ANA_sites/input, which has subfolders for NO3, PT, and
# Q. Each subfolder has data files inside.

# 2. Prepare the input metadata files. There are two: three_ANA_sites.yml and 
# three_ANA_sites/siteInfo.csv. Marcelo has done this for ANA, too, and this 
# RStudio Project has examples.

# 3. Source the file named batch.R. This script loads your input data files,
# runs models, and produces output files.

# See https://www.github.com/USGS-R/loadflexBatch/blob/master/blog.md for
# details on how to prepare data and run loadflexBatch.


#### Quick start for loadflexBatch ####

# Today we will run loadflexBatch for a single site (RONC02800) and solute
# (NO3).

# 1. Browse through the files in three_ANA_sites/input. Open the folder in
# Windows Explorer so that you can open the .csv files in Excel for easy
# viewing.

# 2a. Edit siteInfo.csv so that it only includes the first line (the column 
# names) and the three lines for RONC01200.

# 2b. Edit the 'constituents' line of three_ANA_sites.yml so that only NO3 will
# be modeled. That line should look like this:
constituents: ["NO3"]

# 3. Run loadflexBatch: Open batch.R. First check that line 9 says:
inputs <- yaml::yaml.load_file('three_ANA_sites.yml')
# Then click Source in the top right of the Script Editor pane.

# Wait for the script to finish, then inspect the contents of
# three_ANA_sites/output.

# Later we will look at the output from Marcelo's batch runs, which include more
# sites and solutes.
