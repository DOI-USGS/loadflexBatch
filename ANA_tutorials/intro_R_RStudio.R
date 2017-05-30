#### Exploring RStudio ####

# RStudio is an Integrated Development Environment (IDE). It helps you write R 
# code, run R code, and see the results of the code you run. There are 10 
# RStudio panes that help with these tasks: Script Editor, Console, Environment,
# History, Files, Plots, Packages, Help, and Viewer. They are laid out so you 
# can see many related panes at once.

# 1. Find each of the 10 panes.

# 2. Click to put your cursor in the Console.

# 3. Type your first line of code:
3 * 4
# and press enter to see the result. You should see
[1] 12

# 4. Open a new file in the Script Editor with Ctrl+Shift+N. Type
3 * 4
# in that file and press Enter. Nothing happens. Now place your cursor on the 
# line of code and press Ctrl+Enter. Now the line gets copied to the R console and
# gets run.

# 5. Type a '#' followed by a comment in the Script Editor and press enter.
# Nothing happens, which is good. A comment can be anything you want, as long as
# the line starts with #. See how this whole file is written with comments? They
# don't act like code but let you explain what your code does.

# 6. Can you see everything clearly? You can change the font size, or the text 
# color scheme ('syntax highlighting') using the Tools | Global Options | 
# Appearance menu.

#### Managing files ####

# RStudio includes a file manager. R code files are just text files with a .R 
# extension.

# 1. Bring the Files pane to the front (Ctrl+5) to navigate to files and folders
# on your computer.

# 2. Navigate to a place where you want to create a new folder (I will use the
# Desktop). Click 'New Folder' to create a new folder. Name it 'introR'.

# 3. Save your open R script file in the new introR folder with File | Save or 
# Ctrl+S. Name it 'explore.R'. You should now see the explore.R file in the
# introR folder in your Files pane.

# 4. In the Files pane, click the More drop-down menu and select 'Show Folder in
# New Window' to open to a Windows Explorer window in the same file location.


#### Working directories and RStudio Projects ####

# Like a Windows Explorer window, an R session is always located somewhere in 
# your file system. You can refer to files using file paths relative to that
# location.

# 1. Where is your R session now? Run the getwd function ('get working
# directory') to find out:
getwd()

# 2. List the files in the current working directory:
dir()

# 3. The directory of the Files pane and the directory of your R session are 
# often different. Change your R session working directory so it matches your 
# file directory using the setwd command ('set working directory'). Note that
# you need to use single forward slashes.
setwd('C:/Users/aappling/Desktop/introR')

# 4. Now call
dir()
# again. You should see
[1] "explore.R"

# RStudio uses Projects to do step 3 for you automatically. A Project is a short
# text file with extension '.Rproj' that RStudio creates. When you open that
# file as a Project, RStudio sets both (1) the Files pane and (2) the R session
# working directory to equal the parent directory of the .Rproj file.

# 6. Turn introR into a project. From the menu bar: File | New Project... | 
# Existing Directory | Create Project. The R session will restart. 

# 7. Note that your Files pane is already set to the introR folder. Your working
# directory is set to the same folder. Confirm this by running
getwd()

# 8. Now call
dir()
# again. You should see
[1] "explore.R"    "introR.Rproj"

# 9. Reopen explore.R by clicking on the file name in the Files pane.


#### Running code interactively ####

# When creating an analysis or exploring data, you can use R interactively. This
# means running one command at a time to modify the data in R's memory, pausing 
# to inspect the results and decide what to do next.

# 1. Find the Script editor (Ctrl+1) and Console (Ctrl+2). Also bring the
# Environment pane (Ctrl+8) to the front so you can see it.

# 2. Write your commands in the script editor, where they can be saved. Send the
# commands to the console to run them: place your cursor on the line you want to
# run, then press Ctrl+Enter.

# 3. See a list of what variables are in R's memory in the Environment pane. 
# Explore what is in your Environment by looking at the second column (for 
# simple Values), clicking on the variable name (for Data), or clicking on the 
# blue arrow (for Data or List Values)

#### Objects in R ####

# R stores data in the R session memory. Code mostly works by changing what's in
# memory. Objects in memory are referred to by unquoted names (e.g., mynum 
# below). Data is assigned to objects with the <- (assignment) operator.

# numbers
mynum <- 4.412
mynums <- c(4.412, 5.726)

# text (character strings)
mytext <- "MOGU029000"

# dates
mydate <- as.Date('2017-06-01')
mydatetime <- as.POSIXct('2017-06-01 12:14:01')

# data.frames
mydf <- data.frame(
  date=as.Date(c('2017-06-02','2017-06-03')), 
  discharge=c(9,12)
)

# lists
mylist <- list(
  site=mytext,
  date=mydate,
  mydf=mydf
)

# extracting parts of data.frames and lists
mydf$date
mylist$site

# finding out which parts of data.frames and lists can be extracted
names(mylist)

# assigning parts of data.frames and lists to new variables
mydf_datecolumn <- mydf$date


#### Running an entire script at once ####

# R scripts can also behave more like stand-alone software packages.

# 1. Create a second R script. Save it in the same directory. Call it
# 'script.R'. Copy and paste these commands into the script:
x <- c(1, 2, 3)
y <- c(5, 3, 9)
xydf <- data.frame(x, y)
print(xydf)
# then save the script again.

# 2. Switch back to explore.R. Add this command:
source('script.R')
# then place your cursor on that line of code and use Ctrl+Enter to run the code
# in the Console. You should see the results of script.R appear in the console:
  x y
1 1 5
2 2 3
3 3 9

# SPARROW-R and loadflexBatch both use source() to run more complicated scripts.


#### Running functions from a package

# R packages are collections of R functions. You have already used functions in 
# the base R packages. For example, print() is a function that displays an
# object in the console, and data.frame() is a function that creates a
# data.frame object in the R session's memory. Other functions are only
# available if you load a new package.

# 1. We will be using the loadflex package to estimate river solute fluxes. Load
# the loadflex package with the library() function:
library(loadflex)
# You'll see a lot of package startup messages. These are fine unless you see
# 'Error'.

# 2. The loadflex package has many functions and an example dataset. Load the
# example dataset into memory with the data() function.
data(lamprey_nitrate)
# You can view the contents of lamprey_nitrate by clicking on its name in the
# Environment pane.

# 3. To view a help file for the loadflex package, go to the Help pane (Ctrl+3).
# Type loadflex-package in the search box in the top right and press Enter.

# We will use other functions from loadflex in the next few hours.

install.packages(c('EGRET','loadflex'), repos=union('https://owi.usgs.gov/R', options()$repos))
