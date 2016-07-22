ARQiv
====

ARQiv is a R package for Preparing for and Completing Whole-Organism Screening at High-througput rates.

## Overview

The ARQiv R* package includes functions that fall into two categories - those applied to 'Pre-screening Assay Optimization' 
and 'Compound Analysis'. The functions allow the user to calculate background signal, determine sample size, run quality control tests, perform virtual experiments to simulate compound efficacy - and finally, to perform compound analysis
during iterative drug screen cycles.

The ARQiv package is a part of the published paper[paper name](website).

## Package Installation
* STEP 1:  Install the latest version of R(for Windows (https://cran.r-project.org/bin/windows/base/) or Mac(https://cran.r-project.org/bin/macosx/)) and RStudio(https://www.rstudio.com/products/rstudio/download2/) on your computer.
* STEP 2: To install ARQiv package via Github, the user must have installed [devtools](https://cran.r-project.org/web/packages/devtools/index.html) by running the following commands in RStudio Console window:
```{r}
install.packages("devtools")
library(devtools)
```
* STEP 3: To use the graphical user interface (GUI) in ARQiv package, the user must have first installed RGtk2, GTK2+  with commands in RStudio Console.

##### STEP 3a: for Windows
```{r}
install.packages("RGtk2",depen=T)
library(RGtk2)
```
There will be a notice for missing GTK. Choose "Install GTK+" when prompted - this will take ~1 minute to complete. Afterwards the user needs to restart R/Rstudio.

##### STEP 3b: for Mac
Install [GTK 2.24](http://r.research.att.com/libs/GTK_2.24.17-X11.pkg). And use command
```{r}
install.packages("RGtk2",depen=T)
```
in RStudio console window. Then user need to restart R/Rstudio.

* STEP 4: Complete RGtk2,gWidgets2RGtk2 packages installation by commands in RStudio console window:
```{r}
library(RGtk2)
install.packages("gWidgets2RGtk2")
library(gWidgets2RGtk2)
```

* STEP 5: Install ARQiv package and open GUI window with the following command in RStudio console window:
```{r}
devtools::install_github("mummlab/ARQiv")
library(ARQiv)
GUI()
```

## Further questions
Author:  David T. White, Arife Unal Eroglu, Guohua Wang,Liyun Zhang,Sumitra Sengupta,Ding Ding,Surendra Rajpurohit,Steven L. Walker,Hongkai Ji,Jiang Qian,Jeff S. Mumm

Any questions about the package, please contact
Maintainer: David T. White <dwhite66@jhmi.edu>, Ding Ding <dding6@jhu.edu>
or open an issue at Github.


