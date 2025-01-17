# Title        : Template script for the MLE of haplotype frequencies, MOI, and prevalence
#                from example dataset.
# Objective    : Estimate haplotype frequencies, MOI, and prevalence
#                from genomic/molecular data
# Created by   : Christian Tsoungui Obama
# Created on   : 25.04.22
# Last modified: 26.04.22

# Install the necessary packages if necessary
install.packages('xlsx')   # Comment this line if xlsx installed

# Loading libraries
library(xlsx)

# Import the dataset
DATA <- read.xlsx('/home/janedoe/Documents/example.xlsx', 1, header = TRUE)

# Load external resources
source("/home/janedoe/Documents/SNPModel.R")

# Find the MLEs
est <- mle(DATA, id=TRUE)

# Estimate prevalence
## Unobservable prevalence
unobsprev <- estunobsprev(est)

## Conditional prevalence
condprev <- estcondprev(est)

## Relative prevalence
relprev <- estrelprev(DATA, id=TRUE)
