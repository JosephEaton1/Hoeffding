# Hoeffding's D Test and Measure of Association

Repository contains content relevant to 'Hoeffding's D Test and Measure of Association', written by Joseph W.D. Eaton and William F. Scott.

## Repository Contents

- **Example_dataframes.R** - Includes data frames used in examples 7.1 and 7.2.
- **ExactDnData.RData** - Contains an R environment with complete Dn tables, the data is stored in tabs_clean a list such that tabs_clean[[n-4]] contains a full distribution of Dn. Note, refer to Tables folder for quick quantile tables, this data is quite messy.
- **DnFunction.R** - Contains functions to calculate Dn and the variance of a given sample (using zeta formulae from chapter 6), as well as giving Qn and critical values of Qn (using formulae from chapter 4).
- **QFunctions.R** - An R script to reproduce the limiting distribution, Q, table (using method from chapter 4).
- **ExactDnFunctions.R** - Contains functions used to tabulate the exact discrete distribution of Dn (using Otten's method, see chapter 2).
