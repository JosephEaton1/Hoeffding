# Hoeffding's D Test and Measure of Association

Repository contains content relevant to 'Hoeffding's D Test and Measure of Association', written by Joseph W.D. Eaton and William F. Scott.

## Repository Contents

### Functions
- **DnFunction.R** - Contains functions to calculate Dn and the variance of a given sample (using zeta formulae from Section 6), as well as giving Qn and critical values of Qn (using formulae from Section 4).
- **QFunctions.R** - An R script to reproduce the limiting distribution, Q, table (using method from Section 4).
- **ExactDnFunctions.R** - Contains functions used to tabulate the exact discrete distribution of Dn (using Otten's method, see Section 2).
### Data
- **Example_dataframes.R** - Includes data frames used in Examples 7.1 and 7.2.
- **ExactDnData.RData** - Contains an R environment with complete Dn tables, the data is stored in tabs_clean a list such that tabs_clean[[n-4]] contains a full distribution of Dn. Note, refer to Tables folder for quick quantile tables, this data is quite messy.
### Tables
- **Dn_quantiles.png** - Table of quantiles for Dn from n = 5 to 18
- **Q_distribution.png** - Table of limiting distribution Q
- **Q_quantiles.png** - Table of quantiles of Q
