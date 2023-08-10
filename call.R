################################################
# The function call used for the main analyses #
################################################

## Example based on one dataset, repeated on all 45 datasets ##

# Load the required packages #
library(zoo)
library(boot)
library(fastglm)

# Since the estimation is time consuming parallelization was used #
# Here the clusters are prepared for parallelization #

cl <- makeCluster(8)
clusterExport(cl,"MCsim_func")
clusterEvalQ(cl,library(fastglm))
clusterEvalQ(cl,library(boot))
clusterEvalQ(cl,library(zoo))

# Set seed for reproducibility, different seeds for each data set. Note that due to the parallelization only the #
# point estimates can be reproduced, rerunning the analysis will yield slightly different standard errors.       #

set.seed(220911) 

# Effect estimation and bootstrap SEs #
# Run the estimation function through the boot function to obtain standard errors. The parallel option is used for parallelization. #

time <- proc.time() # Start time

res_int <- boot(data = data_multimp_int1, statistic = MCsim_func,
             mcsim = 200, stype = "i", R = 1000, parallel = "snow", ncpus = 8,cl=cl) 

proc.time() - time # Time consumption

stopCluster(cl)

# The estimated effects are accessed through res_int$t0
# The bootstrap resamples for estimation of SEs are accessed through
# res_int$t. For pooling of the results, see file pooling.
