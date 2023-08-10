#####################################################
# Association models based on multiply imputed data #
#####################################################

# Load the mice-package #
library(mice)


## Unadjusted: ##

# Exposure-outcome: #
expout_unadj <- with(imp_mult_int, glm(DoD~Mid+High, family = binomial))

pool_expout_unadj <- summary(pool(expout_unadj))

# Point estimates and CI-bounds (we take the inverse to get low vs. mid and low vs. high):
round(1/exp(pool_expout_unadj$estimate[2:3]),2)
round(1/exp(pool_expout_unadj$estimate[2:3]-1.96*pool_expout_unadj$std.error[2:3]),2)
round(1/exp(pool_expout_unadj$estimate[2:3]+1.96*pool_expout_unadj$std.error[2:3]),2)

# Exposure-mediator (code only for the first mediator, but repeated for all mediators): #
expmed_unadj <- with(imp_mult_int, glm(KnownSmoker~Mid+High, family = binomial))

pool_expmed_unadj <- summary(pool(expmed_unadj))

# Point estimates and CI-bounds (we take the inverse to get low vs. mid and low vs. high):
round(1/exp(pool_expmed_unadj$estimate[2:3]),2)
round(1/exp(pool_expmed_unadj$estimate[2:3]-1.96*pool_expmed_unadj$std.error[2:3]),2)
round(1/exp(pool_expmed_unadj$estimate[2:3]+1.96*pool_expmed_unadj$std.error[2:3]),2)


# Mediator-outcome (code only for the first mediator, but repeated for all mediators) #
medout_unadj <- with(imp_mult_int, glm(DoD~KnownSmoker, family = binomial))

pool_medout_unadj <- summary(pool(medout_unadj))

# Point estimates and CI-bounds
round(exp(pool_medout_unadj$estimate[2]),2)
round(exp(pool_medout_unadj$estimate[2]-1.96*pool_medout_unadj$std.error[2]),2)
round(exp(pool_medout_unadj$estimate[2]+1.96*pool_medout_unadj$std.error[2]),2)


## Adjusted for confounders ##

# Exposure-outcome: #
expout_adj <- with(imp_mult_int, glm(DoD~Mid+High+sex+Age+ I(Age^2), family = binomial))

pool_expout_adj <- summary(pool(expout_adj))

# Point estimates and CI-bounds (we take the inverse to get low vs. mid and low vs. high):
round(1/exp(pool_expout_adj$estimate[2:3]),2)
round(1/exp(pool_expout_adj$estimate[2:3]-1.96*pool_expout_adj$std.error[2:3]),2)
round(1/exp(pool_expout_adj$estimate[2:3]+1.96*pool_expout_adj$std.error[2:3]),2)

# Exposure-mediator (code only for the first mediator, but repeated for all mediators): #
expmed_adj <- with(imp_mult_int, glm(KnownSmoker~Mid+High+sex+Age+ I(Age^2), family = binomial))

pool_expmed_adj <- summary(pool(expmed_adj))

# Point estimates and CI-bounds (we take the inverse to get low vs. mid and low vs. high):
round(1/exp(pool_expmed_adj$estimate[2:3]),2)
round(1/exp(pool_expmed_adj$estimate[2:3]-1.96*pool_expmed_adj$std.error[2:3]),2)
round(1/exp(pool_expmed_adj$estimate[2:3]+1.96*pool_expmed_adj$std.error[2:3]),2)

# Mediator-outcome (code only for the first mediator, but repeated for all mediators) #
medout_adj <- with(imp_mult_int, glm(DoD~KnownSmoker+Mid+High+sex  + Age+ I(Age^2), family = binomial))

pool_medout_adj <- summary(pool(medout_adj))

# Point estimates and CI-bounds
round(exp(pool_medout_adj$estimate[2]),2)
round(exp(pool_medout_adj$estimate[2]-1.96*pool_medout_adj$std.error[2]),2)
round(exp(pool_medout_adj$estimate[2]+1.96*pool_medout_adj$std.error[2]),2)


## Fully adjusted mediator-outcome model ##

fully_adj <- with(imp_mult_int, glm(DoD~Mid + High + KnownSmoker + diabetes + AF_bl + prevstroke + Antihypertensives + Statins + Antithrombotics +
                                         Anticoagulants +
                                         +NIHSSover5 + reperfusion +
                                         StrokeUnit +sex  + Age+ I(Age^2), family = binomial))

pool_fully_adj <- summary(pool(fully_adj))
 
# Exposure-outcome point estimates and CIs. Take the inverse to get low vs. mid and low vs. high #
round(1/exp(pool_fully_adj$estimate[2:3]),2)
round(1/exp(pool_fully_adj$estimate[2:3]-1.96*pool_fully_adj$std.error[2:3]),2)
round(1/exp(pool_fully_adj$estimate[2:3]+1.96*pool_fully_adj$std.error[2:3]),2)

# Mediator-outcome point estimates and CIs. #
round(exp(pool_fully_adj$estimate[4:14]),2)
round(exp(pool_fully_adj$estimate[4:14]-1.96*pool_fully_adj$std.error[4:14]),2)
round(exp(pool_fully_adj$estimate[4:14]+1.96*pool_fully_adj$std.error[4:14]),2)
