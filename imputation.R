#########################################################
### Multiple imputation of the binary NIHSSA variable ###
#########################################################

## Input:
# data_orig: the original data with missing values. Contains the following variables:
# NIH_NA: Indicator taking the value 1 if NIHSS is missing and 0 otherwise
# NIH_2cat: The dichotomized NIHSS variable, 1 = NIHSS<=5, 0 = NIHSS>5.
# DoD: The outcome dead or dependent in ADL at 3 months (1=yes, 0=no)
# Mid: Dummy variable for mid SES
# High: Dummy variable for high SES 
# KnownSmoker: 1 = Known smoker, 0 = non-smoker or unknown
# diabetes: 1 = yes, 0 = no.
# AF_bl: atrial fibrillation at baseline (1=yes, 0=no)
# prevstroke: previous stroke (1=yes, 0=no)
# BlodtryckssankVIA: Antihypertensives at baseline (1=yes, 0=no)
# LipidVIA: Statins at baseline (1=yes, 0=no)
# TrombocytHammareVIA: Antiplatelets at baseline (1=yes, 0=no)
# AntikoagulantiaVIA: Anticoagulants at baseline (1=yes, 0=no)
# reperfusion: reperfusion therapy received (1=yes, 0=no)
# StrokeUnit: Treated at stroke unit
# sex: 1=Male, 0=Female
# Alder: Age at time of stroke
# Ambulans_3cat: Arrived to hospital in ambulance (yes/no/no information)
# Tid_form: time from stroke onset to hospital arrival (<3h, 3–<4.5h, 4.5–<6h, 6–24h, >24h, no information)
# ConsR: : level of consciousness based on the Reaction Level Scale (RLS =1 vs. RLS >1)


###########################################################################################################################################################
##                                             STEP 1: Prepare dataset for imputations                                                                   ##
## The goal of this step is a dataset including the variables to be used for the imputation, including variables for interactions and higher order terms ##
## This dataset contains only complete variables except for the binary NIHSS variable NIH_2cat that is to be imputed                                     ##
###########################################################################################################################################################

# Model for missing NIHSS values given the analysis variables (including interaction and higher order terms for compatibility with analysis model(s)) and auxiliary variables. #
# The purpose of this model is to obtain a model matrix including relevant variables to be used for the imputations #
model_vars <- glm(NIH_NA~(DoD + Mid + High + KnownSmoker + diabetes + AF_bl + prevstroke + BlodtryckssankVIA + LipidVIA + TrombocytHammareVIA + AntikoagulantiaVIA +
                      reperfusion + StrokeUnit )^2 + sex  + Alder+ I(Alder^2) + Ambulans_3cat + Tid_form + ConsR - Mid:High - AntikoagulantiaVIA:reperfusion, data = data_orig, family = binomial)

model_matrix_vars <- model.matrix(model_vars) # Model matrix with variables for imputation including interaction terms and excluding missing values on variables other than NIHSS
model_matrix_vars <- model_matrix_vars[,-1] # Remove the column of 1:s corresponding to the intercept in model_vars

# Isolate the dichotomized NIHSS variable from the original dataset, removing patients missing values on variables other than NIHSS #
data_complcase <- subset(data_orig,(!is.na(data_orig$DoD) &!is.na(data_orig$Mid) &!is.na(data_orig$High) & !is.na(data_orig$KnownSmoker) & !is.na(data_orig$diabetes) & !is.na(data_orig$AF_bl) & 
                                         !is.na(data_orig$prevstroke) & !is.na(data_orig$BlodtryckssankVIA) & !is.na(data_orig$LipidVIA) &
                                         !is.na(data_orig$TrombocytHammareVIA) & !is.na(data_orig$AntikoagulantiaVIA) &  
                                         !is.na(data_orig$reperfusion) & !is.na(data_orig$StrokeUnit) & !is.na(data_orig$ConsR) ))  

# Put together the dichotomized NIHSS variable that is to be imputed with the variables to be used for the imputations:
data_for_imputations <- cbind(NIH_2cat=data_complcase$NIH_2cat,model_matrix_vars)
data_for_imputations <- data.frame(data_for_imputations)

##########################################################
### STEP 2: Multiple imputation using the mice package ###
##########################################################

library(mice)

imp_mult_int <- mice(data_for_imputations, seed = 591025, m = 45, method = "logreg") # Create 45 imputed datasets

###########################
### STEP 3: Diagnostics ###
###########################

# Check how well observed and imputed data match
densityplot(imp_mult_int, ~NIH_2cat, xlab="NIHSS>5")

# Check convergence
plot(imp_mult_int, c("NIH_2cat"), layout = c(2,1))

#############################################################
### STEP 4: Extract imputed datasets for further analyses ###
#############################################################

# Example with the first dataset, this is repeated for the next 44 datasets #
data_multimp_int1 <- complete(imp_mult_int,action = 1)
data_multimp_int1 <- data_multimp_int1[,c(1:17)] # Select only the analysis variables

# For simplicity, the outcome and mediator names are changed for further use in the analyses: 
# M2: The dichotomized NIHSS variable, 1 = NIHSS<=5, 0 = NIHSS>5.
# Y: The outcome dead or dependent in ADL at 3 months (1=yes, 0=no)
# M11, M12, M13, M14, M15, M16, M17, M18: Risk factors and prescribed medications at baseline (Known smoker, diabetes, atrial fibrillation, previous stroke, antihypertensives,
# statins, antiplatelets, anticoagulants).
# M3: Reperfusion therapy
# M4: Stroke unit
# AlderSq: Age squared

names(data_multimp_int1) <- c("M2","Y","Mid","High","M11","M12","M13",
                              "M14","M15","M16","M17","M18","M3",
                              "M4","sex","Alder","AlderSq") 



