# Study_Mediation_SES_Strokeoutcome_2023
Code to perform the analyses in "Mediation analyses of the mechanisms by which socioeconomic status, comorbidity, stroke severity, and acute care influence stroke outcome"

Contains the following code files:
- "imputations.R" with the code to perform the multiple imputation of the NIHSS>5 variable
- "associations.R" with the code to fit the association models for the results in Table 2 and eTable 1.
- "MCsim_func.R" with the code for the function to estimate the interventional direct and indirect effects.
- "model_indices.R" list of indexes and corresponding variables used in the MCsim_func function
- "call.R" the function call used to estimate the interventional direct and indirect effects and bootstrap standard errors
- "pooling.R" the code used to pool the interventional direct and indirect effects and bootstrap standard errors to obtain the results in Table 3.
