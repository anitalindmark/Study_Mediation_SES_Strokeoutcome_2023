###############################################################
## Function for Monte Carlo simulation to estimate the       ##
##           interventional disparity effects                ##
##                                                           ##
## Note that the function is written specifically for the    ##
## application in the article and needs to be adapted if it  ##
## is to be used in other applications.                      ##
##                                                           ##
## Utilizes the fastglm package for modeling and prediction. ##
## The models are specified based on matrices rather than    ##
## formulas, for a list of indices used see file             ##
## model_indices.                                            ##
##                                                           ##
###############################################################

# Input parameters: 
# dat = data frame, 
# ind = indices for bootstrap, 
# mcsim = number of monte carlo simulations,
# method = fitting method (see documentation fastglm function in package fastglm).

MCsim_func <- function(dat, ind = 1:nrow(dat), mcsim = 200, method = 0) {
  
  #Take bootstrap sample
  data <- dat[ind, ] # 
  
  # Create matrices to be used by fastglm
  X <- model.matrix(Y ~ (Mid + High + M11 + M12 + M13 + M14 + M15 + M16 + M17 + M18 +
                           M2 + M3 + M4)^2 +   
                           sex  + Alder  + AlderSq - Mid:High , data = data) # 
  Y <- as.matrix(data$Y) 
  
  M11 <- as.matrix(data$M11) # Known smoker
  M12 <- as.matrix(data$M12) # Diabetes
  M13 <- as.matrix(data$M13) # Atrial fibrillation
  M14 <- as.matrix(data$M14) # Previous stroke
  M15 <- as.matrix(data$M15) # Antihypertensives
  M16 <- as.matrix(data$M16) # Statins
  M17 <- as.matrix(data$M17) # Antithrombotics
  M18 <- as.matrix(data$M18) # Anticoagulants
  
  M2 <- as.matrix(data$M2) # NIHSS>5
  
  M3 <- as.matrix(data$M3) # Reperfusion therapy
 
  M4 <- as.matrix(data$M4) # Stroke unit care
  
  
  # Set flag to capture bootstrap samples to reject
  flag <- FALSE 
  
  # Model matrices for new predictions. One for the low vs. mid SES contrast, one for the low vs. high contrast.
  X2 <- X[which(data$High!=1),] # Low v mid
  X3 <- X[which(data$Mid!=1),] # Low v high
  n2 <- nrow(X2)
  n3 <- nrow(X3)
  
  # Pre-allocate vectors for effect estimation
  y1_m12340 <- y0_m12340 <- y1_m11_m2340 <- y1_m10_m2340 <- y1_m11_m21_m340 <- y1_m11_m20_m340 <- rep(NA,mcsim)
  y1_m12341 <- y1_m121_m31_m40 <- y1_m121_m30_m40 <- y1_m1231_m41 <- y1_m1231_m40 <- rep(NA,mcsim)
  
  
  ### Model fits ###
  
  ## Marginal models (not given other mediators (mediator groups)) ##
  
  # M1 - given exposure, confounders and previous mediators in the M1 group #
  
  fit_M11_marg <- fastglm(x=X[,c(1:3,15:17)],y=M11, family = binomial(),method=method) #
  if((!fit_M11_marg$converged)|any(is.na(fit_M11_marg$coefficients))) flag <- TRUE
  
  fit_M12_marg <- fastglm(x=X[,c(1:4,15:17,18,29)],y=M12, family = binomial(),method=method) # 
  if((!fit_M12_marg$converged)|any(is.na(fit_M12_marg$coefficients))) flag <- TRUE
  
  fit_M13_marg <- fastglm(x=X[,c(1:5,15:17,18:19,29:30,40)],y=M13, family = binomial(),method=method) # 
  if((!fit_M13_marg$converged)|any(is.na(fit_M13_marg$coefficients))) flag <- TRUE
  
  fit_M14_marg <- fastglm(x=X[,c(1:6,15:17,18:20,29:31,40:41, 50)],y=M14, family = binomial(),method=method) # 
  if((!fit_M14_marg$converged)|any(is.na(fit_M14_marg$coefficients))) flag <- TRUE
  
  fit_M15_marg <- fastglm(x=X[,c(1:7,15:17,18:21,29:32,40:42, 50:51,59)],y=M15, family = binomial(),method=method) # 
  if((!fit_M15_marg$converged)|any(is.na(fit_M15_marg$coefficients))) flag <- TRUE
  
  fit_M16_marg <- fastglm(x=X[,c(1:8,15:17,18:22,29:33,40:43, 50:52,59:60,67)],y=M16, family = binomial(),method=method) # 
  if((!fit_M16_marg$converged)|any(is.na(fit_M16_marg$coefficients))) flag <- TRUE
  
  fit_M17_marg <- fastglm(x=X[,c(1:9,15:17,18:23,29:34,40:44, 50:53,59:61,67:68,74)],y=M17, family = binomial(),method=method) # 
  if((!fit_M17_marg$converged)|any(is.na(fit_M17_marg$coefficients))) flag <- TRUE
  
  fit_M18_marg <- fastglm(x=X[,c(1:10,15:17,18:24,29:35,40:45, 50:54,59:62,67:69,74:75,80)],y=M18, family = binomial(),method=method) # 
  if((!fit_M18_marg$converged)|any(is.na(fit_M18_marg$coefficients))) flag <- TRUE
  
  
  # M2 - given exposure and confounders #
  
  fit_M2_marg <- fastglm(x=X[,c(1:3,15:17)],y=M2,family = binomial(),method=method) # 
  if((!fit_M2_marg$converged)|any(is.na(fit_M2_marg$coefficients))) flag <- TRUE
  
  # M3 - given exposure and confounders #
  
  fit_M3_marg <- fastglm(x=X[,c(1:3,15:17)],y=M3, family = binomial(),method=method) # 
  if((!fit_M3_marg$converged)|any(is.na(fit_M3_marg$coefficients))) flag <- TRUE
  
  # M4 - given exposure and confounders #
  
  fit_M4_marg <- fastglm(x=X[,c(1:3,15:17)],y=M4, family = binomial(),method=method) # 
  if((!fit_M4_marg$converged)|any(is.na(fit_M4_marg$coefficients))) flag <- TRUE
  
   
  
  ### Conditional models ###
  
  ## Conditional on immediately preceding mediator (mediator group) ##
  
  # M2 - given exposure and confounders and mediators in the M1 group
  fit_M2_cond <- fastglm(x=X[,c(1:11,15:17,18:25,29:36,40:46,50:55,59:63,67:70,74:76,80:81,85)],y=M2,family = binomial(),method=method) # 
  if((!fit_M2_cond$converged)|any(is.na(fit_M2_cond$coefficients))) flag <- TRUE
  
  # M3 - given exposure and confounders and M2
  fit_M3_condM2 <- fastglm(x=X[,c(1:3,12,15:17,26,37)],y=M3, family = binomial(),method=method) # 
  if((!fit_M3_condM2$converged)|any(is.na(fit_M3_condM2$coefficients))) flag <- TRUE
  
  # M4 - given exposure and confounders and M3
  fit_M4_condM3 <- fastglm(x=X[,c(1:3,13,15:17,27,38)],y=M4, family = binomial(),method=method) # 
  if((!fit_M4_condM3$converged)|any(is.na(fit_M4_condM3$coefficients))) flag <- TRUE
  
  
  ## Conditional on the two immediately preceding mediators (mediator group) ##
  
  # M3 - given exposure, confounders, mediators in the M1 group and M2
  fit_M3_condM1M2 <- fastglm(x=X[,c(1:11,12,15:17,18:26,29:37,40:47,50:56,59:64,67:71,74:77,80:82,85:86,89)],y=M3, family = binomial(),method=method) # 
  if((!fit_M3_condM1M2$converged)|any(is.na(fit_M3_condM1M2$coefficients))) flag <- TRUE
  
  #M4 - given exposure, confounders, M2 and M3
  fit_M4_condM2M3 <- fastglm(x=X[,c(1:3,12:13,15:17,26:27,37:38,92)],y=M4, family = binomial(),method=method) # 
  if((!fit_M4_condM2M3$converged)|any(is.na(fit_M4_condM2M3$coefficients))) flag <- TRUE
  
  
  ## Conditional on the three immediately preceding mediators (mediator group) ##
  
  # M4 - given exposure, confounders, mediators in the M1 group, M2 and M3
  fit_M4_condM1M2M3 <- fastglm(x=X[,c(1:13,15:17,18:27,29:38,40:48,50:57,59:65,67:72,74:78,80:83,85:87,89,92)],y=M4, family = binomial(),method=method) # 
  if((!fit_M4_condM1M2M3$converged)|any(is.na(fit_M4_condM1M2M3$coefficients))) flag <- TRUE
  
  
  
  ## Outcome model ##
  
  # Y - given exposure, confounders, and all mediators #
  
  fit_Y <- fastglm(x=X,y=Y,family = binomial(),method=method)
  if((!fit_Y$converged)|any(is.na(fit_Y$coefficients))) flag<-TRUE
  
  ### Predictions ###
  
  ## Predicted values that only depend on exposure and confounder values ##
  ##      and therefore do not need to be updated with new MC-draws      ##
  
  X2[,2] <- 1 # Set Mid = 1, which corresponds to unexposed here
  M11_preds0_marg <- predict(fit_M11_marg, newdata = X2[,c(1:3,15:17)], type = "response")
  M2_preds0_marg <- predict(fit_M2_marg, newdata = X2[,c(1:3,15:17)], type = "response")
  M3_preds0_marg <- predict(fit_M3_marg, newdata = X2[,c(1:3,15:17)], type = "response")
  M4_preds0_marg <- predict(fit_M4_marg, newdata = X2[,c(1:3,15:17)], type = "response")
  X2[,2]<-0 # Set Mid = 0, which corresponds to exposed (low ses) here
  M11_preds1_marg <- predict(fit_M11_marg, newdata = X2[,c(1:3,15:17)], type = "response") 
  M2_preds1_marg <- predict(fit_M2_marg, newdata = X2[,c(1:3,15:17)], type = "response")
  M3_preds1_marg <- predict(fit_M3_marg, newdata = X2[,c(1:3,15:17)], type = "response")
  M4_preds1_marg <- predict(fit_M4_marg, newdata = X2[,c(1:3,15:17)], type = "response")
  
  X3[,3] <- 1 # Set High = 1, which corresponds to unexposed here
  M11_preds0_marg3 <- predict(fit_M11_marg, newdata = X3[,c(1:3,15:17)], type = "response")
  M2_preds0_marg3 <- predict(fit_M2_marg, newdata = X3[,c(1:3,15:17)], type = "response")
  M3_preds0_marg3 <- predict(fit_M3_marg, newdata = X3[,c(1:3,15:17)], type = "response")
  M4_preds0_marg3 <- predict(fit_M4_marg, newdata = X3[,c(1:3,15:17)], type = "response")
  X3[,3]<-0 # Set High = 0, which corresponds to exposed (low ses) here
  M11_preds1_marg3 <- predict(fit_M11_marg, newdata = X3[,c(1:3,15:17)], type = "response") 
  M2_preds1_marg3 <- predict(fit_M2_marg, newdata = X3[,c(1:3,15:17)], type = "response")
  M3_preds1_marg3 <- predict(fit_M3_marg, newdata = X3[,c(1:3,15:17)], type = "response")
  M4_preds1_marg3 <- predict(fit_M4_marg, newdata = X3[,c(1:3,15:17)], type = "response")

 
  ## Loop for Monte Carlo simulations for the low vs. mid contrast ##
  
  for(i in 1:mcsim){ 
    
    ### Simulate from mediator distributions ###

    ## Marginal distributions ##

    # For the unexposed condition (mid SES) #
    
    X2[,2] <- 1 # Set the unexposed condition
    
    M11_0_marg <- rbinom(n2, 1, M11_preds0_marg)
    X2[,4] <- M11_0_marg # Update the M11 variable with simulated values
    X2[,18] <- M11_0_marg # Update the Mid:M11-interaction with simulated values
    
    M12_0_marg <- rbinom(n2, 1, predict(fit_M12_marg, newdata = X2[,c(1:4,15:17,18,29)], type = "response"))
    X2[,5] <- M12_0_marg
    X2[,19] <-  M12_0_marg # Mid:M12
    X2[,40] <- M11_0_marg* M12_0_marg # M11:M12
    
    M13_0_marg <- rbinom(n2, 1, predict(fit_M13_marg, newdata = X2[,c(1:5,15:17,18:19,29:30,40)], type = "response"))
    X2[,6] <- M13_0_marg
    X2[,20] <- M13_0_marg # Mid:M13
    X2[,41] <- M11_0_marg*M13_0_marg # M11:M13
    X2[,50] <- M12_0_marg*M13_0_marg # M12:M13
    
    M14_0_marg <- rbinom(n2, 1, predict(fit_M14_marg, newdata = X2[,c(1:6,15:17,18:20,29:31,40:41, 50)], type = "response"))
    X2[,7] <- M14_0_marg
    X2[,21] <- M14_0_marg # Mid:M14
    X2[,42] <- M11_0_marg*M14_0_marg # M11:M14
    X2[,51] <- M12_0_marg*M14_0_marg # M12:M14
    X2[,59] <- M13_0_marg*M14_0_marg # M13:M14
    
    M15_0_marg <- rbinom(n2, 1, predict(fit_M15_marg, newdata = X2[,c(1:7,15:17,18:21,29:32,40:42, 50:51,59)], type = "response"))
    X2[,8] <- M15_0_marg
    X2[,22] <- M15_0_marg # Mid:M15
    X2[,43] <- M11_0_marg*M15_0_marg # M11:M15
    X2[,52] <- M12_0_marg*M15_0_marg # M12:M15
    X2[,60] <- M13_0_marg*M15_0_marg #M13:M15
    X2[,67] <- M14_0_marg*M15_0_marg # M14:M15
    
    M16_0_marg <- rbinom(n2, 1, predict(fit_M16_marg, newdata = X2[,c(1:8,15:17,18:22,29:33,40:43, 50:52,59:60,67)], type = "response"))
    X2[,9] <- M16_0_marg
    X2[,23] <- M16_0_marg # Mid:M16
    X2[,44] <- M11_0_marg*M16_0_marg # M11:M16
    X2[,53] <- M12_0_marg*M16_0_marg # M12:M16
    X2[,61] <- M13_0_marg*M16_0_marg # M13:M16
    X2[,68] <- M14_0_marg*M16_0_marg # M14:M16
    X2[,74] <- M15_0_marg*M16_0_marg # M15:M16
    
    M17_0_marg <- rbinom(n2, 1, predict(fit_M17_marg, newdata = X2[,c(1:9,15:17,18:23,29:34,40:44, 50:53,59:61,67:68,74)], type = "response"))
    X2[,10] <- M17_0_marg
    X2[,24] <- M17_0_marg # Mid:M17
    X2[,45] <- M11_0_marg*M17_0_marg # M11:M17
    X2[,54] <- M12_0_marg*M17_0_marg # M12:M17
    X2[,62] <- M13_0_marg*M17_0_marg # M13:M17
    X2[,69] <- M14_0_marg*M17_0_marg # M14:M17
    X2[,75] <- M15_0_marg*M17_0_marg # M15:M17
    X2[,80] <- M16_0_marg*M17_0_marg # M16:M17
    
    M18_0_marg <- rbinom(n2, 1, predict(fit_M18_marg, newdata = X2[,c(1:10,15:17,18:24,29:35,40:45, 50:54,59:62,67:69,74:75,80)], type = "response"))

    M2_0_marg <- rbinom(n2, 1, M2_preds0_marg)
    M3_0_marg <- rbinom(n2, 1, M3_preds0_marg)
    M4_0_marg <- rbinom(n2, 1, M4_preds0_marg)

    # For the exposed condition (low SES) #
    
    X2[,2]<-0 # Set the exposed condition
    
    M11_1_marg <- rbinom(n2, 1, M11_preds1_marg)
    X2[,4] <- M11_1_marg # Update the M11 variable with simulated values
    X2[,18] <- 0 # Update the Mid:M11-interaction with simulated values
    
    M12_1_marg <- rbinom(n2, 1, predict(fit_M12_marg, newdata = X2[,c(1:4,15:17,18,29)], type = "response"))
    X2[,5] <- M12_1_marg
    X2[,19] <- 0 # Mid:M12
    X2[,40] <- M11_1_marg* M12_1_marg # M11:M12
    
    M13_1_marg <- rbinom(n2, 1, predict(fit_M13_marg, newdata = X2[,c(1:5,15:17,18:19,29:30,40)], type = "response"))
    X2[,6] <- M13_1_marg
    X2[,20] <- 0  # Mid:M13
    X2[,41] <- M11_1_marg*M13_1_marg # M11:M13
    X2[,50] <- M12_1_marg*M13_1_marg # M12:M13
    
    M14_1_marg <- rbinom(n2, 1, predict(fit_M14_marg, newdata = X2[,c(1:6,15:17,18:20,29:31,40:41, 50)], type = "response"))
    X2[,7] <- M14_1_marg
    X2[,21] <- 0 # Mid:M14
    X2[,42] <- M11_1_marg*M14_1_marg # M11:M14
    X2[,51] <- M12_1_marg*M14_1_marg # M12:M14
    X2[,59] <- M13_1_marg*M14_1_marg # M13:M14
    
    M15_1_marg <- rbinom(n2, 1, predict(fit_M15_marg, newdata = X2[,c(1:7,15:17,18:21,29:32,40:42, 50:51,59)], type = "response"))
    X2[,8] <- M15_1_marg
    X2[,22] <- 0 # Mid:M15
    X2[,43] <- M11_1_marg*M15_1_marg # M11:M15
    X2[,52] <- M12_1_marg*M15_1_marg # M12:M15
    X2[,60] <- M13_1_marg*M15_1_marg #M13:M15
    X2[,67] <- M14_1_marg*M15_1_marg # M14:M15
    
    M16_1_marg <- rbinom(n2, 1, predict(fit_M16_marg, newdata = X2[,c(1:8,15:17,18:22,29:33,40:43, 50:52,59:60,67)], type = "response"))
    X2[,9] <- M16_1_marg
    X2[,23] <- 0 # Mid:M16
    X2[,44] <- M11_1_marg*M16_1_marg # M11:M16
    X2[,53] <- M12_1_marg*M16_1_marg # M12:M16
    X2[,61] <- M13_1_marg*M16_1_marg # M13:M16
    X2[,68] <- M14_1_marg*M16_1_marg # M14:M16
    X2[,74] <- M15_1_marg*M16_1_marg # M15:M16
    
    M17_1_marg <- rbinom(n2, 1, predict(fit_M17_marg, newdata = X2[,c(1:9,15:17,18:23,29:34,40:44, 50:53,59:61,67:68,74)], type = "response"))
    X2[,10] <- M17_1_marg
    X2[,24] <- 0 # Mid:M17
    X2[,45] <- M11_1_marg*M17_1_marg # M11:M17
    X2[,54] <- M12_1_marg*M17_1_marg # M12:M17
    X2[,62] <- M13_1_marg*M17_1_marg # M13:M17
    X2[,69] <- M14_1_marg*M17_1_marg # M14:M17
    X2[,75] <- M15_1_marg*M17_1_marg # M15:M17
    X2[,80] <- M16_1_marg*M17_1_marg # M16:M17
    
    M18_1_marg <- rbinom(n2, 1, predict(fit_M18_marg, newdata = X2[,c(1:10,15:17,18:24,29:35,40:45, 50:54,59:62,67:69,74:75,80)], type = "response"))

    M2_1_marg <- rbinom(n2, 1, M2_preds1_marg)
    M3_1_marg <- rbinom(n2, 1, M3_preds1_marg)
    M4_1_marg <- rbinom(n2, 1, M4_preds1_marg)

    ## Conditional distributions ##

    # M2 - given the M1 group #

    # For the unexposed condition (mid SES) #
    
    X2[,2] <- 1
    
    X2[,4] <- M11_0_marg
    X2[,5] <- M12_0_marg
    X2[,6] <- M13_0_marg
    X2[,7] <- M14_0_marg
    X2[,8] <- M15_0_marg
    X2[,9] <- M16_0_marg
    
    X2[,10] <- M17_0_marg
    X2[,11] <- M18_0_marg 
    X2[,18] <- M11_0_marg # Mid:M11
    X2[,19] <- M12_0_marg
    X2[,20] <- M13_0_marg
    X2[,21] <- M14_0_marg
    X2[,22] <- M15_0_marg
    X2[,23] <- M16_0_marg
    X2[,24] <- M17_0_marg
    X2[,25] <- M18_0_marg # Mid:M18
    
    X2[,40] <- M11_0_marg* M12_0_marg # M11:M12
    X2[,41] <- M11_0_marg*M13_0_marg # M11:M13
    X2[,42] <- M11_0_marg*M14_0_marg # M11:M14
    X2[,43] <- M11_0_marg*M15_0_marg # M11:M15
    X2[,44] <- M11_0_marg*M16_0_marg # M11:M16
    X2[,45] <- M11_0_marg*M17_0_marg # M11:M17
    X2[,46] <- M11_0_marg*M18_0_marg # M11:M18
    X2[,50] <- M12_0_marg*M13_0_marg # M12:M13
    X2[,51] <- M12_0_marg*M14_0_marg # M12:M14
    X2[,52] <- M12_0_marg*M15_0_marg # M12:M15
    X2[,53] <- M12_0_marg*M16_0_marg # M12:M16
    X2[,54] <- M12_0_marg*M17_0_marg # M12:M17
    X2[,55] <- M12_0_marg*M18_0_marg # M12:M18
    X2[,59] <- M13_0_marg*M14_0_marg # M13:M14
    X2[,60] <- M13_0_marg*M15_0_marg # M13:M15
    X2[,61] <- M13_0_marg*M16_0_marg # M13:M16
    X2[,62] <- M13_0_marg*M17_0_marg # M13:M17
    X2[,63] <- M13_0_marg*M18_0_marg # M13:M18
    X2[,67] <- M14_0_marg*M15_0_marg # M14:M15
    X2[,68] <- M14_0_marg*M16_0_marg # M14:M16
    X2[,69] <- M14_0_marg*M17_0_marg # M14:M17
    X2[,70] <- M14_0_marg*M18_0_marg # M14:M18
    X2[,74] <- M15_0_marg*M16_0_marg # M15:M16
    X2[,75] <- M15_0_marg*M17_0_marg # M15:M17
    X2[,76] <- M15_0_marg*M18_0_marg # M15:M18
    X2[,80] <- M16_0_marg*M17_0_marg # M16:M17
    X2[,81] <- M16_0_marg*M18_0_marg # M16:M18
    X2[,85] <- M17_0_marg*M18_0_marg # M17:M18

    M2_0_cond <- rbinom(n2, 1, predict(fit_M2_cond, newdata = X2[,c(1:11,15:17,18:25,29:36,40:46,50:55,59:63,67:70,74:76,80:81,85)], type = "response"))


    # For the exposed condition (low SES) #
    
    X2[,2] <- 0
    
    X2[,4] <- M11_1_marg
    X2[,5] <- M12_1_marg
    X2[,6] <- M13_1_marg
    X2[,7] <- M14_1_marg
    X2[,8] <- M15_1_marg
    X2[,9] <- M16_1_marg
    X2[,10] <- M17_1_marg
    X2[,11] <- M18_1_marg
    
    X2[,18] <- 0 # Mid:M11
    X2[,19] <- 0
    X2[,20] <- 0
    X2[,21] <- 0
    X2[,22] <- 0
    X2[,23] <- 0
    X2[,24] <- 0
    X2[,25] <- 0 # Mid:M18
    
    X2[,40] <- M11_1_marg* M12_1_marg # M11:M12
    X2[,41] <- M11_1_marg*M13_1_marg # M11:M13
    X2[,42] <- M11_1_marg*M14_1_marg # M11:M14
    X2[,43] <- M11_1_marg*M15_1_marg # M11:M15
    X2[,44] <- M11_1_marg*M16_1_marg # M11:M16
    X2[,45] <- M11_1_marg*M17_1_marg # M11:M17
    X2[,46] <- M11_1_marg*M18_1_marg# M11:M18
    X2[,50] <- M12_1_marg*M13_1_marg # M12:M13
    X2[,51] <- M12_1_marg*M14_1_marg # M12:M14
    X2[,52] <- M12_1_marg*M15_1_marg # M12:M15
    X2[,53] <- M12_1_marg*M16_1_marg # M12:M16
    X2[,54] <- M12_1_marg*M17_1_marg # M12:M17
    X2[,55] <- M12_1_marg*M18_1_marg# M12:M18
    X2[,59] <- M13_1_marg*M14_1_marg # M13:M14
    X2[,60] <- M13_1_marg*M15_1_marg # M13:M15
    X2[,61] <- M13_1_marg*M16_1_marg # M13:M16
    X2[,62] <- M13_1_marg*M17_1_marg # M13:M17
    X2[,63] <- M13_1_marg*M18_1_marg# M13:M18
    X2[,67] <- M14_1_marg*M15_1_marg # M14:M15
    X2[,68] <- M14_1_marg*M16_1_marg # M14:M16
    X2[,69] <- M14_1_marg*M17_1_marg # M14:M17
    X2[,70] <- M14_1_marg*M18_1_marg# M14:M18
    X2[,74] <- M15_1_marg*M16_1_marg # M15:M16
    X2[,75] <- M15_1_marg*M17_1_marg # M15:M17
    X2[,76] <- M15_1_marg*M18_1_marg# M15:M18
    X2[,80] <- M16_1_marg*M17_1_marg # M16:M17
    X2[,81] <- M16_1_marg*M18_1_marg# M16:M18
    X2[,85] <- M17_1_marg*M18_1_marg# M17:M18

    M2_1_cond <- rbinom(n2, 1, predict(fit_M2_cond, newdata = X2[,c(1:11,15:17,18:25,29:36,40:46,50:55,59:63,67:70,74:76,80:81,85)], type = "response"))


    # M3 - given M2 #
    
    # For the unexposed condition (mid SES) #
    
    X2[,2] <- 1
    X2[,12] <- M2_0_marg
    X2[,26] <- M2_0_marg # Mid:M2

    M3_0_condM2 <- rbinom(n2, 1, predict(fit_M3_condM2, newdata = X2[,c(1:3,12,15:17,26,37)], type = "response"))

    # M4 - given M3 #
   
     # For the unexposed condition (mid SES) #
    
    X2[,13] <- M3_0_marg
    X2[,27] <- M3_0_marg # Mid:M3

    M4_0_condM3 <- rbinom(n2, 1, predict(fit_M4_condM3, newdata = X2[,c(1:3,13,15:17,27,38)], type = "response"))


    # M3 - given the M1 group and M2
    
    # For the exposed condition (low SES) #
    
    X2[,2] <- 0

    X2[,12] <- M2_1_cond
    X2[,26] <- 0 # Mid:M2

    X2[,47] <- M11_1_marg*M2_1_cond # M11:M2
    X2[,56] <- M12_1_marg*M2_1_cond # M12:M2
    X2[,64] <- M13_1_marg*M2_1_cond # M13:M2
    X2[,71] <- M14_1_marg*M2_1_cond # M14:M2
    X2[,77] <- M15_1_marg*M2_1_cond # M15:M2
    X2[,82] <- M16_1_marg*M2_1_cond # M16:M2
    X2[,86] <- M17_1_marg*M2_1_cond # M17:M2
    X2[,89] <- M18_1_marg*M2_1_cond # M18:M2

    M3_1_condM1M2 <- rbinom(n2, 1, predict(fit_M3_condM1M2, newdata = X2[,c(1:11,12,15:17,18:26,29:37,40:47,50:56,59:64,67:71,74:77,80:82,85:86,89)], type = "response"))

    # For the unexposed condition (mid SES) #
    
    X2[,2] <- 1
    
    X2[,4] <- M11_0_marg
    X2[,5] <- M12_0_marg
    X2[,6] <- M13_0_marg
    X2[,7] <- M14_0_marg
    X2[,8] <- M15_0_marg
    X2[,9] <- M16_0_marg
    X2[,10] <- M17_0_marg
    X2[,11] <- M18_0_marg
    
    X2[,18] <- M11_0_marg # Mid:M11
    X2[,19] <- M12_0_marg
    X2[,20] <- M13_0_marg
    X2[,21] <- M14_0_marg
    X2[,22] <- M15_0_marg
    X2[,23] <- M16_0_marg
    X2[,24] <- M17_0_marg
    X2[,25] <- M18_0_marg # Mid:M18

    X2[,40] <- M11_0_marg* M12_0_marg # M11:M12
    X2[,41] <- M11_0_marg*M13_0_marg # M11:M13
    X2[,42] <- M11_0_marg*M14_0_marg # M11:M14
    X2[,43] <- M11_0_marg*M15_0_marg # M11:M15
    X2[,44] <- M11_0_marg*M16_0_marg # M11:M16
    X2[,45] <- M11_0_marg*M17_0_marg # M11:M17
    X2[,46] <- M11_0_marg*M18_0_marg # M11:M18
    X2[,50] <- M12_0_marg*M13_0_marg # M12:M13
    X2[,51] <- M12_0_marg*M14_0_marg # M12:M14
    X2[,52] <- M12_0_marg*M15_0_marg # M12:M15
    X2[,53] <- M12_0_marg*M16_0_marg # M12:M16
    X2[,54] <- M12_0_marg*M17_0_marg # M12:M17
    X2[,55] <- M12_0_marg*M18_0_marg # M12:M18
    X2[,59] <- M13_0_marg*M14_0_marg # M13:M14
    X2[,60] <- M13_0_marg*M15_0_marg # M13:M15
    X2[,61] <- M13_0_marg*M16_0_marg # M13:M16
    X2[,62] <- M13_0_marg*M17_0_marg # M13:M17
    X2[,63] <- M13_0_marg*M18_0_marg # M13:M18
    X2[,67] <- M14_0_marg*M15_0_marg # M14:M15
    X2[,68] <- M14_0_marg*M16_0_marg # M14:M16
    X2[,69] <- M14_0_marg*M17_0_marg # M14:M17
    X2[,70] <- M14_0_marg*M18_0_marg # M14:M18
    X2[,74] <- M15_0_marg*M16_0_marg # M15:M16
    X2[,75] <- M15_0_marg*M17_0_marg # M15:M17
    X2[,76] <- M15_0_marg*M18_0_marg # M15:M18
    X2[,80] <- M16_0_marg*M17_0_marg # M16:M17
    X2[,81] <- M16_0_marg*M18_0_marg # M16:M18
    X2[,85] <- M17_0_marg*M18_0_marg # M17:M18

    X2[,12] <- M2_0_cond
    X2[,26] <- M2_0_cond # Mid:M2

    X2[,47] <- M11_0_marg*M2_0_cond # M11:M2
    X2[,56] <- M12_0_marg*M2_0_cond # M12:M2
    X2[,64] <- M13_0_marg*M2_0_cond # M13:M2
    X2[,71] <- M14_0_marg*M2_0_cond # M14:M2
    X2[,77] <- M15_0_marg*M2_0_cond # M15:M2
    X2[,82] <- M16_0_marg*M2_0_cond # M16:M2
    X2[,86] <- M17_0_marg*M2_0_cond # M17:M2
    X2[,89] <- M18_0_marg*M2_0_cond # M18:M2

    M3_0_condM1M2 <- rbinom(n2, 1, predict(fit_M3_condM1M2, newdata = X2[,c(1:11,12,15:17,18:26,29:37,40:47,50:56,59:64,67:71,74:77,80:82,85:86,89)], type = "response"))


    # M4 - given M2 and M3 #
    
    # For the unexposed condition (mid SES) #

    X2[,12] <- M2_0_marg
    X2[,13] <- M3_0_condM2

    X2[,26] <- M2_0_marg # Mid:M2
    X2[,27] <- M3_0_condM2 # Mid:M3

    X2[,92] <- M2_0_marg*M3_0_condM2 # M2:M3

    M4_0_condM2M3 <- rbinom(n2, 1, predict(fit_M4_condM2M3, newdata = X2[,c(1:3,12:13,15:17,26:27,37:38,92)], type = "response"))


    # M4 - given the M1 group, M2 and M3 #
    
    # For the unexposed condition (mid SES) #

    X2[,12] <- M2_0_cond

    X2[,13] <- M3_0_condM1M2

    X2[,26] <- M2_0_cond # Mid:M2
    X2[,27] <- M3_0_condM1M2 # Mid:M3

    X2[,48] <- M11_0_marg*M3_0_condM1M2 # M11:M3
    X2[,57] <- M12_0_marg*M3_0_condM1M2 # M12:M3
    X2[,65] <- M13_0_marg*M3_0_condM1M2 # M13:M3
    X2[,72] <- M14_0_marg*M3_0_condM1M2 # M14:M3
    X2[,78] <- M15_0_marg*M3_0_condM1M2 # M15:M3
    X2[,83] <- M16_0_marg*M3_0_condM1M2 # M16:M3
    X2[,87] <- M17_0_marg*M3_0_condM1M2 # M17:M3

    X2[,92] <- M2_0_cond*M3_0_condM1M2# M2:M3

    M4_0_condM1M2M3 <- rbinom(n2, 1, predict(fit_M4_condM1M2M3, newdata = X2[,c(1:13,15:17,18:27,29:38,40:48,50:57,59:65,67:72,74:78,80:83,85:87,89,92)], type = "response"))

    # For the exposed condition (low SES) #
    
    X2[,2] <- 0
    X2[,4] <- M11_1_marg
    X2[,5] <- M12_1_marg
    X2[,6] <- M13_1_marg
    X2[,7] <- M14_1_marg
    X2[,8] <- M15_1_marg
    X2[,9] <- M16_1_marg
    X2[,10] <- M17_1_marg
    X2[,11] <- M18_1_marg

    X2[,12] <- M2_1_cond

    X2[,13] <- M3_1_condM1M2

    X2[,18] <- 0
    X2[,19] <- 0
    X2[,20] <- 0
    X2[,21] <- 0
    X2[,22] <- 0
    X2[,23] <- 0
    X2[,24] <- 0
    X2[,25] <- 0


    X2[,40] <- M11_1_marg* M12_1_marg # M11:M12
    X2[,41] <- M11_1_marg*M13_1_marg # M11:M13
    X2[,42] <- M11_1_marg*M14_1_marg # M11:M14
    X2[,43] <- M11_1_marg*M15_1_marg # M11:M15
    X2[,44] <- M11_1_marg*M16_1_marg # M11:M16
    X2[,45] <- M11_1_marg*M17_1_marg # M11:M17
    X2[,46] <- M11_1_marg*M18_1_marg# M11:M18
    X2[,50] <- M12_1_marg*M13_1_marg # M12:M13
    X2[,51] <- M12_1_marg*M14_1_marg # M12:M14
    X2[,52] <- M12_1_marg*M15_1_marg # M12:M15
    X2[,53] <- M12_1_marg*M16_1_marg # M12:M16
    X2[,54] <- M12_1_marg*M17_1_marg # M12:M17
    X2[,55] <- M12_1_marg*M18_1_marg# M12:M18
    X2[,59] <- M13_1_marg*M14_1_marg # M13:M14
    X2[,60] <- M13_1_marg*M15_1_marg # M13:M15
    X2[,61] <- M13_1_marg*M16_1_marg # M13:M16
    X2[,62] <- M13_1_marg*M17_1_marg # M13:M17
    X2[,63] <- M13_1_marg*M18_1_marg# M13:M18
    X2[,67] <- M14_1_marg*M15_1_marg # M14:M15
    X2[,68] <- M14_1_marg*M16_1_marg # M14:M16
    X2[,69] <- M14_1_marg*M17_1_marg # M14:M17
    X2[,70] <- M14_1_marg*M18_1_marg# M14:M18
    X2[,74] <- M15_1_marg*M16_1_marg # M15:M16
    X2[,75] <- M15_1_marg*M17_1_marg # M15:M17
    X2[,76] <- M15_1_marg*M18_1_marg# M15:M18
    X2[,80] <- M16_1_marg*M17_1_marg # M16:M17
    X2[,81] <- M16_1_marg*M18_1_marg# M16:M18
    X2[,85] <- M17_1_marg*M18_1_marg# M17:M18

    X2[,26] <- 0 # Mid:M2
    X2[,47] <- M11_1_marg*M2_1_cond # M11:M2
    X2[,56] <- M12_1_marg*M2_1_cond # M12:M2
    X2[,64] <- M13_1_marg*M2_1_cond # M13:M2
    X2[,71] <- M14_1_marg*M2_1_cond # M14:M2
    X2[,77] <- M15_1_marg*M2_1_cond # M15:M2
    X2[,82] <- M16_1_marg*M2_1_cond # M16:M2
    X2[,86] <- M17_1_marg*M2_1_cond # M17:M2
    X2[,89] <- M18_1_marg*M2_1_cond # M18:M2

    X2[,27] <- 0 # Mid:M3
    X2[,48] <- M11_1_marg*M3_1_condM1M2 # M11:M3
    X2[,57] <- M12_1_marg*M3_1_condM1M2 # M12:M3
    X2[,65] <- M13_1_marg*M3_1_condM1M2 # M13:M3
    X2[,72] <- M14_1_marg*M3_1_condM1M2 # M14:M3
    X2[,78] <- M15_1_marg*M3_1_condM1M2 # M15:M3
    X2[,83] <- M16_1_marg*M3_1_condM1M2 # M16:M3
    X2[,87] <- M17_1_marg*M3_1_condM1M2 # M17:M3
    
    X2[,92] <- M2_1_cond*M3_1_condM1M2# M2:M3

    M4_1_condM1M2M3 <- rbinom(n2, 1, predict(fit_M4_condM1M2M3, newdata = X2[,c(1:13,15:17,18:27,29:38,40:48,50:57,59:65,67:72,74:78,80:83,85:87,89,92)], type = "response"))


    ### Outcome ###

    ## Estimation of conditional expectations ###


    # With joint distribution for M1 and M2 #
    
    # Exposed condition (low SES) for exposure, M1, M2 and M3, unexposed condition for M4 #
    X2[,13] <- M3_1_marg
    X2[,27] <- 0 # Mid:M3
    X2[,48] <- M11_1_marg*M3_1_marg # M11:M3
    X2[,57] <- M12_1_marg*M3_1_marg # M12:M3
    X2[,65] <- M13_1_marg*M3_1_marg # M13:M3
    X2[,72] <- M14_1_marg*M3_1_marg # M14:M3
    X2[,78] <- M15_1_marg*M3_1_marg # M15:M3
    X2[,83] <- M16_1_marg*M3_1_marg # M16:M3
    X2[,87] <- M17_1_marg*M3_1_marg # M17:M3
    X2[,90] <- M18_1_marg*M3_1_marg # M18:M3

    X2[,92] <- M2_1_cond*M3_1_marg # M2:M3

    X2[,14] <- M4_0_marg #
    X2[,28] <- 0 # Mid:M4
    X2[,49] <- M11_1_marg*M4_0_marg # M11:M4
    X2[,58] <- M12_1_marg*M4_0_marg  # M12:M4
    X2[,66] <- M13_1_marg*M4_0_marg  # M13:M4
    X2[,73] <- M14_1_marg*M4_0_marg  # M14:M4
    X2[,79] <- M15_1_marg*M4_0_marg  # M15:M4
    X2[,84] <- M16_1_marg*M4_0_marg  # M16:M4
    X2[,88] <- M17_1_marg*M4_0_marg  # M17:M4
    X2[,91] <- M18_1_marg*M4_0_marg  # M18:M4

    X2[,93] <- M2_1_cond*M4_0_marg # M2:M4
    X2[,94] <- M3_1_marg*M4_0_marg # M3:M4


    y1_m121_m31_m40[i] <- mean(predict(fit_Y, newdata = X2, type = "response")) 

    
    # Exposed condition (low SES) for exposure, M1, and M2, unexposed condition for M3 and M4 #
    
    X2[,13] <- M3_0_marg #

    X2[,27] <- 0 # Mid:M3
    X2[,48] <- M11_1_marg*M3_0_marg # M11:M3
    X2[,57] <- M12_1_marg*M3_0_marg # M12:M3
    X2[,65] <- M13_1_marg*M3_0_marg # M13:M3
    X2[,72] <- M14_1_marg*M3_0_marg # M14:M3
    X2[,78] <- M15_1_marg*M3_0_marg # M15:M3
    X2[,83] <- M16_1_marg*M3_0_marg # M16:M3
    X2[,87] <- M17_1_marg*M3_0_marg # M17:M3
    X2[,90] <- M18_1_marg*M3_0_marg # M18:M3

    X2[,92] <- M2_1_cond*M3_0_marg # M2:M3

    X2[,94] <- M3_0_marg*M4_0_marg # M3:M4

    y1_m121_m30_m40[i] <- mean(predict(fit_Y, newdata = X2, type = "response")) 

    
    # With joint distribution for M1, M2 and M3 #
    
    # Exposed condition (low SES) for exposure, M1, M2 and M3, unexposed condition for M4 #
    
    X2[,13] <- M3_1_condM1M2 #

    X2[,27] <- 0 # Mid:M3
    X2[,48] <- M11_1_marg*M3_1_condM1M2 # M11:M3
    X2[,57] <- M12_1_marg*M3_1_condM1M2 # M12:M3
    X2[,65] <- M13_1_marg*M3_1_condM1M2 # M13:M3
    X2[,72] <- M14_1_marg*M3_1_condM1M2 # M14:M3
    X2[,78] <- M15_1_marg*M3_1_condM1M2 # M15:M3
    X2[,83] <- M16_1_marg*M3_1_condM1M2 # M16:M3
    X2[,87] <- M17_1_marg*M3_1_condM1M2 # M17:M3
    X2[,90] <- M18_1_marg*M3_1_condM1M2 # M18:M3

    X2[,92] <- M2_1_cond*M3_1_condM1M2 # M2:M3

    X2[,94] <- M3_1_condM1M2*M4_0_marg # M3:M4


    y1_m1231_m40[i] <- mean(predict(fit_Y, newdata = X2, type = "response")) #

    
    # Exposed condition (low SES) for exposure, M1, M2, M3, and M4 #
    
    X2[,14] <- M4_1_marg #

    X2[,28] <- 0 # Mid:M4
    X2[,49] <- M11_1_marg*M4_1_marg # M11:M4
    X2[,58] <- M12_1_marg*M4_1_marg  # M12:M4
    X2[,66] <- M13_1_marg*M4_1_marg  # M13:M4
    X2[,73] <- M14_1_marg*M4_1_marg  # M14:M4
    X2[,79] <- M15_1_marg*M4_1_marg  # M15:M4
    X2[,84] <- M16_1_marg*M4_1_marg  # M16:M4
    X2[,88] <- M17_1_marg*M4_1_marg  # M17:M4
    X2[,91] <- M18_1_marg*M4_1_marg  # M18:M4

    X2[,93] <- M2_1_cond*M4_1_marg # M2:M4

    X2[,94] <- M3_1_condM1M2*M4_1_marg # M3:M4

    y1_m1231_m41[i] <- mean(predict(fit_Y, newdata = X2, type = "response")) #

    
    # With joint distribution for M1, M2, M3 and M4 #
    
    # Exposed condition (low SES) for exposure, M1, M2, M3, and M4 #
    
    X2[,14] <- M4_1_condM1M2M3

    X2[,28] <- 0 # Mid:M4
    X2[,49] <- M11_1_marg*M4_1_condM1M2M3 # M11:M4
    X2[,58] <- M12_1_marg*M4_1_condM1M2M3  # M12:M4
    X2[,66] <- M13_1_marg*M4_1_condM1M2M3  # M13:M4
    X2[,73] <- M14_1_marg*M4_1_condM1M2M3  # M14:M4
    X2[,79] <- M15_1_marg*M4_1_condM1M2M3  # M15:M4
    X2[,84] <- M16_1_marg*M4_1_condM1M2M3  # M16:M4
    X2[,88] <- M17_1_marg*M4_1_condM1M2M3  # M17:M4
    X2[,91] <- M18_1_marg*M4_1_condM1M2M3  # M18:M4

    X2[,93] <- M2_1_cond*M4_1_condM1M2M3 # M2:M4

    X2[,94] <- M3_1_condM1M2*M4_1_condM1M2M3 # M3:M4

    y1_m12341[i] <- mean(predict(fit_Y, newdata = X2, type = "response")) #


    # With joint distribution for M3 and M4 #
    
    # Exposed condition (low SES) for exposure, M1, and M2, unexposed condition for M3 and M4 #
    
    X2[,12] <- M2_1_marg
    X2[,26] <- 0 # Mid:M2
    X2[,47] <- M11_1_marg*M2_1_marg # M11:M2
    X2[,56] <- M12_1_marg*M2_1_marg # M12:M2
    X2[,64] <- M13_1_marg*M2_1_marg # M13:M2
    X2[,71] <- M14_1_marg*M2_1_marg # M14:M2
    X2[,77] <- M15_1_marg*M2_1_marg # M15:M2
    X2[,82] <- M16_1_marg*M2_1_marg # M16:M2
    X2[,86] <- M17_1_marg*M2_1_marg # M17:M2
    X2[,89] <- M18_1_marg*M2_1_marg # M18:M2

    X2[,13] <- M3_0_marg

    X2[,27] <- 0 # Mid:M3
    X2[,48] <- M11_1_marg*M3_0_marg # M11:M3
    X2[,57] <- M12_1_marg*M3_0_marg # M12:M3
    X2[,65] <- M13_1_marg*M3_0_marg # M13:M3
    X2[,72] <- M14_1_marg*M3_0_marg # M14:M3
    X2[,78] <- M15_1_marg*M3_0_marg # M15:M3
    X2[,83] <- M16_1_marg*M3_0_marg # M16:M3
    X2[,87] <- M17_1_marg*M3_0_marg # M17:M3
    X2[,90] <- M18_1_marg*M3_0_marg # M18:M3

    X2[,92] <- M2_1_marg*M3_0_marg # M2:M3

    X2[,14] <- M4_0_condM3

    X2[,28] <- 0 # Mid:M4
    X2[,49] <- M11_1_marg*M4_0_condM3 # M11:M4
    X2[,58] <- M12_1_marg*M4_0_condM3  # M12:M4
    X2[,66] <- M13_1_marg*M4_0_condM3  # M13:M4
    X2[,73] <- M14_1_marg*M4_0_condM3  # M14:M4
    X2[,79] <- M15_1_marg*M4_0_condM3  # M15:M4
    X2[,84] <- M16_1_marg*M4_0_condM3  # M16:M4
    X2[,88] <- M17_1_marg*M4_0_condM3  # M17:M4
    X2[,91] <- M18_1_marg*M4_0_condM3  # M18:M4

    X2[,93] <- M2_1_marg*M4_0_condM3 # M2:M4

    X2[,94] <- M3_0_marg*M4_0_condM3 # M3:M4

    y1_m11_m21_m340[i] <- mean(predict(fit_Y, newdata = X2, type = "response")) #

    
    # Exposed condition (low SES) for exposure, M1, unexposed condition for M2, M3 and M4 #
    
    X2[,12] <- M2_0_marg

    X2[,26] <- 0 # Mid:M2
    X2[,47] <- M11_1_marg*M2_0_marg # M11:M2
    X2[,56] <- M12_1_marg*M2_0_marg # M12:M2
    X2[,64] <- M13_1_marg*M2_0_marg # M13:M2
    X2[,71] <- M14_1_marg*M2_0_marg # M14:M2
    X2[,77] <- M15_1_marg*M2_0_marg # M15:M2
    X2[,82] <- M16_1_marg*M2_0_marg # M16:M2
    X2[,86] <- M17_1_marg*M2_0_marg # M17:M2
    X2[,89] <- M18_1_marg*M2_0_marg # M18:M2

    X2[,92] <- M2_0_marg*M3_0_marg # M2:M3

    X2[,93] <- M2_0_marg*M4_0_condM3 # M2:M4

    y1_m11_m20_m340[i] <- mean(predict(fit_Y, newdata = X2, type = "response")) #

    
    # With joint distribution for M2, M3 and M4 #
    
    # Exposed condition (low SES) for exposure, M1, unexposed condition for M2, M3 and M4 #
    
    X2[,13] <- M3_0_condM2
    X2[,14] <- M4_0_condM2M3

    X2[,27] <- 0 # Mid:M3
    X2[,48] <- M11_1_marg*M3_0_condM2 # M11:M3
    X2[,57] <- M12_1_marg*M3_0_condM2 # M12:M3
    X2[,65] <- M13_1_marg*M3_0_condM2 # M13:M3
    X2[,72] <- M14_1_marg*M3_0_condM2 # M14:M3
    X2[,78] <- M15_1_marg*M3_0_condM2 # M15:M3
    X2[,83] <- M16_1_marg*M3_0_condM2 # M16:M3
    X2[,87] <- M17_1_marg*M3_0_condM2 # M17:M3
    X2[,90] <- M18_1_marg*M3_0_condM2 # M18:M3
    X2[,92] <- M2_0_marg*M3_0_condM2 # M2:M3

    X2[,28] <- 0 # Mid:M4
    X2[,49] <- M11_1_marg*M4_0_condM2M3 # M11:M4
    X2[,58] <- M12_1_marg*M4_0_condM2M3  # M12:M4
    X2[,66] <- M13_1_marg*M4_0_condM2M3  # M13:M4
    X2[,73] <- M14_1_marg*M4_0_condM2M3  # M14:M4
    X2[,79] <- M15_1_marg*M4_0_condM2M3  # M15:M4
    X2[,84] <- M16_1_marg*M4_0_condM2M3  # M16:M4
    X2[,88] <- M17_1_marg*M4_0_condM2M3  # M17:M4
    X2[,91] <- M18_1_marg*M4_0_condM2M3  # M18:M4

    X2[,93] <- M2_0_marg*M4_0_condM2M3 # M2:M4

    X2[,94] <- M3_0_condM2*M4_0_condM2M3 # M3:M4

    y1_m11_m2340[i] <- mean(predict(fit_Y, newdata = X2, type = "response")) #

    
    # Exposed condition (low SES) for exposure, unexposed condition for M1, M2, M3 and M4 #
    
    X2[,4] <- M11_0_marg
    X2[,5] <- M12_0_marg
    X2[,6] <- M13_0_marg
    X2[,7] <- M14_0_marg
    X2[,8] <- M15_0_marg
    X2[,9] <- M16_0_marg
    X2[,10] <- M17_0_marg
    X2[,11] <- M18_0_marg

    X2[,18] <- 0 # Mid:M11
    X2[,19] <- 0
    X2[,20] <- 0
    X2[,21] <- 0
    X2[,22] <- 0
    X2[,23] <- 0
    X2[,24] <- 0
    X2[,25] <- 0
    
    X2[,40] <- M11_0_marg* M12_0_marg # M11:M12
    X2[,41] <- M11_0_marg*M13_0_marg # M11:M13
    X2[,42] <- M11_0_marg*M14_0_marg # M11:M14
    X2[,43] <- M11_0_marg*M15_0_marg # M11:M15
    X2[,44] <- M11_0_marg*M16_0_marg # M11:M16
    X2[,45] <- M11_0_marg*M17_0_marg # M11:M17
    X2[,46] <- M11_0_marg*M18_0_marg# M11:M18
    X2[,50] <- M12_0_marg*M13_0_marg # M12:M13
    X2[,51] <- M12_0_marg*M14_0_marg # M12:M14
    X2[,52] <- M12_0_marg*M15_0_marg # M12:M15
    X2[,53] <- M12_0_marg*M16_0_marg # M12:M16
    X2[,54] <- M12_0_marg*M17_0_marg # M12:M17
    X2[,55] <- M12_0_marg*M18_0_marg# M12:M18
    X2[,59] <- M13_0_marg*M14_0_marg # M13:M14
    X2[,60] <- M13_0_marg*M15_0_marg # M13:M15
    X2[,61] <- M13_0_marg*M16_0_marg # M13:M16
    X2[,62] <- M13_0_marg*M17_0_marg # M13:M17
    X2[,63] <- M13_0_marg*M18_0_marg# M13:M18
    X2[,67] <- M14_0_marg*M15_0_marg # M14:M15
    X2[,68] <- M14_0_marg*M16_0_marg # M14:M16
    X2[,69] <- M14_0_marg*M17_0_marg # M14:M17
    X2[,70] <- M14_0_marg*M18_0_marg# M14:M18
    X2[,74] <- M15_0_marg*M16_0_marg # M15:M16
    X2[,75] <- M15_0_marg*M17_0_marg # M15:M17
    X2[,76] <- M15_0_marg*M18_0_marg# M15:M18
    X2[,80] <- M16_0_marg*M17_0_marg # M16:M17
    X2[,81] <- M16_0_marg*M18_0_marg# M16:M18
    X2[,85] <- M17_0_marg*M18_0_marg# M17:M18
    
    X2[,47] <- M11_0_marg*M2_0_marg # M11:M2
    X2[,56] <- M12_0_marg*M2_0_marg # M12:M2
    X2[,64] <- M13_0_marg*M2_0_marg # M13:M2
    X2[,71] <- M14_0_marg*M2_0_marg # M14:M2
    X2[,77] <- M15_0_marg*M2_0_marg # M15:M2
    X2[,82] <- M16_0_marg*M2_0_marg # M16:M2
    X2[,86] <- M17_0_marg*M2_0_marg # M17:M2
    X2[,89] <- M18_0_marg*M2_0_marg # M18:M2
   
    X2[,48] <- M11_0_marg*M3_0_condM2 # M11:M3
    X2[,57] <- M12_0_marg*M3_0_condM2 # M12:M3
    X2[,65] <- M13_0_marg*M3_0_condM2 # M13:M3
    X2[,72] <- M14_0_marg*M3_0_condM2 # M14:M3
    X2[,78] <- M15_0_marg*M3_0_condM2 # M15:M3
    X2[,83] <- M16_0_marg*M3_0_condM2 # M16:M3
    X2[,87] <- M17_0_marg*M3_0_condM2 # M17:M3
    X2[,90] <- M18_0_marg*M3_0_condM2 # M17:M3
    
    X2[,49] <- M11_0_marg*M4_0_condM2M3 # M11:M4
    X2[,58] <- M12_0_marg*M4_0_condM2M3  # M12:M4
    X2[,66] <- M13_0_marg*M4_0_condM2M3  # M13:M4
    X2[,73] <- M14_0_marg*M4_0_condM2M3  # M14:M4
    X2[,79] <- M15_0_marg*M4_0_condM2M3  # M15:M4
    X2[,84] <- M16_0_marg*M4_0_condM2M3  # M16:M4
    X2[,88] <- M17_0_marg*M4_0_condM2M3  # M17:M4
    X2[,91] <- M18_0_marg*M4_0_condM2M3  # M18:M4

    y1_m10_m2340[i] <- mean(predict(fit_Y, newdata = X2, type = "response")) #

    
    # Exposed condition (low SES) for exposure, unexposed condition for M1, M2, M3 and M4 #
    
    X2[,12] <- M2_0_cond

    X2[,26] <- 0 # Mid:M2
    X2[,47] <- M11_0_marg*M2_0_cond # M11:M2
    X2[,56] <- M12_0_marg*M2_0_cond # M12:M2
    X2[,64] <- M13_0_marg*M2_0_cond # M13:M2
    X2[,71] <- M14_0_marg*M2_0_cond # M14:M2
    X2[,77] <- M15_0_marg*M2_0_cond # M15:M2
    X2[,82] <- M16_0_marg*M2_0_cond # M16:M2
    X2[,86] <- M17_0_marg*M2_0_cond # M17:M2
    X2[,89] <- M18_0_marg*M2_0_cond # M18:M2

    X2[,92] <- M2_0_cond*M3_0_condM1M2 # M2:M3

    X2[,93] <- M2_0_cond*M4_0_condM1M2M3 # M2:M4

    X2[,13] <- M3_0_condM1M2 #

    X2[,27] <- 0 # Mid:M3
    
    X2[,48] <- M11_0_marg*M3_0_condM1M2 # M11:M3
    X2[,57] <- M12_0_marg*M3_0_condM1M2 # M12:M3
    X2[,65] <- M13_0_marg*M3_0_condM1M2 # M13:M3
    X2[,72] <- M14_0_marg*M3_0_condM1M2 # M14:M3
    X2[,78] <- M15_0_marg*M3_0_condM1M2 # M15:M3
    X2[,83] <- M16_0_marg*M3_0_condM1M2 # M16:M3
    X2[,87] <- M17_0_marg*M3_0_condM1M2 # M17:M3
    X2[,90] <- M18_0_marg*M3_0_condM1M2 # M17:M3
    
    X2[,94] <- M3_0_condM1M2*M4_0_condM1M2M3 # M3:M4

    X2[,14] <- M4_0_condM1M2M3 #

    X2[,28] <- 0 # Mid:M4
    
    X2[,49] <- M11_0_marg*M4_0_condM1M2M3 # M11:M4
    X2[,58] <- M12_0_marg*M4_0_condM1M2M3  # M12:M4
    X2[,66] <- M13_0_marg*M4_0_condM1M2M3  # M13:M4
    X2[,73] <- M14_0_marg*M4_0_condM1M2M3  # M14:M4
    X2[,79] <- M15_0_marg*M4_0_condM1M2M3  # M15:M4
    X2[,84] <- M16_0_marg*M4_0_condM1M2M3  # M16:M4
    X2[,88] <- M17_0_marg*M4_0_condM1M2M3  # M17:M4
    X2[,91] <- M18_0_marg*M4_0_condM1M2M3  # M18:M4

    y1_m12340[i] <- mean(predict(fit_Y, newdata = X2, type = "response")) #

    # Unexposed condition (mid SES) for exposure, M1, M2, M3 and M4 #
    
    X2[,2] <- 1
    
    X2[,18] <- M11_0_marg # Mid:M11
    X2[,19] <- M12_0_marg
    X2[,20] <- M13_0_marg
    X2[,21] <- M14_0_marg
    X2[,22] <- M15_0_marg
    X2[,23] <- M16_0_marg
    X2[,24] <- M17_0_marg
    X2[,25] <- M18_0_marg # Mid:M18

    X2[,26] <- M2_0_cond # Mid:M2

    X2[,27] <- M3_0_condM1M2 # Mid:M3

    X2[,28] <- M4_0_condM1M2M3 # Mid:M4

    y0_m12340[i] <- mean(predict(fit_Y, newdata = X2, type = "response")) #

  }
  
  ### Estimation of the interventional disparity effects ###

  ## Adjusted total association ##

  Adj_TA <- mean(y1_m12341 - y0_m12340)


  ## Direct effect

  IDM_DE <- mean(y1_m12340 - y0_m12340)


  ## Indirect effects

  # Via the M1 group #
  IDM_IE1 <- mean(y1_m11_m2340 - y1_m10_m2340)

  # Via M2 (NIHSS>5) #
  IDM_IE2 <- mean(y1_m11_m21_m340 - y1_m11_m20_m340)

  # Via M3 (Reperfusion) #
  IDM_IE3 <- mean(y1_m121_m31_m40 - y1_m121_m30_m40)

  # Via M4 (Stroke unit) #
  IDM_IE4 <- mean(y1_m1231_m41 - y1_m1231_m40)

  # Via the interdependence between the mediators #
  IDM_MD <- Adj_TA - (IDM_DE + IDM_IE1 + IDM_IE2 + IDM_IE3 + IDM_IE4)

  # Via all mediators taken together #
  IDM_IE <- mean(y1_m12341 - y1_m12340)


  ### Collect results

  resMid <- c(Adj_TA, IDM_DE, IDM_IE, IDM_IE1, IDM_IE2, IDM_IE3, IDM_IE4, IDM_MD)



 ### Loop for Monte Carlo simulations for the low vs. high contrast ###
  
  for(i in 1:mcsim){ 
    
    ### Simulate from mediator distributions ###
    
    ## Marginal distributions ##
    
    # For the unexposed condition (high SES) #
    
    X3[,3] <- 1 # Set the unexposed condition
    M11_0_marg <- rbinom(n3, 1, M11_preds0_marg3)
    X3[,4] <- M11_0_marg # Update the M11 variable with simulated values
    X3[,29] <- M11_0_marg # Update the High:M11-interaction with simulated values
    
    M12_0_marg <- rbinom(n3, 1, predict(fit_M12_marg, newdata = X3[,c(1:4,15:17,18,29)], type = "response"))
    X3[,5] <- M12_0_marg
    X3[,30] <- M12_0_marg # High:M12
    X3[,40] <- M11_0_marg* M12_0_marg # M11:M12
    
    M13_0_marg <- rbinom(n3, 1, predict(fit_M13_marg, newdata = X3[,c(1:5,15:17,18:19,29:30,40)], type = "response"))
    X3[,6] <- M13_0_marg
    X3[,31] <- M13_0_marg # High:M13
    X3[,41] <- M11_0_marg*M13_0_marg # M11:M13
    X3[,50] <- M12_0_marg*M13_0_marg # M12:M13
    
    M14_0_marg <- rbinom(n3, 1, predict(fit_M14_marg, newdata = X3[,c(1:6,15:17,18:20,29:31,40:41, 50)], type = "response"))
    X3[,7] <- M14_0_marg
    X3[,32] <- M14_0_marg # High:M14
    X3[,42] <- M11_0_marg*M14_0_marg # M11:M14
    X3[,51] <- M12_0_marg*M14_0_marg # M12:M14
    X3[,59] <- M13_0_marg*M14_0_marg # M13:M14
    
    M15_0_marg <- rbinom(n3, 1, predict(fit_M15_marg, newdata = X3[,c(1:7,15:17,18:21,29:32,40:42, 50:51,59)], type = "response"))
    X3[,8] <- M15_0_marg
    X3[,33] <- M15_0_marg # High:M15
    X3[,43] <- M11_0_marg*M15_0_marg # M11:M15
    X3[,52] <- M12_0_marg*M15_0_marg # M12:M15
    X3[,60] <- M13_0_marg*M15_0_marg #M13:M15
    X3[,67] <- M14_0_marg*M15_0_marg # M14:M15
    
    M16_0_marg <- rbinom(n3, 1, predict(fit_M16_marg, newdata = X3[,c(1:8,15:17,18:22,29:33,40:43, 50:52,59:60,67)], type = "response"))
    X3[,9] <- M16_0_marg
    X3[,34] <- M16_0_marg # High:M16
    X3[,44] <- M11_0_marg*M16_0_marg # M11:M16
    X3[,53] <- M12_0_marg*M16_0_marg # M12:M16
    X3[,61] <- M13_0_marg*M16_0_marg # M13:M16
    X3[,68] <- M14_0_marg*M16_0_marg # M14:M16
    X3[,74] <- M15_0_marg*M16_0_marg # M15:M16
    
    M17_0_marg <- rbinom(n3, 1, predict(fit_M17_marg, newdata = X3[,c(1:9,15:17,18:23,29:34,40:44, 50:53,59:61,67:68,74)], type = "response"))
    X3[,10] <- M17_0_marg
    X3[,35] <- M17_0_marg # High: M17
    X3[,45] <- M11_0_marg*M17_0_marg # M11:M17
    X3[,54] <- M12_0_marg*M17_0_marg # M12:M17
    X3[,62] <- M13_0_marg*M17_0_marg # M13:M17
    X3[,69] <- M14_0_marg*M17_0_marg # M14:M17
    X3[,75] <- M15_0_marg*M17_0_marg # M15:M17
    X3[,80] <- M16_0_marg*M17_0_marg # M16:M17
    
    M18_0_marg <- rbinom(n3, 1, predict(fit_M18_marg, newdata = X3[,c(1:10,15:17,18:24,29:35,40:45, 50:54,59:62,67:69,74:75,80)], type = "response"))
    
    M2_0_marg <- rbinom(n3, 1, M2_preds0_marg3)
    M3_0_marg <- rbinom(n3, 1, M3_preds0_marg3)
    M4_0_marg <- rbinom(n3, 1, M4_preds0_marg3)
    
    # For the exposed condition (low SES) #
    
    X3[,3]<-0 # Set the exposed condition
    
    M11_1_marg <- rbinom(n3, 1, M11_preds1_marg3)
    X3[,4] <- M11_1_marg # Update the M11 variable with simulated values
    X3[,29] <- 0 # Update the High:M11-interaction with simulated values
    
    M12_1_marg <- rbinom(n3, 1, predict(fit_M12_marg, newdata = X3[,c(1:4,15:17,18,29)], type = "response"))
    X3[,5] <- M12_1_marg
    X3[,30] <- 0 # High:M12
    X3[,40] <- M11_1_marg* M12_1_marg # M11:M12
    
    M13_1_marg <- rbinom(n3, 1, predict(fit_M13_marg, newdata = X3[,c(1:5,15:17,18:19,29:30,40)], type = "response"))
    X3[,6] <- M13_1_marg
    X3[,31] <- 0 # High:M13
    X3[,41] <- M11_1_marg*M13_1_marg # M11:M13
    X3[,50] <- M12_1_marg*M13_1_marg # M12:M13
    
    M14_1_marg <- rbinom(n3, 1, predict(fit_M14_marg, newdata = X3[,c(1:6,15:17,18:20,29:31,40:41, 50)], type = "response"))
    X3[,7] <- M14_1_marg
    X3[,32] <- 0 # High:M14
    X3[,42] <- M11_1_marg*M14_1_marg # M11:M14
    X3[,51] <- M12_1_marg*M14_1_marg # M12:M14
    X3[,59] <- M13_1_marg*M14_1_marg # M13:M14
    
    M15_1_marg <- rbinom(n3, 1, predict(fit_M15_marg, newdata = X3[,c(1:7,15:17,18:21,29:32,40:42, 50:51,59)], type = "response"))
    X3[,8] <- M15_1_marg
    X3[,33] <- 0 # High:M15
    X3[,43] <- M11_1_marg*M15_1_marg # M11:M15
    X3[,52] <- M12_1_marg*M15_1_marg # M12:M15
    X3[,60] <- M13_1_marg*M15_1_marg #M13:M15
    X3[,67] <- M14_1_marg*M15_1_marg # M14:M15
    
    M16_1_marg <- rbinom(n3, 1, predict(fit_M16_marg, newdata = X3[,c(1:8,15:17,18:22,29:33,40:43, 50:52,59:60,67)], type = "response"))
    X3[,9] <- M16_1_marg
    X3[,34] <- 0 # High:M16
    X3[,44] <- M11_1_marg*M16_1_marg # M11:M16
    X3[,53] <- M12_1_marg*M16_1_marg # M12:M16
    X3[,61] <- M13_1_marg*M16_1_marg # M13:M16
    X3[,68] <- M14_1_marg*M16_1_marg # M14:M16
    X3[,74] <- M15_1_marg*M16_1_marg # M15:M16
    
    M17_1_marg <- rbinom(n3, 1, predict(fit_M17_marg, newdata = X3[,c(1:9,15:17,18:23,29:34,40:44, 50:53,59:61,67:68,74)], type = "response"))
    X3[,10] <- M17_1_marg
    X3[,35] <- 0 # High: M17
    X3[,45] <- M11_1_marg*M17_1_marg # M11:M17
    X3[,54] <- M12_1_marg*M17_1_marg # M12:M17
    X3[,62] <- M13_1_marg*M17_1_marg # M13:M17
    X3[,69] <- M14_1_marg*M17_1_marg # M14:M17
    X3[,75] <- M15_1_marg*M17_1_marg # M15:M17
    X3[,80] <- M16_1_marg*M17_1_marg # M16:M17
    
    M18_1_marg <- rbinom(n3, 1, predict(fit_M18_marg, newdata = X3[,c(1:10,15:17,18:24,29:35,40:45, 50:54,59:62,67:69,74:75,80)], type = "response"))
    
    M2_1_marg <- rbinom(n3, 1, M2_preds1_marg3)
    M3_1_marg <- rbinom(n3, 1, M3_preds1_marg3)
    M4_1_marg <- rbinom(n3, 1, M4_preds1_marg3)
    
    
    # M2 - given the M1 group #
    
    # For the unexposed condition (high SES) #
    
    X3[,3] <- 1
    X3[,4] <- M11_0_marg
    X3[,5] <- M12_0_marg
    X3[,6] <- M13_0_marg
    X3[,7] <- M14_0_marg
    X3[,8] <- M15_0_marg
    X3[,9] <- M16_0_marg
    X3[,10] <- M17_0_marg
    X3[,11] <- M18_0_marg
    
    X3[,29] <- M11_0_marg # High:M11
    X3[,30] <- M12_0_marg # High:M12
    X3[,31] <- M13_0_marg # High:M13
    X3[,32] <- M14_0_marg # High:M14
    X3[,33] <- M15_0_marg # High:M15
    X3[,34] <- M16_0_marg # High:M16
    X3[,35] <- M17_0_marg # High:M17
    X3[,36] <- M18_0_marg # High:M18
    
    X3[,40] <- M11_0_marg* M12_0_marg # M11:M12
    X3[,41] <- M11_0_marg*M13_0_marg # M11:M13
    X3[,42] <- M11_0_marg*M14_0_marg # M11:M14
    X3[,43] <- M11_0_marg*M15_0_marg # M11:M15
    X3[,44] <- M11_0_marg*M16_0_marg # M11:M16
    X3[,45] <- M11_0_marg*M17_0_marg # M11:M17
    X3[,46] <- M11_0_marg*M18_0_marg # M11:M18
    X3[,50] <- M12_0_marg*M13_0_marg # M12:M13
    X3[,51] <- M12_0_marg*M14_0_marg # M12:M14
    X3[,52] <- M12_0_marg*M15_0_marg # M12:M15
    X3[,53] <- M12_0_marg*M16_0_marg # M12:M16
    X3[,54] <- M12_0_marg*M17_0_marg # M12:M17
    X3[,55] <- M12_0_marg*M18_0_marg # M12:M18
    X3[,59] <- M13_0_marg*M14_0_marg # M13:M14
    X3[,60] <- M13_0_marg*M15_0_marg # M13:M15
    X3[,61] <- M13_0_marg*M16_0_marg # M13:M16
    X3[,62] <- M13_0_marg*M17_0_marg # M13:M17
    X3[,63] <- M13_0_marg*M18_0_marg # M13:M18
    X3[,67] <- M14_0_marg*M15_0_marg # M14:M15
    X3[,68] <- M14_0_marg*M16_0_marg # M14:M16
    X3[,69] <- M14_0_marg*M17_0_marg # M14:M17
    X3[,70] <- M14_0_marg*M18_0_marg # M14:M18
    X3[,74] <- M15_0_marg*M16_0_marg # M15:M16
    X3[,75] <- M15_0_marg*M17_0_marg # M15:M17
    X3[,76] <- M15_0_marg*M18_0_marg # M15:M18
    X3[,80] <- M16_0_marg*M17_0_marg # M16:M17
    X3[,81] <- M16_0_marg*M18_0_marg # M16:M18
    X3[,85] <- M17_0_marg*M18_0_marg # M17:M18
    
    M2_0_cond <- rbinom(n3, 1, predict(fit_M2_cond, newdata = X3[,c(1:11,15:17,18:25,29:36,40:46,50:55,59:63,67:70,74:76,80:81,85)], type = "response"))
    
    # For the exposed condition (low SES) #
    
    X3[,3] <- 0
    X3[,4] <- M11_1_marg
    X3[,5] <- M12_1_marg
    X3[,6] <- M13_1_marg
    X3[,7] <- M14_1_marg
    X3[,8] <- M15_1_marg
    X3[,9] <- M16_1_marg
    X3[,10] <- M17_1_marg
    X3[,11] <- M18_1_marg
    X3[,29] <- 0 # High:M11
    X3[,30] <- 0
    X3[,31] <- 0
    X3[,32] <- 0
    X3[,33] <- 0
    X3[,34] <- 0
    X3[,35] <- 0
    X3[,36] <- 0 # High:M18
    
    X3[,40] <- M11_1_marg* M12_1_marg # M11:M12
    X3[,41] <- M11_1_marg*M13_1_marg # M11:M13
    X3[,42] <- M11_1_marg*M14_1_marg # M11:M14
    X3[,43] <- M11_1_marg*M15_1_marg # M11:M15
    X3[,44] <- M11_1_marg*M16_1_marg # M11:M16
    X3[,45] <- M11_1_marg*M17_1_marg # M11:M17
    X3[,46] <- M11_1_marg*M18_1_marg# M11:M18
    X3[,50] <- M12_1_marg*M13_1_marg # M12:M13
    X3[,51] <- M12_1_marg*M14_1_marg # M12:M14
    X3[,52] <- M12_1_marg*M15_1_marg # M12:M15
    X3[,53] <- M12_1_marg*M16_1_marg # M12:M16
    X3[,54] <- M12_1_marg*M17_1_marg # M12:M17
    X3[,55] <- M12_1_marg*M18_1_marg# M12:M18
    X3[,59] <- M13_1_marg*M14_1_marg # M13:M14
    X3[,60] <- M13_1_marg*M15_1_marg # M13:M15
    X3[,61] <- M13_1_marg*M16_1_marg # M13:M16
    X3[,62] <- M13_1_marg*M17_1_marg # M13:M17
    X3[,63] <- M13_1_marg*M18_1_marg# M13:M18
    X3[,67] <- M14_1_marg*M15_1_marg # M14:M15
    X3[,68] <- M14_1_marg*M16_1_marg # M14:M16
    X3[,69] <- M14_1_marg*M17_1_marg # M14:M17
    X3[,70] <- M14_1_marg*M18_1_marg# M14:M18
    X3[,74] <- M15_1_marg*M16_1_marg # M15:M16
    X3[,75] <- M15_1_marg*M17_1_marg # M15:M17
    X3[,76] <- M15_1_marg*M18_1_marg# M15:M18
    X3[,80] <- M16_1_marg*M17_1_marg # M16:M17
    X3[,81] <- M16_1_marg*M18_1_marg# M16:M18
    X3[,85] <- M17_1_marg*M18_1_marg# M17:M18
    
    M2_1_cond <- rbinom(n3, 1, predict(fit_M2_cond, newdata = X3[,c(1:11,15:17,18:25,29:36,40:46,50:55,59:63,67:70,74:76,80:81,85)], type = "response"))
    
    
    # M3 - given M2 #
    
    # For the unexposed condition (high SES) #
    
    X3[,3] <- 1
    X3[,12] <- M2_0_marg
    X3[,37] <- M2_0_marg # High:M2
    
    M3_0_condM2 <- rbinom(n3, 1, predict(fit_M3_condM2, newdata = X3[,c(1:3,12,15:17,26,37)], type = "response"))
    
    # M4 - given M3 #
    
    # For the unexposed condition (high SES) #
    
    X3[,13] <- M3_0_marg
    X3[,38] <- M3_0_marg # High:M3
    M4_0_condM3 <- rbinom(n3, 1, predict(fit_M4_condM3, newdata = X3[,c(1:3,13,15:17,27,38)], type = "response"))
    
    
    # M3 - given the M1 group and M2
    
    # For the exposed condition (low SES) #
    
    X3[,3] <- 0
    
    X3[,12] <- M2_1_cond
    X3[,37] <- 0 # High:M2
    
    X3[,47] <- M11_1_marg*M2_1_cond # M11:M2
    X3[,56] <- M12_1_marg*M2_1_cond # M12:M2
    X3[,64] <- M13_1_marg*M2_1_cond # M13:M2
    X3[,71] <- M14_1_marg*M2_1_cond # M14:M2
    X3[,77] <- M15_1_marg*M2_1_cond # M15:M2
    X3[,82] <- M16_1_marg*M2_1_cond # M16:M2
    X3[,86] <- M17_1_marg*M2_1_cond # M17:M2
    X3[,89] <- M18_1_marg*M2_1_cond # M18:M2
    
    M3_1_condM1M2 <- rbinom(n3, 1, predict(fit_M3_condM1M2, newdata = X3[,c(1:11,12,15:17,18:26,29:37,40:47,50:56,59:64,67:71,74:77,80:82,85:86,89)], type = "response"))
    
    # For the unexposed condition (high SES) #
    
    X3[,3] <- 1
    
    X3[,4] <- M11_0_marg
    X3[,5] <- M12_0_marg
    X3[,6] <- M13_0_marg
    X3[,7] <- M14_0_marg
    X3[,8] <- M15_0_marg
    X3[,9] <- M16_0_marg
    X3[,10] <- M17_0_marg
    X3[,11] <- M18_0_marg
    
    X3[,29] <- M11_0_marg # High:M11
    X3[,30] <- M12_0_marg # High:M12
    X3[,31] <- M13_0_marg # High:M13
    X3[,32] <- M14_0_marg # High:M14
    X3[,33] <- M15_0_marg # High:M15
    X3[,34] <- M16_0_marg # High:M16
    X3[,35] <- M17_0_marg # High:M17
    X3[,36] <- M18_0_marg # High:M18
    
    X3[,40] <- M11_0_marg* M12_0_marg # M11:M12
    X3[,41] <- M11_0_marg*M13_0_marg # M11:M13
    X3[,42] <- M11_0_marg*M14_0_marg # M11:M14
    X3[,43] <- M11_0_marg*M15_0_marg # M11:M15
    X3[,44] <- M11_0_marg*M16_0_marg # M11:M16
    X3[,45] <- M11_0_marg*M17_0_marg # M11:M17
    X3[,46] <- M11_0_marg*M18_0_marg # M11:M18
    X3[,50] <- M12_0_marg*M13_0_marg # M12:M13
    X3[,51] <- M12_0_marg*M14_0_marg # M12:M14
    X3[,52] <- M12_0_marg*M15_0_marg # M12:M15
    X3[,53] <- M12_0_marg*M16_0_marg # M12:M16
    X3[,54] <- M12_0_marg*M17_0_marg # M12:M17
    X3[,55] <- M12_0_marg*M18_0_marg # M12:M18
    X3[,59] <- M13_0_marg*M14_0_marg # M13:M14
    X3[,60] <- M13_0_marg*M15_0_marg # M13:M15
    X3[,61] <- M13_0_marg*M16_0_marg # M13:M16
    X3[,62] <- M13_0_marg*M17_0_marg # M13:M17
    X3[,63] <- M13_0_marg*M18_0_marg # M13:M18
    X3[,67] <- M14_0_marg*M15_0_marg # M14:M15
    X3[,68] <- M14_0_marg*M16_0_marg # M14:M16
    X3[,69] <- M14_0_marg*M17_0_marg # M14:M17
    X3[,70] <- M14_0_marg*M18_0_marg # M14:M18
    X3[,74] <- M15_0_marg*M16_0_marg # M15:M16
    X3[,75] <- M15_0_marg*M17_0_marg # M15:M17
    X3[,76] <- M15_0_marg*M18_0_marg # M15:M18
    X3[,80] <- M16_0_marg*M17_0_marg # M16:M17
    X3[,81] <- M16_0_marg*M18_0_marg # M16:M18
    X3[,85] <- M17_0_marg*M18_0_marg # M17:M18
    
    X3[,12] <- M2_0_cond
    X3[,37] <- M2_0_cond # High:M2
    
    X3[,47] <- M11_0_marg*M2_0_cond # M11:M2
    X3[,56] <- M12_0_marg*M2_0_cond # M12:M2
    X3[,64] <- M13_0_marg*M2_0_cond # M13:M2
    X3[,71] <- M14_0_marg*M2_0_cond # M14:M2
    X3[,77] <- M15_0_marg*M2_0_cond # M15:M2
    X3[,82] <- M16_0_marg*M2_0_cond # M16:M2
    X3[,86] <- M17_0_marg*M2_0_cond # M17:M2
    X3[,89] <- M18_0_marg*M2_0_cond # M18:M2
    
    M3_0_condM1M2 <- rbinom(n3, 1, predict(fit_M3_condM1M2, newdata = X3[,c(1:11,12,15:17,18:26,29:37,40:47,50:56,59:64,67:71,74:77,80:82,85:86,89)], type = "response"))
    
    
    # M4 - given M2 and M3 #
    
    # For the unexposed condition (high SES) #
    
    X3[,12] <- M2_0_marg
    X3[,13] <- M3_0_condM2
    
    X3[,37] <- M2_0_marg # High:M2
    X3[,38] <- M3_0_condM2 # High:M3
    
    X3[,92] <- M2_0_marg*M3_0_condM2 # M2:M3
    
    M4_0_condM2M3 <- rbinom(n3, 1, predict(fit_M4_condM2M3, newdata = X3[,c(1:3,12:13,15:17,26:27,37:38,92)], type = "response"))
    
    
    # M4 - given the M1 group, M2 and M3 #
    
    # For the unexposed condition (high SES) #
    
    X3[,12] <- M2_0_cond
    
    X3[,13] <- M3_0_condM1M2
    
    
    X3[,37] <- M2_0_cond # High:M2
    X3[,38] <- M3_0_condM1M2 # High:M3
    
    X3[,48] <- M11_0_marg*M3_0_condM1M2 # M11:M3
    X3[,57] <- M12_0_marg*M3_0_condM1M2 # M12:M3
    X3[,65] <- M13_0_marg*M3_0_condM1M2 # M13:M3
    X3[,72] <- M14_0_marg*M3_0_condM1M2 # M14:M3
    X3[,78] <- M15_0_marg*M3_0_condM1M2 # M15:M3
    X3[,83] <- M16_0_marg*M3_0_condM1M2 # M16:M3
    X3[,87] <- M17_0_marg*M3_0_condM1M2 # M17:M3
    
    X3[,92] <- M2_0_cond*M3_0_condM1M2# M2:M3
    
    M4_0_condM1M2M3 <- rbinom(n3, 1, predict(fit_M4_condM1M2M3, newdata = X3[,c(1:13,15:17,18:27,29:38,40:48,50:57,59:65,67:72,74:78,80:83,85:87,89,92)], type = "response"))
    
    # For the exposed condition (low SES) #
    
    X3[,3] <- 0
    
    X3[,4] <- M11_1_marg
    X3[,5] <- M12_1_marg
    X3[,6] <- M13_1_marg
    X3[,7] <- M14_1_marg
    X3[,8] <- M15_1_marg
    X3[,9] <- M16_1_marg
    X3[,10] <- M17_1_marg
    X3[,11] <- M18_1_marg
    
    X3[,12] <- M2_1_cond
    
    X3[,13] <- M3_1_condM1M2
    
    X3[,29] <- 0 # High:M11
    X3[,30] <- 0
    X3[,31] <- 0
    X3[,32] <- 0
    X3[,33] <- 0
    X3[,34] <- 0
    X3[,35] <- 0
    X3[,36] <- 0 # High:M18
    
    
    X3[,40] <- M11_1_marg* M12_1_marg # M11:M12
    X3[,41] <- M11_1_marg*M13_1_marg # M11:M13
    X3[,42] <- M11_1_marg*M14_1_marg # M11:M14
    X3[,43] <- M11_1_marg*M15_1_marg # M11:M15
    X3[,44] <- M11_1_marg*M16_1_marg # M11:M16
    X3[,45] <- M11_1_marg*M17_1_marg # M11:M17
    X3[,46] <- M11_1_marg*M18_1_marg# M11:M18
    X3[,50] <- M12_1_marg*M13_1_marg # M12:M13
    X3[,51] <- M12_1_marg*M14_1_marg # M12:M14
    X3[,52] <- M12_1_marg*M15_1_marg # M12:M15
    X3[,53] <- M12_1_marg*M16_1_marg # M12:M16
    X3[,54] <- M12_1_marg*M17_1_marg # M12:M17
    X3[,55] <- M12_1_marg*M18_1_marg# M12:M18
    X3[,59] <- M13_1_marg*M14_1_marg # M13:M14
    X3[,60] <- M13_1_marg*M15_1_marg # M13:M15
    X3[,61] <- M13_1_marg*M16_1_marg # M13:M16
    X3[,62] <- M13_1_marg*M17_1_marg # M13:M17
    X3[,63] <- M13_1_marg*M18_1_marg# M13:M18
    X3[,67] <- M14_1_marg*M15_1_marg # M14:M15
    X3[,68] <- M14_1_marg*M16_1_marg # M14:M16
    X3[,69] <- M14_1_marg*M17_1_marg # M14:M17
    X3[,70] <- M14_1_marg*M18_1_marg# M14:M18
    X3[,74] <- M15_1_marg*M16_1_marg # M15:M16
    X3[,75] <- M15_1_marg*M17_1_marg # M15:M17
    X3[,76] <- M15_1_marg*M18_1_marg# M15:M18
    X3[,80] <- M16_1_marg*M17_1_marg # M16:M17
    X3[,81] <- M16_1_marg*M18_1_marg# M16:M18
    X3[,85] <- M17_1_marg*M18_1_marg# M17:M18
    
    X3[,37] <- 0 # High:M2
    X3[,47] <- M11_1_marg*M2_1_cond # M11:M2
    X3[,56] <- M12_1_marg*M2_1_cond # M12:M2
    X3[,64] <- M13_1_marg*M2_1_cond # M13:M2
    X3[,71] <- M14_1_marg*M2_1_cond # M14:M2
    X3[,77] <- M15_1_marg*M2_1_cond # M15:M2
    X3[,82] <- M16_1_marg*M2_1_cond # M16:M2
    X3[,86] <- M17_1_marg*M2_1_cond # M17:M2
    X3[,89] <- M18_1_marg*M2_1_cond # M18:M2
    
    X3[,38] <- 0 # High:M3
    X3[,48] <- M11_1_marg*M3_1_condM1M2 # M11:M3
    X3[,57] <- M12_1_marg*M3_1_condM1M2 # M12:M3
    X3[,65] <- M13_1_marg*M3_1_condM1M2 # M13:M3
    X3[,72] <- M14_1_marg*M3_1_condM1M2 # M14:M3
    X3[,78] <- M15_1_marg*M3_1_condM1M2 # M15:M3
    X3[,83] <- M16_1_marg*M3_1_condM1M2 # M16:M3
    X3[,87] <- M17_1_marg*M3_1_condM1M2 # M17:M3
    
    X3[,92] <- M2_1_cond*M3_1_condM1M2# M2:M3
    
    M4_1_condM1M2M3 <- rbinom(n3, 1, predict(fit_M4_condM1M2M3, newdata = X3[,c(1:13,15:17,18:27,29:38,40:48,50:57,59:65,67:72,74:78,80:83,85:87,89,92)], type = "response"))
    
    ### Outcome ###
    
    ## Estimation of conditional expectations ###
    
    
    # With joint distribution for M1 and M2 #
    
    # Exposed condition (low SES) for exposure, M1, M2 and M3, unexposed condition for M4 #
    
    X3[,13] <- M3_1_marg
    X3[,38] <- 0 # High:M3
    X3[,48] <- M11_1_marg*M3_1_marg # M11:M3
    X3[,57] <- M12_1_marg*M3_1_marg # M12:M3
    X3[,65] <- M13_1_marg*M3_1_marg # M13:M3
    X3[,72] <- M14_1_marg*M3_1_marg # M14:M3
    X3[,78] <- M15_1_marg*M3_1_marg # M15:M3
    X3[,83] <- M16_1_marg*M3_1_marg # M16:M3
    X3[,87] <- M17_1_marg*M3_1_marg # M17:M3
    X3[,90] <- M18_1_marg*M3_1_marg # M18:M3
    
    X3[,92] <- M2_1_cond*M3_1_marg # M2:M3
    
    X3[,14] <- M4_0_marg #
    X3[,39] <- 0 # High:M4
    X3[,49] <- M11_1_marg*M4_0_marg # M11:M4
    X3[,58] <- M12_1_marg*M4_0_marg  # M12:M4
    X3[,66] <- M13_1_marg*M4_0_marg  # M13:M4
    X3[,73] <- M14_1_marg*M4_0_marg  # M14:M4
    X3[,79] <- M15_1_marg*M4_0_marg  # M15:M4
    X3[,84] <- M16_1_marg*M4_0_marg  # M16:M4
    X3[,88] <- M17_1_marg*M4_0_marg  # M17:M4
    X3[,91] <- M18_1_marg*M4_0_marg  # M18:M4
    
    X3[,93] <- M2_1_cond*M4_0_marg # M2:M4
    X3[,94] <- M3_1_marg*M4_0_marg # M3:M4
    
    
    y1_m121_m31_m40[i] <- mean(predict(fit_Y, newdata = X3, type = "response")) # rätt
    
    
    # Exposed condition (low SES) for exposure, M1, and M2, unexposed condition for M3 and M4 #
    
    X3[,13] <- M3_0_marg #
    
    X3[,38] <- 0 # High:M3
    X3[,48] <- M11_1_marg*M3_0_marg # M11:M3
    X3[,57] <- M12_1_marg*M3_0_marg # M12:M3
    X3[,65] <- M13_1_marg*M3_0_marg # M13:M3
    X3[,72] <- M14_1_marg*M3_0_marg # M14:M3
    X3[,78] <- M15_1_marg*M3_0_marg # M15:M3
    X3[,83] <- M16_1_marg*M3_0_marg # M16:M3
    X3[,87] <- M17_1_marg*M3_0_marg # M17:M3
    X3[,90] <- M18_1_marg*M3_0_marg # M18:M3
    
    X3[,92] <- M2_1_cond*M3_0_marg # M2:M3
    
    X3[,94] <- M3_0_marg*M4_0_marg # M3:M4
    
    y1_m121_m30_m40[i] <- mean(predict(fit_Y, newdata = X3, type = "response")) # rätt
    
    
    
    # With joint distribution for M1, M2 and M3 #
    
    # Exposed condition (low SES) for exposure, M1, M2 and M3, unexposed condition for M4 #
    
    X3[,13] <- M3_1_condM1M2 #
    
    X3[,38] <- 0 # High:M3
    X3[,48] <- M11_1_marg*M3_1_condM1M2 # M11:M3
    X3[,57] <- M12_1_marg*M3_1_condM1M2 # M12:M3
    X3[,65] <- M13_1_marg*M3_1_condM1M2 # M13:M3
    X3[,72] <- M14_1_marg*M3_1_condM1M2 # M14:M3
    X3[,78] <- M15_1_marg*M3_1_condM1M2 # M15:M3
    X3[,83] <- M16_1_marg*M3_1_condM1M2 # M16:M3
    X3[,87] <- M17_1_marg*M3_1_condM1M2 # M17:M3
    X3[,90] <- M18_1_marg*M3_1_condM1M2 # M18:M3
    
    X3[,92] <- M2_1_cond*M3_1_condM1M2 # M2:M3
    
    X3[,94] <- M3_1_condM1M2*M4_0_marg # M3:M4
    
    
    y1_m1231_m40[i] <- mean(predict(fit_Y, newdata = X3, type = "response")) # rätt
    
    
    # Exposed condition (low SES) for exposure, M1, M2, M3, and M4 #
    
    X3[,14] <- M4_1_marg #
    
    X3[,39] <- 0 # High:M4
    X3[,49] <- M11_1_marg*M4_1_marg # M11:M4
    X3[,58] <- M12_1_marg*M4_1_marg  # M12:M4
    X3[,66] <- M13_1_marg*M4_1_marg  # M13:M4
    X3[,73] <- M14_1_marg*M4_1_marg  # M14:M4
    X3[,79] <- M15_1_marg*M4_1_marg  # M15:M4
    X3[,84] <- M16_1_marg*M4_1_marg  # M16:M4
    X3[,88] <- M17_1_marg*M4_1_marg  # M17:M4
    X3[,91] <- M18_1_marg*M4_1_marg  # M18:M4
    
    X3[,93] <- M2_1_cond*M4_1_marg # M2:M4
    
    X3[,94] <- M3_1_condM1M2*M4_1_marg # M3:M4
    
    y1_m1231_m41[i] <- mean(predict(fit_Y, newdata = X3, type = "response")) # fel
    
    
    # With joint distribution for M1, M2, M3 and M4 #
    
    # Exposed condition (low SES) for exposure, M1, M2, M3, and M4 #
    
    X3[,14] <- M4_1_condM1M2M3
    
    X3[,39] <- 0 # High:M4
    X3[,49] <- M11_1_marg*M4_1_condM1M2M3 # M11:M4
    X3[,58] <- M12_1_marg*M4_1_condM1M2M3  # M12:M4
    X3[,66] <- M13_1_marg*M4_1_condM1M2M3  # M13:M4
    X3[,73] <- M14_1_marg*M4_1_condM1M2M3  # M14:M4
    X3[,79] <- M15_1_marg*M4_1_condM1M2M3  # M15:M4
    X3[,84] <- M16_1_marg*M4_1_condM1M2M3  # M16:M4
    X3[,88] <- M17_1_marg*M4_1_condM1M2M3  # M17:M4
    X3[,91] <- M18_1_marg*M4_1_condM1M2M3  # M18:M4
    
    X3[,93] <- M2_1_cond*M4_1_condM1M2M3 # M2:M4
    
    X3[,94] <- M3_1_condM1M2*M4_1_condM1M2M3 # M3:M4
    
    y1_m12341[i] <- mean(predict(fit_Y, newdata = X3, type = "response")) #
    
    
    # With joint distribution for M3 and M4 #
    
    # Exposed condition (low SES) for exposure, M1, and M2, unexposed condition for M3 and M4 #
    
    X3[,12] <- M2_1_marg
    X3[,37] <- 0 # High:M2
    
    X3[,47] <- M11_1_marg*M2_1_marg # M11:M2
    X3[,56] <- M12_1_marg*M2_1_marg # M12:M2
    X3[,64] <- M13_1_marg*M2_1_marg # M13:M2
    X3[,71] <- M14_1_marg*M2_1_marg # M14:M2
    X3[,77] <- M15_1_marg*M2_1_marg # M15:M2
    X3[,82] <- M16_1_marg*M2_1_marg # M16:M2
    X3[,86] <- M17_1_marg*M2_1_marg # M17:M2
    X3[,89] <- M18_1_marg*M2_1_marg # M18:M2
    
    X3[,13] <- M3_0_marg
    
    X3[,38] <- 0 # High:M3
    X3[,48] <- M11_1_marg*M3_0_marg # M11:M3
    X3[,57] <- M12_1_marg*M3_0_marg # M12:M3
    X3[,65] <- M13_1_marg*M3_0_marg # M13:M3
    X3[,72] <- M14_1_marg*M3_0_marg # M14:M3
    X3[,78] <- M15_1_marg*M3_0_marg # M15:M3
    X3[,83] <- M16_1_marg*M3_0_marg # M16:M3
    X3[,87] <- M17_1_marg*M3_0_marg # M17:M3
    X3[,90] <- M18_1_marg*M3_0_marg # M18:M3
    
    X3[,92] <- M2_1_marg*M3_0_marg # M2:M3
    
    X3[,14] <- M4_0_condM3
    
    X3[,39] <- 0 # High:M4
    X3[,49] <- M11_1_marg*M4_0_condM3 # M11:M4
    X3[,58] <- M12_1_marg*M4_0_condM3  # M12:M4
    X3[,66] <- M13_1_marg*M4_0_condM3  # M13:M4
    X3[,73] <- M14_1_marg*M4_0_condM3  # M14:M4
    X3[,79] <- M15_1_marg*M4_0_condM3  # M15:M4
    X3[,84] <- M16_1_marg*M4_0_condM3  # M16:M4
    X3[,88] <- M17_1_marg*M4_0_condM3  # M17:M4
    X3[,91] <- M18_1_marg*M4_0_condM3  # M18:M4
    
    X3[,93] <- M2_1_marg*M4_0_condM3 # M2:M4
    
    X3[,94] <- M3_0_marg*M4_0_condM3 # M3:M4
    
    y1_m11_m21_m340[i] <- mean(predict(fit_Y, newdata = X3, type = "response")) # rätt
    
    # Exposed condition (low SES) for exposure, M1, unexposed condition for M2, M3 and M4 #
    
    X3[,12] <- M2_0_marg
    
    X3[,37] <- 0 # High:M2
    X3[,47] <- M11_1_marg*M2_0_marg # M11:M2
    X3[,56] <- M12_1_marg*M2_0_marg # M12:M2
    X3[,64] <- M13_1_marg*M2_0_marg # M13:M2
    X3[,71] <- M14_1_marg*M2_0_marg # M14:M2
    X3[,77] <- M15_1_marg*M2_0_marg # M15:M2
    X3[,82] <- M16_1_marg*M2_0_marg # M16:M2
    X3[,86] <- M17_1_marg*M2_0_marg # M17:M2
    X3[,89] <- M18_1_marg*M2_0_marg # M18:M2
    
    X3[,92] <- M2_0_marg*M3_0_marg # M2:M3
    
    X3[,93] <- M2_0_marg*M4_0_condM3 # M2:M4
    
    y1_m11_m20_m340[i] <- mean(predict(fit_Y, newdata = X3, type = "response")) # rätt
    
    
    # With joint distribution for M2, M3 and M4 #
    
    # Exposed condition (low SES) for exposure, M1, unexposed condition for M2, M3 and M4 #
    
    X3[,13] <- M3_0_condM2
    X3[,14] <- M4_0_condM2M3
    
    X3[,38] <- 0 # High:M3
    X3[,48] <- M11_1_marg*M3_0_condM2 # M11:M3
    X3[,57] <- M12_1_marg*M3_0_condM2 # M12:M3
    X3[,65] <- M13_1_marg*M3_0_condM2 # M13:M3
    X3[,72] <- M14_1_marg*M3_0_condM2 # M14:M3
    X3[,78] <- M15_1_marg*M3_0_condM2 # M15:M3
    X3[,83] <- M16_1_marg*M3_0_condM2 # M16:M3
    X3[,87] <- M17_1_marg*M3_0_condM2 # M17:M3
    X3[,90] <- M18_1_marg*M3_0_condM2 # M18:M3
    X3[,92] <- M2_0_marg*M3_0_condM2 # M2:M3
    
    X3[,39] <- 0 # High:M4
    X3[,49] <- M11_1_marg*M4_0_condM2M3 # M11:M4
    X3[,58] <- M12_1_marg*M4_0_condM2M3  # M12:M4
    X3[,66] <- M13_1_marg*M4_0_condM2M3  # M13:M4
    X3[,73] <- M14_1_marg*M4_0_condM2M3  # M14:M4
    X3[,79] <- M15_1_marg*M4_0_condM2M3  # M15:M4
    X3[,84] <- M16_1_marg*M4_0_condM2M3  # M16:M4
    X3[,88] <- M17_1_marg*M4_0_condM2M3  # M17:M4
    X3[,91] <- M18_1_marg*M4_0_condM2M3  # M18:M4
    
    X3[,93] <- M2_0_marg*M4_0_condM2M3 # M2:M4
    
    X3[,94] <- M3_0_condM2*M4_0_condM2M3 # M3:M4
    
    y1_m11_m2340[i] <- mean(predict(fit_Y, newdata = X3, type = "response")) # rätt
    
    
    # Exposed condition (low SES) for exposure, unexposed condition for M1, M2, M3 and M4 #
    
    X3[,4] <- M11_0_marg
    X3[,5] <- M12_0_marg
    X3[,6] <- M13_0_marg
    X3[,7] <- M14_0_marg
    X3[,8] <- M15_0_marg
    X3[,9] <- M16_0_marg
    X3[,10] <- M17_0_marg
    X3[,11] <- M18_0_marg
    
    X3[,29] <- 0 # High:M11
    X3[,30] <- 0
    X3[,31] <- 0
    X3[,32] <- 0
    X3[,33] <- 0
    X3[,34] <- 0
    X3[,35] <- 0
    X3[,36] <- 0 # High:M18
    
    X3[,40] <- M11_0_marg* M12_0_marg # M11:M12
    X3[,41] <- M11_0_marg*M13_0_marg # M11:M13
    X3[,42] <- M11_0_marg*M14_0_marg # M11:M14
    X3[,43] <- M11_0_marg*M15_0_marg # M11:M15
    X3[,44] <- M11_0_marg*M16_0_marg # M11:M16
    X3[,45] <- M11_0_marg*M17_0_marg # M11:M17
    X3[,46] <- M11_0_marg*M18_0_marg# M11:M18
    X3[,50] <- M12_0_marg*M13_0_marg # M12:M13
    X3[,51] <- M12_0_marg*M14_0_marg # M12:M14
    X3[,52] <- M12_0_marg*M15_0_marg # M12:M15
    X3[,53] <- M12_0_marg*M16_0_marg # M12:M16
    X3[,54] <- M12_0_marg*M17_0_marg # M12:M17
    X3[,55] <- M12_0_marg*M18_0_marg# M12:M18
    X3[,59] <- M13_0_marg*M14_0_marg # M13:M14
    X3[,60] <- M13_0_marg*M15_0_marg # M13:M15
    X3[,61] <- M13_0_marg*M16_0_marg # M13:M16
    X3[,62] <- M13_0_marg*M17_0_marg # M13:M17
    X3[,63] <- M13_0_marg*M18_0_marg# M13:M18
    X3[,67] <- M14_0_marg*M15_0_marg # M14:M15
    X3[,68] <- M14_0_marg*M16_0_marg # M14:M16
    X3[,69] <- M14_0_marg*M17_0_marg # M14:M17
    X3[,70] <- M14_0_marg*M18_0_marg# M14:M18
    X3[,74] <- M15_0_marg*M16_0_marg # M15:M16
    X3[,75] <- M15_0_marg*M17_0_marg # M15:M17
    X3[,76] <- M15_0_marg*M18_0_marg# M15:M18
    X3[,80] <- M16_0_marg*M17_0_marg # M16:M17
    X3[,81] <- M16_0_marg*M18_0_marg# M16:M18
    X3[,85] <- M17_0_marg*M18_0_marg# M17:M18
    
    X3[,47] <- M11_0_marg*M2_0_marg # M11:M2
    X3[,56] <- M12_0_marg*M2_0_marg # M12:M2
    X3[,64] <- M13_0_marg*M2_0_marg # M13:M2
    X3[,71] <- M14_0_marg*M2_0_marg # M14:M2
    X3[,77] <- M15_0_marg*M2_0_marg # M15:M2
    X3[,82] <- M16_0_marg*M2_0_marg # M16:M2
    X3[,86] <- M17_0_marg*M2_0_marg # M17:M2
    X3[,89] <- M18_0_marg*M2_0_marg # M18:M2
    
    X3[,48] <- M11_0_marg*M3_0_condM2 # M11:M3
    X3[,57] <- M12_0_marg*M3_0_condM2 # M12:M3
    X3[,65] <- M13_0_marg*M3_0_condM2 # M13:M3
    X3[,72] <- M14_0_marg*M3_0_condM2 # M14:M3
    X3[,78] <- M15_0_marg*M3_0_condM2 # M15:M3
    X3[,83] <- M16_0_marg*M3_0_condM2 # M16:M3
    X3[,87] <- M17_0_marg*M3_0_condM2 # M17:M3
    X3[,90] <- M18_0_marg*M3_0_condM2 # M17:M3
    
    X3[,49] <- M11_0_marg*M4_0_condM2M3 # M11:M4
    X3[,58] <- M12_0_marg*M4_0_condM2M3  # M12:M4
    X3[,66] <- M13_0_marg*M4_0_condM2M3  # M13:M4
    X3[,73] <- M14_0_marg*M4_0_condM2M3  # M14:M4
    X3[,79] <- M15_0_marg*M4_0_condM2M3  # M15:M4
    X3[,84] <- M16_0_marg*M4_0_condM2M3  # M16:M4
    X3[,88] <- M17_0_marg*M4_0_condM2M3  # M17:M4
    X3[,91] <- M18_0_marg*M4_0_condM2M3  # M18:M4
    
    y1_m10_m2340[i] <- mean(predict(fit_Y, newdata = X3, type = "response")) # rätt
    
    
    # Exposed condition (low SES) for exposure, unexposed condition for M1, M2, M3 and M4 #
    
    X3[,12] <- M2_0_cond
    
    X3[,37] <- 0 # High:M2
    X3[,47] <- M11_0_marg*M2_0_cond # M11:M2
    X3[,56] <- M12_0_marg*M2_0_cond # M12:M2
    X3[,64] <- M13_0_marg*M2_0_cond # M13:M2
    X3[,71] <- M14_0_marg*M2_0_cond # M14:M2
    X3[,77] <- M15_0_marg*M2_0_cond # M15:M2
    X3[,82] <- M16_0_marg*M2_0_cond # M16:M2
    X3[,86] <- M17_0_marg*M2_0_cond # M17:M2
    X3[,89] <- M18_0_marg*M2_0_cond # M18:M2
    
    X3[,92] <- M2_0_cond*M3_0_condM1M2 # M2:M3
    
    X3[,93] <- M2_0_cond*M4_0_condM1M2M3 # M2:M4
    
    X3[,13] <- M3_0_condM1M2 #
    
    X3[,38] <- 0 # High:M3
    
    X3[,48] <- M11_0_marg*M3_0_condM1M2 # M11:M3
    X3[,57] <- M12_0_marg*M3_0_condM1M2 # M12:M3
    X3[,65] <- M13_0_marg*M3_0_condM1M2 # M13:M3
    X3[,72] <- M14_0_marg*M3_0_condM1M2 # M14:M3
    X3[,78] <- M15_0_marg*M3_0_condM1M2 # M15:M3
    X3[,83] <- M16_0_marg*M3_0_condM1M2 # M16:M3
    X3[,87] <- M17_0_marg*M3_0_condM1M2 # M17:M3
    X3[,90] <- M18_0_marg*M3_0_condM1M2 # M17:M3
    
    X3[,94] <- M3_0_condM1M2*M4_0_condM1M2M3 # M3:M4
    
    X3[,14] <- M4_0_condM1M2M3 #
    
    X3[,39] <- 0 # High:M4
    
    X3[,49] <- M11_0_marg*M4_0_condM1M2M3 # M11:M4
    X3[,58] <- M12_0_marg*M4_0_condM1M2M3  # M12:M4
    X3[,66] <- M13_0_marg*M4_0_condM1M2M3  # M13:M4
    X3[,73] <- M14_0_marg*M4_0_condM1M2M3  # M14:M4
    X3[,79] <- M15_0_marg*M4_0_condM1M2M3  # M15:M4
    X3[,84] <- M16_0_marg*M4_0_condM1M2M3  # M16:M4
    X3[,88] <- M17_0_marg*M4_0_condM1M2M3  # M17:M4
    X3[,91] <- M18_0_marg*M4_0_condM1M2M3  # M18:M4
    
    y1_m12340[i] <- mean(predict(fit_Y, newdata = X3, type = "response")) # rätt
    
    
    # Unexposed condition (high SES) for exposure, M1, M2, M3 and M4 #
    
    X3[,3] <- 1
    
    X3[,29] <- M11_0_marg # High:M11
    X3[,30] <- M12_0_marg # High:M12
    X3[,31] <- M13_0_marg # High:M13
    X3[,32] <- M14_0_marg # High:M14
    X3[,33] <- M15_0_marg # High:M15
    X3[,34] <- M16_0_marg # High:M16
    X3[,35] <- M17_0_marg # High:M17
    X3[,36] <- M18_0_marg # High:M18
    
    X3[,37] <- M2_0_cond # High:M2
    
    X3[,38] <- M3_0_condM1M2 # High:M3
    
    X3[,39] <- M4_0_condM1M2M3 # High:M4
    
    y0_m12340[i] <- mean(predict(fit_Y, newdata = X3, type = "response")) # rätt
    
    
  }
  
  ### Estimation of the interventional disparity effects ###
  
  ## Adjusted total association ##
  
  Adj_TA <- mean(y1_m12341 - y0_m12340)
  
  
  ## Direct effect
  
  IDM_DE <- mean(y1_m12340 - y0_m12340)
  
  
  ## Indirect effects
  
  # Via the M1 group #
  IDM_IE1 <- mean(y1_m11_m2340 - y1_m10_m2340)
  
  # Via M2 (NIHSS>5) #
  IDM_IE2 <- mean(y1_m11_m21_m340 - y1_m11_m20_m340)
  
  # Via M3 (Reperfusion) #
  IDM_IE3 <- mean(y1_m121_m31_m40 - y1_m121_m30_m40)
  
  # Via M4 (Stroke unit) #
  IDM_IE4 <- mean(y1_m1231_m41 - y1_m1231_m40)
  
  # Via the interdependence between the mediators #
  IDM_MD <- Adj_TA - (IDM_DE + IDM_IE1 + IDM_IE2 + IDM_IE3 + IDM_IE4)
  
  # Via all mediators taken together #
  IDM_IE <- mean(y1_m12341 - y1_m12340)
  
  
  ### Collect and return results
  
  resHigh <- c(Adj_TA, IDM_DE, IDM_IE, IDM_IE1, IDM_IE2, IDM_IE3, IDM_IE4, IDM_MD)
  
  res <- c(resMid, resHigh)
  
  if(!flag) return(res) else
    return(rep(NA,length(res)))

  
}