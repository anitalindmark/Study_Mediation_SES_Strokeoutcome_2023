#######################################################################
# Pooling of the estimates from the main analyses using Rubin's rules #
#######################################################################

# Create a matrix with the point estimates from the 45 imputed datasets #
#      Each column corresponds to one effect, each row to one dataset   # 
Estimates <- rbind(res_int$t0, res_int2$t0, res_int3$t0, res_int4$t0,
                  res_int5$t0, res_int6$t0, res_int7$t0, res_int8$t0,
                  res_int9$t0, res_int10$t0, res_int11$t0, res_int12$t0,
                  res_int13$t0, res_int14$t0, res_int15$t0, res_int16$t0,
                  res_int17$t0, res_int18$t0, res_int19$t0, res_int20$t0,
                  res_int21$t0, res_int22$t0, res_int23$t0, res_int24$t0,
                  res_int25$t0, res_int26$t0, res_int27$t0, res_int28$t0,
                  res_int29$t0, res_int30$t0, res_int31$t0, res_int32$t0,
                  res_int33$t0, res_int34$t0, res_int35$t0, res_int36$t0,
                  res_int37$t0, res_int38$t0, res_int39$t0, res_int40$t0,
                  res_int41$t0, res_int42$t0, res_int43$t0, res_int44$t0, res_int45$t0)


# Create a matrix with the bootstrap variances from the 45 imputed datasets #
#      Each column corresponds to one effect, each row to one dataset       # 

Variances <- rbind(apply(res_int$t, 2, var, na.rm=T), apply(res_int2$t, 2, var, na.rm=T), apply(res_int3$t, 2, var, na.rm=T), apply(res_int4$t, 2, var, na.rm=T),
                  apply(res_int5$t, 2, var, na.rm=T), apply(res_int6$t, 2, var, na.rm=T), apply(res_int7$t, 2, var, na.rm=T), apply(res_int8$t, 2, var, na.rm=T),
                  apply(res_int9$t, 2, var, na.rm=T), apply(res_int10$t, 2, var, na.rm=T), apply(res_int11$t, 2, var, na.rm=T), apply(res_int12$t, 2, var, na.rm=T),
                  apply(res_int13$t, 2, var, na.rm=T), apply(res_int14$t, 2, var, na.rm=T), apply(res_int15$t, 2, var, na.rm=T), apply(res_int16$t, 2, var, na.rm=T),
                  apply(res_int17$t, 2, var, na.rm=T), apply(res_int18$t, 2, var, na.rm=T), apply(res_int19$t, 2, var, na.rm=T), apply(res_int20$t, 2, var, na.rm=T),
                  apply(res_int21$t, 2, var, na.rm=T), apply(res_int22$t, 2, var, na.rm=T), apply(res_int23$t, 2, var, na.rm=T), apply(res_int24$t, 2, var, na.rm=T),
                  apply(res_int25$t, 2, var, na.rm=T), apply(res_int26$t, 2, var, na.rm=T), apply(res_int27$t, 2, var, na.rm=T), apply(res_int28$t, 2, var, na.rm=T),
                  apply(res_int29$t, 2, var, na.rm=T), apply(res_int30$t, 2, var, na.rm=T), apply(res_int31$t, 2, var, na.rm=T), apply(res_int32$t, 2, var, na.rm=T),
                  apply(res_int33$t, 2, var, na.rm=T), apply(res_int34$t, 2, var, na.rm=T), apply(res_int35$t, 2, var, na.rm=T), apply(res_int36$t, 2, var, na.rm=T),
                  apply(res_int37$t, 2, var, na.rm=T), apply(res_int38$t, 2, var, na.rm=T), apply(res_int39$t, 2, var, na.rm=T), apply(res_int40$t, 2, var, na.rm=T),
                  apply(res_int41$t, 2, var, na.rm=T), apply(res_int42$t, 2, var, na.rm=T), apply(res_int43$t, 2, var, na.rm=T), apply(res_int44$t, 2, var, na.rm=T),
                  apply(res_int45$t, 2, var, na.rm=T))
  

# Pooling of the results using Rubin's rules #

pooled_estimates <- colMeans(Estimates) # Vector with the pooled estimated effects 

within_imp_var <- colMeans(Variances) # Vector with the within imputation variances 

between_imp_var <- apply(Estimates,2,var) # Vector with the between imputation variances 

total_variance <- within_imp_var + (1+1/45)*between_imp_var # Vector with the total variances 

dfs <- (45-1)*(1+within_imp_var/((1+1/45)*between_imp_var))^2 # Vector with the degrees of freedom

P_vals <- 2*pt(abs(pooled_estimates/sqrt(total_variance)),df=dfs, lower.tail=F) # Vector with p-values based on the t-distribution

Prop_med <- round(c(100*pooled_estimates[2:8]/pooled_estimates[1],100*pooled_estimates[10:16]/pooled_estimates[9]),1) # Vector with the proportions of the adjusted total association

# Lower and upper confidence interval limits #
CIlow <- pooled_estimates - qt(0.975,dfs)*sqrt(total_variance)
CIupp <- pooled_estimates + qt(0.975,dfs)*sqrt(total_variance)


