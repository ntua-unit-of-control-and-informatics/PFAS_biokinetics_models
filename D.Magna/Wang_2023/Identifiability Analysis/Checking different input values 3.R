library(parallel)
library(deSolve)
library(nloptr)
#setwd("C:/Users/ptsir/Documents/GitHub/PFAS_biokinetics_models/D.Magna/Wang_2023")
setwd("/Users/elenistrompoula/Documents/GitHub/PFAS_biokinetics_models/D.Magna/Wang_2023")

rmse <- function(observed, predicted){
  sqrt(mean((observed-predicted)^2))
}

#=====================================#
#  Weighted Sum of Squared Residuals  #
#=====================================#

WSSR <- function(observed, predicted, weights, comp.names =NULL){
  # Check if the user provided the correct input format
  if (!is.list(observed) || !is.list(predicted) || !is.list(weights)){
    stop(" The observations, predictions and weights must be lists")
  }
  # Check if the user provided equal length lists
  if (length(observed) != length(predicted) || length(observed) != length(weights)){
    stop(" The observations, predictions and weights must have the same compartments")
  }

  # Define the number of observed outputs
  N_outputs <- length(predicted)
  # Define the number of observations per output
  N_obs <- rep(NA, N_outputs)

  # A vector to store the values of the weighted squared sum residuals of each compartment
  outputs_res <- c()
  for (i in 1:N_outputs) { # loop over the observed outputs
    N_obs[i] <- length(observed[[i]])

    # Check that all observed, predicted and weights vectors have the same length
    if(N_obs[i] != length(predicted[[i]]) || N_obs[i] != length(weights[[i]])){
      stop(paste0("Compartment ",i," had different length in the observations and predictions"))
    }
    # The number of observations for output i
    N <- N_obs[i]

    # Initiate a variable to estimate the sum of squared residuals for output j
    sq_weighted_res_j <- 0
    for (j in 1:N) { #loop over the experimental points i compartment i
      sq_weighted_res_j <- sq_weighted_res_j + ((observed[[i]][j] - predicted[[i]][j]) / weights[[i]][j])^2
    }
    outputs_res[i] <- sq_weighted_res_j
  }

  WSSR_results <- sum(outputs_res)

  # Name the list of compartment discrepancy indices
  if ( !is.null(comp.names)){
    names(WSSR_results) <- comp.names
  }else if (!is.null(names(observed))){
    names(WSSR_results) <- names(observed)
  } else if (!is.null(names(predicted)) && is.null(comp.names) ){
    names(WSSR_results) <- names(predicted)
  } else if (!is.null(names(weights)) && is.null(comp.names) ){
    names(WSSR_results) <- names(weights)
  }

  return(WSSR_results)
}

#=================#
# Model functions #
#=================#

# Function for estimating length of D. magna based on age (from Betini et al. (2019))
# Input: age [days], temperature [oC], food["low"/"high"]/ Output: length [mm]
# Considers female D. magna
Size_estimation <<- function(age, temperature = 22, food="high"){

  # T = 15 o C
  a_low_15 <-0.354
  b_low_15 <- 0.527
  a_high_15 <- 0.105
  b_high_15 <- 0.953

  # T = 25 o C
  a_low_25 <- 0.811
  b_low_25 <- 0.355
  a_high_25 <- 0.698
  b_high_25 <- 0.83

  if(food == "low"){
    if(temperature <= 15){
      a <- a_low_15
      b <- b_low_15
    }else if(temperature >= 25){
      a <- a_low_25
      b <- b_low_25
    }else{
      a <- approx(c(15,25), c(a_low_15, a_low_25), temperature)$y
      b <- approx(c(15,25), c(b_low_15, b_low_25), temperature)$y
    }
  }else if (food == "high"){
    if(temperature <= 15){
      a <- a_high_15
      b <- b_high_15
    }else if(temperature >= 25){
      a <- a_high_25
      b <- b_high_25
    }else{
      a <- approx(c(15,25), c(a_high_15, a_high_25), temperature)$y
      b <- approx(c(15,25), c(b_high_15, b_high_25), temperature)$y
    }
  }else{
    stop('food must be either "low" or "high" ')
  }
  return(a + b * log(age))
}

# Dumont et al. (1975)
# Input: length [mm]/ Output: dry weight[mg]
dry_weight_estimation <<- function(L){

  w1 = (1.89e-06*(L*1000)^2.25)/1000 #Donka Lake
  w2 = (4.88e-05*(L*1000)^1.80)/1000 #River Sambre
  # Selected w1 after validation with Martin-Creuzburg et al. (2018
  return(w1)
}

#==============
# ODEs System
#==============

ode_func <- function(time, inits, params){
  with(as.list(c(inits, params)),{
    # Units explanation:
    # ke: 1/day
    # ku: Lwater/day
    # Cw: ng PFAS/L water
    # C_daphnia_unbound/C_daphnia_bound: mol PFAS/L D.magna
    # C_daphnia_unbound_unmol/C_daphnia_bound_unmol: ng PFAS/g D.magna
    # C_prot_un: mol prot/L

    age <- init_age + time
    #size in mm
    size <- Size_estimation(age, temperature = 23, food="high")


    # dry weight mg
    DW <- dry_weight_estimation(size)
    # Convert DW to WW
    WW <- 15.5 * DW  # Conversion rate 11-20 from DW to WW (Garner et al., 2018)
    # Another post discussing DW to WW can be accessed through:
    #https://www.madsci.org/posts/archives/2005-09/1127049424.Zo.r.html

    # Water concentration in ng/L
    dCw <- 0
    
    ku <- 10^ku
    ke <- 10^ke
    kon <- (10^kon)*60*60*24 #convert to mol/l/d
    Ka <- 10^Ka 
    koff <- kon/Ka
    
    # Reported values for Kon and koff range between 1e2-1e04 and 1e-3-1e-1 in L/mol/s
    # and s-^-1 respectively. Our concentrations are in ng/g. Assuming density = 1000g/L
    # then the concentration in ng/g is multyplied by 1000 to make it ng/L and then by
    # multiply by 1e-9 to make it grams and divide by MW. We do this directly to kon
    # and koff to make the concentration mol/L so that we can match the literature
    # values of kon and koff with ours
    
    
    dC_daphnia_unbound <-  ku*(Cw*1e-09/MW)/(WW/1000)  - kon*C_prot_un*C_daphnia_unbound +
                           koff*C_daphnia_bound - ke*C_daphnia_unbound
    dC_daphnia_bound <- kon*C_prot_un*C_daphnia_unbound - koff*C_daphnia_bound
    dC_prot_un <-   koff*C_daphnia_bound -  kon*C_prot_un*C_daphnia_unbound
    
    # Multiply MW by 1e09 to convert g/mol to ng/mol and then 
    # divide by the density of water to make L to g wet weight
    C_daphnia_unbound_unmol <- C_daphnia_unbound*(MW*1e9)/1000
    C_daphnia_bound_unmol <- C_daphnia_bound*(MW*1e9)/1000
    C_tot <- C_daphnia_unbound_unmol + C_daphnia_bound_unmol
    
    return(list(c("dCw" = dCw,   "dC_daphnia_unbound" = dC_daphnia_unbound,
                  "dC_daphnia_bound" = dC_daphnia_bound, "dC_prot_un" = dC_prot_un),
                "WW" = WW, "C_tot" = C_tot))
  })
}


obj_f <- function(x, params_names, constant_theta, constant_theta_names,
                  constant_params=NULL,data_df, errors_df,
                  PFAS_name, Cwater, age, temperatures, MW){

  # Assign the values of the x vector to the corresponding parameters
  if(!is.null(constant_theta)){
    if(length(constant_theta_names) != length(constant_theta)){
      stop("The constant_theta_names vector must be of equal length with the constant_theta vector")
    }
    for (j in 1:length(constant_theta)){
      assign(constant_theta_names[j], constant_theta[[j]])
    }
  }
  # Assign the values of the x vector to the corresponding parameters
  if(length(x) != length(params_names)){
    stop("The params_names must be of equal length with the x vector")
  }
  for (k in 1:length(x)) {
    assign(params_names[k], x[k])
  }

  if(!is.null(constant_params)){
    for (k in 1:length(constant_params)) {
      assign(names(constant_params)[k], constant_params[[k]])
    }
  }

  # Indexes of body burden and exposure time in data frame
  BB_index <- c(2,4,6)
  ExpTime_index <- c(1,3,5)
  Errors_index <- c(2,4,6)
  # Age of D.magna at beginning of exposure
  init_age <- age
  # Initiate a vactor to sae the score values per temperature
  score <- rep(NA, length(temperatures))
  # Load PFAS data
  df <- data_df

  # Iterate over number of distinct temperature used in the experiment
  for (temp_iter in 1:length(temperatures)){
    # Initial water concentration of PFAS at selected temperature
    C_water <-  Cwater[[temp_iter]]
    # Temperature of experiment
    Temp <- temperatures[temp_iter]
    # Time of measurement of selected PFAS at selected temperature
    exp_time <- round(df[!is.na(df[,ExpTime_index[temp_iter]]),ExpTime_index[temp_iter]],1)
    # Body burden of selected PFAS at selected temperature
    BodyBurden <- df[!is.na(df[,BB_index[temp_iter]]),BB_index[temp_iter]]
    # Time used by numerical solver that integrates the system of ODE
    sol_times <- seq(0,15, 0.1 )
    # Fitted parameters


    inits <- c( "Cw" = C_water,  "C_daphnia_unbound" = 0,
                "C_daphnia_bound" = 0, "C_prot_un" = C_prot_init)

    params <- c("init_age"=age, "Temp" = Temp, "ku"= ku,
                "kon" = kon, "Ka" = Ka, "ke"= ke, "MW" = MW)
    solution <- data.frame(deSolve::ode(times = sol_times,  func = ode_func,
                                        y = inits,
                                        parms = params,
                                        method="lsodes",
                                        rtol = 1e-6, atol = 1e-6))

    if(sum(round(solution$time,2) %in% exp_time) == length(exp_time)){
      results <- solution[which(round(solution$time,2) %in% exp_time), 'C_tot']
      score[temp_iter] <- rmse(BodyBurden, results)
    }else{
     # stop(print("Length of predictions is not equal to the length of data"))
      score[temp_iter]=50000 
      }

    #score[temp_iter] <- rmse(BodyBurden, results)
  }

  # Take the average score of all PFAS and temperatures
  final_score <- mean(score)
  return(final_score)
}
################################################################################


# This is a wrapper for the optimization process in order to work in parallel
wrapper_opt <- function(X){

  PFAS_names <- c("PFBA", "F-53B", "GenX", "PFBS", "PFDA", "PFDoA", "PFHpA",
                  "PFHxA", "PFNA", "PFOA", "PFOS", "PFPeA", "PFUnA")
  Molecular_weights <- list("PFBA" = 214, "F-53B" = 570, "GenX" = 330, "PFBS" = 300, "PFDA" = 514,
                            "PFDoA" = 614, "PFHpA" = 364, "PFHxA" = 314,"PFNA" = 464, "PFOA" = 414,  "PFOS" = 500,
                            "PFPeA" = 364, "PFUnA" = 564)
  data_ls <- list()
  data_plot <- list()

  for(sheet_name in PFAS_names){
    data_ls[[sheet_name]] <- openxlsx::read.xlsx ('Wang_data_reduced2.xlsx', sheet = sheet_name)
    data_plot[[sheet_name]] <- openxlsx::read.xlsx ('Wang_data.xlsx', sheet = sheet_name)
  }

  # Water concentration in ng/mL
  Cwater = matrix(c(1.44, 4.05, 9.56, 19.40, 18.5, 22.7, 22.5, 22.9, 22.6, 23.2,22.9,17.3, 22.5,
                    2.05, 4.73, 10.4, 20.4, 19.6, 23.2, 23.7, 23.6, 23.4, 22.7, 23.2, 19.2, 23.5,
                    2.31, 5.3, 12.1, 21.5, 21.9, 25.0, 25.2, 25.8, 26.9, 25.7, 24.6, 21.9, 24.8), ncol = 3)
  colnames(Cwater) <-  c("16oC", "20oC", "24oC")
  rownames(Cwater) <- PFAS_names
  # Convert water concentration in ng/L
  Cwater = Cwater*1000
  age = 7+7 # age of D.magna at the beginning of exposure in days
  temperatures <- c(16, 20, 24) #experiment temperature in oC

  #x, constant_theta, constant_theta_names, params_names,
  # constant_params=NULL,data_df, error_df
  # x0 must be given in log10-scale
  # x0 contains the initial values of the Ka and ke
  x0 <- c(7,3)
  params_names <- c("kon", "ke")
  constant_theta = X
  constant_theta_names =  c("ku", "Ka", "C_prot_init")
  constant_params <- NULL

  # the selected settings for the optimizer
  opts <- list( "algorithm" = "NLOPT_LN_SBPLX", #"NLOPT_LN_NELDERMEAD" ,#"NLOPT_LN_SBPLX", #"NLOPT_LN_BOBYQA" #"NLOPT_LN_COBYLA"
                "xtol_rel" = 1e-10,
                "ftol_rel" = 1e-10,
                "ftol_abs" = 0,
                "xtol_abs" = 0 ,
                "maxeval" = 3000,
                "print_level" = 1)

  opt_params_per_substance <- list()
  score_per_substance <- c()

  for (i in 1:length(PFAS_names)) {
    optimization <- nloptr::nloptr(x0 = x0,
                                   eval_f = obj_f,
                                   lb	=  c(2,-4),
                                   ub =   c(10,5),
                                   constant_theta = constant_theta,
                                   constant_theta_names = constant_theta_names,
                                   params_names = params_names,
                                   constant_params=NULL,
                                   data_df = data_ls[[PFAS_names[i]]],
                                   errors_df = NULL,
                                   age=age,
                                   PFAS_name=PFAS_names[i],
                                   Cwater=Cwater[PFAS_names[i],],
                                   temperatures=temperatures,
                                   MW=Molecular_weights[[PFAS_names[i]]],
                                   opts = opts
                                   )
    score_per_substance[i] <- optimization$objective
    names(score_per_substance)[i] <- PFAS_names[i]
    opt_params_per_substance[[i]] <- optimization$solution
    names(opt_params_per_substance)[i] <- PFAS_names[i]
  }
  return(list("Fixed_params_used"=X,
              "Overall_score"=mean(score_per_substance),
              "opt_params_per_substance"=opt_params_per_substance))
}


################################################################################

# In this section, various values will be employed for the parameters Ku,
# C_prot_init, and kon. Following each parameter combination, the remaining
# parameters will undergo re-optimization. Once optimization is completed,
# the subsequent step involves applying the identifiability analysis function
# to derive conclusions.

# Here are the values of the parameters that will be tested
ku_values <- log10(c( 5e-3, 1e-2, 5e-2))
C_prot_init_values <- log10(c( 5e-4, 1e-3, 5e-3))
ka_values = log10(c(1e3, 5e3, 1e4))


# Generate all possible combinations of the parameters for each PFAS substance
# Create a list to save all the possible triplets of parameters
fixed_params <- list()
for (i in 1:length(ku_values)) {
  for (j in 1:length(C_prot_init_values)) {
    for (k in 1:length(ka_values)) {
      triplet <- c(ku = ku_values[i], Ka = ka_values[k], C_prot_init = C_prot_init_values[j])
      fixed_params[[length(fixed_params) + 1]] <- triplet
    }
  }
}


input_ls = fixed_params
start.time <- Sys.time()
# Set up the cluster.
cluster <- makeCluster(detectCores()-2)
# Export to the cluster any function or parameter that the obj_f needs to work.
clusterExport(cl=cluster,c("obj_f", "rmse", "ode_func",
                           "Size_estimation", "dry_weight_estimation"))
output <- parLapply(cluster, input_ls, wrapper_opt)
# Terminate the cluster.
stopCluster(cluster)
total.duration <- Sys.time() - start.time


scores <- c()
for (i in 1:length(output)) {
  scores[i] <- output[[i]]$Overall_score
}
print(output[which.min(scores)])
