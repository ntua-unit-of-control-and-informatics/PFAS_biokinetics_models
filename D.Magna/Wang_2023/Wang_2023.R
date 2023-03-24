# Working directory

#setwd("/Users/elenistrompoula/Documents/GitHub/PFAS_biokinetics_models/D.Magna")
setwd("C:/Users/ptsir/Documents/GitHub/PFAS_biokinetics_models/D.Magna/Wang_2023")

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


# Filtering rate of Daphnia magna is calculated based on Burns et al. 1969 or Preuss.
# Input: length [mm],Temperature [oC]/ Output: filtration rate[mL/h]
Filtration_rate_estimation <<- function(length, temperature = 22, method = "Preuss"){
  if(method == "Burns"){
    F_rate_15 <- 0.153 * length^2.16
    F_rate_20 <- 0.208 * length^2.80
    F_rate_25 <- 0.202 * length^2.38
    
    if(temperature <= 15){
      F_rate <- F_rate_15
    }else if(temperature >= 25){  
      F_rate <- F_rate_25
    }else{ 
      F_rate <- approx(c(15,20,25), c(F_rate_15, F_rate_20, F_rate_25), temperature)$y
    }
  }else if(method == "Preuss"){
    F_rate <- 0.5*length^2.45
  }else{
    stop("Please select a valid estimation method; either 'Burns' or 'Preuss' ")
  }
  return(F_rate)
}


# Dumont et al. (1975)
# Input: length [mm]/ Output: dry weight[mg]
dry_weight_estimation <<- function(L){
  
  w1 = (1.89e-06*(L*1000)^2.25)/1000 #Donka Lake
  w2 = (4.88e-05*(L*1000)^1.80)/1000 #River Sambre
  # Selected w1 after validation with Martin-Creuzburg et al. (2018
  return(w1)
}


################################################################################
# Load the data for PFAS concentration

sheet_names <- c("PFBA", "F-53B", "GenX", "PFBS", "PFDA", "PFDoA", "PFHpA", 
                 "PFHxA", "PFNA", "PFOA", "PFOS", "PFPeA", "PFUnA")
data_ls <- list()
for(sheet_name in sheet_names){
  data_ls[[sheet_name]] <- openxlsx::read.xlsx ('Wang_data.xlsx', sheet = sheet_name)
}

# *** metrics ***

# The metric used for the optimization
mse_custom <- function(observed, predicted){
  mean((observed - predicted)^2)
}

mape <- function(observed, predicted){
  mean(abs(observed-predicted)*100/observed)
}

rmse <- function(observed, predicted){
  sqrt(mean((observed-predicted)^2)) 
}

AAFE <- function(observations, predictions, times=NULL){
  y_obs <- unlist(observations)
  y_pred <- unlist(predictions)
  # Total number of observations
  N<- length(y_obs)
  log_ratio <- rep(NA, N) 
  for ( i in 1:N){
    log_ratio[i] <- abs(log((y_pred[i]/y_obs[i]), base = 10))
  }
  aafe <- 10^(sum(log_ratio)/N) 
  return(aafe)
}

PBKOF <- function(observed, predicted, comp.names =NULL){
  # Check if the user provided the correct input format
  if (!is.list(observed) || !is.list(predicted)){
    stop(" The observations and predictions must be lists")
  }
  # Check if the user provided equal length lists
  if (length(observed) != length(predicted)){
    stop(" The observations and predictions must have the same compartments")
  }
  Ncomp <- length(observed) # Number of compartments
  I <- rep(NA, Ncomp) # Compartment discrepancy index
  N_obs <- rep(NA, Ncomp) #Number of observations per compartment
  #loop over the compartments
  for (i in 1:Ncomp){
    Et <- 0 #relative error with observations
    St <- 0  #relative error with simulations
    N <- length(observed[[i]]) # number of observations for compartment i
    # Check if observations and predictions have equal length
    if(N != length(predicted[[i]])){
      stop(paste0("Compartment ",i," had different length in the observations and predictions"))
    }
    N_obs[i] <- N # populate the N_obs vector
    for (j in 1:N){
      # sum of relative squared errors (error = observed - predicted)
      Et <- Et + ( abs(observed[[i]][j] - predicted[[i]][j])  / observed[[i]][j] )  ^2
      St <- St + ( abs(observed[[i]][j] - predicted[[i]][j])  / predicted[[i]][j] )  ^2
    }
    
    # root mean of the square of observed values
    RMEt <- sqrt(Et/N)
    # root mean of the square of simulated values
    RMSt <- sqrt( St/N)
    
    I[i] <- (RMEt + RMSt)/2   
  }
  # Total number of observations
  Ntot <- sum(N_obs)
  # Initialise the consolidated discrepancy index
  Ic <-0
  for (i in 1:Ncomp){
    # Give weight to compartments with more observations (more information)
    Ic <- Ic +  I[i]* N_obs[i]/Ntot
  }
  # Name the list of compartment discrepancy indices
  if ( !is.null(comp.names)){
    names(I) <- comp.names
  }else if (!is.null(names(observed))){
    names(I) <- names(observed)
  } else if (!is.null(names(predicted)) && is.null(comp.names) ){
    names(I) <- names(predicted)
  }
  return(Ic)
  #return(list(Total_index = Ic, Compartment_index= I))
}

#=====================================#
# Functions used for the optimization #
#=====================================#

# ode_func(): the differential equation system that describes the model

ode_func <- function(time, inits, params){
  with(as.list(c(inits, params)),{
    
    # Units explanation:
    # C_daphnia: ng PFAS/g D.magna WW
    # ke: 1/h
    # Cw: ng/L
    # Fsorption:  L/d/g D.magna WW
    
    age <- init_age + time
    #size in mm
    size <- Size_estimation(age, temperature = Temp, food="high")
    
    # Filtation rate in mL/h
    Frate <- Filtration_rate_estimation(size, temperature = Temp, method = "Preuss")
    # Filtation rate in L/day
    Frate <- Frate * 24/1000
    
    # dry weight mg
    DW <- dry_weight_estimation(size)
    # Convert DW to WW
    WW <- 15.5 * DW  # Conversion rate 11-20 from DW to WW (Garner et al., 2018)
    # Another post discussing DW to WW can be accessed through:
    #https://www.madsci.org/posts/archives/2005-09/1127049424.Zo.r.html
    
    
    # Water concentration in ng/L
    dCw <- 0
    # D.magna concentration in ng/g
    dC_daphnia <- Frate*Cw/WW +  Fsorption*Cw/WW  - ke*C_daphnia  
    
    return(list(c("dCw" = dCw, "dC_daphnia" = dC_daphnia), "Frate" = Frate,
                "WW" = WW))
  })
}



obj_func <- function(x, PFAS_data, PFAS_names, Cwater, age, temperatures, metric){
     
    # Indexes of body burden and exposure time in data frame
    BB_index <- c(2,4,6)
    ExpTime_index <- c(1,3,5)
    # Age of D.magna at beginning of exposure
    init_age <- age
    # Create a counter to mark the position of fitted parameters of x
    # that corresponds to a specific combination of PFAS and temperature
    counter <- 1
    score <- matrix(rep(NA, length(PFAS_names)*length(temperatures)), ncol = length(temperatures))
    # Iterate over PFAS names 
    for (PFAS in PFAS_names){
      # Load PFAS data
      df <- PFAS_data[[PFAS]]
      # Iterate over number of distinct temperature used in the experiment
      for (temp_iter in 1:length(temperatures)){
         # Initial water concentration of PFAS at selected temperature
         C_water <-  Cwater[PFAS,temp_iter]
         # Temperature of experiment
         Temp <- temperatures[temp_iter]
         # Time of measurement of selected PFAS at selected temperature
         exp_time <- df[!is.na(df[,ExpTime_index[temp_iter]]),ExpTime_index[temp_iter]]
         # Body burden of selected PFAS at selected temperature
         BodyBurden <- df[!is.na(df[,BB_index[temp_iter]]),BB_index[temp_iter]]
         # Time used by numerical solver that integrates the system of ODE
         sol_times <- seq(0,29, 0.01 )
         # Fitted parameters
         Fsorption <- x[counter]
         ke <- x[counter+1]
         counter <- counter + 2
         inits <- c( "Cw" = C_water, 'C_daphnia'= 0)
         params <- c("init_age"=init_age, "Temp" = Temp, "Fsorption"= Fsorption, "ke"  = ke)
         solution <- data.frame(deSolve::ode(times = sol_times,  func = ode_func,
                                             y = inits,
                                             parms = params,
                                             method="lsodes",
                                             rtol = 1e-5, atol = 1e-5))
                                
        if(sum(solution$time %in% exp_time) == length(exp_time)){
              results <- solution[which(solution$time %in% exp_time), 'C_daphnia']
        }else{
              stop(print("Length of predictions is not equal to the length of data"))
        }
        
        # Find the position of the current PFAS in the PFAS_names vector
         PFAS_position <- match(PFAS, PFAS_names)
        if(metric == "AAFE"){
          score[PFAS_position, temp_iter] <- AAFE(BodyBurden, results) 
        }else if (metric =="rmse"){
          score[PFAS_position, temp_iter] <- rmse(BodyBurden, results)
        }else if(metric == "PBKOF"){
          score[PFAS_position, temp_iter] <- PBKOF(list(BodyBurden), list(results))
       }       
      } 
    }
    # Take the average score of all PFAS and temperatures 
    final_score <- mean(score)  
  return(final_score)
}

plot_func <- function(params,PFAS_data, PFAS_name, Cwater, age, temperatures, metric){
  library(ggplot2)
  # User defined parameters
  init_age <- age
  # Fitted parameters
  Fsorption <- x[1]
  ke <- x[2]
  # Water concentration
  Cw <- Cwater * 1e06 #ng/L
  exp_time <- PFAS$Time
  sol_times <- seq(0,round(max(PFAS$Time))+1, 0.01 )
  inits <- c( "Cw" = Cw, 'C_daphnia'= 0)
  params <- c("init_age"=init_age, "Fsorption"= Fsorption, "ke"  = ke)
  solution <- data.frame(deSolve::ode(times = sol_times,  func = ode_func,
                                      y = inits,
                                      parms = params,
                                      method="lsodes",
                                      rtol = 1e-5, atol = 1e-5))
  keep_predictions <- data.frame(matrix(NA, nrow = length(sol_times), ncol =2))
  keep_predictions[,1] <- sol_times
  colnames(keep_predictions) <- c('Time', 'BodyBurden')
  keep_predictions[,2] <- solution[,"C_daphnia"]

    
      ggplot()+
      geom_line(data = keep_predictions, aes(x=Time, y=BodyBurden), size=1.7)+
      geom_point(data = PFAS, aes(x=Time, y=Concentration), size=5)+
      labs(title = paste(substitute(PFAS),"body burden in D.magna", sep = " "),
           y = "Body burden (ng/g WW)", x = "Time (days)")+
      theme(plot.title = element_text(hjust = 0.5,size=30), 
            axis.title.y =element_text(hjust = 0.5,size=25,face="bold"),
            axis.text.y=element_text(size=22),
            axis.title.x =element_text(hjust = 0.5,size=25,face="bold"),
            axis.text.x=element_text(size=22),
            legend.title=element_text(hjust = 0.5,size=25), 
            legend.text=element_text(size=22)) + 
      
      theme(legend.key.size = unit(1.5, 'cm'),  
            legend.title = element_text(size=14),
            legend.text = element_text(size=14),
            axis.text = element_text(size = 14))
}
  
##############################################################################

opts <- list( "algorithm" = "NLOPT_LN_SBPLX", #"NLOPT_LN_NEWUOA", #"NLOPT_LN_SBPLX" , #"NLOPT_LN_BOBYQA" #"NLOPT_LN_COBYLA"
              "xtol_rel" = 1e-07, 
              "ftol_rel" = 1e-07,
              "ftol_abs" = 0.0,
              "xtol_abs" = 0.0 ,
              "maxeval" = 1000,
              "print_level" = 1)

#---------------------------------------------
#.         PFOS
#---------------------------------------------
# Input preparation
PFAS_names <- c("PFDoA","PFUnA", "PFDA", "PFNA", "PFOS", "PFOA",   "PFHpA", 
                "PFHxA", "PFPeA",  "PFBS", "PFBA", "F-53B", "GenX")
# Water concentration in ng/mL
Cwater = matrix(c(1.44, 4.05, 9.56, 19.40, 18.5, 22.7, 22.5, 22.9, 22.6, 23.2,22.9,17.3, 22.5,
                  2.05, 4.73, 10.4, 20.4, 19.6, 23.2, 23.7, 23.6, 23.4, 22.7, 23.2, 19.2, 23.5,
                  2.31, 5.3, 12.1, 21.5, 21.9, 25.0, 25.2, 25.8, 26.9, 25.7, 24.6, 21.9, 24.8), ncol = 3)
colnames(Cwater) <-  c("16oC", "20oC", "24oC")
rownames(Cwater) <- PFAS_names
# Convert water concentration in ng/L
Cwater = Cwater*1000
age = 7 + 2*7 # age of D.magna at the beginning of exposure in days
temperatures <- c(16, 20, 24) #experiment temperature in oC 
# Define initial values of fitted parameters to provide to the optimization routine
# For each PFAS and temperature combination we have two parameters
x0 <- rep(c(1e-04, 0.05),length(PFAS_names) * length(temperatures))
optimization<- nloptr::nloptr(x0 = x0,
                               eval_f = obj_func,
                               lb	= c(0,0),
                               ub = c(100,100),
                               opts = opts,
                               PFAS_data = data_ls,
                               PFAS_names = PFAS_names, 
                               Cwater = Cwater, 
                               age = age , 
                               temperatures = temperatures,
                               metric = "AAFE")

PFAS_params <- array(optimization$solution, dim = c(2, 3, 13),
                     dimnames = list(c("Fsoprtion", "ke"),
                                     c("16oC", "20oC", "24oC"),
                                     PFAS_names))


plot_func(optimization_pfos, PFAS = PFOS, Cwater = 0.005, age = 1)

# 
# # Fitted parameters
# Fsorption <- PFOS_params["Fsorption"]
# ke <- PFOS_params["ke"]
# # Water concentration
# Cw <- 0.005 * 1e06 #ng/L
# exp_time <- PFOS$Time
# sol_times <- seq(0,round(max(PFOS$Time))+1, 0.01 )
# inits <- c('C_daphnia'= 0, "Cw" = Cw)
# params <- c("init_age"=1, "Fsorption"= Fsorption, "ke"  = ke)
# solution <- data.frame(deSolve::ode(times = sol_times,  func = ode_func,
#                                     y = inits,
#                                     parms = params,
#                                     method="lsodes",
#                                     rtol = 1e-5, atol = 1e-5))

