# Working directory

#setwd("/Users/elenistrompoula/Documents/GitHub/PFAS_biokinetics_models/D.Magna")
#setwd("C:/Users/ptsir/Documents/GitHub/PFAS_biokinetics_models/D.Magna/Wang_2023")
setwd("C:/Users/user/Documents/GitHub/PFAS_biokinetics_models/D.Magna/Wang_2023")

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
    # Fsorption:  L/day/ gWW
    
    age <- init_age + time
    #size in mm
    size <- Size_estimation(age, temperature = 23, food="high")
    
    # Filtation rate in mL/h
    Frate <- Filtration_rate_estimation(size, temperature = 23, method = "Preuss")
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
    # D.magna concentration in lumenconcentration in ng/g
    dC_lumen <- alpha*Frate*Cw/WW  -  ku*C_lumen + ke*C_daphnia_unbound- alpha*Frate*C_lumen
    dC_daphnia_unbound <-  ku*C_lumen - ke*C_daphnia_unbound - kon*C_daphnia_unbound +  koff*C_daphnia_bound
    dC_daphnia_bound <- kon*C_daphnia_unbound - koff*C_daphnia_bound
    C_tot <- C_daphnia_unbound + C_daphnia_bound
    return(list(c("dCw" = dCw, "dC_lumen" = dC_lumen, "dC_daphnia_unbound" = dC_daphnia_unbound,
                  "dC_daphnia_bound" = dC_daphnia_bound), "Frate" = Frate,
                "WW" = WW, "C_tot" = C_tot))
  })
}



obj_func <- function(x, PFAS_data, PFAS_name, Cwater, age, temperatures, metric){
  
  # Indexes of body burden and exposure time in data frame
  BB_index <- c(2,4,6)
  ExpTime_index <- c(1,3,5)
  # Age of D.magna at beginning of exposure
  init_age <- age
  # Create a counter to mark the position of fitted parameters of x
  # that corresponds to a specific combination of PFAS and temperature
  counter <- 1
  score <- rep(NA, length(temperatures))
  # Iterate over PFAS names 
  # Load PFAS data
  df <- PFAS_data[[PFAS_name]]
  # Iterate over number of distinct temperature used in the experiment
  ku <- x[1]
  ke <- x[2]
  kon <- x[3]
  koff <- x[4]
  for (temp_iter in 1:length(temperatures)){
    # Initial water concentration of PFAS at selected temperature
    C_water <-  Cwater[PFAS_name,temp_iter]
    # Temperature of experiment
    Temp <- temperatures[temp_iter]
    # Time of measurement of selected PFAS at selected temperature
    exp_time <- round(df[!is.na(df[,ExpTime_index[temp_iter]]),ExpTime_index[temp_iter]],2)
    # Body burden of selected PFAS at selected temperature
    BodyBurden <- df[!is.na(df[,BB_index[temp_iter]]),BB_index[temp_iter]]
    # Time used by numerical solver that integrates the system of ODE
    sol_times <- seq(0,29, 0.01 )
    # Fitted parameters
    alpha <- x[4+temp_iter]
    
    
    inits <- c( "Cw" = C_water, "C_lumen" = 0, "C_daphnia_unbound" = 0,
                "C_daphnia_bound" = 0)
    params <- c("init_age"=age, "Temp" = Temp, "ku"= ku, 
                "ke"  = ke,  "kon" = kon, "koff" = koff, "alpha" = alpha)
    solution <- data.frame(deSolve::ode(times = sol_times,  func = ode_func,
                                        y = inits,
                                        parms = params,
                                        method="lsodes",
                                        rtol = 1e-5, atol = 1e-5))
    
    if(sum(round(solution$time,2) %in% exp_time) == length(exp_time)){
      results <- solution[which(round(solution$time,2) %in% exp_time), 'C_tot']
    }else{
      stop(print("Length of predictions is not equal to the length of data"))
    }
    
    # Find the position of the current PFAS in the PFAS_names vector
    if(metric == "AAFE"){
      score[temp_iter] <- AAFE(BodyBurden, results) 
    }else if (metric =="rmse"){
      score[temp_iter] <- rmse(BodyBurden, results)
    }else if(metric == "PBKOF"){
      score[temp_iter] <- PBKOF(list(BodyBurden), list(results))
    }       
  } 
  
  # Take the average score of all PFAS and temperatures 
  final_score <- mean(score)  
  return(final_score)
}

plot_func <- function(params,PFAS_data, PFAS_name, Cwater, age, temperatures){
  library(ggplot2)
  
  # Age of D.magna at beginning of exposure
  init_age <- age
  # Create a counter to mark the position of fitted parameters of x
  # that corresponds to a specific combination of PFAS and temperature
  
  # Load PFAS data
  df <- PFAS_data[[PFAS_name]]
  # Load parameters
  parameters <- params
  # Time used by numerical solver that integrates the system of ODE
  sol_times <- seq(0,15, 0.01 )
  # Data frame to store predictions for each temperature
  predictions <- data.frame("time" = sol_times, "BB_16" = rep(NA, length(sol_times)),
                            "BB_20" = rep(NA, length(sol_times)), "BB_24" = rep(NA, length(sol_times)))
  ku <- parameters[1]
  ke <-  parameters[2]
  kon <-  parameters[3]
  koff <-  parameters[4]
  # Iterate over number of distinct temperature used in the experiment
  for (temp_iter in 1:length(temperatures)){
    # Initial water concentration of PFAS at selected temperature
    C_water <-  Cwater[PFAS_name,temp_iter]
    # Temperature of experiment
    Temp <- temperatures[temp_iter]
    # Fitted parameters
    alpha <- parameters[4+temp_iter]
    
    inits <- c( "Cw" = C_water, "C_lumen" = 0, "C_daphnia_unbound" = 0,
                "C_daphnia_bound" = 0)
    params <- c("init_age"=age, "Temp" = Temp, "ku"= ku, 
                "ke"  = ke,  "kon" = kon, "koff" = koff, "alpha" = alpha )
    solution <- data.frame(deSolve::ode(times = sol_times,  func = ode_func,
                                        y = inits,
                                        parms = params,
                                        method="lsodes",
                                        rtol = 1e-5, atol = 1e-5))
    
    predictions[,temp_iter+1] <- solution$C_tot
  }
  ggplot()+
    geom_line(data = predictions, aes(x=time, y=BB_16,  colour = "16oC"), size=1.7)+
    geom_line(data = predictions, aes(x=time, y=BB_20, colour = "20oC"), size=1.7)+
    geom_line(data = predictions, aes(x=time, y=BB_24, colour = "24oC"), size=1.7)+
    geom_point(data = df, aes(x=time_16  , y=BB_16 ,  colour = "16oC"), size=5)+
    geom_point(data = df, aes(x=time_20, y=BB_20,  colour = "20oC"), size=5)+
    geom_point(data = df, aes(x=time_24, y=BB_24,  colour = "24oC"), size=5)+
    
    labs(title = paste(PFAS_name,"body burden in D.magna", sep = " "),
         y = "Body burden (ng/g WW)", x = "Time (days)")+
    scale_colour_manual("Temperature", 
                        breaks = c("16oC", "20oC", "24oC"),
                        values = c("red", "green", "blue")) +
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
  ggsave(paste0(PFAS_name,".png"))
}

##############################################################################
################################################################################
# Load the data for PFAS concentration

sheet_names <- c("PFBA", "F-53B", "GenX", "PFBS", "PFDA", "PFDoA", "PFHpA", 
                 "PFHxA", "PFNA", "PFOA", "PFOS", "PFPeA", "PFUnA")
PFAS_names <- sheet_names

data_ls <- list()
data_plot <- list()

for(sheet_name in sheet_names){
  data_ls[[sheet_name]] <- openxlsx::read.xlsx ('Wang_data_reduced.xlsx', sheet = sheet_name)
  data_plot[[sheet_name]] <- openxlsx::read.xlsx ('Wang_data.xlsx', sheet = sheet_name)
}

opts <- list( "algorithm" = "NLOPT_LN_SBPLX",#"NLOPT_LN_SBPLX", #"NLOPT_LN_NEWUOA", #"NLOPT_LN_SBPLX" , #"NLOPT_LN_BOBYQA" #"NLOPT_LN_COBYLA"
              "xtol_rel" = 1e-06, 
              "ftol_rel" = 1e-06,
              "ftol_abs" = 0.0,
              "xtol_abs" = 0.0 ,
              "maxeval" = 1000,
              "print_level" = 1)

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
# List to store parameters
parameters <- list()
# List to store solution for best parameter set
solutions <- list()
# List to store optimization results
optimizations <- list()
for (i in 1:length(PFAS_names)){
  #short chain PFAS behave differently and follow biexponential
  
  # Define initial values of fitted parameters to provide to the optimization routine
  # For each PFAS and temperature combination we have two parameters
  x0 <- c(1, 1, 1,1,rep(1,length(temperatures)))
  optimization<- nloptr::nloptr(x0 = x0,
                                eval_f = obj_func,
                                lb	=  c(0,0,0,0,0.2,0.2,0.2),
                                ub = c(1000,1000,1000,1000,5,5,5),
                                opts = opts,
                                PFAS_data = data_ls,
                                PFAS_name = PFAS_names[i],
                                Cwater = Cwater,
                                age = age ,
                                temperatures = temperatures,
                                metric = "PBKOF")
  optimizations[[PFAS_names[i]]] <- optimization
  parameters[[PFAS_names[i]]] <- optimization$solution
  names(parameters[[PFAS_names[i]]]) = c("ku", "ke", "kon","koff", "alpha_16", "alpha_20", "alpha_24")
  
  sol_times <- seq(0,15, 0.01 )
  # Iterate over number of distinct temperature used in the experiment
  temp_iter <- 2
  # Initial water concentration of PFAS at selected temperature
  C_water <-  Cwater[PFAS_names[i],temp_iter]
  # Temperature of experiment
  Temp <- temperatures[temp_iter]
  # Fitted parameters
  ku <- parameters[[PFAS_names[i]]][1]
  ke <- parameters[[PFAS_names[i]]][2]
  kon <-  parameters[[PFAS_names[i]]][3]
  koff <-  parameters[[PFAS_names[i]]][4]
  alpha  <-  parameters[[PFAS_names[i]]][4+temp_iter]
  
  inits <- c( "Cw" = C_water, "C_lumen" = 0, "C_daphnia_unbound" = 0,
              "C_daphnia_bound" = 0)
  params <- c("init_age"=age, "Temp" = Temp, "ku"= ku, 
              "ke"  = ke,  "kon" = kon, "koff" = koff)
  solutions[[PFAS_names[i]]] <- data.frame(deSolve::ode(times = sol_times,  func = ode_func,
                                                        y = inits,
                                                        parms = params,
                                                        method="lsodes",
                                                        rtol = 1e-5, atol = 1e-5))
  
  
  plot_func(params = parameters[[PFAS_names[i]]], PFAS_data = data_plot, PFAS_name  = PFAS_names[i], 
            Cwater = Cwater, age = age,  temperatures = temperatures )
}

