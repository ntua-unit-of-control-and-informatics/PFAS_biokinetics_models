# Working directory

#setwd("/Users/elenistrompoula/Documents/GitHub/PFAS_biokinetics_models/D.Magna")
setwd("C:/Users/ptsir/Documents/GitHub/PFAS_biokinetics_models/D.Magna")

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

PFOA <- read.csv('PFOA.csv')
PFOA <- PFOA[2:dim(PFOA)[1],]
PFOA[,1] <-  round(PFOA[,1])

PFOS <- read.csv('PFOS.csv')
PFOS <- PFOS[2:dim(PFOS)[1],]
PFOS[,1] <-  round(PFOS[,1])

PFNA <- read.csv('PFNA.csv')
PFNA <- PFNA[2:dim(PFNA)[1],]
PFNA[,1] <-  round(PFNA[,1])

PFDA <- read.csv('PFDA.csv')
PFDA <- PFDA[2:dim(PFDA)[1],]
PFDA[,1] <-  round(PFDA[,1])

PFUnA <- read.csv('PFUnA.csv')
PFUnA <- PFUnA[2:dim(PFUnA)[1],]
PFUnA[,1] <-  round(PFUnA[,1])

PFDoA <- read.csv('PFDoA.csv')
PFDoA <- PFDoA[2:dim(PFDoA)[1],]
PFDoA[,1] <-  round(PFDoA[,1])



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

AAFE <- function(predictions, observations, times=NULL){
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
    size <- Size_estimation(age, temperature = 20, food="high")
    
    # Filtation rate in mL/h
    Frate <- Filtration_rate_estimation(size, temperature = 20, method = "Preuss")
    # Filtation rate in L/day
    Frate <- Frate * 24/1000
    
    # dry weight mg
    DW <- dry_weight_estimation(size)
    # Convert DW to WW
    WW <- 15 * DW  # Conversion rate 11-20 from DW to WW (Garner et al., 2018)
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



obj_func <- function(x, PFAS, Cwater, age, metric){
     
    # User defined parameters
    init_age <- age
    # Fitted parameters
    Fsorption <- x[1]
    ke <- x[2]
    # Water concentration
    Cw <- Cwater * 1e06 #ng/L
    
    BodyBurden <- PFAS$Concentration
    exp_time <- PFAS$Time
    sol_times <- seq(0,round(max(PFAS$Time))+1, 0.01 )
    inits <- c('C_daphnia'= 0, "Cw" = Cw)
    params <- c("init_age"=init_age, "Fsorption"= Fsorption, "ke"  = ke)
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
      if(metric == "AAFE"){
        score <- AAFE(results,BodyBurden) 
      }else if (metric =="rmse"){
        score <- rmse(BodyBurden, results)
      }else if(metric == "PBKOF"){
        score <- PBKOF(list(BodyBurden), list(results))
    }
  return(score)
}

plot_func <- function(optimization, PFAS, Cwater, age){
  library(ggplot2)
  x <- optimization$solution
  # User defined parameters
  init_age <- age
  # Fitted parameters
  Fsorption <- x[1]
  ke <- x[2]
  # Water concentration
  Cw <- Cwater * 1e06 #ng/L
  exp_time <- PFAS$Time
  sol_times <- seq(0,round(max(PFAS$Time))+1, 0.01 )
  inits <- c('C_daphnia'= 0, "Cw" = Cw)
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

    
  draw_plot <- ggplot()+
      geom_line(data = keep_predictions, aes(x=Time, y=BodyBurden), size=1.7)+
      geom_point(data = PFAS, aes(x=Time, y=Concentration), size=5)+
      labs(title = "PFOS body burden in D.magna",
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
    print(draw_plot)
}
  
##############################################################################

x0 <- c(1e-04, 0.05)
opts <- list( "algorithm" = "NLOPT_LN_SBPLX", #"NLOPT_LN_NEWUOA", #"NLOPT_LN_SBPLX" , #"NLOPT_LN_BOBYQA" #"NLOPT_LN_COBYLA"
              "xtol_rel" = 1e-07, 
              "ftol_rel" = 1e-07,
              "ftol_abs" = 0.0,
              "xtol_abs" = 0.0 ,
              "maxeval" = 1000,
              "print_level" = 1)
optimization <- nloptr::nloptr(x0 = x0,
                               eval_f = obj_func,
                               lb	= c(0,0),
                               ub = c(100,100),
                               opts = opts,
                               metric = "AAFE",
                               PFAS = PFOS,
                               Cwater = 0.005, #mg/L
                               age = 1)

PFOS_params <- c("Fsorption" = optimization$solution[1],  
                 "ke" = optimization$solution[2])


plot_func(optimization, PFAS = PFOS, Cwater = 0.005, age = 1)


# Fitted parameters
Fsorption <- PFOS_params["Fsorption"]
ke <- PFOS_params["ke"]
# Water concentration
Cw <- 0.005 * 1e06 #ng/L
exp_time <- PFOS$Time
sol_times <- seq(0,round(max(PFOS$Time))+1, 0.01 )
inits <- c('C_daphnia'= 0, "Cw" = Cw)
params <- c("init_age"=1, "Fsorption"= Fsorption, "ke"  = ke)
solution <- data.frame(deSolve::ode(times = sol_times,  func = ode_func,
                                    y = inits,
                                    parms = params,
                                    method="lsodes",
                                    rtol = 1e-5, atol = 1e-5))



