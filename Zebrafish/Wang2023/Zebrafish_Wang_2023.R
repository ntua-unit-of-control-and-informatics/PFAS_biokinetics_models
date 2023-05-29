# Working directory

#etwd("/Users/elenistrompoula/Documents/GitHub/PFAS_biokinetics_models/Zebrafish/Wang2023")
setwd("C:/Users/ptsir/Documents/GitHub/PFAS_biokinetics_models/Zebrafish/Wang2023")


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
    # Zebrafish: mol PFAS/L zebrafish
    # ke: 1/h
    # Cw: mol/L
    
    
    # Water concentration in ng/L
    dCw <- 0
    # Zebrafish concentration in lumenconcentration in ng/g
    ku <- 10^ku
    ke <- 10^ke
    kon <- 10^kon
    Ka <- 10^Ka
    koff <- kon/Ka
    # Reported values for Kon and koff range between 1e2-1e04 and 1e-3-1e-1 in L/mol/s 
    # and s-^-1 respectively. Our concentrations are in ng/g. Assuming density = 1000g/L
    # then the concentration in ng/g is multyplied by 1000 to make it ng/L and then by
    # multiply by 1e-9 to make it grams and divide by MW. We do this directly to kon 
    # and koff to make the concentration mol/L so that we can match the literature 
    # values of kon and koff with ours
    C_zebrafish_unbound_unmol <- C_zebrafish_unbound*MW/(1000*1e-09)
    C_zebrafish_bound_unmol <- C_zebrafish_bound*MW/(1000*1e-09)
    dC_zebrafish_unbound <-  ku*(Cw*1e-09/MW)  - kon*C_prot_un*C_zebrafish_unbound +   koff*C_zebrafish_bound - ke*C_zebrafish_unbound
    dC_zebrafish_bound <- kon*C_prot_un*C_zebrafish_unbound - koff*C_zebrafish_bound
    dC_prot_un <-   koff*C_zebrafish_bound -  kon*C_prot_un*C_zebrafish_unbound
    C_tot <- C_zebrafish_unbound_unmol + C_zebrafish_bound_unmol
    return(list(c("dCw" = dCw,   "dC_zebrafish_unbound" = dC_zebrafish_unbound,
                  "dC_zebrafish_bound" = dC_zebrafish_bound, "dC_prot_un" = dC_prot_un), 
             "C_tot" = C_tot))
  })
}



obj_func <- function(x, PFAS_data, PFAS_name, Cwater, age, temperatures, MW, metric){
  
  # Indexes of body burden and exposure time in data frame
  BB_index <- c(2,4,6)
  ExpTime_index <- c(1,3,5)
  # Age of zebrafish at beginning of exposure
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
  kon <- x[2]
  Ka <- x[3]
  ke <- x[4]
  C_prot_init <- x[5]
  for (temp_iter in 1:length(temperatures)){
    # Initial water concentration of PFAS at selected temperature
    C_water <-  Cwater[PFAS_name,temp_iter]
    # Temperature of experiment
    Temp <- temperatures[temp_iter]
    # Time of measurement of selected PFAS at selected temperature
    exp_time <- round(df[!is.na(df[,ExpTime_index[temp_iter]]),ExpTime_index[temp_iter]],1)
    # Body burden of selected PFAS at selected temperature
    BodyBurden <- df[!is.na(df[,BB_index[temp_iter]]),BB_index[temp_iter]]
    # Time used by numerical solver that integrates the system of ODE
    sol_times <- seq(0,29, 0.1 )
    # Fitted parameters
    
    
    inits <- c( "Cw" = C_water,  "C_zebrafish_unbound" = 0,
                "C_zebrafish_bound" = 0, "C_prot_un" = C_prot_init)
    
    params <- c("init_age"=age, "Temp" = Temp, "ku"= ku, 
                "kon" = kon, "Ka" = Ka, "ke"= ke, "MW" = MW)
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



plot_func <- function(params,PFAS_data, PFAS_name, Cwater, age, temperatures,MW){
  library(ggplot2)
  #setwd("/Users/elenistrompoula/Documents/GitHub/PFAS_biokinetics_models/Zebrafish/Wang2023/Results")
  setwd("C:/Users/ptsir/Documents/GitHub/PFAS_biokinetics_models/Zebrafish/Wang2023/results")
        
  # Age of D.magna at beginning of exposure
  init_age <- age
  # Create a counter to mark the position of fitted parameters of x
  # that corresponds to a specific combination of PFAS and temperature
  
  # Load PFAS data
  df <- PFAS_data[[PFAS_name]]
  # Load parameters
  parameters <- params
  # Time used by numerical solver that integrates the system of ODE
  sol_times <- seq(0,29, 0.1 )
  # Data frame to store predictions for each temperature
  predictions <- data.frame("time" = sol_times, "BB_16" = rep(NA, length(sol_times)),
                            "BB_20" = rep(NA, length(sol_times)), "BB_24" = rep(NA, length(sol_times)))
  ku <- parameters[1]
  kon <-  parameters[2]
  Ka <-  parameters[3]
  ke <- parameters[4]
  C_prot_init <- unname(parameters[5])
  # Iterate over number of distinct temperature used in the experiment
  for (temp_iter in 1:length(temperatures)){
    # Initial water concentration of PFAS at selected temperature
    C_water <-  Cwater[PFAS_name,temp_iter]
    # Temperature of experiment
    Temp <- temperatures[temp_iter]
    # Fitted parameters
    
    inits <- c( "Cw" = C_water,  "C_zebrafish_unbound" = 0,
                "C_zebrafish_bound" = 0, "C_prot_un" = C_prot_init)
    
    params <- c("init_age"=age, "Temp" = Temp, "ku"= ku, 
                "kon" = kon, "Ka" = Ka, "ke"= ke, "MW" = MW)
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
    geom_point(data = df, aes(x=TIME16  , y=BB16 ,  colour = "16oC"), size=5)+
    geom_point(data = df, aes(x=TIME20, y=BB20,  colour = "20oC"), size=5)+
    geom_point(data = df, aes(x=TIME24, y=BB24,  colour = "24oC"), size=5)+
    
    labs(title = paste(PFAS_name,"body burden in zebrafish", sep = " "),
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
  ggsave(paste0(PFAS_name,".png"), width=15, height=10)
}

##############################################################################
################################################################################
# Load the data for PFAS concentration

sheet_names <- c("PFBA", "F-53B", "GenX", "PFBS", "PFDA", "PFDoA", "PFHpA", 
                 "PFHxA", "PFNA", "PFOA", "PFOS", "PFPeA", "PFUnA")
PFAS_names <- c("PFBA", "F-53B", "GenX", "PFBS", "PFDA", "PFDoA", "PFHpA", 
                "PFHxA", "PFNA", "PFOA", "PFOS", "PFPeA", "PFUnA")
Molecular_weights <- list("PFBA" = 214, "F-53B" = 570, "GenX" = 330, "PFBS" = 300, "PFDA" = 514, 
                          "PFDoA" = 614, "PFHpA" = 364, "PFHxA" = 314,"PFNA" = 464, "PFOA" = 414,  "PFOS" = 500,  
                          "PFPeA" = 364, "PFUnA" = 564)
data_ls <- list()
data_plot <- list()

for(sheet_name in sheet_names){
  data_ls[[sheet_name]] <- openxlsx::read.xlsx ('Wang_Data.xlsx', sheet = sheet_name)
  data_plot[[sheet_name]] <- openxlsx::read.xlsx ('Wang_Data.xlsx', sheet = sheet_name)
}

opts <- list( "algorithm" = "NLOPT_LN_SBPLX",#"NLOPT_LN_SBPLX", #"NLOPT_LN_NEWUOA", #"NLOPT_LN_SBPLX" , #"NLOPT_LN_BOBYQA" #"NLOPT_LN_COBYLA"
              "xtol_rel" = 1e-08, 
              "ftol_rel" = 1e-08,
              "ftol_abs" = 0.0,
              "xtol_abs" = 0.0 ,
              "maxeval" = 3000,
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
  MW <- Molecular_weights[[PFAS_names[i]]]
  # Define initial values of fitted parameters to provide to the optimization routine
  # For each PFAS and temperature combination we have two parameters
  x0 <- c(8, 5, 5, 7, 1e-04)
  #x0 <- c(7, 6, 6, 8, 5e-05)# For PFBS
  set.seed(12312)
  optimization<- nloptr::nloptr(x0 = x0,
                                eval_f = obj_func,
                              lb	=  c(-5,-5,-5,-5, 1e-08),
                               ub =   c(15, 15, 15, 15,  1e-03),
                                opts = opts,
                                PFAS_data = data_ls,
                                PFAS_name = PFAS_names[i],
                                Cwater = Cwater,
                                age = age ,
                                temperatures = temperatures,
                                MW = MW,
                                metric = "PBKOF")
  optimizations[[PFAS_names[i]]] <- optimization
  parameters[[PFAS_names[i]]] <- optimization$solution
  names(parameters[[PFAS_names[i]]]) = c("ku",  "kon","Ka", "ke", "C_prot_init")
  
  sol_times <- seq(0,29, 0.01 )
  # Iterate over number of distinct temperature used in the experiment
  temp_iter <- 2
  # Initial water concentration of PFAS at selected temperature
  C_water <-  Cwater[PFAS_names[i],temp_iter]
  # Temperature of experiment
  Temp <- temperatures[temp_iter]
  # Fitted parameters
  ku <- parameters[[PFAS_names[i]]][1]
  kon <-  parameters[[PFAS_names[i]]][2]
  Ka <-  parameters[[PFAS_names[i]]][3]
  ke <- parameters[[PFAS_names[i]]][4]
  C_prot_init <- unname(parameters[[PFAS_names[i]]][5])
  
  inits <- c( "Cw" = C_water,  "C_zebrafish_unbound" = 0,
                              "C_zebrafish_bound" = 0, "C_prot_un" = C_prot_init)
  params <- c("init_age"=age, "Temp" = Temp, "ku"= ku, 
              "kon" = kon, "Ka" = Ka, "ke"= ke, "MW" = MW)
  solutions[[PFAS_names[i]]] <- data.frame(deSolve::ode(times = sol_times,  func = ode_func,
                                                        y = inits,
                                                        parms = params,
                                                        method="lsodes",
                                                        rtol = 1e-5, atol = 1e-5))
  
  
  plot_func(params = parameters[[PFAS_names[i]]], PFAS_data = data_plot, PFAS_name  = PFAS_names[i], 
            Cwater = Cwater, age = age,  temperatures = temperatures, MW = MW )
}

