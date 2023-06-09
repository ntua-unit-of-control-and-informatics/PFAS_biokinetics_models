# Working directory
setwd('C:/Users/ptsir/Desktop/projects/NanoSolveIt/polystyrene')

################################################################################
# Load the data for concentration = 0.1 mg/l
mode1 <- read.csv('mode1.csv')
mode2 <- read.csv('mode2.csv')
mode3 <- read.csv('mode3.csv')
mode4 <- read.csv('mode4.csv')

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

# ode_func(): the differential equation system that deiscribes the model

ode_func <- function(time, inits, params){
  with(as.list(c(inits, params)),{
    
    # Units explanation:
    # A_daphnia: daphnia fluorescence intensity
    # ku: intensity/h
    # ke: 1/h
    
    # For 1 hour uptake experiment
    if (time<1){
        dA_daphnia <- ku - ke*A_daphnia  
    }else{
      dA_daphnia <- -ke*A_daphnia  
    }
    
    return(list(c( dA_daphnia)))
  })
}



obj_func <- function(x, list_of_experiments, metric){
  score_per_experiment <- c()
  lag <- x[1]
  for (k in 0:(length(list_of_experiments)-1)) { # loop over the experiments from different papers
    params_to_fit <- x[(k*4+2):(k*4+5)]
    exp_data <- list_of_experiments[[k+1]]
    exp_time <- exp_data[,1]
    score_per_type <- c()

    #loop over PS50 and PS500
    for (j in 1:2){
      ku <- params_to_fit[2*j-1];ke = params_to_fit[2*j]
      intensity <- exp_data[,1+j]
      sol_times <- seq(0,7, 0.01)
      inits <- c('A_daphnia'= intensity[1])
      params <- c("ku"=ku, "ke"=ke, "lag"= lag)
      solution <- data.frame(deSolve::ode(times = sol_times,  func = ode_func,
                                            y = inits,
                                            parms = params,
                                            method="lsodes",
                                            rtol = 1e-5, atol = 1e-5))
        
        if(sum(solution$time %in% exp_time) == dim(exp_data)[1]){
          results <- solution[which(solution$time %in% exp_time), 'A_daphnia']
        } else{
          stop(print("Length of predictions is not equal to the length of data"))
        }
        if(metric == "AAFE"){
          score_per_type[j] <- AAFE(results, exp_data[,j+1]) 
        }else if (metric =="rmse"){
          score_per_type[j] <- rmse(exp_data[,j+1], results)
        }else if(metric == "PBKOF")
          score_per_type[j] <- PBKOF(list(exp_data[,j+1]), list(results))
      }
    score_per_experiment[k+1] <- mean(score_per_type)
  }
  
  return(mean(score_per_experiment))
  
}

plot_func <- function(optimization, list_of_experiments){
  library(ggplot2)
  x <- optimization$solution
  plots <- list()
  
  for (k in 1:length(list_of_experiments)-1) { # loop over the experiments from different papers
    params_to_fit <- x[(k*4+2):(k*4+5)]
    exp_data <- list_of_experiments[[k+1]]
    exp_time <- exp_data[,1]
    sol_times <- seq(0,7, 0.01)
    
    keep_predictions <- data.frame(matrix(NA, nrow = length(sol_times), ncol = 3))
    keep_predictions[,1] <- sol_times
    colnames(keep_predictions) <- c('Time', 'PS50', 'PS500')
      #loop over PS50 and PS500
      for (j in 1:2){
        ku <- params_to_fit[2*j-1];ke = params_to_fit[2*j]
        intensity <- exp_data[,1+j]
        inits <- c('A_daphnia'= intensity[1])
        params <- c("ku"=ku, "ke"=ke, "lag"= lag)
        solution <- data.frame(deSolve::ode(times = sol_times,  func = ode_func,
                                            y = inits,
                                            parms = params,
                                            method="lsodes",
                                            rtol = 1e-5, atol = 1e-5))
        keep_predictions[,j+1] <- solution$A_daphnia
        
      }
      
    
    strings <- c("PS50", "PS500")
    color_codes <- scales::hue_pal()(2) # to return 3 color codes 
    cls <- c()  
    for (i in 1:length(strings)) {
      cls[i] <- color_codes[i]
      names(cls)[i] <- strings[i]
    }
     
      draw_plot <- ggplot()+
        geom_line(data = keep_predictions, aes(x=Time, y=PS50, color=strings[1]), size=1.7)+
        geom_line(data = keep_predictions, aes(x=Time, y=PS500, color=strings[2]), size=1.7)+

        geom_point(data = exp_data, aes(x=Time, y=PS50, color=strings[1]), size=5)+
        geom_point(data = exp_data, aes(x=Time, y=PS500, color=strings[2]), size=5)+
        #scale_y_log10()+
        
        
        labs(title = paste0("Fluorescence intensity in D.magna of mode ", k+1),
             y = "Fluorescence intensity", x = "Time (hours)")+
        theme(plot.title = element_text(hjust = 0.5,size=30), 
              axis.title.y =element_text(hjust = 0.5,size=25,face="bold"),
              axis.text.y=element_text(size=22),
              axis.title.x =element_text(hjust = 0.5,size=25,face="bold"),
              axis.text.x=element_text(size=22),
              legend.title=element_text(hjust = 0.5,size=25), 
              legend.text=element_text(size=22)) + 
        
        scale_color_manual("PS type", values=cls)+
        theme(legend.key.size = unit(1.5, 'cm'),  
              legend.title = element_text(size=14),
              legend.text = element_text(size=14),
              axis.text = element_text(size = 14))
      print(draw_plot)
      
      
    }
  
  
  
}
##############################################################################


list_of_experiments <- list("mode1" = mode1,"mode2" = mode2,
                             "mode3" = mode3, "mode4" = mode4)

x0 <- c(0.5, rep(c(5000,0.1), 8))
opts <- list( "algorithm" = "NLOPT_LN_SBPLX", #"NLOPT_LN_NEWUOA", #"NLOPT_LN_SBPLX" , #"NLOPT_LN_BOBYQA" #"NLOPT_LN_COBYLA"
              "xtol_rel" = 1e-07, 
              "ftol_rel" = 1e-07,
              "ftol_abs" = 0.0,
              "xtol_abs" = 0.0 ,
              "maxeval" = 3000,
              "print_level" = 1)
metric <- "PBKOF"
optimization <- nloptr::nloptr(x0 = x0,
                               eval_f = obj_func,
                               lb	= c(0,rep(0,16)),
                               ub = c(1,rep(20000,16)),
                               opts = opts,
                               list_of_experiments=list_of_experiments,
                               metric = metric)

fitted_params <- optimization$solution
lag <- fitted_params[1]
params_mode1 <- fitted_params[2:5]
params_mode2 <- fitted_params[6:9]
params_mode3 <- fitted_params[10:13]
params_mode4 <- fitted_params[14:17]
parameters <- data.frame(rbind(params_mode1,params_mode2,
                               params_mode3, params_mode4))
colnames(parameters) <- c("ku_PS50", "ke_PS50", "ku_PS500", "ke_PS500")
rownames(parameters) <- c("mode1", "mode2", "mode3", "mode4")


plot_func(optimization, list_of_experiments)

