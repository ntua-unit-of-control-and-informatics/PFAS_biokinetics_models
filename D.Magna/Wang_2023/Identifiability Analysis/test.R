  PFAS_data = data_plot;
  PFAS_name  = "PFBA"; 
  Cwater = Cwater; age = age;  temperatures = temperatures; MW = MW 
  
  ku <-  1
  kon <- 5
  Ka <- 2
  ke <- 3
  C_prot_init <- -1
  pars <-c(ku, kon, Ka, ke, C_prot_init)
  
  
  # Age of D.magna at beginning of exposure
  init_age <- age
  # Create a counter to mark the position of fitted parameters of x
  # that corresponds to a specific combination of PFAS and temperature
  
  # Load PFAS data
  df <- PFAS_data[[PFAS_name]]
  # Load parameters
  params <- pars
  # Time used by numerical solver that integrates the system of ODE
  sol_times <- seq(0,15, 0.1 )
  # Data frame to store predictions for each temperature
  predictions <- data.frame("time" = sol_times, "BB_16" = rep(NA, length(sol_times)),
                            "BB_20" = rep(NA, length(sol_times)), "BB_24" = rep(NA, length(sol_times)))
  
  "BB_20" = rep(NA, length(sol_times)), "BB_24" = rep(NA, length(sol_times)))
  ku <- unname(params[1])
  kon <-  unname(params[2])
  Ka <-  unname(params[3])
  ke <- unname(params[4])
  C_prot_init <- 10^unname(params[5])
  # Iterate over number of distinct temperature used in the experiment
  for (temp_iter in 1:length(temperatures)){
    # Initial water concentration of PFAS at selected temperature
    C_water <-  Cwater[PFAS_name,temp_iter]
    # Temperature of experiment
    Temp <- temperatures[temp_iter]
    # Fitted parameters
    inits <- c( "Cw" = C_water,  "C_daphnia_unbound" = 0,
                "C_daphnia_bound" = 0, "C_prot_un" = C_prot_init)
    
    params <- c("init_age"=age, "Temp" = Temp, "ku"= ku, 
                "kon" = kon, "Ka" = Ka, "ke"= ke, "MW" = MW)
    solution <- data.frame(deSolve::ode(times = sol_times,  func = ode_func,
                                        y = inits,
                                        parms = params,
                                        method="lsodes",
                                        rtol = 1e-7, atol = 1e-7))
    
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
  
  
  
