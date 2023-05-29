i=1

x0 <- c(8, 5, 5, 7, 1e-04)

ku <- 1
kon <-  9
Ka <- 1
ke <- 
C_prot_init <- 1e-08

inits <- c( "Cw" = C_water,  "C_zebrafish_unbound" = 0,
            "C_zebrafish_bound" = 0, "C_prot_un" = C_prot_init)
params <- c("init_age"=age, "Temp" = Temp, "ku"= ku, 



PFAS_data = data_plot; PFAS_name  = PFAS_names[i]; 
          Cwater = Cwater; age = age;  temperatures = temperatures; MW = MW




df <- PFAS_data[[PFAS_name]]
# Load parameters
# Time used by numerical solver that integrates the system of ODE
sol_times <- seq(0,29, 0.1 )
# Data frame to store predictions for each temperature
predictions <- data.frame("time" = sol_times, "BB_16" = rep(NA, length(sol_times)),
                          "BB_20" = rep(NA, length(sol_times)), "BB_24" = rep(NA, length(sol_times)))

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