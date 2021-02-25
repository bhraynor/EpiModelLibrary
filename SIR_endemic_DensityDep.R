#Brinkley Raynor
#March 25, 2020
#Endemic SIR model, frequency dependent

###############################################################################################################
#Notes
###############################################################################################################
#Frequancy dependent assumes that contact rate (c) remeains the same despite density
#Examples: Very large popns were mixing between all individuals not possible, STDs, Tasmanian devel tumors
#Implications: Culling (decreasing density) is not useful

#Endemic models: easist to model closed populations were birth rate = overall death rate (sum derivatives)

##Compartmental model schematic

#         |
#         |m*N + a*I
#         V
#        _____            _____            _____
#       |     |   BI     |     |    d     |     |
#       |  S  | ------>  |  I  | ------>  |  R  |
#       |_____|          |_____|          |_____|
#          |              |   |              |
#          |m             |m  |a             |m
#          V              V   V              V

#Definitions
#L= force of infection (incidence) = c*p*r = B*I, where: 
#c= contact rate (# of hosts individuals encounter/ time unit)
#r= probability of contact being w/ an infectious host
#p= porbability that transmission occurs given contact w/ infected host
#B= transmission constant
#a= per capita rate of deat attributable to infection
#d= per capita rate of recovery (reciprocal of time spent in box I)
#m= per capita normal death rate


###############################################################################################################
#Load libraries
###############################################################################################################
library(deSolve)
library(ggplot2)

###############################################################################################################
#Parameterize
###############################################################################################################
init <- c(S=99, I=1, R=0)  # popn initial values
parameters <- c(B= 0.003,
                a=0.0001,
                d=0.1,
                m=0.001
)

times <- seq(0, 100, by = 1) #unit=days, seq(start day, end day, step)

###############################################################################################################
#Model
###############################################################################################################
#Model function
#model 1: mu simple
disease_dynamics <- function(time, state, parameters) {
  with(
    as.list(c(state, parameters)), {
      #Equations
      dS <- - B*S*I - m*S + m*(S+I+R) + a*I
      dI <- B*S*I - a*I - d*I - m*I
      dR <- d*I - m*R
      
      return(list(c(dS, dI, dR)))
    }
  )
}


#Run model
out <- as.data.frame(ode(y = init, times = times, func = disease_dynamics, parms = parameters))
RESULTS<-data.frame(out$S, out$I, out$R) #results

###############################################################################################################
#Plot
###############################################################################################################
#Real plot
p <- ggplot()+
  theme_classic()+ #white background w/ out grid
  geom_line(data=RESULTS, aes(x=times, y=RESULTS$out.S), color="gold")+ 
  geom_line(data=RESULTS, aes(x=times, y=RESULTS$out.I), color="red")+ 
  geom_line(data=RESULTS, aes(x=times, y=RESULTS$out.R), color="blue")+ 
  ggtitle("Disease dynamics") + #title label
  labs(y= "Number of hosts", x = "Time (days)") #axis labels
p

