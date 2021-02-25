#Brinkley Raynor
#March 25, 2020
#Endemic SEIR model, density dependent

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
#        _____            _____          _____           _____
#       |     |   BI     |     |    g   |     |    d    |     |
#       |  S  | ------>  |  E  | ------>|  I  | ------> |  R  |
#       |_____|          |_____|        |_____|         |_____|
#          |              |   |            |               |
#          |m             |m  |a           |m              |m
#          V              V   V            V               V

#Definitions
#L= force of infection (incidence) = c*p*r = B*I, where: 
#c= contact rate (# of hosts individuals encounter/ time unit)
#r= probability of contact being w/ an infectious host
#p= porbability that transmission occurs given contact w/ infected host
#B= transmission constant
#a= instantaneous, per capita rate of deat attributable to infection
#d= instantaneous, per capita rate of recovery (reciprocal of time spent in box I)
#m= instantaneous, per capita normal death rate
#g= instantaneous, per capita rate exposed individuals become infectious

###############################################################################################################
#Load libraries
###############################################################################################################
library(deSolve)
library(ggplot2)

###############################################################################################################
#Parameterize
###############################################################################################################
init <- c(S=99, E=0, I=1, R=0)  # popn initial values
parameters <- c(B= 0.003,
                a=0.0001,
                d=0.1,
                m=0.0001,
                g= 0.2
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
      dS <- m*(S+E+I+R)+a*I- B*I*S - m*S
      dE <- B*S*I - g*E - m*E
      dI <- g*E -a*I - d*I - m*I
      dR <- d*I - m*R
      
      return(list(c(dS, dE, dI, dR)))
    }
  )
}


#Run model
out <- as.data.frame(ode(y = init, times = times, func = disease_dynamics, parms = parameters))
RESULTS<-data.frame(out$S, out$E, out$I, out$R) #results

###############################################################################################################
#Plot
###############################################################################################################
#Real plot
p <- ggplot()+
  theme_classic()+ #white background w/ out grid
  geom_line(data=RESULTS, aes(x=times, y=RESULTS$out.S), color="gold")+ 
  geom_line(data=RESULTS, aes(x=times, y=RESULTS$out.E), color="pink")+ 
  geom_line(data=RESULTS, aes(x=times, y=RESULTS$out.I), color="red")+ 
  geom_line(data=RESULTS, aes(x=times, y=RESULTS$out.R), color="blue")+ 
  ggtitle("Disease dynamics") + #title label
  labs(y= "Number of hosts", x = "Time (days)") #axis labels
p

