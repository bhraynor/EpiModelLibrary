#Brinkley Raynor
#March 25, 2020
#Epidemic SEIR model, density dependent

###############################################################################################################
#Notes
###############################################################################################################
#Density dependent assumes that contact rate (c) increases w/ density
#This is a good assumption for smaller populations where as more individuals are added, more contacts occur
#Examples: Small herds, small populations
#Implications: Culling (decreasing density) animal population is a useful control strategy

##Compartmental model schematic

#        _____            _____          _____           _____
#       |     | BI/N     |     |    g   |     |    d    |     |
#       |  S  | ------>  |  E  | ------>|  I  | ------> |  R  |
#       |_____|          |_____|        |_____|         |_____|
#                                          |
#                                          |a
#                                          V  
#Definitions
#B= transmission constant
#a= per capita rate of deat attributable to infection
#g= per capita, instantaneous incubation rate
#d= per capita rate of recovery (reciprocal of time spent in box I)
#N= Total population number (S+I+R)


###############################################################################################################
#Load libraries
###############################################################################################################
library(deSolve)
library(ggplot2)

###############################################################################################################
#Parameterize
###############################################################################################################
init <- c(S=99, I=1, R=0)  # popn initial values
parameters <- c(B= 0.3,
                a=0.0001,
                d=0.1
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
      dS <- - B*S*I/(S+I+R)
      dI <- B*S*I/(S+I+R) - a*I - d*I
      dR <- d*I
      
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

