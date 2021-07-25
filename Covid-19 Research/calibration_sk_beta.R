#packages
require(deSolve)
require(reshape2)
require(ggplot2)

setwd("/Users/brandonpae/Desktop/Covid-19 Research") #set working directory

reported_data <- read.csv("SK Case Data.csv") #reported COVID-19 cases in South Korea

sir_model <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    
    #these parameters are in units of days^-1 (i.e. per day)
    
    b = 1 / (83.21*365) #birth rate
    mu = 1 / (83.21*365) #mortality rate (1 / average lifespan in days)
    
    N <- S+I+R #susceptible, infected, removed people
    
    lambda <- beta * I/N #beta is the rate of infection; lambda = force of infection
    
    #the basic differential equations: how each compartment changes
    
    dS <- b*N - lambda*S - mu*S
    
    dI <- lambda*S - gamma*I - mu*I
    
    dR <- gamma*I - mu*R
    
    return(list(c(dS, dI, dR))) 
  })
}

likelihood_calibration <- function(func, parameters, dat) {
  beta <- parameters[1] #first parameter value
  gamma <- parameters[2] #second parameter value
  result = as.data.frame(ode(y = init_values, times = times, func = func, parms = c(beta = beta, gamma = gamma)))
  
  #calculates how well a certain SIR model fits the South Korea data
  LL <- sum(dpois(x = dat$number_reported, lambda = 1 * result$I[result$time %in% dat$time], log = TRUE))
  return(LL)
}

#with a population of 51 million in South Korea, these are the starting values for each compartment
init_values <- c(S = 5.1*10^7 - 19, I = 19, R = 0)
times = seq(from = 1, to = 294, by = 1)

#calculate the optimal values for beta and gamma
optimised <- optim(par = c(0.5, 0.5),
                   fn  = likelihood_calibration,
                   dat = reported_data,
                   control = list(fnscale=-1),
                   func = sir_model)

optimised

parameters = c(beta = 2.533772, gamma = 2.503224) #optimal values

#ode() to solve the differential equations at each day (find the state values)
opt_mod1 = as.data.frame(ode(y = init_values, times = times, func = sir_model, parms = parameters))

plot = ggplot() +
  geom_line(data = opt_mod1, aes(x = time, y = I), colour = "blue") + #our model
  geom_point(data = reported_data, aes(x = time, y = number_reported), colour = "red", shape = "x") + #data
  xlab("Time (Days)") +
  ylab("Number of Infected (People)") +
  labs(color = "Compartments") +
  labs(title = "Prevalence of Infected People in South Korea") + #add a title
  theme(legend.position = "bottom")

plot #show plot
