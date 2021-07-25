#packages
require(deSolve)
require(reshape2)
require(ggplot2)

N = 5.1 * 10^7 #population of 51 million in South Korea
p = 0.70 #proportion of population that is vaccinated

#at t = 250, S = 48271805, I = 499, R = 2727696
init_values = c(S = (1-p)*48271805, I = 499, R = p*48271805 + 2727696) #initial state values

parameters = c(beta = 0.5520629, gamma = 0.4721594)

times = seq(from = 1, to = 44, by = 1) #from 250 to 294 days

sir_model = function(time, states, parameters) {
  with(as.list(c(states, parameters)), {
    
    #these parameters are in units of days^-1 (i.e. per day)
    
    b = 1 / (83.21*365) #birth rate
    mu = 1 / (83.21*365) #mortality rate (1 / average lifespan in days)
    
    N = S+I+R #total population size
    
    lambda = beta * I/N
    
    dS = b*N - lambda*S - mu*S
    
    dI = lambda*S - gamma*I - mu*I
    
    dR = gamma*I - mu*R
    
    return(list(c(dS, dI, dR)))
  })
}

output = as.data.frame(ode(y = init_values, times = times, func = sir_model, parms = parameters))

plot = ggplot(data = output, aes(x = time, y = I)) +
  geom_line(color = "blue") +
  xlab("Time (Days)") +
  ylab("Number of Infected (People)") +
  labs(color = "Compartments") +
  labs(title = "70% Vaccination Coverage After 250 Days") +
  theme(legend.position = "bottom")

plot
