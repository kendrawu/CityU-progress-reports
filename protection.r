library(deSolve)
library(EpiEstim)

# Parameters
vaccine1 = 0.5
vaccine2 = 0.7
vaccine3 = 0.9

# Initialization
N=1000
init       <- c(S1 = N/6-1, S2 = N/6, S3 = N/6, S4 = N/6, S5 = N/6, S6 = N/6, E = 0, I = 1, R = 0.0)
parameters <- c(beta1 = 0.003, beta2 = 0.008, gamma = 0.03, sigma = 0.01)
t = 1000
times      <- seq(0, t, by = 1)


sir <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    
    dS1 <- -beta1 * S1 * I * (1-vaccine1) * 0.5
    dS2 <- -beta1 * S2 * I * (1-vaccine2) * 0.5
    dS3 <- -beta1 * S3 * I * (1-vaccine3) * 0.5
    dS4 <- -beta2 * S4 * I * (1-vaccine1) * 0.5
    dS5 <- -beta2 * S5 * I * (1-vaccine2) * 0.5
    dS6 <- -beta2 * S6 * I * (1-vaccine3) * 0.5
    dE  <- beta1*S1*I*(1-vaccine1)/2+beta1*S2*I*(1-vaccine2)/2+beta1*S3*I*(1-vaccine3)/2+beta2*S4*I*(1-vaccine1)/2+beta2*S5*I*(1-vaccine2)/2+beta2*S6*I*(1-vaccine3)/2- sigma * E
    dI  <- sigma * E - gamma * I
    dR  <- gamma * I
    
    return(list(c(dS1, dS2, dS3, dS4, dS5, dS6, dE, dI, dR)))
  })
}

## Solve using ode (General Solver for Ordinary Differential Equations)
out <- ode(y = init, times = times, func = sir, parms = parameters)
## change to data frame
out <- as.data.frame(out)
## Delete time variable
out$time <- NULL
## Show data
head(out, 10)

## Plot
matplot(x = times, y = out, type = "l",
        xlab = "Time", ylab = "Number of people", main = "SEIR Model",
        lwd = 1, lty = 1, bty = "l", col = 2:12)

## Add legend
legend(40, 0.7, c("Susceptible", "Exposed", "Infected", "Recovered"), pch = 1, col = 2:6, bty = "n")

## Estimate reproductive number based on incidence data
result <- estimate_R(incid=out[50:(t+1),5],method="parametric_si",config=make_config(list(mean_si=5.4,std_si=1.5)))
plot(result$R[,3], xlab="Days", ylab="Reproductive number")