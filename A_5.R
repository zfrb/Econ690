#==============================================================
# Assignment 5: Random Preference Model

# by Bahar Zafer
#==============================================================

# simulation methods
# expected values same, variance increases.
# use dat_choices

rm(list=ls())

### Writing crra utility function 

crra = function(c, theta) 
      {
        ifelse(abs(theta-1)<0.01, log(c), c^(1-theta)/(1-theta))
      }

### matrix of Table 1 

table1 = rbind(c(48,48,40,64), c(40,64,32,80), c(32,80,24,96), c(24,96,16,112), c(16,112,8,120),
               c(48,48,42,66), c(42,66,36,84), c(36,84,30,102), c(30,102,24,120), c(24,120,16,128),
               c(48,48,38,62), c(38,62,28,76), c(28,76,18,90), c(18,90,8,104), c(8,104,0,112),
               c(42,42,36,60), c(36,60,30,78), c(30,78,24,96), c(24,96,18,114), c(18,114,10,122),
               c(54,54,44,68), c(44,68,34,82), c(34,82,24,96), c(24,96,14,110), c(14,110,6,118))

colnames(table1) = c("L1P1","L1P2","L2P1","L2P2")

# =========================================================================#
# Exercise 1 Part 1 
# =========================================================================#

### Define functions to be used in the simulation step
value = function(param, theta)
      {
        out = 0.5*crra(param[1], theta) + 0.5*crra(param[2], theta)
        return(out)
      }

delta = function(lottery, theta)
      {
        C1 = 1:2
        C2 = 3:4
        value(table1[lottery,C1], theta) - value(table1[lottery,C2], theta)
      }

frequency_simulator = function(theta_f, lottery_f)
      {
      choice_prob = vector(length = ncol(theta_f))  
      num = ncol(theta_f)
      
      for(j in 1:num)
          {
            simulated_prob = vector()
      
            for(i in 1:steps)
                  {
                    t = theta_f[i,j]
                    if(delta(lottery = lottery_f, theta = t)>=0) {simulated_prob[i]=1}
                    else {simulated_prob[i]=0}
                  }
                  choice_prob[j] = mean(simulated_prob)
          }
          return(choice_prob)  
      }

### Plot the relationship between choice probability and risk aversion

plot_lottery = function(lottery, grid_risk, st_dev)
      {
        # Simulation Step
        thetas = sapply(X=grid_risk, FUN=rnorm, n=steps, sd=st_dev)
        pr = frequency_simulator(theta_f = thetas, lottery_f = lottery)
        
        plot(grid_risk, pr, xlab = "mean risk aversion", ylab = "choice probability")
      }

# 1. Assume that the scale parameter is known 
  # for a grid of risk aversion parameters report the relationship between the choice probability and the risk aversion.

# Initialization Step
sig = 0.5                       # standard deviation of risk aversion
grid_theta = seq(0.1,2,0.01)    # initial mean values of risk aversion
steps = 1000                    # number of simulations
select_lottery = 10             # draw for all lotteries, {1, 2, ..., 25}

plot_lottery(lottery = select_lottery, grid_risk = grid_theta, st_dev = sig) # draws the probability of Choice 1 (C1)


# 2. Implement the frequency simulator and estimate the risk aversion for individuals 120, 280 and 1200.

setwd("C:\\Users\\bahar\\Desktop\\Econ690_Computational_Economics\\midterm")
library(nloptr)
library(haven)
dat_choices = read_dta("dat_choices.dta") # rows are individuals. columns are choices.
ind120 = dat_choices[120,]
ind280 = dat_choices[280,]
ind1200 = dat_choices[1200,]


# Initialization Step
sig = 0.5                       # initial value of standard deviation
grid_theta = seq(0.1,2,0.1)     # initial mean values of risk aversion
steps = 1000

# Likelihood Step
Llike = function(p) 
{
  r = p[1]
  
  w = 20 # wealth
  nc = 4
  pr   = NULL;
  for (iX in 1:nrow(table1))
        {
          expU    = crra(w+table1[iX,],rep(r,nc))
          # Simulation Step
          thetas = sapply(X=grid_theta, FUN=rnorm, n=steps, sd=sig)
          pr = frequency_simulator(theta_f = thetas, lottery_f = iX)
        }
  pr[pr<0.00001] = 0.00001;
  pr[pr>0.99999] = 0.99999;
  logLik = ydata * log(pr) + (1-ydata)*log(1-pr);
  return(-sum(logLik));  # because algorithms do minimization and we want to maximization.
}

# individual 120
ydata = ind120
B120_theta = bobyqa(c(0.1), Llike, lower = c(0), upper = c(5))
L120_theta = lbfgs(c(0.1), Llike, lower = c(0), upper = c(5))

# individual 280
ydata = ind280
B280_theta = bobyqa(c(0.1), Llike, lower = c(0), upper = c(5))
L280_theta = lbfgs(c(0.1), Llike, lower = c(0), upper = c(5))

# individual 1200  
ydata = ind1200
B1200_theta = bobyqa(c(0.1), Llike, lower = c(0), upper = c(5))
L1200_theta = lbfgs(c(0.1), Llike, lower = c(0), upper = c(5))

