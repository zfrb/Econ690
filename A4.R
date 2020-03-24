#====== Assignment 3: Nonlinear Equations and Numerical Optimization ======#
# =========================================================================#
# Bahar Zafer
# 3/3/2020

# =========================================================================#
# Exercise 1 Part 1 
# =========================================================================#
rm(list=ls())
### Writing crra utility function 

#//comments
# perfect

crra = function(c, theta) 
{
  ifelse(abs(theta-1)<0.01, log(c), c^(1-theta)/(1-theta))
}

### The shape of the utility function for c=10 and theta (risk aversion) ranging from -2 to 2. 

theta_grid = seq(-2,2,0.01)
utility = crra(10, theta_grid)

plot(theta_grid,utility)

# =========================================================================#
# Exercise 2 Part 2 
# =========================================================================#

### Crete a matrix that records all choices in Table 1.

table1 = rbind(c(48,48,40,64),
                  c(40,64,32,80),
                  c(32,80,24,96),
                  c(24,96,16,112),
                  c(16,112,8,120),
                  
                  c(48,48,42,66),
                  c(42,66,36,84),
                  c(36,84,30,102),
                  c(30,102,24,120),
                  c(24,120,16,128),
                  
                  c(48,48,38,62),
                  c(38,62,28,76),
                  c(28,76,18,90),
                  c(18,90,8,104),
                  c(8,104,0,112),
                  
                  c(42,42,36,60),
                  c(36,60,30,78),
                  c(30,78,24,96),
                  c(24,96,18,114),
                  c(18,114,10,122),
                  
                  c(54,54,44,68),
                  c(44,68,34,82),
                  c(34,82,24,96),
                  c(24,96,14,110),
                  c(14,110,6,118))

colnames(table1) = c("L1P1","L1P2","L2P1","L2P2")

### Find the level of risk aversion that makes an individual indifferent 
# between each of the 25 choices
# using bisection techniques

# defining a function that gives the expected utility of each choice.
# P1 = P2 = 0.5 for both lottary 1 and 2.

# expected utility function taking a choice of gamble, risk aversion as inputs
# THIS IS VERY CONFUSING NAMING!! 
exp_crra = function(param, theta)
{
  out = 0.5*crra(param[1], theta) + 0.5*crra(param[2], theta)
  return(out)
}

# function of difference between expected utilities
utility_diff = function(gamble, theta)
{
Lot1 = 1:2
Lot2 = 3:4
exp_crra(table1[gamble,Lot1], theta) - exp_crra(table1[gamble,Lot2], theta)
}

# bisection function, given the objective function
my_bisec = function(a, b, gamble, objFunc = utility_diff)
{
  f = function(theta) {utility_diff(gamble, theta)}
  c = (a+b)/2
  tol = 0.00000001
  
  while(f(c)!=0 && b-a>tol)
  {
    if(f(c)==0) {c}
    if(f(a)*f(c)<0) {b=c}
    else if(f(b)*f(c)<0)  {a=c}
    {
      c=(a+b)/2
    }
  }
  return(c) 
}

# vector of risk aversion parameters that make individuals indifferent 
indifference_thetas = vector()

for(num in 1:25)    # there are 25 gambles.
{
  indifference_thetas[num] = my_bisec(-2, 10, gamble = num)
}

indifference_thetas

# THis works well. Well done. 


### Use turning points to construct identified sets for each set of lotteries.

### Use these identified sets along with dat_choice 
# to construct the distribution of risk aversion in the population for each list of questions.

## MISSING


# =========================================================================#
# Exercise 3 Part 3 
# =========================================================================#

# 1. Likelihood function associated to least risky lottery
Llike3 = function(ydata, theta, sig = 1)
{
  #assuming that sigma = 1
  w = 20
  nc = 4
  pr   = NULL;
  for (iX in 1:nrow(table1))
  {
    expU    = crra(w+table1[iX,],rep(theta,nc))
    # THIS IS MY CODE.. YOU have defined exp_CRRA higher
    pr[iX]  = pnorm((expU[1]+expU[2])/2-(expU[3]+expU[4])/2, sig)
  }
  pr[pr<0.00001] = 0.00001;
  pr[pr>0.99999] = 0.99999;
  logLik = ydata * log(pr) + (1-ydata)*log(1-pr);
  return(sum(logLik));  
}

# 2.
grid1 = seq(-2,2,0.01)
library(haven)
setwd("Dropbox/Teaching/2020/Methods/Assignment/A4")
dat_choices = read_dta("dat_choices.dta")   # rows are individuals and columns are choices.

dat_choices = read_dta("dat_choices (1).dta")   # rows are individuals and columns are choices.

# creating the vectors of choices for individuals 900 and 115.
y900 = as.vector(dat_choices[900,])
y115 = as.vector(dat_choices[115,])

# vectors of likelihoods 
L900 = vector()
L115 = vector()
i=1
for(theta in grid1)
{
  L900[i] = Llike3(y900, theta)
  L115[i] = Llike3(y115, theta)
  i = i + 1
}

# this works, but this is very bad coding - something better
out900 = sapply(X=grid1, FUN = Llike3,ydata=dat_choices[900,])

plot(grid1,L900)
plot(grid1,L115)

# finding theta value that gives max of L900 and L115
m = cbind(grid1, L900, L115)
max900 = m[m[,2]==max(m[,2]),-3]  # theta* = 2
max115 = m[m[,3]==max(m[,3]),-2]  # theta* = 1.68


# =========================================================================#
# Exercise 4 Part 4
# =========================================================================#

beale = function(param) 
{
  x = param[1]
  y = param[2]
  out = (1.5 - x + x*y)^2 + 
    (2.25 - x + x*y^2)^2 + 
    (2.625 - x + x*y^3)^2
  return(out)
}

values_x_y = expand.grid(seq(-5, 5, 0.01), seq(-5,5,0.01))

z = beale(values_x_y)
# Seems to produce something.. but what is it? You take a function that takes as input a vector and you supply a matrix.. 
# I have no idea what R is doing!!!! 
# cbind(values_x_y,apply(values_x_y,1,beale))

beale_values = cbind(values_x_y, z)

min_beale = beale_values[which(beale_values[,3]==min(beale_values[,3])), ] # selecting the row where z has its minimum.
min_beale         # minimum of function equals 0 at (3, 0.5). 

x0_y0 = vector()
x0_y0[1] = min_beale[,1]  
x0_y0[2] = min_beale[,2]   # (x0, y0) = (3, 0.5), function equals 0 at (3, 0.5) 

Dbeale = function(param)
{
  eps = 0.00001
  out = NULL
  for(ix in length(param))
  {
    z1     = z2 = param
    z1[ix] = z1[ix] + eps
    z2[ix] = z2[ix] - eps
    out[ix] = beale(z1) - beale(z2)
  }
  return(out/2*eps)
}

Dbeale(x0_y0) # returns NA and 3.824999e-19


# =========================================================================#
# Exercise 5 Part 5
# =========================================================================#

Banana = function(vec)
{
  x = vec[1]
  y = vec[2]
  
  (1-x)^2 + 5*(y-x^2)^2
}

gradient = function(vec)
{
  eps = 0.0001;
  out = NULL
  for (iX in 1:length(vec))
  {
    z1      = z2 = vec;
    z1[iX]  = z1[iX] + eps;
    z2[iX]  = z2[iX] - eps;
    out[iX] = Banana(z1) - Banana(z2)
  }
  return(out/(2*eps));
}

# custom steepest descent function to find the minimum of Banana
B_steep = function(x0, tol = 0.000001, maxIter = 1000)
{
  diff = 1      # Initialize the difference
  iter = 1      # Initialize the iteration
  alpha_vec = c(0.0001, 0.001, 0.01, 0.1, 1);
  while (diff>tol & iter<maxIter) 
  {
    fx = Banana(x0)
    dx = gradient(x0)
    
    # find the optimal scaling parameter
    fPrime = NULL
    xPrime = NULL
    # in the spirit this is wrong.. 
    # what you want to do is 
    # evaluate various x1 for different values of alpha.. just like here 
    #val = sapply(X = alphaVec, FUN = line_search, vals=x0)
    # then select the scaling parameter that gives the max step increase
    #alpha = alphaVec[which(val==min(val))]
    # then update 
    # start = start - alpha*deriv
    for (iX in 1:length(alpha_vec))
    {
    x1 = x0 - alpha_vec[iX]*dx
    dx = gradient(x1);
    
    diff = Banana(x1) - Banana(x0)
    
    iter = iter + 1
    print(iter) 
    x0 = x1
    }
  }
 xprime = x1
 fPrime = Banana(x1)
    
a = array()
a$parameters = xprime
a$value = fPrime

return(a)
}

# =========================================================================#
# Exercise 6 Part 6
# this is exactly what I expected... 
# =========================================================================#
library(sjmisc)
library(AER)

data("SmokeBan")
mydat = SmokeBan

# converting categorical variables to numeric
mydat[,"smoker"] = as.numeric(as.numeric(mydat[,"smoker"])>1)
mydat[,"ban"] = as.numeric(as.numeric(mydat[,"ban"])>1)
mydat[,"gender"] = as.numeric(as.numeric(mydat[,"gender"])>1)  #female = 1
mydat[,"hispanic"] = as.numeric(as.numeric(mydat[,"hispanic"])>1)
mydat[,"afam"] = as.numeric(as.numeric(mydat[,"afam"])>1)
binary_education = to_dummy(mydat[,"education"])
colnames(binary_education) = c("hs drop out", "hs", "some college", "college", "master")
mydat = cbind(mydat, binary_education)

# 1. Probit model of smoking using optim()
xx = as.matrix(mydat[,-c(1,4,8)]) 
xx = cbind(matrix(1,10000,1), xx)   
xx = as.matrix(xx)
xx = xx/100
yy = as.matrix(mydat[,"smoker"])

betas = vector(length = 10L)

LogLike_Probit = function(y,x,b) {
  P = pnorm(x %*% b)
  
  f = sum(y*log(P) + (1-y)*log(1-P))
  return(-f)    
}

Gradient_LogLike_Probit = function (y,x,b) {
  P = pnorm(x %*% b) 
  D = dnorm(x %*% b) 
  
  n = length(y)     
  coef_num = length(b)
  
  g = t(matrix(rep(D/P,coef_num),nrow=n)*x) %*% y - 
    t(matrix(rep(D/(1-P),coef_num),nrow=n)*x) %*% (1-y)
  
  return(-g)
} 

probit_smoke = optim(par = betas, LogLike_Probit, y = yy, x = xx, gr = Gradient_LogLike_Probit,
                     method = c("BFGS"), hessian = TRUE)

# 2. Using Fisher Information to derive the standard error of the estimates.

Fish = solve(probit_smoke$hessian)
Std_errs = sqrt(diag(Fish))

# =========================================================================#
# Exercise 7 Part 7
# =========================================================================#
library(nloptr)
# 1. Estimating risk aversion and sigma for individual 5
ydata = dat_choices[5,]

Llike7 = function(p) 
{
  r = p[1]
  s = p[2]
  # using y5 as data
  w = 20
  nc = 4
  pr   = NULL;
  for (iX in 1:nrow(table1))
  {
    expU    = crra(w+table1[iX,],rep(r,nc))
    pr[iX]  = pnorm((expU[1]+expU[2])/2-(expU[3]+expU[4])/2, sd = s)
  }
  pr[pr<0.00001] = 0.00001;
  pr[pr>0.99999] = 0.99999;
  logLik = ydata * log(pr) + (1-ydata)*log(1-pr);
  return(-sum(logLik));  # because algorithms do minimization and we want to maximization.
}

# this take forever on my computer.. You need to set things like max iteration.. 
I = isres(c(2.1, 0.8), Llike7, lower = c(0, 0), upper = c(15, 15))  
L = lbfgs(c(2.1, 0.8), Llike7, lower = c(0, 0), upper = c(15, 15))
B = bobyqa(c(2.1, 0.8), Llike7, lower = c(0, 0), upper = c(15, 15))
NM = neldermead(c(2.1, 0.8), Llike7, lower = c(0, 0), upper = c(15, 15))

# Bobyqa algorithm is chosen for the next part.

# 2. Estimating risk aversion and sigma for all individuals using LBFGS

r7 = vector()
sig7 = vector()
hist_sig7 = hist(sig7, breaks = seq(min(sig7)-0.00001, max(sig7)+0.00001, 0.0000000003))

for(i in 1:nrow(dat_choices)){
  ydata = dat_choices[i,]
  B = bobyqa(c(2.1, 0.8), Llike7, lower = c(0, 0), upper = c(15, 15))
  r7[i] = B$par[1]
  sig7[i] = L$par[2] #typo probably
  print(i)
}

# 3. Representing the distribution of risk aversion across the sample

hist_r7 = hist(r7, breaks = seq(min(r7)-0.1, max(r7)+0.1, 0.03), 
               main = "Distribution of risk aversion parameter in the sample (estimated with BOBYQA)")


# =========================================================================#
# Exercise 8 Part 8
#==========================================================================#

#1.

grid_r = seq(0.01, 2, by = 0.01)

V_crra = function(c, w, theta) 
{
  ifelse(abs(theta-1)<0.01, log(c+w), (c+w)^(1-theta)/(1-theta))
}

# expected value difference between lottery 1 and 2
exp_value_diff = function(gamble, w, theta)
{
  Lot11 = 1
  Lot12 = 2
  Lot21 = 3
  Lot22 = 4
  
  V1_V2 = 0.5*V_crra(table1[gamble,Lot11], w, theta) + 0.5*V_crra(table1[gamble,Lot12], w, theta)-
    (0.5*V_crra(table1[gamble,Lot21], w, theta) + 0.5*V_crra(table1[gamble,Lot22], w, theta))
  
  return(V1_V2)
}

gambles = 1:25
probs = matrix(ncol=length(gambles), nrow=length(grid_r))  

for(g in gambles)
{
  w = 20
  r_num = 1
  for(r in grid_r)
  {
    probs[r_num, g] = pnorm(exp_value_diff(g, w, r), sd = 0.5)
    # I am really not a big fan of sneaking iterators in a loop.. Double looping will be more efficient in any matrix based language
    r_num = r_num +1
  }
}
# plots for the choices 1, 2, 5, 19 and 25
plot(grid_r, probs[,1])   
plot(grid_r, probs[,2])   
plot(grid_r, probs[,5])
plot(grid_r, probs[,19])
plot(grid_r, probs[,25])

#2. Monte Carlo Study to assess risk aversion

#a.# set values for risk aversion and standard deviation
r = 0.8
std = 0.5

#b.# simulate y_sim given r and std

#y_sim = 1 if utility difference is positive, then choose lottery 1.

exp_value_diff = function(gamble, theta = 0.8)
{
  w = 20
  # simulate error terms for both choice of lotteries
  err1 = rnorm(length(gamble), sd = std)
  err2 = rnorm(length(gamble), sd = std)
  
  Lot11 = 1
  Lot12 = 2
  Lot21 = 3
  Lot22 = 4
  
  V1_V2 = err1 + 0.5*V_crra(table1[gamble,Lot11], w, theta) + 0.5*V_crra(table1[gamble,Lot12], w, theta)-
    (err2 + 0.5*V_crra(table1[gamble,Lot21], w, theta) + 0.5*V_crra(table1[gamble,Lot22], w, theta))
  
  return(V1_V2)
}

y_sim = vector()

for(g in gambles)
{
  if(exp_value_diff(g, r) > 0) {y_sim[g] = 1}
  else {y_sim[g] = 0}
}

#c.# Estimate risk aversion and standard deviation using y_sim
# using likelihood function at Part 7
ydata = y_sim

I_8c = isres(c(2.1, 0.8), Llike7, lower = c(0, 0), upper = c(15, 15))  
L_8c = lbfgs(c(2.1, 0.8), Llike7, lower = c(0, 0), upper = c(15, 15))
B_8c = bobyqa(c(2.1, 0.8), Llike7, lower = c(0, 0), upper = c(15, 15))
NM_8c = neldermead(c(2.1, 0.8), Llike7, lower = c(0, 0), upper = c(15, 15))

# recording estimates 

MonteCarlo_r_sig = rbind(I_8c$par, L_8c$par, B_8c$par, NM_8c$par)
row.names(MonteCarlo_r_sig) = c("ISRES", "LBFGS", "BOBYQA", "Nelder-Mead Simplex")
colnames(MonteCarlo_r_sig) = c("risk aversion", "standard deviation")

#d.# Redo Monte Carlo study above 100 times and construct confidence interval

# simulating y 100 times 
sim = 100       # number of simulations
y = matrix(nrow = 100, ncol = 25)

for(i in 1:sim) 
{
  for(g in gambles)
  {
    if(exp_value_diff(g, r) > 0) {y[i,g] = 1}
    else {y[i,g] = 0}
  }
}

r_8d = vector()
sig_8d = vector()

for(i in 1:nrow(y)){
  ydata = y[i,]
  L = lbfgs(c(2.1, 0.8), Llike7, lower = c(0, 0), upper = c(15, 15))
  r_8d[i] = L$par[1]
  sig_8d[i] = L$par[2]
  print(i)    # to track the progress
}

# simulating standard errors to construct 95% confidence intervals
error_r = qnorm(0.975) * sd(r_8d)/sqrt(length(r_8d))
error_sig = qnorm(0.975) * sd(sig_8d)/sqrt(length(sig_8d))

CI_r = c(mean(r_8d)-error_r, mean(r_8d)+error_r)
CI_sig = c(mean(sig_8d)-error_sig, mean(sig_8d)+error_sig)


#3. Do the monte carlo study using the grid of parameters

grid_par = rbind(c(0.5, 0.1),c(0.5, 0.9),c(0.8, 0.1),c(0.8, 0.9),
                 c(1.5, 0.1),c(1.5, 0.9),c(2.8, 0.1),c(2.8, 0.9))

# simulate y
y_sim2 = vector()
r = 0.8

for(g in gambles)
{
  if(exp_value_diff(g, r) > 0) {y_sim2[g] = 1}
  else {y_sim2[g] = 0}
}


r_8.3 = vector()
sig_8.3 = vector()
ydata = y_sim2

for(i in 1:nrow(grid_par))
{
  x0 = as.vector(grid_par[i,])

  B = bobyqa(x0, Llike7, lower = c(0, 0), upper = c(5, 5))
  
  r_8.3[i] = B$par[1]
  sig_8.3[i] = B$par[2]
}

r_8.3
sig_8.3

# overall this part is really hard to follow...
# starting a piece code with a general feel - WHat I am doing here would have helped..
# ALL our estimates under bobyqa hit the bound.. Which is a sign that there is computing problem.. 
# I would have organized this into a single function.. that calls
  # simulateData = function(t, s, choicerow, w=wealth) {
  # myLikelihood = function(params, ydat, choiceSet=choice, w=wealth) {
  # estimateLikelihood = function(person, start=myParams, choiceSet = choice, w=wealth) {


#==========================================================================#
# Exercise 9 Part 9
#==========================================================================#

w = 20 # by assumption

# Table 2: Time Preference Lotteries 

Table2 = rbind(c(75,1,75.31,31), c(75,1,75.63,31), c(75,1,76.25,31),
               c(75,1,78.13,31), c(75,1,81.25,31), c(75,1,87.50,31),
               
               c(75,8,75.31,31), c(75,8,75.63,31), c(75,8,76.25,31),
               c(75,8,78.13,31), c(75,8,81.25,31), c(75,8,87.50,31),
               
               c(75,31,75.31,61), c(75,31,75.63,61), c(75,31,76.25,61),
               c(75,31,78.13,61), c(75,31,81.25,61), c(75,31,87.50,61),
               
               c(75,91,75.31,121), c(75,91,75.63,121), c(75,91,76.25,121),
               c(75,91,78.13,121), c(75,91,81.25,121), c(75,91,87.50,121),
               
               c(75,1,75.31,361), c(75,1,75.63,361), c(75,1,76.25,361),
               c(75,1,78.13,361), c(75,1,81.25,361), c(75,1,87.50,361),
               
               c(75,8,75.31,368), c(75,8,75.63,368), c(75,8,76.25,368),
               c(75,8,78.13,368), c(75,8,81.25,368), c(75,8,87.50,368),
               
               c(75,31,75.31,391), c(75,31,75.63,391), c(75,31,76.25,391),
               c(75,31,78.13,391), c(75,31,81.25,391), c(75,31,87.50,391),
               
               c(75,91,75.31,451), c(75,91,75.63,451), c(75,91,76.25,451),
               c(75,91,78.13,451), c(75,91,81.25,451), c(75,91,87.50,451))

colnames(Table2) = c("P1","Time1", "P2", "Time2")
# utility function is crra 
V_crra = function(c, w, theta) 
{
  ifelse(abs(theta-1)<0.01, log(c+w), (c+w)^(1-theta)/(1-theta))
}

# utility at period 0 accounting for time discounting 

Utility_t0 = function(row, preference, theta, param)
{
  Time1 = Table2[row, "Time1"]
  Time2 = Table2[row, "Time2"]
  P1    = Table2[row, "P1"]
  P2    = Table2[row, "P2"]
  
  beta  = param[1]
  delta = param[2]
  
  w = 20
  
  if(preference == 1) {U = beta*delta^Time1*V_crra(P1, w, theta)}
  else if (preference == 2) {U = beta*delta^Time2*V_crra(P2, w, theta)}
  
  return(U)
}

# Utility difference between earlier and later choice assuming normally dist. errors

Diff_time_pref = function(row, theta, params)
{
  Time1 = Table2[row, "Time1"]
  Time2 = Table2[row, "Time2"]
  P1    = Table2[row, "P1"]
  P2    = Table2[row, "P2"]
  
  beta  = params[1]
  delta = params[2]
  sig   = params[3]
  
  w = 20
  
  err1 = rnorm(length(row), sd = sig)
  err2 = rnorm(length(row), sd = sig)
  
  Out = err1 + Utility_t0(row, 1, theta, params[1:2]) - err2 - Utility_t0(row, 2, theta, params[1:2])

  return(Out)
}
# Log of Likelihood of the choice problem (y=1 if P1 is chosen.)

LL = function(params) 
{
  beta  = params[1]
  delta = params[2]
  sig   = params[3]
  
  w = 20
  
  pr = NULL;
  for(ix in 1:nrow(Table2))
  {
    pr[ix] = pnorm(Diff_time_pref(ix, theta, params))
  }
  pr[pr<0.00001] = 0.00001;
  pr[pr>0.99999] = 0.99999;
  LL = y*log(pr) + (1-y)*log(1-pr)
  return(-sum(LL))
}

# 2. 
library(haven)
dat_time = read_dta("dat_time.dta")

betas9  = vector()
deltas9 = vector()
sigs9   = vector()

for(i in 1:nrow(dat_time)){
  y     = dat_time[i,]  # selecting individual data
  theta = r7[i]         # selecting theta estimated at Part7
  B = bobyqa(c(0, 0,0), LL, lower = c(0,0,0), upper = c(10,10,10))
  betas9[i]  = B$par[1]
  deltas9[i] = B$par[2]
  sigs9[i]   = B$par[3]
  print(i)
}

# Derive the standard deviation of these estimates.

sd_betas = sd(betas9)
sd_deltas = sd(deltas9)
sd_sigs = sd(sigs9)

# 3. Evaluate the correlation between theta (risk aversion), and beta and delta.

corr1 = cor(r7, betas9)
corr2 = cor(r7, deltas9)


#========================================================================

Output <- file.path(getwd())
e <- environment()
save(file = file.path(Output, "A4_B.RData"),
     list=ls(), env=e)
