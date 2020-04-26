### FINAL ASSIGNMANT ### 
# 04/25/2020

# Bahar Zafer 
########################
rm(list=ls())

library(Rcpp)
library(AER)
library(nloptr)

data("GSOEP9402", package = "AER")
dat = GSOEP9402
w = dat$income

# Exercise 1  Kernel Density Estimation

# 1. Bandwidth
band = function(wage_vec)
{
  sig = sd(wage_vec)
  n = length(wage_vec)
  
  h = (4*sig^5/3*n)^(1/5) 
  return(h)
}

# 2. Gaussian kernel
gaussian = function(z) 
{
  m = as.matrix(z)
  K = (sqrt(2*pi))^(-1) * (exp((t(m)%*%m)/2))^(-1)
  K = as.vector(K)
  return(K)
}

# 3. Kernel density estimator at x0
kdens = function(wage_vec, x0)
{
  n = length(wage_vec)
  h = band(wage_vec)
  nh = n*h
  
  vec = vector()
  for(j in 1:n)
  {
    input = (x0 - wage_vec[j])/h
    vec[j] = gaussian(input)
  }
  s = sum(vec)/nh
  return(s)
}

# Exercise 2 Empirical CDF

# 1. Empirical cdf at point x0
ecdf = function(vecX, x0)
{ 
  n = length(vecX)
  vec = vector()
  for(i in 1:n)
  { 
    if(vecX[i] <= x0) {vec[i] = 1} 
    else {vec[i] = 0} 
  } 
  s = sum(vec)
  return(s/n)
}

####################################################################
# to see how these work so far
# w = dat$income

ecdf(w, median(w))  # this returns 0.5007407
ecdf(w, max(w))     # this returns 1

s = Sys.time()

k_density = vector()
for(i in 1:length(w))
{
  k_density[i] = kdens(w, w[i])
}

ecdd = vector()
for(i in 1:length(w))
{
  ecdd[i] = ecdf(w, w[i])
}

a = Sys.time()
a-s   # Time difference of 5.349357 secs

par(mfrow=c(1,3))
plot(w, k_density)
hist(w)
plot(w, ecdd)
####################################################################

# Exercise 3  Numerical Integration

# 1. Trapezoidal integral rule

# a general trapezoidal integral function 
trapzf = function(f, c, d, nn, X)
{
  hh = (d-c)/nn
  
  t1 = (f(wage_vec=X, c) + f(wage_vec=X, d))*hh/2
  
  vecc = vector()
  for(i in 1:nn)
  {
    if((c+i*hh) <= d) {vecc[i] = f(wage_vec=X, (c+i*hh))}
    else break
  }
  t2 = hh*sum(vecc)
  
  T = t1+t2
  L = list(result=T, points=(length(vecc)+2))
  
  return(L)
}

# trapezoidal integral function customized for kdens function and income data
trapzf2 = function(c, d)
{
  w = dat$income
  nn = length(w)/3
  hh = (d-c)/nn
  
  t1 = (kdens(wage_vec=w, c) + kdens(wage_vec=w, d))*hh/2
  
  vecc = vector()
  for(i in 1:nn)
  {
    if((c+i*hh) <= d) {vecc[i] = kdens(wage_vec=w, c+i*hh)}
    else break
  }
  t2 = hh*sum(vecc)
  
  T = t1+t2
  L = list(result=T, points=(length(vecc)+2))
  
  return(L)
}

# w = dat$income
max(w) # 258341.4
min(w)  # 1247.59
B = trapzf2(1247.59, 258341.4) # B$result = 0.6219899 & B$points = 227

# trapezoid function is tried with known functions and it works
#==============================================================
trapz_integ = function(func, c, d)
{
  nn = 1000000  
  hh = (d-c)/nn
  
  t1 = (func(c) + func(d))*hh/2
  
  vecc = vector()
  for(i in 1:nn)
  {
    if((c+i*hh) <= d) {vecc[i] = func(c+i*hh)}
    else break
  }
  t2 = hh*sum(vecc)
  
  T = t1+t2
  
  L = list(result=T, points=(length(vecc)+2))
  return(L)
}

m = function(x)
{
  r = x + 90
  return(r)
}

A = trapz_integ(m,-1,1)
C = trapz_integ(func = m, -2, 1)

####################################

# Exercise 4 Reservation Wages

# w = dat$income
# w_bar = max(w)

# lambda_u = 2
# lambda_e = 0.5
# delt = 0.2
# r = 0.01

# 1. 
fcdf = function(wages, x0)
{
  delt = 0.2
  lambda_e = 0.5
  G = ecdf(vecX=wages, x0=x0)
  
  out = (delt/G - delt)/lambda_e
   
  return(1-out)
}


fpdf = function(wages, x0)
{
  delt = 0.2
  lambda_e = 0.5
  G = kdens(wage_vec = wages, x0 = x0)
  
  F_bar = (delt/G - delt)/lambda_e
    
  return(1-F_bar)
}


# 2. 
reserve_wage = function(phi)
{
  w = dat$income
  w_bar = max(w)
  
  lambda_e = 0.5
  lambda_u = 2
  dif = lambda_u-lambda_e
  # modify trapezoidal function to use in reserve wage function 
          trapzf4 = function(c, d)
            {
              w = dat$income      
              
              f = function(x)
                  {
                    lambda_e = 0.5
                    delt = 0.02
                    r = 0.01
                    # F_bar = 1 - fcdf 
                    R = (1-fcdf(wages = w, x0 = x))/
                      (r+delt+lambda_e*(1-fcdf(wages = w, x0 = x)))
                    return(R)
                   }
              func = f
              
              nn = length(w)/10     # lower nn can be chosen to decrease the computation time
              hh = (d-c)/nn
              
              t1 = (func(c) + func(d))*hh/2
              
              vecc = vector()
              for(i in 1:nn)
              {
                if((c+i*hh) < d) {vecc[i] = func(c+i*hh)}
                else break
              }
              t2 = hh*sum(vecc)
              
              T = t1+t2
              # L = list(result=T, points=(length(vecc)+2))
              
              return(T)
            }
  Obj_funct = leisure - phi + dif*trapzf4(c = phi, d = w_bar)  
  return(Obj_funct)
}

leisure = max(w) # = 258341.4
U = uniroot(reserve_wage, lower = 2000, upper = 260000)
phi = U$root     # = 258341.4

# Exercise 5 Indirect Inference
#===============================
# 1. 
library(sjmisc)

marital = to_dummy(dat[,"marital"])
colnames(marital) = c("married", "single" , "widowed" ,"divorced", "separated")
dat = cbind(dat, marital)

reg = glm(income ~ kids + meducation + married + single + widowed + divorced, data = dat)
reg_out = summary(reg)
err = residuals(reg)
coeffs = coefficients(reg)
sd_e = sd(err)

theta_data = c(coeffs, sd_e=sd_e)   # set of parameters 

# 2. 
# set initial values 
theta_0 = c(intercept = -10000, kids=3500, meducation=5000, married=37000, 
            single=4000, widowed=24000, divorced=12000, sd_e=2500)
length(dat$marital) # there are 675 individuals
n = 675
# estimate b_0 for each individuals in the dataset.
errors = rnorm(675, sd=theta_0["sd_e"])
b_0 = vector()
for (i in 1:n)
{
  b_0[i] = theta_0["intercept"] + theta_0["kids"]*dat$kids[i] + 
           theta_0["meducation"]*dat$meducation[i] + theta_0["married"]*dat$married[i] + 
           theta_0["single"]*dat$single[i] + theta_0["widowed"]*dat$widowed[i] +
           theta_0["divorced"]*dat$divorced[i] + errors[i]
}  

# recover phi(b_0)
phi_0 = vector()

for(i in 1:n)
{
  leisure = b_0[i]
  U = uniroot(reserve_wage, lower = 2000, upper = 260000)
  phi_0[i] = U$root
  print(i)
}
# plot(b_0, phi_0)

# simulation 
steps = 10
phi = vector()
for (j in 1:n)  # n = 675 , number of indivduals
{
e = rnorm(steps, sd=theta_0["sd_e"])
b0 = vector()
phi0 = vector()
for (i in 1:steps) 
  {
    b0[i] = theta_0["intercept"] + theta_0["kids"]*dat$kids[j] + 
             theta_0["meducation"]*dat$meducation[j] + theta_0["married"]*dat$married[j] + 
             theta_0["single"]*dat$single[j] + theta_0["widowed"]*dat$widowed[j] +
             theta_0["divorced"]*dat$divorced[j] + e[i]
   
    leisure = b0[i]
    U = uniroot(reserve_wage, lower = 2000, upper = 260000)
    
    phi0[i] = U$root
}
phi[j] = mean(phi0)
print(j)   # to track the progress, takes 6-7 seconds for each individual.
}

# regression of phi(b0) on observed attributes 
reg_sim = glm(phi ~ kids + meducation + married + single + widowed + divorced, data = dat)

reg_sim_out = summary(reg_sim)
err_sim = residuals(reg_sim)
coeffs_sim = coefficients(reg_sim)
sd_e_sim = sd(err_sim)

theta_model = c(coeffs_sim, sd_e_sim=sd_e_sim)

# 3. Minimization problem
n = 675     # number of indivduals
steps = 10  # number of simulations
obj_F = function(t0)
{
  phi = vector()
  
  for (j in 1:n) 
  {
    e = rnorm(steps, sd=t0["sd_e"])
    b0 = vector()
    phi0 = vector()
    for (i in 1:steps) 
    {
      b0[i] = t0["intercept"] + t0["kids"]*dat$kids[j] + 
        t0["meducation"]*dat$meducation[j] + t0["married"]*dat$married[j] + 
        t0["single"]*dat$single[j] + t0["widowed"]*dat$widowed[j] +
        t0["divorced"]*dat$divorced[j] + e[i]
      
      leisure = b0[i]
      U = uniroot(reserve_wage, lower = 2000, upper = 260000)
      
      phi0[i] = U$root
    }
    phi[j] = mean(phi0)
  }
  
  reg_sim = glm(phi ~ kids + meducation + married + single + widowed + divorced, data = dat[1:n,])
  
  err_sim = residuals(reg_sim)
  coeffs_sim = coefficients(reg_sim)
  sd_e_sim = sd(err_sim)
  
  theta_model = c(coeffs_sim, sd_e_sim=sd_e_sim)
  theta_model = as.matrix(theta_model)
  theta_data = as.matrix(theta_data)
  
  return(t(theta_model-theta_data) %*% (theta_model-theta_data))
}

optm = bobyqa(theta_0, obj_F)

optimized_pars = optm$par

# 4.  
boot = function(samplesize, bootNum, simNum) 
{ 
  M = matrix(nrow = bootNum, ncol = 8)    # 8 is the number of parameters 
  for(i in 1:bootNum)
  {
    #1  estimate theta_data for samples 
    indx = sample(1:nrow(dat), samplesize, replace = T)
    datB = dat[indx,]
    Model = glm(income ~ kids + meducation + married + single + widowed + divorced, data = datB)
    theta_data = c(coef(Model), sd(residuals(Model)))
    theta_data = as.matrix(theta_data)
    
    #2  estimate theta_model for samples 
    Model_Obj = function(t0){
      n = samplesize
      steps = simNum
      dat = datB
      obj_F(t0)
    }
    OPT = bobyqa(theta_0, Model_Obj)
    theta_model = OPT$par
    M[i,] = theta_model 
  }
return(M)
}

# estimate bootstrap
M = boot(100, 500, 10)

SD_Estimator = function(parChoice)  
{
  parMean = mean(M[,parChoice]) # colmean 
  bootNum = nrow(M)
  vec = vector()
  for(i in 1:bootNum)
  {
    vec[i] = (M[i, parChoice] - parMean)^2
  }
  Variance = sum(vec)/(bootNum-1)
  SD = sqrt(Variance)
  return(SD)
}

SD_intercept = SD_Estimator[1] 
SD_kids      = SD_Estimator[2] 
SD_education = SD_Estimator[3] 
SD_married   = SD_Estimator[4] 
SD_single    = SD_Estimator[5] 
SD_widowed   = SD_Estimator[6] 
SD_divorced  = SD_Estimator[7] 
SD_sd_e      = SD_Estimator[8] 

