
?`Rcpp-package`

rm(list=ls())

library(Rcpp)
library(AER)
data("GSOEP9402", package = "AER")
dat = GSOEP9402
#  kids, marital, meducation, income
levels(dat$marital)  # "married" "single" "widowed" "divorced" "separated" (might need to transform these to numberical)

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
w = dat$income

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

trapzf = function(f, c, d, nn, X)
{
  hh = (d-c)/nn
  
  t1 = (f(wage_vec=X, c) + f(wage_vec=X, d))*hh/2
  
  vecc = vector()
  for(i in 1:nn)
  {
    while ((c+i*hh) < d) {vecc[i] = f(wage_vec=X, (c+i*hh))}
  }
  t2 = hh*sum(vecc)
  
  T = t1+t2
  return(list(T, vecc))
}

######### Alternative Trapezoidal integral function############

trapzf2 = function(c, d)
{
  w = dat$income
  nn = length(w)/2  
  hh = (d-c)/nn
  
  t1 = (kdens(wage_vec=w, x0 = c) + kdens(wage_vec=w, x0 = d))*hh/2
  
  vecc = vector()
  for(i in 1:nn)
  {
    if((c+i*hh) <= d) {vecc[i] = kdens(wage_vec=w, x0 = (c+i*hh))}
    else break
  }
  t2 = hh*sum(vecc)
  
  T = t1+t2
  return(list(T, vecc))
}
max(w) # 258341.4
min(w)  # 1247.59
B = trapzf2(1247.59, 258341.4) # B[[1]] = 0.6211489 & length(B[[2]]) = 337
####################################

# Exercise 4 Reservation Wages

w = dat$income
w_bar = max(w)

lambda_u = 2
lambda_e = 0.5
delt = 0.2
r = 0.01

# 1. 

fcdf = function(wages, x0)
{
  delt = 0.2
  lambda_e = 0.5
  
  out = (delt/ecdf(vecX=wages, x0=x0) - delt)/lambda_e
  
  return(1-out)
}


fpdf = function(wages, x0)
{
  delt = 0.2
  lambda_e = 0.5
  
  F_bar = (delt/kdens(wage_vec = wages, x0 = x0) - delt)/lambda_e
    
  return(1-F_bar)
}

offer_cdf = vector()
for(i in 1:length(w))
{
  offer_cdf[i] = fcdf(w, w[i])
}

offer_pdf = vector()
for(i in 1:length(w))
{
  offer_pdf[i] = fpdf(w, w[i])
}

par(mfrow=c(1,3))
plot(w, offer_cdf, ylim = range(0,1))
plot(w, offer_cdf)
plot(w, offer_pdf)    # all values are negative here
##########################################################
# 2. 
reserv_wage = function(phi)
{
  w = dat$income
  w_bar = max(w)
  select_wage = which(w>=phi)
  
  nn = length(w)/3
  hh = (w_bar - phi)/nn
  
  lambda_e = 0.5
  lambda_u = 2
  diff = lambda_u - lambda_e
  r = 0.01
  delt = 0.2
  
  f1 = leisure - phi
  f2 = function(x)
    {
     F_bar = (delt/ecdf(vecX=select_wage, x0=x) - delt)/lambda_e
     R = F_bar/(r+delt+lambda_e*F_bar)
     return(R)
    }
  
  vecc = vector()
  for(i in nn)
  {
    x = phi + i*hh
    if(x < w_bar) {vecc[i] = f2(x)}
    else break    
  }
  
  Trpz_integral = hh*((f2(w_bar) + f2(phi))/2 + sum(vecc))
    
  Obj_funct = f1 + diff * Trpz_integral
  return(Obj_funct)
}

maxx = max(w)
minn = min(w)

leisure = 4000
U = uniroot(reserv_wage, lower = minn, upper =  maxx)
reserv_wage(U$root)


# Exercise 5 Indirect Inference

# 1. 
library(sjmisc)

marital = to_dummy(dat[,"marital"])
colnames(marital) = c("married", "single" , "widowed" ,"divorced", "separated")
dat = cbind(dat, marital)

reg = glm(income ~ kids + meducation + married + single + widowed +  divorced, data = dat)
summary(reg)
err = residuals(reg)
coeffs = coefficients(reg)
sig = sd(err)

theta_data = c(coeffs, sig)

# 2. 
