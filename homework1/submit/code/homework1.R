# ==================== UTILS ====================
check_num = function(a, b, fun){
  if (abs(a[1] - b[1]) / max(abs(a[1]), 1) > 0.05){
    #print(c(a[1], b[1]))
    #print(c(fun(a), fun(b)))
    #print('check numerical error! [a]')
  }
  if (abs(a[2] - b[2]) / max(abs(a[2]), 1) > 0.05){
    #print('check numerical error! [b]')
  }
}

cmp = function(par1, par2, fun){
  if (fun(par1) < fun(par2))
    return(par1)
  else
    return(par2)
}

len = function(x){
  return(length(x))
}

# ==================== METHOD ====================
# (vertical) Least square estimate
VLSE = function(x, y){
  sx = sum(x)
  sx2 = sum(x * x)
  sxy = sum(x * y)
  sy = sum(y)
  delta = n * sx2 - sx^2
  a = (n * sxy - sx * sy) / delta
  b = (sx2 * sy - sxy * sx) / delta
  lsev = c(a, b)
  return(lsev)
}

# (horizontal) Lesat square estimate
HLSE = function(x, y){
  # analytical solution
  u = x - mean(x)
  v = y - mean(y)
  a = sum(v * v) / sum(v * u)
  b = mean(y) - a * mean(x)
  lseh = c(a, b)
  # numerical solution
  f_lseh = function(par){
    a = par[1]
    b = par[2]
    hd = x - (y - b) / a
    return(sum(hd * hd) / len(x))
  }
  gr_lseh = function(par){
    a = par[1]
    b = par[2]
    hd = (x - (y - b) / a)
    gr_a = 2 * sum(hd * (y - b) / a^2) / len(x)
    gr_b = 2 * sum(hd / a) / len(x)
    return(c(gr_a, gr_b))
  }
  num_lseh = optim(VLSE(x, y), f_lseh, gr_lseh)$par
  check_num(lseh, num_lseh, f_lseh)
  return(cmp(lseh, num_lseh, f_lseh))
}

# (perpendicular) Least square estimate
PLSE = function(x, y){
  # analytical solution
  u = x - mean(x)
  v = y - mean(y)
  A = sum(u * v)
  B = sum(u * u) - sum(v * v)
  C = -sum(u * v)
  delta = B^2 - 4 * A * C
  a1 = (-B + sqrt(delta)) / (2 * A)
  b1 = mean(y) - a1 * mean(x)
  a2 = (-B - sqrt(delta)) / (2 * A)
  b2 = mean(y) - a2 * mean(x)
  f = function(a, b){
    return(sum((y - a * x - b)^2) / (a^2 + 1))
  }
  if (f(a1, b1) < f(a2, b2)) a = a1 else a = a2
  b = mean(y) - a * mean(x)
  lsep = c(a, b)
  
  # numerical solution
  f_lsep = function(par){
    a = par[1]
    b = par[2]
    hd = y - a * x - b
    return(sum(hd * hd) / (a^2 + 1) / len(x))
  }
  gr_lsep = function(par){
    a = par[1]
    b = par[2]
    hd = y - a * x - b
    ihd = -x - a * y - b
    gr_a = sum(hd * ihd) * 2 / (a^2 + 1)^2 / len(x)
    gr_b = -2 / (a^2 + 1) * sum(hd) / len(x)
    return(c(gr_a, gr_b))
  }
  lsep_num = optim(VLSE(x, y), f_lsep, gr_lsep)$par
  check_num(lsep, lsep_num, f_lsep)
  return(cmp(lsep, lsep_num, f_lsep))
}

# (vertical) Least Absolute Deviation
VLAD = function(x, y){
  f_ladv = function(par){
    a = par[1]
    b = par[2]
    hd = (y - a * x - b)
    return(sum(abs(hd)) / len(x))
  }
  gr_ladv = function(par){
    a = par[1]
    b = par[2]
    hd = (y - a * x - b)
    gr_a = (((hd > 0) * (-x)) - ((hd < 0) * (-x))) / len(x)
    gr_b = (((hd > 0) * (-1)) - ((hd < 0) * (-1))) / len(x)
    return(c(gr_a, gr_b))
  }
  ladv = optim(VLSE(x, y), f_ladv, gr_ladv)
  return(ladv$par)
}

# (vertical) Least Maximum Estimate
VLME = function(x, y){
  f_lmev = function(par){
    a = par[1]
    b = par[2]
    hd = y - a * x - b
    return(max(abs(hd)))
  }
  gr_lmev = function(par){
    a = par[1]
    b = par[2]
    hd = y - a * x - b
    i = which.max(abs(hd))
    gr_a = (y[i] - a * x[i] - b) * (-x[i])
    gr_b = (y[i] - a * x[i] - b) * (-1)
    return(c(gr_a, gr_b))
  }
  lmev = optim(VLSE(x, y), f_lmev, gr_lmev)
  return(lmev$par)
}

# ==================== EVALUATION ====================
# environment:
require(ggplot2)
require(xtable)
set.seed(123469)

# model:
#   y_i = a x_i + b + e_i

# set true parameter
a = 1
b = 2

# set number of repeats
T = 100000
vlad = vlse = vlme = hlse = plse = matrix(rep(0, 2 * T), nrow = T, ncol = 2)

# set number of observations
n = 30

for (t in seq(1, T)){
  # generate x data
  # x = runif(n, min = 0, max = 1)
  x = seq(0, 1, length.out = n)
    
  # set noise
  # noise = rnorm(n, mean = 0, sd = 2)
  noise = rcauchy(n, location = 0, scale = 0.1)
  # noise = runif(n, min = -0.2, max = 0.2)
  y = a*x + b + noise
  
  vlad[t,1] = VLAD(x, y)[1]; vlad[t,2] = VLAD(x, y)[2]
  vlse[t,1] = VLSE(x, y)[1]; vlse[t,2] = VLSE(x, y)[2]
  vlme[t,1] = VLME(x, y)[1]; vlme[t,2] = VLME(x, y)[2] 
  hlse[t,1] = HLSE(x, y)[1]; hlse[t,2] = HLSE(x, y)[2] 
  plse[t,1] = PLSE(x, y)[1]; plse[t,2] = PLSE(x, y)[2]
}

# ==================== summary ====================
am = c(mean(vlad[,1]), mean(vlse[,1]), mean(vlme[,1]), mean(hlse[,1]), mean(plse[,1]))-a
as = c(sd(vlad[,1]), sd(vlse[,1]), sd(vlme[,1]), sd(hlse[,1]), sd(plse[,1]))
bm = c(mean(vlad[,2]), mean(vlse[,2]), mean(vlme[,2]), mean(hlse[,2]), mean(plse[,2]))-b
bs = c(sd(vlad[,2]), sd(vlse[,2]), sd(vlme[,2]), sd(hlse[,2]), sd(plse[,2]))
summary = data.frame(estimator = c('VLAD', 'VLSE', 'VLME', 'HLSE', 'PLSE'), a.bias = am, a.sd = as, b.bias = bm, b.sd = bs)
xtable(summary,digits=3)

# ==================== 5 method plot ==================== 
estimate = c(vlad[,1], vlad[,2], vlse[,1], vlse[,2], vlme[,1], vlme[,2],
          hlse[,1], hlse[,2], plse[,1], plse[,2])
parameter = rep(c(rep('a', T), rep('b', T)), times=5)
method = rep(c('VLAD','VLSE','VLME','HLSE','PLSE'), each=2*T)
id = rep(seq(1,1000), 10)
result = data.frame(idx=id, estimate=estimate, parameter=parameter, method=method)
ggplot(result, aes(x = method, y = estimate, col = parameter)) + geom_boxplot()

# ==================== 3 method plot ==================== 
estimate = c(vlad[,1], vlad[,2], vlse[,1], vlse[,2], vlme[,1], vlme[,2])
parameter = rep(c(rep('a', T), rep('b', T)), times=3)
method = rep(c('VLAD','VLSE','VLME'), each=2*T)
id = rep(seq(1,1000), 6)
result = data.frame(idx=id, estimate=estimate, parameter=parameter, method=method)
ggplot(result, aes(x = method, y = estimate, col = parameter)) + geom_boxplot()

# ==================== fit line visualization ====================
# set.seed(123469)
set.seed(1234679)
# set true parameter
a = 1
b = 2

n = 100

# get x data
x = runif(n, min = 0, max = 1)
# x = seq(0, 1, length.out = n)

# noise
# noise = rnorm(n, mean = 0, sd = 1)
# noise = rcauchy(n, location = 0, scale = 0.2)
noise = runif(n, min = -1, max = 1)

y = x * a + b + noise

scat = data.frame(x = x, y = y)
g = ggplot(scat, aes(x = x, y = y)) + geom_point()
g = g + geom_abline(slope=a, intercept=b, show.legend=TRUE)
g = g + geom_abline(slope=VLSE(x,y)[1], intercept=VLSE(x,y)[2], col='red')
g = g + geom_abline(slope=VLAD(x,y)[1], intercept=VLAD(x,y)[2], col='orange')
g = g + geom_abline(slope=VLME(x,y)[1], intercept=VLME(x,y)[2], col='green')
g = g + geom_abline(slope=HLSE(x,y)[1], intercept=HLSE(x,y)[2], col='blue')
g = g + geom_abline(slope=PLSE(x,y)[1], intercept=PLSE(x,y)[2], col='purple')
g



