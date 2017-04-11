# ==================== EVALUATION ====================
evaluate = function(l, r, delta, estimator, value, is.plot=F){
  x = seq(l + delta, r, by = delta)
  y = value(x)
  n = length(x)
  y.hat = rep(0, n)
  id = 0
  for (i in x){
    id = id + 1
    y.hat[id] = estimator(x[id])
  }
  if (is.plot){
    lines(x, y, col = 'black')
    lines(x, y.hat, col = 'red')
  }
  l2.error = sum((y.hat - y)^2) * delta / (r - l)
  return(l2.error)
}

# environment:
require(ggplot2)
require(xtable)

# model:
#   y_i = m(x_i) + e_i

# set true parameter
m = function(x){
  return(log(1+1/x))
  #return(x * sin(x))
}

# set number of repeats
T = 1000

# set number of observations
n = 1000

#para = seq(0.01, 3.0, 0.01)
para = c(1.98)
mean.l2.error = para
id = 0
for (C in para){
  id = id + 1
  l2.error = rep(0, T)
  set.seed(123469)
  for (t in seq(1, T)){
    # generate x data
    x = runif(n, min = 0, max = 1)
    #x = seq(1.0/n, 1, length.out = n)
    x = sort(x)
  
    # set noise
    #noise = rnorm(n, mean = 0, sd = 0.1)
    noise = rt(n=n, df=1)
    
    # generate data
    y = m(x) + noise
  
    # estimator parameter
    h = as.integer(C * n ^ (4/5))
  
    # estimator
    m.hat = function(qx){
      # :para qx: float
      # :return float: estimate
      idx = h
      idx = min(h, sum(x <= qx))
      idx = min(idx, sum(x >= qx))
      idx = idx * 2 + 1
      if (idx > n)
        idx = n
      x.dist = abs(qx - x)
      x.dist = sort(x.dist)
      bound = x.dist[idx]
      check = abs(qx - x) <= bound
      ans = sum(check * y) / sum(check)
      return(ans)
    }
    if (t == T)
      plot(x, y)
    l2.error[t] = evaluate(0, 1, delta=0.01, m.hat, m, is.plot=(t == T))
  }
  mean.l2.error[id] = mean(l2.error)
}
