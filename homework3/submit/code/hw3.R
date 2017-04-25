# library
library('nortest')
library('ggplot2')
library('reshape2')

# set significance level
alpha = 0.05

# calculate type 1 error
T = 100000

set.seed(123469)
vec = rep(0, 5)
results = data.frame(n.30 = vec, n.100 = vec, n.1000 = vec)
rownames(results) = c('ad', 'cvm', 'lillie', 'pearson', 'sf')
colnames(results) = c('n=30', 'n=100', 'n=1000')

n.vec = c(30, 100, 1000)

for (k in seq(1, length(n.vec))){
  n = n.vec[k]
  ad.reject = 0
  cvm.reject = 0
  lillie.reject = 0
  pearson.reject = 0
  sf.reject = 0
  for (i in seq(1, T)){
    x = rnorm(n, 0, 1)
    ad.reject = ad.reject + (ad.test(x)$p.value < alpha)
    cvm.reject = cvm.reject + (cvm.test(x)$p.value < alpha)
    lillie.reject = lillie.reject + (lillie.test(x)$p.value < alpha)
    pearson.reject = pearson.reject + (pearson.test(x)$p.value < alpha)
    sf.reject = sf.reject + (sf.test(x)$p.value < alpha)
  }
  results[1, k] = (ad.reject + 0.0) / T
  results[2, k] = (cvm.reject + 0.0) / T
  results[3, k] = (lillie.reject + 0.0) / T
  results[4, k] = (pearson.reject + 0.0) / T
  results[5, k] = (sf.reject + 0.0) / T
}
print(results)
xtable(results, digits=4)


# calculate type 2 error - trival

x = seq(-5, 5, by=0.01)
density = dnorm(x, 0, 1)
plot(x, density, type='l', ylim = c(0, 1.0))
lines(x, dunif(x*sqrt(1/12)+0.5, 0, 1)*sqrt(1/12), col = 'red')
lines(x, dcauchy(x, 0, 1), col = 'orange')
m = mean(rgamma(100000, 5, 3))
s = sd(rgamma(100000, 5, 3))
lines(x, dgamma(x*s + m, 5, 3) * s, col = 'green')
m = mean(rexp(100000, 1))
s = sd(rexp(100000, 1))
lines(x, dexp(x*s + m, 1) * s, col = 'blue')
m = mean(rchisq(100000, 10))
s = sd(rchisq(100000, 10))
lines(x, dchisq(x*s + m, 10) * s, col = 'purple')
m = mean(rt(100000, 10))
s = sd(rt(100000, 10))
lines(x, dt(x*s + m, 10) * s, col = 'brown')

T = 10000

set.seed(123469)
vec = rep(0, 5)
results = data.frame(n.30 = vec, n.100 = vec, n.1000 = vec)
rownames(results) = c('ad', 'cvm', 'lillie', 'pearson', 'sf')
colnames(results) = c('n=30', 'n=100', 'n=1000')

n.vec = c(30, 100, 1000)

for (k in seq(1, length(n.vec))){
  n = n.vec[k]
  ad.reject = 0
  cvm.reject = 0
  lillie.reject = 0
  pearson.reject = 0
  sf.reject = 0
  for (i in seq(1, T)){
    # x = runif(n, 0, 1)
    # x = rcauchy(n, 0, 1)
    # x = rgamma(n, 5, 3)
    # x = rexp(n, 1)
    # x = rchisq(n, 10)
    x = rt(n, 10)
    ad.reject = ad.reject + (ad.test(x)$p.value < alpha)
    cvm.reject = cvm.reject + (cvm.test(x)$p.value < alpha)
    lillie.reject = lillie.reject + (lillie.test(x)$p.value < alpha)
    pearson.reject = pearson.reject + (pearson.test(x)$p.value < alpha)
    sf.reject = sf.reject + (sf.test(x)$p.value < alpha)
  }
  results[1, k] = 1 - (ad.reject + 0.0) / T
  results[2, k] = 1 - (cvm.reject + 0.0) / T
  results[3, k] = 1 - (lillie.reject + 0.0) / T
  results[4, k] = 1 - (pearson.reject + 0.0) / T
  results[5, k] = 1 - (sf.reject + 0.0) / T
}
print(results)
xtable(results, digits=4)

# Type II Error: large sample property

T = 10000

set.seed(123469)
vec = rep(0, 100)
results = data.frame(vec, vec, vec, vec, vec, vec)
colnames(results) = c('n', 'ad', 'cvm', 'lillie', 'pearson', 'sf')

for (para in seq(1, 200) * 10){
  print(para)
  ad.reject = 0
  cvm.reject = 0
  lillie.reject = 0
  pearson.reject = 0
  sf.reject = 0
  for (i in seq(1, T)){
    n = para
    # x = runif(n, 0, 1)
    # x = rcauchy(n, 0, 1)
    # x = rgamma(n, 5, 3)
    # x = rexp(n, 1)
    # x = rchisq(n, 10)
    x = rt(n, 10)
    ad.reject = ad.reject + (ad.test(x)$p.value < alpha)
    cvm.reject = cvm.reject + (cvm.test(x)$p.value < alpha)
    lillie.reject = lillie.reject + (lillie.test(x)$p.value < alpha)
    pearson.reject = pearson.reject + (pearson.test(x)$p.value < alpha)
    sf.reject = sf.reject + (sf.test(x)$p.value < alpha)
  }
  results[para/10, 1] = para
  results[para/10, 2] = 1 - (ad.reject + 0.0) / T
  results[para/10, 3] = 1 - (cvm.reject + 0.0) / T
  results[para/10, 4] = 1 - (lillie.reject + 0.0) / T
  results[para/10, 5] = 1 - (pearson.reject + 0.0) / T
  results[para/10, 6] = 1 - (sf.reject + 0.0) / T
}
results.trans = melt(results[1:200,], id="n", value.name='T2.Error')
g = ggplot(data=results.trans, aes(x=n, y=T2.Error, colour=variable))
g = g + geom_line()

# calculate power

# mixture
T = 10000

set.seed(123469)
vec = rep(0, 51)
results = data.frame(vec, vec, vec, vec, vec, vec)
colnames(results) = c('pi', 'ad', 'cvm', 'lillie', 'pearson', 'sf')

n = 1000
for (para in seq(1, 51)){
  print(para)
  ad.reject = 0
  cvm.reject = 0
  lillie.reject = 0
  pearson.reject = 0
  sf.reject = 0
  for (i in seq(1, T)){
    z = rbinom(n, 1, (para - 1.0) / 100)
    x1 = rnorm(n, 0, 1)
    x2 = rnorm(n, 5, 1)
    x = x1
    x[z == 1] = x2[1:sum(z)]
    ad.reject = ad.reject + (ad.test(x)$p.value < alpha)
    cvm.reject = cvm.reject + (cvm.test(x)$p.value < alpha)
    lillie.reject = lillie.reject + (lillie.test(x)$p.value < alpha)
    pearson.reject = pearson.reject + (pearson.test(x)$p.value < alpha)
    sf.reject = sf.reject + (sf.test(x)$p.value < alpha)
  }
  results[para, 1] = (para - 1) / 100
  results[para, 2] = (ad.reject + 0.0) / T
  results[para, 3] = (cvm.reject + 0.0) / T
  results[para, 4] = (lillie.reject + 0.0) / T
  results[para, 5] = (pearson.reject + 0.0) / T
  results[para, 6] = (sf.reject + 0.0) / T
}
results.trans = melt(results[1:51,], id="pi", value.name='power')
g = ggplot(data=results.trans, aes(x=pi, y=power, colour=variable))
g = g + geom_line()
g


# student's t approximation
T = 10000

set.seed(123469)
vec = rep(0, 50)
results = data.frame(vec, vec, vec, vec, vec, vec)
colnames(results) = c('df', 'ad', 'cvm', 'lillie', 'pearson', 'sf')

n = 1000
for (para in seq(1, 200)){
  print(para)
  ad.reject = 0
  cvm.reject = 0
  lillie.reject = 0
  pearson.reject = 0
  sf.reject = 0
  for (i in seq(1, T)){
    x = rt(n, para)
    ad.reject = ad.reject + (ad.test(x)$p.value < alpha)
    cvm.reject = cvm.reject + (cvm.test(x)$p.value < alpha)
    lillie.reject = lillie.reject + (lillie.test(x)$p.value < alpha)
    pearson.reject = pearson.reject + (pearson.test(x)$p.value < alpha)
    sf.reject = sf.reject + (sf.test(x)$p.value < alpha)
  }
  results[para, 1] = para
  results[para, 2] = (ad.reject + 0.0) / T
  results[para, 3] = (cvm.reject + 0.0) / T
  results[para, 4] = (lillie.reject + 0.0) / T
  results[para, 5] = (pearson.reject + 0.0) / T
  results[para, 6] = (sf.reject + 0.0) / T
}
results.trans = melt(results[1:200,], id="df", value.name='power')
g = ggplot(data=results.trans, aes(x=df, y=power, colour=variable))
g = g + geom_line()
g


# noise
T = 10000

set.seed(123469)
vec = rep(0, 100)
results = data.frame(vec, vec, vec, vec, vec, vec)
colnames(results) = c('epsilon', 'ad', 'cvm', 'lillie', 'pearson', 'sf')

n = 30
for (para in seq(1, 101)){
  print(para)
  ad.reject = 0
  cvm.reject = 0
  lillie.reject = 0
  pearson.reject = 0
  sf.reject = 0
  for (i in seq(1, T)){
    x = rnorm(n, 0, 1)
    x = x + (para - 1) / 100 * sin(x)
    ad.reject = ad.reject + (ad.test(x)$p.value < alpha)
    cvm.reject = cvm.reject + (cvm.test(x)$p.value < alpha)
    lillie.reject = lillie.reject + (lillie.test(x)$p.value < alpha)
    pearson.reject = pearson.reject + (pearson.test(x)$p.value < alpha)
    sf.reject = sf.reject + (sf.test(x)$p.value < alpha)
  }
  results[para, 1] = (para - 1) / 100
  results[para, 2] = (ad.reject + 0.0) / T
  results[para, 3] = (cvm.reject + 0.0) / T
  results[para, 4] = (lillie.reject + 0.0) / T
  results[para, 5] = (pearson.reject + 0.0) / T
  results[para, 6] = (sf.reject + 0.0) / T
}
results.trans = melt(results[1:101,], id="epsilon", value.name='power')
g = ggplot(data=results.trans, aes(x=epsilon, y=power, colour=variable))
g = g + geom_line()
g
