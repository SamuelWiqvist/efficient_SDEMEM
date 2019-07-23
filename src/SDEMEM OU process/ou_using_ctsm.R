# load packages 
library(ctsmr)
library(sde)

# set correct dir 
setwd("~/Documents/projects/corrPMMH/code")

# simulate data 

t1_true = 0.5 # ground-truth parameter values 
t2_true = 4 
t3_true = 0.1
sigma_true = 0.2
dt = 0.1

# simulate measurment process
set.seed(123)
x_path = sde.sim(t0 = 0, 
                 T = 20.000000, 
                 X0 = 0, 
                 N = 199, 
                 delta = 0.1, 
                 drift = expression(t1_true*(t2_true-x)), 
                 sigma = expression(t3_true), 
                  sigma.x= expression(0)) 

y_path = x_path + rnorm(200, mean = 0, sd = sigma_true)

t = 0:199/10

# plot latent process and obs process
par(mfrow=c(1,2))
plot(t, x_path,"l")
plot(t, y_path,"l")

# path to data data/SDEMEM OU/kalman/data_ctms.csv
data <- read.csv(file="data/SDEMEM OU/data_julia_small_theta_3.csv", header=TRUE, sep=",")
names(data) = c("t", "Y")

# save data 
data = data.frame(t, y_path)
names(data) = c("t", "Y")
write.csv(data, file = "data/SDEMEM OU/data_ctms_small_t3.csv",row.names=FALSE)

# initlize model 
model <- ctsm()

# set model equations 
model$addSystem(dX ~ t1*(t2-X)*dt + t3*dw1)
model$addObs(Y ~ X)
model$setVariance(Y ~ sigma^2)

# initilize parameters and set parameter restrictions 
model$setParameter(t1 = c(init = 1, lower = 0, upper = 100))
model$setParameter(t2 = c(init = 1, lower = -10, upper = 100))
model$setParameter(t3 = c(init = 1, lower = 0, upper = 100))
model$setParameter(sigma = c(init = 1, lower = 0, upper = 100))

# initlize data and set restrictions 
model$setParameter( X = c(init=0, lower=-10, upper=1E3) )

plot(data$t,data$Y,"l")

# run estimation
fit <- model$estimate(data=data)

fit

# print parameter estimations 
summary(fit, extended=TRUE)
