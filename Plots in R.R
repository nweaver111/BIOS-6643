#### This R Script will create nice visuals for the mean estimate over time for each model in my project ####

##Set directory to access the data:
setwd("C:\\Users\\weavenic\\Documents\\Analysis of Longitudinal\\SASUniversityedition\\myfolders\\Project")
dat <- read.csv("Theoph.csv", header = T, as.is = T)

#Create a sequence so that I can plot a pseudo curve representing mean effect for a model
x <- seq(from = 0, to = 24, length.out = 1000)

#Create a function that represents model 1 (piecwise cubic regression with knot at mean maximum concentration value (1.5)). Numbers come from SAS output.
mod1 <- function(x, s, d, w) {
  reg1 = -5.3046 + (11.0276)*(x) + 0.7159*(x)*d + -0.07769*x*x*d + 0.002052*x*x*x*d -7.3608*x*x + 1.2584*x*x*x
  reg2 = -2.7889*max((x - s),0) - 3.6779*max((x - s),0)*max((x - s),0) + -1.2680*max((x - s),0)*max((x - s),0)*max((x - s),0) 
  y = reg1 + reg2 + 0.1368*d + 0.06685*w
  return(y)
}

#Create a function that represents model 2 (gamma-like nonlinear mixed model). Numbers come from SAS output.
mod2 <- function(x, d, w) {
  y = 6.9664*x^(0.5732)*exp(-x*0.2098) + 0.2807*d + -0.00336*w
  return(y)
}

#Create a function that represents model 3 (one-compartment model). Numbers come from SAS output.
mod3 <- function(x,d){
  y = 0.08499*d*49.9690/(0.08499-1.9500)*(exp(-1.9500*x)-exp(-.08499*x))
  return(y)
}

#I needed to run model 1 in a loop function to ensure that values are calculated correctly (this was my best work-around for an error I kept getting when trying
# to run a vector into the function)
for (i in 1:length(x)){
  tmp <- mod1(x[i], 3, d = mean(dat$Dose), w = mean(dat$Wt))
  y1[i] <- tmp
}

#Create mean concentration for each time period (X) for models 2 and 3:
y2 <- mod2(x, d = mean(dat$Dose), w = mean(dat$Wt))
y3 <- mod3(x, d = mean(dat$Dose))

#In order to save as a data base with numeric values I use the following work-around where I write the data frame to a csv file and the read it as-is back into R.
dat2 <- as.data.frame(cbind(x,y1,y2, y3))
write.csv(dat2, file = "forR.csv", row.names = F)
dat3 <- read.csv("forR.csv", header = T, as.is = T)

#Now we plot the mean effect for all three models:
library(ggplot2)

ggplot(dat, aes(x = IdenticalTime, y = conc)) + geom_point(aes(color = as.factor(Subject))) + geom_point(data = dat3, aes(x = x, y = y1)) + geom_line(aes(group = as.factor(Subject), color = as.factor(Subject))) +
  labs(title = "piecewise cubic regression fit", y = "Concentration", x = "Time") 

#It is clear I picked the wrong knot and that cubic terms are NOT needed. Try 1.5 instead with piecwise squared regression!
mod1_updated <- function(x, s, d, w) {
  reg1 = -7.1527 + (11.8403)*(x) + 0.2285*(x)*d + -0.00916*x*x*d -4.9043*x*x
  reg2 = 1.4406*max((x - s),0) + 4.95549*max((x - s),0)*max((x - s),0) 
  y = reg1 + reg2 + 0.5570*d + 0.06685*w
  return(y)
}

### Update everything:
for (i in 1:length(x)){
  tmp <- mod1_updated(x[i], 1.5, d = mean(dat$Dose), w = mean(dat$Wt))
  y1[i] <- tmp
}

dat2 <- as.data.frame(cbind(x,y1,y2, y3))
write.csv(dat2, file = "forR.csv", row.names = F)
dat3 <- read.csv("forR.csv", header = T, as.is = T)

#Plot again:
ggplot(dat, aes(x = IdenticalTime, y = conc)) + geom_point(aes(color = as.factor(Subject))) + geom_point(data = dat3, aes(x = x, y = y1)) + geom_line(aes(group = as.factor(Subject), color = as.factor(Subject))) +
  labs(title = "piecewise cubic regression fit", y = "Concentration", x = "Time") 

#Looks much better!

#Plot model 2:
ggplot(dat, aes(x = IdenticalTime, y = conc)) + geom_point(aes(color = as.factor(Subject))) + geom_point(data = dat3, aes(x = x, y = y2)) + geom_line(aes(group = as.factor(Subject), color = as.factor(Subject))) +
  labs(title = "Gamma regression fit", y = "Concentration", x = "Time")

#Plot model 3:
ggplot(dat, aes(x = IdenticalTime, y = conc)) + geom_point(aes(color = as.factor(Subject))) + geom_point(data = dat3, aes(x = x, y = y3)) + geom_line(aes(group = as.factor(Subject), color = as.factor(Subject))) +
  labs(title = "Pharmacokinetic regression fit", y = "Concentration", x = "Time")

#It is possible to create subject-specific curves for the models by using EB estimates for the RE's and plotting each new curve. It is done below.

rand_eff_mod_1 <- c(1.2982, -0.05665, 0.1507, 0.02642, 0.3503, -1.2183, -0.9392, -0.5519, 0.4847, 0.6273, -0.4150, 0.2435)
rand_eff_mod_2 <- c(1.4651214438, -0.061894453, 0.1430467816, 0.0455273814, 0.4255095098, -1.136198573, -1.060775654, -0.603032455, 0.3452632094, 0.6354732729, -0.490586701, 0.231517603)
rand_eff_mod_3 <- matrix(c(-0.016515737, -0.250623685, 21.672420047,
                           0.0065457044, 0.1940563111, 1.26674436648,
                           0.0003827421, 0.5360525223, 11.01950098,
                           -0.001008228, -0.495800256, -10.26806722,
                           0.0006978232, -0.327151983, -12.10110863,
                           0.0013575673, -0.457285375, -18.39708688,
                           0.0016108084, -1.149712357, -33.10276611,
                           0.0007987044, -0.437484296, -16.25759362,
                           0.0145171866, 2.0585717551, 64.716552065,
                           -0.010065732, -1.107117062, -24.69119488,
                           0.0086186938, 1.335576046, 10.287443298,
                           0.0057016642, -0.998102913, -25.49423263), nrow = 3, ncol = 12)

#Create a function that represents model 1 and allows for RE for time (piecwise cubic regression with knot at mean maximum concentration value (1.5)). Numbers come from SAS output.
mod1_re <- function(x, s, d, w, re) {
  reg1 = -7.1527 + re + (11.8403)*(x) + 0.2285*(x)*d + -0.00916*x*x*d -4.9043*x*x
  reg2 = 1.4406*max((x - s),0) + 4.95549*max((x - s),0)*max((x - s),0) 
  y = reg1 + reg2 + 0.5570*d + 0.06685*w
  return(y)
}

#Create a matrix to store subject specific trend lines:
y1_re <- matrix(nrow = length(x), ncol = 12)

#Find subject-specific trend lines and assign to matrix:
for (j in 1:12){
  for (i in 1:length(x)){
    tmp <- mod1_re(x[i], 1.5, d = mean(dat$Dose), w = mean(dat$Wt), re = rand_eff_mod_1[j])
    y1_re[i,j] <- tmp
  }
}

#Append x values for plotting purposes
y1_re <- cbind(x, y1_re)

#Work-around to get numeric elements in the Dataframe:
write.csv(y1_re, file = "forR_re.csv", row.names = F)
dat4 <- read.csv("forR_re.csv", header = T, as.is = T)

#Plot:
ggplot(dat, aes(x = IdenticalTime, y = conc)) + geom_point(aes(color = as.factor(Subject))) + geom_line(data = dat4, aes(x = x, y = dat4[,2], color = "1")) + 
  geom_line(data = dat4, aes(x = x, y = dat4[,3], color = "2")) +
  geom_line(data = dat4, aes(x = x, y = dat4[,4], color = "3")) +
  geom_line(data = dat4, aes(x = x, y = dat4[,5], color = "4")) +
  geom_line(data = dat4, aes(x = x, y = dat4[,6], color = "5")) +
  geom_line(data = dat4, aes(x = x, y = dat4[,7], color = "6")) +
  geom_line(data = dat4, aes(x = x, y = dat4[,8], color = "7")) +
  geom_line(data = dat4, aes(x = x, y = dat4[,9], color = "8")) +
  geom_line(data = dat4, aes(x = x, y = dat4[,10], color = "9")) +
  geom_line(data = dat4, aes(x = x, y = dat4[,11], color = "10")) +
  geom_line(data = dat4, aes(x = x, y = dat4[,12], color = "11")) +
  geom_line(data = dat4, aes(x = x, y = dat4[,13], color = "12")) +
  geom_point(aes(group = as.factor(Subject), color = as.factor(Subject))) +
  labs(title = "piecewise cubic regression fit", y = "Concentration", x = "Time") 

#It is clear by this plot that mixed models have the effect of shrinking towards the mean estimate (for subject-specific slopes). Also,
# it appears that future work may benefit from messing around with the RE that is used. This could get somewhat involved.

### Mod 2 subject-specific:
mod2_re <- function(x, d, w, re) {
  y = 6.9664*x^(0.5732)*exp(-x*0.2098) + 0.2807*d + -0.00336*w + re
  return(y)
}

y2_re <- matrix(nrow = length(x), ncol = 12)
for(j in 1:12){
  y2_re[,j] <- mod2_re(x, d = mean(dat$Dose), w = mean(dat$Wt), re = rand_eff_mod_2[j])
}
y2_re <- cbind(x, y2_re)

#Work-around to get numeric elements in the Dataframe:
write.csv(y2_re, file = "forR_re2.csv", row.names = F)
dat5 <- read.csv("forR_re2.csv", header = T, as.is = T)

ggplot(dat, aes(x = IdenticalTime, y = conc)) + geom_point(aes(color = as.factor(Subject))) + geom_line(data = dat5, aes(x = x, y = dat5[,2], color = "1")) + 
  geom_line(data = dat5, aes(x = x, y = dat5[,3], color = "2")) +
  geom_line(data = dat5, aes(x = x, y = dat5[,4], color = "3")) +
  geom_line(data = dat5, aes(x = x, y = dat5[,5], color = "4")) +
  geom_line(data = dat5, aes(x = x, y = dat5[,6], color = "5")) +
  geom_line(data = dat5, aes(x = x, y = dat5[,7], color = "6")) +
  geom_line(data = dat5, aes(x = x, y = dat5[,8], color = "7")) +
  geom_line(data = dat5, aes(x = x, y = dat5[,9], color = "8")) +
  geom_line(data = dat5, aes(x = x, y = dat5[,10], color = "9")) +
  geom_line(data = dat5, aes(x = x, y = dat5[,11], color = "10")) +
  geom_line(data = dat5, aes(x = x, y = dat5[,12], color = "11")) +
  geom_line(data = dat5, aes(x = x, y = dat5[,13], color = "12")) +
  geom_point(aes(group = as.factor(Subject), color = as.factor(Subject))) +
  labs(title = "Gamma Like Regression", y = "Concentration", x = "Time") 

### Mod 3 subject-specific:
mod3_re <- function(x,d,re){
  y = (0.08499 + re[1])*d*(49.9690+ re[3])/((0.08499 + re[1])-(1.9500 + re[2]))*(exp(-(1.9500+ re[2])*x)-exp(-(.08499+ re[1])*x))
  return(y)
}

y3_re <- matrix(nrow = length(x), ncol = 12)
for(j in 1:12){
  y3_re[,j] <- mod3_re(x, d = mean(dat$Dose), re = rand_eff_mod_3[,j])
}
y3_re <- cbind(x, y3_re)

#Work-around to get numeric elements in the Dataframe:
write.csv(y3_re, file = "forR_re3.csv", row.names = F)
dat6 <- read.csv("forR_re3.csv", header = T, as.is = T)

ggplot(dat, aes(x = IdenticalTime, y = conc)) + geom_point(aes(color = as.factor(Subject))) + geom_line(data = dat6, aes(x = x, y = dat6[,2], color = "1")) + 
  geom_line(data = dat6, aes(x = x, y = dat6[,3], color = "2")) +
  geom_line(data = dat6, aes(x = x, y = dat6[,4], color = "3")) +
  geom_line(data = dat6, aes(x = x, y = dat6[,5], color = "4")) +
  geom_line(data = dat6, aes(x = x, y = dat6[,6], color = "5")) +
  geom_line(data = dat6, aes(x = x, y = dat6[,7], color = "6")) +
  geom_line(data = dat6, aes(x = x, y = dat6[,8], color = "7")) +
  geom_line(data = dat6, aes(x = x, y = dat6[,9], color = "8")) +
  geom_line(data = dat6, aes(x = x, y = dat6[,10], color = "9")) +
  geom_line(data = dat6, aes(x = x, y = dat6[,11], color = "10")) +
  geom_line(data = dat6, aes(x = x, y = dat6[,12], color = "11")) +
  geom_line(data = dat6, aes(x = x, y = dat6[,13], color = "12")) +
  geom_point(aes(group = as.factor(Subject), color = as.factor(Subject))) +
  labs(title = "One-Compartment Regression", y = "Concentration", x = "Time") 

#we finally get very different subject-specific curves. The trend lines don't line up with the observed points too well though. Maybe there is room for model improvement?