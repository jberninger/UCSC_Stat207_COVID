## Stat 207
## Final Take Home
## Covid Analysis
####################################################################################################################################
library(readr)
library(dplyr)
library(truncnorm)
library(MASS)
library(tmvtnorm)
library(invgamma)
library(mvtnorm)
library(janitor)
library(LearnBayes)
set.seed(11)
library(latex2exp)
library(tmvtnorm)
library(LearnBayes)
library(plyr)
library(dplyr)
library(janitor)
library(readr)
library(ggmap)
library(ggplot2)
library(gridExtra)
library(ggmap)
library(maps)
library(mapdata)
library(geoR)
setwd("~/Desktop/Stat 207/Final/")
# dataset
data<-read_csv("CountyData.csv")%>%clean_names()
####################################################################################################################################
# name the fields
covid<-data%>%dplyr::mutate(ny=total_cases,
                     y=log(total_cases),
                     cy=population_original,
                     county=tolower(county))%>%
  rename(c(group="census.group"))
ny<-covid$total_cases
ny[ny==0]<-0.01
y<-log(ny) # has some -Inf values that must be handled
#y[y==-Inf]<--9999
cy<-covid$population_original
d<-covid$density
group<-covid$census.group
N<-nrow(covid)
####################################################################################################################################
# EDA 
## Perform an exploratory data analysis involving 
## yij (log cases), the five different groups, the population, and the population density. 
## Discuss possible associations and clusters.
# data prep for the map
{
states <- map_data("state")
cali <- subset(states, region == "california")
counties <- map_data("county")
ca_county <- subset(counties, region == "california") %>% 
  inner_join(covid, by = c("subregion" = "county")) %>%
  mutate(death.rate = 100*round(deaths/population_original, 10))

ditch_the_axes <- theme(
  axis.text = element_blank(),
  axis.line = element_blank(),
  axis.ticks = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.title = element_blank()
)

ca_base <- ggplot(data = cali, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(1.3) + 
  geom_polygon(color = "black", fill = "gray")
}
# maps
{
# Death rate
ca_base + 
  geom_polygon(data = ca_county, aes(fill =ordered(round(100*death.rate, 2))), color = "white") +
  geom_polygon(color = "black", fill = NA) +
  ggtitle("Death Rates by County") +
  theme_bw() +
  labs(fill = "Deaths / Infected (%)") +
  ditch_the_axes
# Regions
ca_base + 
  geom_polygon(data = ca_county, aes(fill =factor(census.group)), color = 'white') +
  geom_polygon(color = "black", fill = NA) +
  ggtitle("Census Groups in CA") +
  theme_bw() +
  labs(fill = "Census Group") +
  ditch_the_axes
# Population Density
ca_base + 
  geom_polygon(data = ca_county, aes(fill =ordered(round(density/100, 0)))) +
  geom_polygon(color = "black", fill = NA) +
  ggtitle("Population Density by County") +
  theme_bw() +
  labs(fill = "Scaled Density (Density/100)") +
  ditch_the_axes
}
# plots
p1<-ggplot(data=covid)+geom_point(aes(x=population_original,y=density,color=factor(census.group)),size=2) +
  ggtitle("Density vs Population") +
  labs(color = "Census Group") + 
  xlab("Population")
p2<-ggplot(data=covid)+geom_point(aes(x=log(population_original),y=log(density),color=factor(census.group)),size=2)+
  ggtitle("log(Density) vs log(Population)") +
  labs(color = "Census Group") 
xlab("log(Population)")
grid.arrange(p1, p2, nrow=2)

# population vs log cases
plot(cy,y,col=group,type="p",pch=4)
legend(group)
plot(log(cy),y,col=group,type="p",pch=4)

ggplot(data=covid)+geom_point(aes(x=log(population_original),y=log(total_cases),color=factor(census.group)),size=2) +
  ggtitle("Response vs log(population)") +
  labs(color = "Census Group") + 
  xlab("log(population)") +
  ylab("log(cases)")

ggplot(data=covid)+geom_point(aes(x=log(density),y=log(total_cases),color=factor(census.group)),size=2) +
  ggtitle("Response vs log(density)") +
  labs(color = "Census Group") + 
  xlab("log(density)") +
  ylab("log(cases)")

# population density vs log cases
plot(d,y,col=group,type="p",pch=4)
plot(log(d),y,col=group,type="p")

ggplot(data=covid,aes(x=census.group,y=log(density),group=census.group,fill=factor(census.group)))+geom_boxplot()+
  theme(legend.position = "none") +
  ggtitle("Distribution of log(Density) by Group")+
  xlab("Census group")+
  ylab("log(Density)")
ggplot(data=covid,aes(x=census.group,y=log(population_original),group=census.group,fill=factor(census.group)))+geom_boxplot()+
  theme(legend.position = "none") +
  ggtitle("Distribution of log(Population) by Group")+
  xlab("Census group")+
  ylab("log(Population)")
####################################################################################################################################
# model 1 : y = \mu + e
X=matrix(rep(1,N),ncol=1)
Q=diag(10^3/cy)
beta.hat<-solve(t(X)%*%solve(Q)%*%X)%*%t(X)%*%solve(Q)%*%y
V.beta<-solve(t(X)%*%solve(Q)%*%X)
# mu
beta.hat.sample<-rmvnorm(n=1000,mean=beta.hat,sigma=V.beta)
hist(beta.hat.sample,main=TeX('Model 1 Posterior Distribution for $\\mu$'),xlab=TeX("$\\mu$"))
abline(v=mean(beta.hat.sample),col="red")
abline(v=quantile(beta.hat.sample,0.05),col="blue")
abline(v=quantile(beta.hat.sample,0.95),col="blue")
# sigma sq
k<-1
s.sq<-(1/(N-k))*t(y-X%*%beta.hat)%*%solve(Q)%*%(y-X%*%beta.hat)
sigma.sq.sample<-geoR::rinvchisq(n=1000,df=N-k,scale=s.sq)
hist(sigma.sq.sample,main=TeX('Model 1 Posterior Distribution for $\\sigma^2$'),xlab=TeX("$\\sigma^2$"))
abline(v=mean(sigma.sq.sample),col="red")
abline(v=quantile(sigma.sq.sample,0.05),col="blue")
abline(v=quantile(sigma.sq.sample,0.95),col="blue")

# g prior
L=chol(Q)
W=solve(L)%*%X
Z=solve(L)%*%y
beta.hat<-solve(t(W)%*%W)%*%t(W)%*%Z
V.beta<-solve(t(W)%*%W)
# mu
beta.hat.sample<-rmvnorm(n=10000,mean=beta.hat,sigma=V.beta)
hist(beta.hat.sample[,1])
abline(v=mean(beta.hat.sample[,1]),col="red")
# R squared
Yhat= X %*% beta.hat
# Compute Residuals
Resid=y-Yhat
SSE=t(Resid) %*% Resid
R.sq.1 <- 1-SSE/sum((y-mean(y))^2)
R.sq.1 # it does worse than the mean, it can be negative because it does not have a centering parameter

plot(log(covid$population_original),Resid,xlab="log(population)",main="Model 1 Residuals vs log(Population)") 
# clear trend in the residuals, this model is garbage, as expected need to include population
mu<-scale(W[,1])
X<-scale(W[,-1])
Z<-scale(Z) # doesnt work if Z isnt scaled, but it def should be
G<-N
# Double check this in office hour
beta.hat.g<-solve(t(X)%*%X)%*%t(X)%*%Z 
# is this the right beta.hat? think so
# they have shrunk from the other ones
sigma.sq.sample.g<-rinvgamma(n=1000,N/2,0.5*var(Z)+(1/(2*G+1))*t(beta.hat.g)%*%t(X)%*%X%*%beta.hat.g)
# it might not be correct that I had to sample sigma.squared in order to sample beta

mu.samples<-rnorm(n=1000,mean=mean(y),sd=sd(y)/sqrt(N))

hist(mu.samples,main=TeX('Model 1 Posterior Distribution for $\\mu$'),xlab=TeX("$\\mu$"))
abline(v=mean(mu.samples),col="red")
abline(v=quantile(mu.samples,0.05),col="blue")
abline(v=quantile(mu.samples,0.95),col="blue")
####################################################################################################################################
# model 2 : y = \mu + B*d + e
X=matrix(c(rep(1,N),d),ncol=2)
Q=diag(10^3/cy)
beta.hat<-solve(t(X)%*%solve(Q)%*%X)%*%t(X)%*%solve(Q)%*%y
V.beta<-solve(t(X)%*%solve(Q)%*%X)
# mu
beta.hat.sample<-rmvnorm(n=1000,mean=beta.hat,sigma=V.beta)
hist(beta.hat.sample[,1])
abline(v=lm(y~1+d,weights=1/sqrt(10^3/cy))$coef[1],col="red")
# beta
hist(beta.hat.sample[,2])
abline(v=lm(y~1+d,weights=1/sqrt(10^3/cy))$coef[2],col="red")
# sigma sq
k<-2
s.sq<-(1/(N-k))*t(y-X%*%beta.hat)%*%solve(Q)%*%(y-X%*%beta.hat)
sigma.sq.sample<-geoR::rinvchisq(n=1000,df=N-k,scale=s.sq)
hist(sigma.sq.sample)
# unequal var
# g prior
L=chol(Q)
W=solve(L)%*%X
Z=solve(L)%*%y
beta.hat<-solve(t(W)%*%W)%*%t(W)%*%Z
V.beta<-solve(t(W)%*%W)
# mu
beta.hat.sample<-rmvnorm(n=10000,mean=beta.hat,sigma=V.beta)
hist(beta.hat.sample[,1])
abline(v=mean(beta.hat.sample[,1]),col="red")
# G-prior
{
  # G Prior
  # design matrix needs to have mu removed and colmeans of 0
  # should the mu be scaled?
  mu<-scale(W[,1])
  X<-scale(W[,-1])
  Z<-scale(Z) # doesnt work if Z isnt scaled, but it def should be
  G<-N
  # Double check this in office hour
  beta.hat.g<-solve(t(X)%*%X)%*%t(X)%*%Z 
  # is this the right beta.hat? think so
  # they have shrunk from the other ones
  sigma.sq.sample.g<-rinvgamma(n=1000,N/2,0.5*var(Z)+(1/(2*G+1))*t(beta.hat.g)%*%t(X)%*%X%*%beta.hat.g)
  # it might not be correct that I had to sample sigma.squared in order to sample beta
  beta.samples.g<-matrix(NA,1000,length(beta.hat.g))
  for(i in 1:1000){
    beta.samples.g[i,]<-rmvnorm(n=1,mean=(G/(1+G))*beta.hat.g,sigma=(G/(1+G))*sigma.sq.sample.g[i]*solve(t(X)%*%X))
  }
  mu.samples<-rnorm(n=1000,mean=mean(y),sd=sd(y)/sqrt(N))

  
  hist(beta.samples.g[,1],main=TeX('Model 2 Posterior Distribution for $\\beta$'),xlab=TeX("$\\beta$"))
  abline(v=quantile(beta.samples.g[,1],0.05),col="blue")
  abline(v=quantile(beta.samples.g[,1],0.95),col="blue")
  abline(v=mean(beta.samples.g[,1]),col="red")

  hist(mu.samples,main=TeX('Model 2 Posterior Distribution for $\\mu$'),xlab=TeX("$\\mu$"))
  abline(v=mean(mu.samples),col="red")
  abline(v=quantile(mu.samples,0.05),col="blue")
  abline(v=quantile(mu.samples,0.95),col="blue")

  
  # G Prior R^2
  Yhat= X %*% beta.hat.g
  # Compute Residuals
  Resid=Z-Yhat
  SSE=t(Resid) %*% Resid
  R.sq.2 <- 1-SSE/sum((Z-mean(Z))^2)
}
R.sq.2
plot(log(covid$population_original),Resid,xlab="log(population)",main="Model 2 Residuals vs log(Population)")
# residuals increasing wrt log(pop) and are heteroskedastic
plot(Yhat,Resid,xlab=TeX("$\\hat{y}$"),main="Model 2 Residuals vs Fitted Values")
####################################################################################################################################
# model 3 : y = \mu + eta + e
S=matrix(0,5,4)
diag(S)=rep(1,4)
S[5,]=rep(-1,4)
k=5
# Define Design Matrix
X=matrix(0,N,k)
X[,1]=1
X[,2]=ifelse(group==1,1,0)
X[,3]=ifelse(group==2,1,0)
X[,4]=ifelse(group==3,1,0)
X[,5]=ifelse(group==4,1,0)
for(i in 1:N){
  if(group[i]==5){X[i,2:5]=-1}
}

Q=diag(10^3/cy)
beta.hat<-solve(t(X)%*%solve(Q)%*%X)%*%t(X)%*%solve(Q)%*%y
V.beta<-solve(t(X)%*%solve(Q)%*%X)
# mu
beta.hat.sample<-rmvnorm(n=10000,mean=beta.hat,sigma=V.beta)
hist(beta.hat.sample[,1])
abline(v=mean(beta.hat.sample[,1]),col="red")

hist(beta.hat.sample[,2])
abline(v=mean(beta.hat.sample[,2]),col="red")
# abline(v=lm(y~1+factor(group),weights=1/sqrt(10^3/cy))$coef[2],col="red")
# regression model puts the group 1 factor into the intercept term
# beta
# unequal var
L=chol(Q)
W=solve(L)%*%X
Z=solve(L)%*%y
beta.hat<-solve(t(W)%*%W)%*%t(W)%*%Z
V.beta<-solve(t(W)%*%W)
# mu
beta.hat.sample<-rmvnorm(n=10000,mean=beta.hat,sigma=V.beta)
hist(beta.hat.sample[,1])
abline(v=mean(beta.hat.sample[,1]),col="red")
# G prior
{
# model: L^-1y = L^-1XB + v, v \sim N(0,I)
# L = cholesky demcomp of the error covariance matrix
# beta hat is given by
# 3 steps from slides:
# (1) cholesky decomp of V (error covaraince)
# (2) solve LW = X and LZ=y
# (3) compute LSE
# G Prior
# design matrix needs to have mu removed and colmeans of 0
# should the mu be scaled?
mu<-scale(W[,1])
X<-scale(W[,-1])
Z<-scale(Z) # doesnt work if Z isnt scaled, but it def should be
G<-22
G<-N
# Double check this in office hour
beta.hat.g<-solve(t(X)%*%X)%*%t(X)%*%Z 
# is this the right beta.hat? think so
# they have shrunk from the other ones
sigma.sq.sample.g<-rinvgamma(n=1000,N/2,0.5*var(Z)+(1/(2*G+1))*t(beta.hat.g)%*%t(X)%*%X%*%beta.hat.g)
# it might not be correct that I had to sample sigma.squared in order to sample beta
beta.samples.g<-matrix(NA,1000,length(beta.hat.g))
for(i in 1:1000){
  beta.samples.g[i,]<-rmvnorm(n=1,mean=(G/(1+G))*beta.hat.g,sigma=(G/(1+G))*sigma.sq.sample.g[i]*solve(t(X)%*%X))
}

par(mfrow=c(2,2))
for(i in 1:6){
  hist(beta.samples.g[,i])
  abline(v=mean(beta.samples.g[,i]),col="red")
  abline(v=quantile(beta.samples.g[,i],0.05),col="blue")
  abline(v=quantile(beta.samples.g[,i],0.95),col="blue")
}
dev.off()
etas<-cbind(beta.samples.g,-rowSums(beta.samples.g))
boxplot(etas,use.cols=TRUE,
        main=TeX("Model 3 Posterior Distribution of Group Coefficients"),xlab=TeX("$\\eta _i$"))
abline(h=0,col="red")

mu.samples<-rnorm(n=1000,mean=mean(y),sd=sd(y)/sqrt(N))
hist(mu.samples,main=TeX('Model 2 Posterior Distribution for $\\mu$'),xlab=TeX("$\\mu$"))
abline(v=mean(mu.samples),col="red")
abline(v=quantile(mu.samples,0.05),col="blue")
abline(v=quantile(mu.samples,0.95),col="blue")
abline(v=mean(mu.samples[,1]),col="red")
# distribution for mu is not changing?

# G Prior R^2
Yhat= X %*% beta.hat.g
# Compute Residuals
Resid=Z-Yhat
SSE=t(Resid) %*% Resid
R.sq.3 <- 1-SSE/sum((Z-mean(Z))^2)
}
R.sq.3

plot(log(covid$population_original),Resid,xlab="log(population)",main="Model 3 Residuals vs log(Population)") 
# residuals clearly have a trend and are heteroskedastic
plot(Yhat,Resid,xlab=TeX("$\\hat{y}$"),main="Model 3 Residuals vs Fitted Values")
####################################################################################################################################
# model 4 : y = \mu + B*d + eta_j + e
k=6
# Define Design Matrix
X=matrix(0,N,k)
X[,1]=1
X[,2]=d
X[,3]=ifelse(group==1,1,0)
X[,4]=ifelse(group==2,1,0)
X[,5]=ifelse(group==3,1,0)
X[,6]=ifelse(group==4,1,0)
for(i in 1:N){
  if(group[i]==5){X[i,3:6]=-1}
}

Q=diag(10^3/cy)
beta.hat<-solve(t(X)%*%solve(Q)%*%X)%*%t(X)%*%solve(Q)%*%y
V.beta<-solve(t(X)%*%solve(Q)%*%X)
# mu
beta.hat.sample<-rmvnorm(n=10000,mean=beta.hat,sigma=V.beta)
hist(beta.hat.sample[,1])
abline(v=mean(beta.hat.sample[,1]),col="red")

par(mfrow=c(2,3))
for(i in 1:5){
  hist(beta.hat.sample[,1+i])
  abline(v=mean(beta.hat.sample[,1+i]),col="red")
}
# one of the params post distributions is centered on zero.
# G prior
{
# model: L^-1y = L^-1XB + v, v \sim N(0,I)
# L = cholesky demcomp of the error covariance matrix
# beta hat is given by
# 3 steps from slides:
# (1) cholesky decomp of V (error covaraince)
# (2) solve LW = X and LZ=y
# (3) compute LSE
L=chol(Q)
W=solve(L)%*%X
Z=solve(L)%*%y
# ok this seems to work
# confirmed it re-captures beta.hat
# G Prior
# design matrix needs to have mu removed and colmeans of 0
# should the mu be scaled?
mu<-scale(W[,1])
X<-scale(W[,-1])
Z<-scale(Z) # doesnt work if Z isnt scaled, but it def should be
G<-32
# Double check this in office hour
beta.hat.g<-solve(t(X)%*%X)%*%t(X)%*%Z 
# is this the right beta.hat? think so
# they have shrunk from the other ones
sigma.sq.sample.g<-rinvgamma(n=1000,N/2,0.5*var(Z)+(1/(2*G+1))*t(beta.hat.g)%*%t(X)%*%X%*%beta.hat.g)
# it might not be correct that I had to sample sigma.squared in order to sample beta
beta.samples.g<-matrix(NA,1000,length(beta.hat.g))
for(i in 1:1000){
  beta.samples.g[i,]<-rmvnorm(n=1,mean=(G/(1+G))*beta.hat.g,sigma=(G/(1+G))*sigma.sq.sample.g[i]*solve(t(X)%*%X))
}
par(mfrow=c(2,3))
for(i in 1:6){
  hist(beta.samples.g[,i])
  abline(v=mean(beta.samples.g[,i]),col="red")
}
# G Prior R^2
Yhat= X %*% beta.hat.g
# Compute Residuals
Resid=Z-Yhat
SSE=t(Resid) %*% Resid
R.sq.4 <- 1-SSE/sum((Z-mean(Z))^2)
}
R.sq.4
# plots
dev.off()
etas<-cbind(beta.samples.g[,2:5],-rowSums(beta.samples.g[,2:5]))
boxplot(etas,use.cols=TRUE,
        main=TeX("Model 4 Posterior Distribution of Group Coefficients"),xlab=TeX("$\\eta _i$"))
abline(h=0,col="red")

hist(beta.samples.g[,2],main=TeX('Model 4 Posterior Distribution for $\\beta$'),xlab=TeX("$\\beta$"))
abline(v=mean(beta.samples.g[,2]),col="red")
abline(v=quantile(beta.samples.g[,2],0.05),col="blue")
abline(v=quantile(beta.samples.g[,2],0.95),col="blue")

plot(log(covid$population_original),Resid,xlab="log(population)",main="Model 4 Residuals vs log(Population)") 
# errors have a trend
plot(Yhat,Resid,xlab=TeX("$\\hat{y}$"),main="Model 4 Residuals vs Fitted Values")
####################################################################################################################################
# model 5 : y = \mu + B*d +C*log(d) + e
X=matrix(c(rep(1,N),d,log(d)),ncol=3)
Q=diag(10^3/cy)
beta.hat<-solve(t(X)%*%solve(Q)%*%X)%*%t(X)%*%solve(Q)%*%y
V.beta<-solve(t(X)%*%solve(Q)%*%X)
# mu
beta.hat.sample<-rmvnorm(n=1000,mean=beta.hat,sigma=V.beta)
hist(beta.hat.sample[,1])
abline(v=lm(y~1+d,weights=1/sqrt(10^3/cy))$coef[1],col="red")
# beta
hist(beta.hat.sample[,2])
# C
hist(beta.hat.sample[,3])
# sigma sq
k<-3
s.sq<-(1/(N-k))*t(y-X%*%beta.hat)%*%solve(Q)%*%(y-X%*%beta.hat)
sigma.sq.sample<-geoR::rinvchisq(n=1000,df=N-k,scale=s.sq)
hist(sigma.sq.sample)
# unequal var
# g prior
L=chol(Q)
W=solve(L)%*%X
Z=solve(L)%*%y
beta.hat<-solve(t(W)%*%W)%*%t(W)%*%Z
V.beta<-solve(t(W)%*%W)
# mu
beta.hat.sample<-rmvnorm(n=10000,mean=beta.hat,sigma=V.beta)
hist(beta.hat.sample[,1])
abline(v=mean(beta.hat.sample[,1]),col="red")
# G-prior
{
  # G Prior
  # design matrix needs to have mu removed and colmeans of 0
  # should the mu be scaled?
  mu<-scale(W[,1])
  X<-scale(W[,-1])
  Z<-scale(Z) # doesnt work if Z isnt scaled, but it def should be
  G<-643
  # Double check this in office hour
  beta.hat.g<-solve(t(X)%*%X)%*%t(X)%*%Z 
  # is this the right beta.hat? think so
  # they have shrunk from the other ones
  sigma.sq.sample.g<-rinvgamma(n=1000,N/2,0.5*var(Z)+(1/(2*G+1))*t(beta.hat.g)%*%t(X)%*%X%*%beta.hat.g)
  # it might not be correct that I had to sample sigma.squared in order to sample beta
  beta.samples.g<-matrix(NA,1000,length(beta.hat.g))
  for(i in 1:1000){
    beta.samples.g[i,]<-rmvnorm(n=1,mean=(G/(1+G))*beta.hat.g,sigma=(G/(1+G))*sigma.sq.sample.g[i]*solve(t(X)%*%X))
  }
  mu.samples<-rnorm(n=1000,mean=mean(y),sd=sd(y)/sqrt(N))
  
  
  hist(beta.samples.g[,1],main=TeX('Model 5 Posterior Distribution for $\\beta$'),xlab=TeX("$\\beta$"))
  abline(v=quantile(beta.samples.g[,1],0.05),col="blue")
  abline(v=quantile(beta.samples.g[,1],0.95),col="blue")
  abline(v=mean(beta.samples.g[,1]),col="red")
  # quadratic term
  hist(beta.samples.g[,2],main=TeX('Model 5 Posterior Distribution for $\\gamma$'),xlab=TeX("$\\gamma$"))
  abline(v=quantile(beta.samples.g[,2],0.05),col="blue")
  abline(v=quantile(beta.samples.g[,2],0.95),col="blue")
  abline(v=mean(beta.samples.g[,2]),col="red")
  
  hist(mu.samples,main=TeX('Model 2 Posterior Distribution for $\\mu$'),xlab=TeX("$\\mu$"))
  abline(v=mean(mu.samples),col="red")
  abline(v=quantile(mu.samples,0.05),col="blue")
  abline(v=quantile(mu.samples,0.95),col="blue")
  
  
  # G Prior R^2
  Yhat= X %*% beta.hat.g
  # Compute Residuals
  Resid=Z-Yhat
  SSE=t(Resid) %*% Resid
  R.sq.5 <- 1-SSE/sum((Z-mean(Z))^2)
}
R.sq.5
plot(log(covid$population_original),Resid,xlab="log(population)",main="Model 5 Residuals vs log(Population)")
plot(Yhat,Resid,xlab=TeX("$\\hat{y}$"),main="Model 5 Residuals vs Fitted Values")
####################################################################################################################################
# model 6 : y = \mu + Beta*d +gamma*log(d) + eta_j + e
k=7
# Define Design Matrix
X=matrix(0,N,k)
X[,1]=1
X[,2]=d
X[,3]=log(d)
X[,4]=ifelse(group==1,1,0)
X[,5]=ifelse(group==2,1,0)
X[,6]=ifelse(group==3,1,0)
X[,7]=ifelse(group==4,1,0)
for(i in 1:N){
  if(group[i]==5){X[i,3:6]=-1}
}

Q=diag(10^3/cy)
beta.hat<-solve(t(X)%*%solve(Q)%*%X)%*%t(X)%*%solve(Q)%*%y
V.beta<-solve(t(X)%*%solve(Q)%*%X)
# mu
beta.hat.sample<-rmvnorm(n=10000,mean=beta.hat,sigma=V.beta)
hist(beta.hat.sample[,1])
abline(v=mean(beta.hat.sample[,1]),col="red")

par(mfrow=c(2,3))
for(i in 1:5){
  hist(beta.hat.sample[,1+i])
  abline(v=mean(beta.hat.sample[,1+i]),col="red")
}
# one of the params post distributions is centered on zero.
# G prior
{
  # model: L^-1y = L^-1XB + v, v \sim N(0,I)
  # L = cholesky demcomp of the error covariance matrix
  # beta hat is given by
  # 3 steps from slides:
  # (1) cholesky decomp of V (error covaraince)
  # (2) solve LW = X and LZ=y
  # (3) compute LSE
  L=chol(Q)
  W=solve(L)%*%X
  Z=solve(L)%*%y
  # ok this seems to work
  # confirmed it re-captures beta.hat
  # G Prior
  # design matrix needs to have mu removed and colmeans of 0
  # should the mu be scaled?
  mu<-scale(W[,1])
  X<-scale(W[,-1])
  Z<-scale(Z) # doesnt work if Z isnt scaled, but it def should be
  G<-273
  # Double check this in office hour
  beta.hat.g<-solve(t(X)%*%X)%*%t(X)%*%Z 
  # is this the right beta.hat? think so
  # they have shrunk from the other ones
  sigma.sq.sample.g<-rinvgamma(n=1000,N/2,0.5*var(Z)+(1/(2*G+1))*t(beta.hat.g)%*%t(X)%*%X%*%beta.hat.g)
  # it might not be correct that I had to sample sigma.squared in order to sample beta
  beta.samples.g<-matrix(NA,1000,length(beta.hat.g))
  for(i in 1:1000){
    beta.samples.g[i,]<-rmvnorm(n=1,mean=(G/(1+G))*beta.hat.g,sigma=(G/(1+G))*sigma.sq.sample.g[i]*solve(t(X)%*%X))
  }
  par(mfrow=c(2,3))
  for(i in 1:6){
    hist(beta.samples.g[,i])
    abline(v=mean(beta.samples.g[,i]),col="red")
  }
  # G Prior R^2
  Yhat= X %*% beta.hat.g
  # Compute Residuals
  Resid=Z-Yhat
  SSE=t(Resid) %*% Resid
  R.sq.6 <- 1-SSE/sum((Z-mean(Z))^2)
}
R.sq.6
# plots
dev.off()

hist(beta.samples.g[,1],main=TeX('Model 6 Posterior Distribution for $\\beta$'),xlab=TeX("$\\beta$"))
abline(v=quantile(beta.samples.g[,1],0.05),col="blue")
abline(v=quantile(beta.samples.g[,1],0.95),col="blue")
abline(v=mean(beta.samples.g[,1]),col="red")
# quadratic term
hist(beta.samples.g[,2],main=TeX('Model 6 Posterior Distribution for $\\gamma$'),xlab=TeX("$\\gamma$"))
abline(v=quantile(beta.samples.g[,2],0.05),col="blue")
abline(v=quantile(beta.samples.g[,2],0.95),col="blue")
abline(v=mean(beta.samples.g[,2]),col="red")

etas<-cbind(beta.samples.g[,3:6],-rowSums(beta.samples.g[,3:6]))
boxplot(etas,use.cols=TRUE,
        main=TeX("Model 6 Posterior Distribution of Group Coefficients"),xlab=TeX("$\\eta _i$"))
abline(h=0,col="red")

hist(beta.samples.g[,2],main=TeX('Model 6 Posterior Distribution for $\\gamma$'),xlab=TeX("$\\gamma$"))
abline(v=mean(beta.samples.g[,2]),col="red")
abline(v=quantile(beta.samples.g[,2],0.05),col="blue")
abline(v=quantile(beta.samples.g[,2],0.95),col="blue")

plot(log(covid$population_original),Resid,xlab="log(population)",main="Model 6 Residuals vs log(Population)") 
# errors have a trend
plot(Yhat,Resid,xlab=TeX("$\\hat{y}$"),main="Model 6 Residuals vs Fitted Values")
####################################################################################################################################
# Bayes Factor for G-prior model comparison
# change G
# compare models with the reference prior?
# F-test for nested models?
# posterior prediction?
# LaTeX?
p1<-1
p2<-2
p3<-5
p4<-6
p5<-3
p6<-7
BF.2.0<-(1+G)^((N-p2-1)/2)/(1+G*(1-R.sq.2))^((N-1)/2)
BF.3.0<-(1+G)^((N-p3-1)/2)/(1+G*(1-R.sq.3))^((N-1)/2)
BF.4.0<-(1+G)^((N-p4-1)/2)/(1+G*(1-R.sq.4))^((N-1)/2)
BF.5.0<-(1+G)^((N-p4-1)/2)/(1+G*(1-R.sq.5))^((N-1)/2)
# calling model 1 as Model 0 here, since it is the null model w/ intercept only
# formulas on the G-prior slides
BF.2.3<-(1+G)^((p3-p2)/2)*((1+G*(1-R.sq.3))/(1+G*(1-R.sq.2)))^((N-1)/2)
BF.2.4<-(1+G)^((p4-p2)/2)*((1+G*(1-R.sq.4))/(1+G*(1-R.sq.2)))^((N-1)/2)
BF.2.5<-(1+G)^((p5-p2)/2)*((1+G*(1-R.sq.5))/(1+G*(1-R.sq.2)))^((N-1)/2)
BF.2.6<-(1+G)^((p6-p2)/2)*((1+G*(1-R.sq.6))/(1+G*(1-R.sq.2)))^((N-1)/2)

BF.3.4<-(1+G)^((p4-p3)/2)*((1+G*(1-R.sq.4))/(1+G*(1-R.sq.3)))^((N-1)/2)
BF.3.5<-(1+G)^((p5-p3)/2)*((1+G*(1-R.sq.5))/(1+G*(1-R.sq.3)))^((N-1)/2)
BF.3.6<-(1+G)^((p6-p3)/2)*((1+G*(1-R.sq.6))/(1+G*(1-R.sq.3)))^((N-1)/2)

BF.4.5<-(1+G)^((p5-p4)/2)*((1+G*(1-R.sq.5))/(1+G*(1-R.sq.4)))^((N-1)/2)
BF.4.6<-(1+G)^((p6-p4)/2)*((1+G*(1-R.sq.6))/(1+G*(1-R.sq.4)))^((N-1)/2)

BF.5.6<-(1+G)^((p6-p5)/2)*((1+G*(1-R.sq.6))/(1+G*(1-R.sq.5)))^((N-1)/2)

BF.2.0
BF.3.0
BF.4.0
BF.5.0
BF.2.3
BF.2.4
BF.2.5
BF.2.6
BF.3.4
BF.3.5
BF.3.6
BF.4.5
BF.4.6
BF.5.6

# 5 is the best now
R.sq.2
R.sq.3
R.sq.4
R.sq.5
R.sq.6
####################################################################################################################################
# Conclusions
# are there other models worth considering?
# is the grouping significant? 
# Is the population density significant?
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Residual plots for all the models
# Bayes Factor for the uninformative reference prior models?
# Posterior distributions of beta - do the credible intervals include 0?
# bayes factor does not include the prior information
# marginal does include the prior information
####################################################################################################################################
# marginal2
m2<-((gamma((N-1)/0.5)/(pi^((N+1)/2)*sqrt(N)))*norm(Z-mean(Z),type="2")^(1-N))*(((1+G)^((N-1-p2))/2)/(1+G*(1-R.sq.2))^((N-1)/2))
m3<-((gamma((N-1)/0.5)/(pi^((N+1)/2)*sqrt(N)))*norm(Z-mean(Z),type="2")^(1-N))*(((1+G)^((N-1-p3))/2)/(1+G*(1-R.sq.3))^((N-1)/2))
m4<-((gamma((N-1)/0.5)/(pi^((N+1)/2)*sqrt(N)))*norm(Z-mean(Z),type="2")^(1-N))*(((1+G)^((N-1-p4))/2)/(1+G*(1-R.sq.4))^((N-1)/2))
m5<-((gamma((N-1)/0.5)/(pi^((N+1)/2)*sqrt(N)))*norm(Z-mean(Z),type="2")^(1-N))*(((1+G)^((N-1-p5))/2)/(1+G*(1-R.sq.5))^((N-1)/2))
m6<-((gamma((N-1)/0.5)/(pi^((N+1)/2)*sqrt(N)))*norm(Z-mean(Z),type="2")^(1-N))*(((1+G)^((N-1-p6))/2)/(1+G*(1-R.sq.6))^((N-1)/2))
m2
m3
m4
m5
m6
# marginal for 5 is the biggest?
####################################################################################################################################
# G chosen from empirical bayes
g.eb.2<-max((R.sq.2/p2)/((1-R.sq.2)/(N-1-p2))-1,0)
g.eb.3<-max((R.sq.3/p3)/((1-R.sq.3)/(N-1-p3))-1,0)
g.eb.4<-max((R.sq.4/p4)/((1-R.sq.4)/(N-1-p4))-1,0)
g.eb.5<-max((R.sq.5/p5)/((1-R.sq.5)/(N-1-p5))-1,0)
g.eb.6<-max((R.sq.6/p6)/((1-R.sq.6)/(N-1-p6))-1,0)

BF.eb.2.0<-(1+g.eb.2)^((N-p2-1)/2)/(1+g.eb.2*(1-R.sq.2))^((N-1)/2)
BF.eb.3.0<-(1+g.eb.3)^((N-p3-1)/2)/(1+g.eb.3*(1-R.sq.3))^((N-1)/2)
BF.eb.4.0<-(1+g.eb.4)^((N-p4-1)/2)/(1+g.eb.4*(1-R.sq.4))^((N-1)/2)
BF.eb.5.0<-(1+g.eb.5)^((N-p5-1)/2)/(1+g.eb.5*(1-R.sq.5))^((N-1)/2)
BF.eb.6.0<-(1+g.eb.6)^((N-p6-1)/2)/(1+g.eb.6*(1-R.sq.6))^((N-1)/2)

BF.eb.2.3<-BF.eb.2.0/BF.eb.3.0
BF.eb.2.4<-BF.eb.2.0/BF.eb.4.0
BF.eb.2.5<-BF.eb.2.0/BF.eb.5.0
BF.eb.2.6<-BF.eb.2.0/BF.eb.6.0
BF.eb.3.4<-BF.eb.3.0/BF.eb.4.0
BF.eb.3.5<-BF.eb.3.0/BF.eb.5.0
BF.eb.3.6<-BF.eb.3.0/BF.eb.6.0
BF.eb.4.5<-BF.eb.4.0/BF.eb.5.0
BF.eb.4.6<-BF.eb.4.0/BF.eb.6.0
BF.eb.5.6<-BF.eb.5.0/BF.eb.6.0
BF.eb.2.3
BF.eb.2.4
BF.eb.2.5
BF.eb.2.6
BF.eb.3.4
BF.eb.3.5
BF.eb.3.6
BF.eb.4.5
BF.eb.4.6
BF.eb.5.6
# model 5 is best by the G-prior chosen by empirical bayes
####################################################################################################################################
# hyper g-prior
a<-4 # how to choose a?
((gamma((N-1)/0.5)/(pi^((N+1)/2)*sqrt(N)))*norm(Z-mean(Z),type="2")^(1-N))*((a-2)/(p2+a-2))*hypergeo((N-1)/2,1,(p2+a)/2,R.sq.2)

((gamma((N-1)/0.5)/(pi^((N+1)/2)*sqrt(N)))*norm(Z-mean(Z),type="2")^(1-N))*((a-2)/(p3+a-2))*hypergeo((N-1)/2,1,(p3+a)/2,R.sq.3)
# marginal for hyper-g
