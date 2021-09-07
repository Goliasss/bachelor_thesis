library(extraDistr)
library(stats)
set.seed(1)

# US Army 
dataset <- read.csv(file = 'ANSURIIMALEPublic.csv')
army_height = dataset$stature/10 # To have centimeters

# Data exploration 
hist(army_height, prob=TRUE, xlab = "Height [cm]", main = "", col = "bisque", 
     breaks = seq(from=145, to=210, by=5))
summary(army_height)
length(army_height)
var(army_height)
sd(army_height)

# Setting up truncation boundaries
a = 147
b = 180

# Truncating and exploring the truncated data
truncated_height = army_height[army_height<b & army_height>a]
hist(truncated_height, prob = TRUE,xlab = "Height [cm]", main = "", 
     col = "bisque", breaks = seq(from=145, to=210, by=5))
summary(truncated_height)
var(truncated_height)
length(truncated_height)

# Manual fitting
log_likelihood<-function(theta,y){
  mu<-theta[1]
  sigma2<-theta[2]
  n<-length(y)
  logl<- -0.5*n*log(2*pi) -0.5*n*log(sigma2) - n*log(pnorm((b-mu)/sqrt(sigma2)) - pnorm((a-mu)/sqrt(sigma2))) -
    (1/(2*sigma2))*sum((y-mu)**2)
  return(logl)
}
(manual_estimator = optim(c(mean(truncated_height),var(truncated_height)),
                          log_likelihood,y=truncated_height, method="BFGS", 
                          control = list(fnscale=-1)))

manual_estimator$par[1] # fitted mean
sqrt(manual_estimator$par[2]) # fitted standard deviation

# Plotting fitted truncated normal density over histogram
hist(truncated_height, ylim =c(0,0.1),prob = TRUE,xlab = "Height [cm]", main = "", col = "bisque", breaks = seq(from=145, to=210, by=5))
curve(dtnorm(x, manual_estimator$par[1], sqrt(manual_estimator$par[2]),a,b),col = "red", add = TRUE, lwd=2)

# Confidence region -- Fisher information matrix
(mu = manual_estimator$par[1])
(sigma = sqrt(manual_estimator$par[2]))
n = length(army_height)
chis = c(qchisq(0.99, df=2),qchisq(0.95, df=2),qchisq(0.9, df=2))
beta = (b-mu)/sigma
alfa = (a-mu)/sigma

J_11 = (1/sigma^2)*( 1 + (alfa*dnorm(alfa) - beta*dnorm(beta))/(pnorm(beta)-pnorm(alfa)) -
                      ((dnorm(alfa)-dnorm(beta)) / 
                                    (pnorm(beta)-pnorm(alfa)))^2  
                    )
J_22 = (1/(2*sigma^4))*( 1 + 
                        (-dnorm(beta)*beta+dnorm(alfa)*alfa)/(pnorm(beta)-pnorm(alfa)) -
                        (1/2)*((-dnorm(beta)*beta+dnorm(alfa)*alfa)/(pnorm(beta)-pnorm(alfa)))^2 - 
                        (1/2)*((dnorm(beta)*(3*beta+beta^3) - dnorm(alfa)*(3*alfa+alfa^3)) /
                                 (pnorm(beta) - pnorm(alfa)))
                        )
J_12 = (1/(sigma^3))*( (-dnorm(beta)+dnorm(alfa))/(pnorm(beta)-pnorm(alfa)) +
                         (1/2)*((dnorm(beta)*(1+beta^2) - dnorm(alfa)*(1+alfa^2)) /
                                  (pnorm(beta) - pnorm(alfa))) -
                         (1/2)*(((dnorm(alfa) - dnorm(beta))*(alfa*dnorm(alfa)-beta*dnorm(beta))) /
                            ((pnorm(beta)-pnorm(alfa))^2))
                       )

# Confidence region -- plotting
x<-seq(174,176,length=1000)
y<-seq(30,55,length=1000)
z<-outer(x,y,function(x,y) n*(J_11*x^2 + J_22*y^2 + 2*J_12*x*y + J_11*mu^2 - 
          2*J_11*x*mu - 2*J_12*y*mu + J_22*sigma^4 -2*J_12*x*sigma^2 -
          2*J_22*y*sigma^2 + 2*J_12*mu*sigma^2) )
contour(x,y,z,levels=chis, labels = c("0.99","0.95","0.9"), drawlabels = TRUE, 
        col = "red", xlab=expression(mu), ylab=expression(sigma^2))
points(mu, sigma^2,pch=10, col="red")
abline(v=mu,col="brown",lty=2)
abline(h=sigma^2,col="brown",lty=2)
