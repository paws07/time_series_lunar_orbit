## Clear console
rm(list=ls())

##Check added packages 

library(features)
library(TSA)


##Loading the Lunar data into R.
luna<-read.csv("luna.csv")

##Setting the numbers up as a vector.
luna.nrs<-as.vector(na.omit(luna[,2]))

##Creating a time series out of the given data.
luna.ts<-ts(luna.nrs, start=c(1951,1,1),freq=365)

##Checking and plotting time series 
luna.ts
plot(luna.ts)

##Truncating a window in the time series for better analysis and plotting it.

luna.tr<-window(luna.ts, end=c(1953,12))
plot(luna.tr)
luna.tr2<-window(luna.ts, end=c(1954,12)) ## For the future. 
plot(luna.tr2)

##Getting the time part of luna.tr
t<-time(luna.tr)

##Generating a spectrum of the time series to observe relevant frequencies and get the corresponding frequencies and getting a summary of it. (See text)  
spec<-spectrum(luna.ts)
spec<-spectrum(luna.tr)

##Using the features package to extract all the critcial points from the raw peridogram of the series. 
luna.feat<-features(spec$freq, spec$spec)
luna.fget<-fget(luna.feat)


##For critical points; luna.fget$crit.pts[i]

###### Regression fit

##Attempt at regression:(frequencies through luna.fget)
luna.fget$crit.pts[1]


## Fitting the data

##best fit around the first critical point, corresponds with the period of the moon.
# for (i in 1200:1500){
#   w2=i/100
#   w=w2*2*pi
#   luna.sinemod <- lm(luna.tr ~  cos(t*w)+sin(t*w))
#   if ((summary(luna.sinemod)$r.squared) > .9799){
#     print(w2)
#     print(summary(luna.sinemod)$r.squared)
# 
#   }
# }

w=13.25*pi*2

##Looking for the best second frequency [Long computation, do not execute]
# for (i in 1:36000){
#   w2=i/100
#   w=w2*2*pi
#   luna.sinemod <- lm(luna.tr ~  cos(t*w1)+sin(t*w1)+ cos(t*w) +sin(t*w) + cos(2*t*w) + sin(2*t*w) + cos(3*t*w) + sin(3*t*w) + cos(4*t*w) + sin(4*t*w) + cos(5*t*w) + sin(5*t*w) +   cos(6*t*w) + sin(6*t*w) + cos(7*t*w) + sin(7*t*w) +
#                        cos(8*t*w) + sin(8*t*w))
#   if ((summary(luna.sinemod)$r.squared) > .9799){
#     print(w2)
#     print(summary(luna.sinemod)$r.squared)
# 
#   }
# }

##Choosing 1.91 and the harmonics with the highest P-values (Directly or using a for loop with [if (summary(model1)$coef[,"Pr(>|t|)"][i]<.01) ]

w1=1.91*pi*2

##Further analysis shows a third frequency 
w2=3.53*pi*2

luna.sinemod <- lm(luna.tr ~  cos(t*w)+sin(t*w)+ cos(6*t*w1) + sin(6*t*w1) + sin(7*t*w1)+cos(3*t*w2) + sin(4*t*w2) +  cos(7*t*w2) + sin(7*t*w2))
summary(luna.sinemod)

## R^2=0.99

plot(t, luna.sinemod$fitted.values, type="l")


## Getting a line
luna.sinemodline <- luna.sinemodline <- luna.sinemod$coefficients[1]+luna.sinemod$coefficients[2]*cos(t*w)+luna.sinemod$coefficients[3]*sin(t*w)+luna.sinemod$coefficients[4]*cos(6*t*w1)+luna.sinemod$coefficients[5]*sin(6*t*w1)+ luna.sinemod$coefficients[6]*sin(7*t*w1)+luna.sinemod$coefficients[7]*cos(3*t*w2) + luna.sinemod$coefficients[8]*sin(4*t*w2) +  luna.sinemod$coefficients[9]*cos(7*t*w2) + luna.sinemod$coefficients[10]*sin(7*t*w2)

##Checking regression 
plot(luna.sinemodline)
plot(luna.tr)
lines(luna.sinemodline, col=2)

##Future Prediction for a different window (The year of 1953)

##Plotting with increased window size

plot(luna.tr, xlim=c(1951, 1954))

##LM model:

lines(luna.sinemodline, col=4)

luna.predts <- window(luna.ts, start=c(1953,12), end=c(1954,12))
t2 <- time(luna.predts)


luna.predline <-  luna.sinemod$coefficients[1]+luna.sinemod$coefficients[2]*cos(t2*w)+luna.sinemod$coefficients[3]*sin(t2*w)+luna.sinemod$coefficients[4]*cos(6*t2*w1)+luna.sinemod$coefficients[5]*sin(6*t2*w1)+ luna.sinemod$coefficients[6]*sin(7*t2*w1)+luna.sinemod$coefficients[7]*cos(3*t2*w2) + luna.sinemod$coefficients[8]*sin(4*t2*w2) +  luna.sinemod$coefficients[9]*cos(7*t2*w2) + luna.sinemod$coefficients[10]*sin(7*t2*w2)
lines(luna.predline, col = 2)


##Residual 
luna.res = luna.tr-luna.sinemodline

plot(luna.res)
plot(decompose(luna.res))

acf(luna.res)


##ARIMA

get.best.arima <- function(x.ts, maxord = c(1,1,1,1,1,1))
{
  best.aic <- 1e8
  n <- length(x.ts)
  
  for(p in 0:maxord[1])
    for (d in 0:maxord[2])
      for(q in 0:maxord[3])
        for(P in 0:maxord[4])
          for(D in 0:maxord[5])
            for(Q in 0:maxord[6])
            {
              try( 
                {
                  fit <- arima(x.ts, order=c(p,d,q), seas = list(order=c(P,D,Q), frequency(x.ts) ), method = "CSS") 
                  fit.aic <- -2*fit$loglik + (log(n)+1) * length(fit$coef) #suppressing code after the + gives better fits (using loglik as the only criterion that way)
                  if (fit.aic < best.aic)
                  {
                    best.aic <- fit.aic
                    best.fit <- fit
                    best.model <- c(p,d,q,P,D,Q)
                  } #end if
                } # end first argument of try
                , FALSE) #end try
              
              print(c(p,d,q,P,D,Q))
              flush.console()
            } # end for	
  
  
  
  dev.new()
  
  
  over <- paste("Process Fit with ARIMA(", toString(best.model[1]), ",", toString(best.model[2]), ",", toString(best.model[3]), ") Process. \n Coefficients:", toString(round(best.fit$coef, digits=3))) 
  under <- paste("Periodic Coefficients:", toString(best.model[4]), ",", toString(best.model[5]), ",", toString(best.model[6]))
  
  acf(na.omit(best.fit$resid), lag.max=100, main=over, xlab=under)
  
  
  
  list(akaike=best.aic, data=best.fit, orders=best.model)
} # end get.best.arima



luna.arima <- get.best.arima(luna.res, c(10,2,10,0,0,0) )

## Fit with ARIMA (7,0,10)

luna.arma <- luna.res - luna.arima$data$residuals
##Checking ARIMA progression plot(luna.arma, xlim=c(1951, 1954))

luna.final <- (luna.sinemodline + luna.arma)

plot(luna.tr)
lines(luna.final, col=4)

##Divide the window into two parts
par(mfrow=c(2,1))

## Plot over the line luna.tr
plot(luna.tr, type="l") 
lines(luna.final, col =2)

##plot over the points in luna.tr
plot(luna.tr, type="p") 
lines(luna.final, col =2)

##Reset the window ratio
par(mfrow=c(1,1))


plot(luna.tr2)
lines(luna.final, col =4)

##Adding the ARIMA approximations 

luna.predarima <- predict(luna.arima$data, n.ahead = length(t2))

lines(luna.predline+luna.predarima$pred, col=2)

 plot(luna.predts)

 lunapredlm <- lm(luna.predts ~ luna.predline)
 summary(lunapredlm)


##ACF residual 
acf(luna.final-luna.tr)
pacf((luna.final-luna.tr))
acf(luna.predline+luna.predarima$pred-luna.predts)


##############################

plot(luna.tr2)
##Orcale is the last period of the data (a larger underlying period)
luna.oracle <- window(luna.ts, start=c(1953.43), end=c(1954.03))
plot(luna.oracle)

plot(luna.tr)
##Almanac would be the year 1952
luna.almanac <- window(luna.ts, start=c(1952.3), end=c(1952.9))

plot(luna.almanac)

##Seeing them side by side
par(mfrow=c(2,1))

plot(luna.oracle)
plot(luna.almanac)

par(mfrow=c(1,1))

##Function 

a=c(1,2,3,4,5)
b=c(3,4,5,6,7)

quality <- function(x, y) {
  x <- as.vector(x)
  y <- as.vector(y)
  i=0
  j=0
    for(n in 1:length(x)){
      i=i+abs(x[n]-y[n])
      j=j+abs(y[n])
    }
  print((i/j)*100)
}

quality(a, b)


## Testing linear significance 
tx <- time(luna.ts)
lunatslm <- lm (luna.ts ~ cos(tx*2*w)+sin(tx*2*w))
summary(lunatslm)
lunatslm$coef[1]
# 
# 
# model1 <-  lm(luna.ts~cos(tx*2*w)+sin(tx*2*w))
# model2 <-  lm(luna.ts~ (0.002573513+cos(tx*2*w)+sin(tx*2*w)))
# anova(model1,model2)

## The temporal pseudoreplication command (ANOVA) did not work out in terms of my code
print(mean(luna.ts)-lunatslm$coef[1])

##The intercept (average distance to the moon at 0 AD (x=0)?) is too small to consider)

# lunalmtry <- lm(luna.ts ~ cos(tx*w)+sin(tx*w)+ cos(6*tx*w1) + sin(6*tx*w1) + sin(7*tx*w1)+cos(3*tx*w2) + sin(4*tx*w2) +  cos(7*tx*w2) + sin(7*tx*w2))
# summary(lunalmtry)
