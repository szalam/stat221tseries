setwd('D:/PhD/Coursework/Course/statistics/time_series/Final Project')
library(PerformanceAnalytics)
library(xts)
library(lubridate)
library(astsa)

#read data
df=read.csv('Traffic_collision_dataset.csv')
head(df)
df=df[,1:7]
head(df)
df$Date.Occurred=as.Date(as.character(df$Date.Occurred),format="%m/%d/%Y")
df=cbind(df,month(df[,1]),year(df[,1]))
head(df)
colnames(df)=c('Date','Collision','Time','Sex','Male','Female','Area','Month','Year')


# barplot
barplot(prop.table(table(df$Area)))
library(ggplot2)
ggplot(data.frame(df),aes(x=Area))+geom_bar()+theme_bw()+
   theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()
# counts
ggp + geom_histogram(fill="lightgreen")

ggp + geom_histogram(fill="lightblue",aes(y=..count../sum(..count..)))
#---------------------
collision_all=df[,c(1,2,8,9)]
head(collision_all)
collision_daily=aggregate(data=collision_all,Collision~(Date),FUN=sum)
head(collision_daily)
plot(collision_daily[,1],collision_daily[,2],ylab='Total Collision')
lines(collision_daily[,1],collision_daily[,2])




collision_monthly=aggregate(data=collision_all,Collision~month(Date)+year(Date),FUN=sum)
monthly_date=seq(as.Date("2010-01-15"), as.Date("2017-12-15"), by="month")
head(collision_monthly)
collision_monthly=cbind(monthly_date,collision_monthly)
plot(collision_monthly[,1],collision_monthly[,4],ylab='Total Collision')
lines(collision_monthly[,1],collision_monthly[,4])
abline(reg=lm(collision_monthly[,4]~collision_monthly[,1]))
class(collision_monthly[,1])
monthly_ts=ts(collision_monthly[,4], frequency=12, start=c(2010,1))
plot.ts(monthly_ts)


collision_yearly=aggregate(data=collision_all,Collision~year(Date),FUN=sum)
yearly_date=seq(as.Date("2010-01-15"), as.Date("2017-12-15"), by="year")
head(collision_yearly)
collision_yearly=cbind(yearly_date,collision_yearly)
plot(collision_yearly[,2],collision_yearly[,3])
lines(collision_yearly[,2],collision_yearly[,3])

#Auto correlation-ACF and PACF
acf(collision_daily[,2], lag.max=800, plot=T) 

acf(collision_monthly[,2],lag.max = 100, plot=T,main='ACF of original data') 
pacf(collision_monthly[,2],lag.max = 100, plot=T,main='PACF of original data') 

#Decomposing data
library(TTR)
kingstimeseriesSMA8 <- SMA(monthly_ts,n=8)
plot.ts(kingstimeseriesSMA8)

birthstimeseries=monthly_ts
birthstimeseriescomponents <- decompose(monthly_ts)
plot(birthstimeseriescomponents)# it shows treand seasonal and random component
# removing seasonal component
birthstimeseriesseasonallyadjusted <- birthstimeseries - birthstimeseriescomponents$seasonal
plot(birthstimeseriesseasonallyadjusted,ylab='Collisions_adjusted',main='Data after seasonality removed')


acf(birthstimeseriesseasonallyadjusted,lag.max = 1000,main='Seasonality removed')
pacf(birthstimeseriesseasonallyadjusted,lag.max = 1000,main='Seasonality removed')

# making dataset stationary by differencing n times
plot.ts(monthly_ts)
skirtsseriesdiff1 <- diff(monthly_ts, differences=1)
#By taking the time series of first differences, we have removed the trend component of the time series
#and are left with an irregular component. We can now examine whether there are correlations between successive terms of this irregular component; if so, this could help us to make a predictive model 
#for the ages at death of the kings
plot.ts(skirtsseriesdiff1,main='Detrend data',ylab='Collisions')
acf(skirtsseriesdiff1,main='Detrend data')
pacf(skirtsseriesdiff1,main='Detrend data')

#----------------------AIC BIC thing
spaic = spec.ar(monthly_ts, log="no") # min AIC spec
spaic = spec.ar(skirtsseriesdiff1, log="no") # min AIC spec
abline(v=frequency(skirtsseriesdiff1)*1/21, lty=3) # El Nino peak
(sun.ar = ar(skirtsseriesdiff1, order.max=30)) # estimates and AICs
dev.new()
plot(1:30, sun.ar$aic[-1], type="o")

temp_data=monthly_ts
n = length(temp_data)
AIC = rep(0, 30) -> AICc -> BIC
for (k in 1:30){
  sigma2 = ar(temp_data, order=k, aic=FALSE)$var.pred
  BIC[k] = log(sigma2) + (k*log(n)/n)
  AICc[k] = log(sigma2) + ((n+k)/(n-k-2))
  AIC[k] = log(sigma2) + ((n+2*k)/n) }
IC = cbind(AIC, BIC+1)
ts.plot(IC, type="o", xlab="p", ylab="AIC / BIC")


#Periodogram

pgram <- spec.pgram(as.vector(monthly_ts))     
abline(v=1/12, col="green")              
abline(v=2/12, col="green")            
abline(v=3/12, col="green")            
abline(v=4/12, col="green")            
abline(v=1/(5*12), col="blue", lty=4)   
# Reducing noise
pgram <- spec.pgram(as.vector(monthly_ts), detrend=TRUE, spans=c(2,7))
abline(v=(1/12)*(1:4), col="blue",lty=4)         
abline(v=(1/60), col="red", lty=3)         

#Periodogram another try
mvspec(monthly_ts, log="no")
specvalues = mvspec(monthly_ts, 4, log="no")


# non parametric and parametric estimate
plot.ts(monthly_ts)  
mvspec(monthly_ts, spans=c(5,5), plot=TRUE, taper=.1, log="no") # nonparametric spectral estimate                           
spec.ar(monthly_ts, log="no")                     # parametric spectral estimate
arma.spec(ar = c(1,-.9), log="no")       # model spectral density 
#----------------------------------------------------------
#Arima
temp=0
for(i in 0:1)
{
  for(j in 0:3)
  {
    for(k in 0:3)
    {
      a=(arima(monthly_ts, order=c(i,j,k)))
      print(paste0('ARIMA','(',i,',',j,',',k,')','=',a$aic))
      temp=rbind(temp,a$aic)
      sarima(monthly_ts,i,j,k)
    }
  }
  
}
sarima(monthly_ts,0,2,2)

temp=temp[-1]
head(temp)
min(temp)

library(forecast)
accuracy(a)
plot(a,axes=F)
arima(monthly_ts, order=c(0,1,1))

# Sarima
for(i in 0:1)
{
  for(j in 0:3)
  {
    for(k in 0:3)
    {
sarima(monthly_ts, i,j,k, 0,0,0,12)
sarima(monthly_ts, i,j,k, 1,0,0,12)
sarima(monthly_ts, i,j,k, 1,1,0,12)
sarima(monthly_ts, i,j,k, 1,1,1,12)
sarima(monthly_ts, i,j,k, 0,1,1,12)
sarima(monthly_ts, i,j,k, 1,1,1,12)
sarima(monthly_ts, i,j,k, 1,0,1,12)
sarima(monthly_ts, i,j,k, 0,1,0,12)
    }
  }
}
sarima(monthly_ts,0,3,1,0,1,1,12)
sarima(monthly_ts,0,1,4,1,0,1,12)
sarima(monthly_ts,0,1,2,1,0,0,12)
sarima(monthly_ts,0,1,2,1,0,0,12)
sarima(monthly_ts,0,2,3,1,1,1,12)
sarima(monthly_ts,0,2,3,1,0,0,12)
sarima(monthly_ts,0,2,2,1,0,0,12)

sarima(monthly_ts,0,0,2)
#final
sarima(monthly_ts,1,1,4,1,0,1,12)
#Prediction
sarima.for(monthly_ts, 10, 1,1,4,1,0,1,12) 
#Male Female statistics
summary(df)



temp1=df[,c(1,2)]
head(temp1)
x <- xts(temp1$Traffic.Collision,temp1$Date.Occurred)

apply.annually(x, colMeans)
    