cor(project_data$volact,project_data$involact
    )
project_data$volact
lm.r<-lm(project_data$volact~project_data$race+project_data$fire+project_data$theft+project_data$age+project_data$involact+project_data$income,data=project_data)
summary(lm.r)
plot(lm.r)
project_data$volact
residlm.r<-resid(lm.r)
residlm.r
fittedlm.r<-fitted(lm.r)
fittedlm.r
# observed vs fitted plot<- our model is good
plot(fittedlm.r,project_data$volact,ylab="fitted",xlab="observed")
abline(a=0,b=1,col="purple")
plot(project_data$volact)
lines(fittedlm.r,col="green")
summary(lm.r)
?lm
?summary(lm)
RSE<-1.928
df<-40
RSS<-(RSE^2)*df
RSS
#acf and pacf plot<- there is no indiaction of dependence among the residuals
acf(residlm.r,lag.max=NULL,type="correlation",plot=TRUE) 
pacf(residlm.r,lag.max=NULL,plot=TRUE)
#plotting square of the residuals vs fitted values<- indiactes that variance function is an increasing function of the mean
residsqlm.r<-1:47
for(i in 1:47){residsqlm.r[i]<-residlm.r[i]*residlm.r[i]}
residsqlm.r
plot(fittedlm.r,residsqlm.r,ylab="residual sqaured",xlab="fitted")
#installing rmarkdown
install.packages("rmarkdown")
