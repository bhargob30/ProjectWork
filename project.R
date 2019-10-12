cor(project_data$volact,project_data$involact
    )
project_data$volact
lm.r<-lm(project_data$volact~project_data$race+project_data$fire+project_data$theft+project_data$age+project_data$involact+project_data$income,data=project_data)
summary(lm.r)
plot(lm.r)
project_data_WO <- project_data[c(-24, -7),] 
model2 = lm(project_data_WO$volact ~ 1 + project_data_WO$race + project_data_WO$fire +
              project_data_WO$theft + project_data_WO$age + project_data_WO$involact + 
              project_data_WO$income, data = project_data_WO)
plot(model2)
project_data$volact
residlm.r<-resid(lm.r)
residlm.r
fittedlm.r<-fitted(lm.r)
fittedlm.r
plot(fittedlm.r,project_data$volact,ylab="fitted",xlab="observed")
abline(a=0,b=1,col="purple")
plot(project_data$volact)
lines(fittedlm.r,col="green")
RSS(lm.r)
summary(lm.r)
?lm
?summary(lm)
RSE<-1.928
df<-40
RSS<-(RSE^2)*df
RSS
acf(residlm.r,lag.max=NULL,type="correlation",plot=TRUE) 
warnings()
pacf(residlm.r,lag.max=NULL,plot=TRUE)
residsqlm.r<-1:47
for(i in 1:47){residsqlm.r[i]<-residlm.r[i]*residlm.r[i]}
residsqlm.r
plot(fittedlm.r,residsqlm.r,ylab="residual sqaured",xlab="fitted")
install.packages("rmarkdown")
