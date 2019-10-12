#importing the dataset
library(readxl)
project_data <- read_excel("E:/study/project/project_data.xlsx")
View(project_data)

#initial model
model1 <- lm(project_data$volact ~ 1 + project_data$race + project_data$fire + project_data$theft +
               project_data$age + project_data$involact + project_data$income, data = project_data)

summary(model1)
plot(model1) #explanation of the plots needed..!!

#remove two outliers blindly
project_data_WO <- project_data[c(-24, -7),] #data without outliers

model2 = lm(project_data_WO$volact ~ 1 + project_data_WO$race + project_data_WO$fire +
              project_data_WO$theft + project_data_WO$age + project_data_WO$involact + 
              project_data_WO$income, data = project_data_WO)
plot(model2)

#residuals vs y_hat graph showed some sort of pattern
b = (model1$residuals^2) / (1 - hatvalues(model1)) 
plot(model1$fitted.values, b, xlab = "Fitted Values") #this gives a better picture for non-constant variance
#this plot matches seberly pg:283 fig (a)


#10/10/19
plot(model2$residuals)

e_star <- model2$residuals + model2$coefficients[2] * project_data_WO[,2]
added_variable <- data.frame(project_data_WO[,2], e_star)
plot(added_variable)

par(mfrow = c(2,3))
for (i in 2:8) {
  if (i == 6)
    next
  e_star <- model2$residuals + model2$coefficients[i] * project_data_WO[,i]
  added_variable <- data.frame(project_data_WO[,i], e_star)
  plot(added_variable)
}

12/10/19
# observed vs fitted plot<- our model is good
plot(project_data$fitted.values,project_data$volact,ylab="fitted",xlab="observed")
abline(a=0,b=1,col="purple")
plot(project_data$volact)
lines(project_data$fitted.values,col="green")
summary(lm.r)
?lm
?summary(lm)

#acf and pacf plot<- there is no indiaction of dependence among the residuals
acf(model1$residuals,lag.max=NULL,type="correlation",plot=TRUE) 
pacf(model1$residuals,lag.max=NULL,plot=TRUE)
#plotting square of the residuals vs fitted values<- indiactes that variance function is an increasing function of the mean
residsqlm.r<-1:47
for(i in 1:47){residsqlm.r[i]<-model1_residuals[i]*model1_residuals[i]}
residsqlm.r
plot(model1_fitted.values,residsqlm.r,ylab="residual sqaured",xlab="fitted")
#installing rmarkdown
install.packages("rmarkdown")
