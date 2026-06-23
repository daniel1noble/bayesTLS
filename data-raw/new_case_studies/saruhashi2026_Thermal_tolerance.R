## Survival model 28-06-2024
# Stefane Saruhashi
#Empty list ----
if(TRUE){
  rm(list = ls())
  
  #get and set work directory ----
  getwd()
  setwd("C:/Users/Stefane/OneDrive - The University of Western Ontario/PhD/Zebrafish/R folder")
  
  #load packages ----
  library(dplyr)
  library(lme4)
  library(multcomp)
  library(ggplot2)
  library(sjPlot)
  library(sjmisc)
  library(sjlabelled)
  library(boot)
  library(AICcmodavg)
  library(emmeans)
  library(plotly)
}

#load and restructure data ----
if(TRUE){
  #data<- read.csv("C:/Users/Stefane/OneDrive - The University of Western Ontario/PhD/Zebrafish/R folder/Complete dataset_Zb.csv", sep=";") 
  data <- read.csv("C:/Users/Stefane/OneDrive - Cornell University/Side projects/Zebrafish/R folder/Survival_UTL_clean_SS.csv")
  if(FALSE){
    #Show data
    str(data)
  }
  
  #Set as factor
  data$Cohort <- as.factor(data$Cohort)
  data$Ploidy <-as.factor(data$Ploidy)
  data$Treatment<-as.factor(data$Treatment)
  
  #rescale data
  data$temp_sc <- scale(data$T) #add column temp_sc to data
  data$temperature_sc <- scale(data$temperature) #add column temp_sc to data
  data$time_sc <- scale(log10(data$time)) #take log10 before rescaling, otherwise NaS produced! (log of negative not possible)
  data$oxygen_sc <- scale(data$oxygen) #add column temp_sc to data
}

#Subsetting 38C and making graphs ----


if(TRUE){
  #subset
  data_38C<- subset(data, data$T==38) 
  
  #bind survival data for graph of 38C
  Y38<-cbind(data_38C$survival,data_38C$total-data_38C$survival) # number of alive and dead larvae
  
  #Make a graph for all data points of 38C
  #normoxia,3n, t 38, grey
  #plot(jitter(log10(data_38C$time[which(data_38C$Treatment=="normoxia"& data_38C$Ploidy=="3")]),a=.025),Y38[which(data_38C$Treatment=="normoxia"&data_38C$Ploidy=="3")]/rowSums(Y38[which(data_38C$Treatment=="normoxia"&data_38C$Ploidy=="3"),]), xlab = expression("Time (log"[10]*" [minutes])"), pch= 2, ylab = "Survival", col="grey",ylim = c(0,1),xlim=c(0,3.5)) 
  
  
  plot(
    jitter(log10(data_38C$time[which(data_38C$Treatment == "normoxia" & data_38C$Ploidy == "3")]), a = .025),
    Y38[which(data_38C$Treatment == "normoxia" & data_38C$Ploidy == "3")] / rowSums(Y38[which(data_38C$Treatment == "normoxia" & data_38C$Ploidy == "3"),]),
    xlab = expression("Time (log"[10]*" [minutes])"),
    ylab = "Survival",
    pch = 2,
    col = "grey",
    ylim = c(0, 1),
    xlim = c(0, 3.5),
    ylwd = 2 # Adjust the width of the plot line
  )
  
  #normoxia,2n, t 38,grey
  points(jitter(log10(data_38C$time[which(data_38C$Treatment=="normoxia"&data_38C$Ploidy=="2")]),a=.025),Y38[which(data_38C$Treatment=="normoxia"& data_38C$Ploidy=="2")]/rowSums(Y38[which(data_38C$Treatment=="normoxia"&data_38C$Ploidy=="2"),]), xlab = expression("Time (log"[10]*" [minutes])"), pch= 21, ylab = "Survival", col="grey",ylim = c(0,1),xlim=c(0,3.5)) 
  
  #hypoxia,3n, t 38, blue
  points(jitter(log10(data_38C$time[which(data_38C$Treatment=="hypoxia"& data_38C$Ploidy=="3" )]),a=.025),Y38[which(data_38C$Treatment=="hypoxia"& data_38C$Ploidy=="3")]/rowSums(Y38[which(data_38C$Treatment=="hypoxia"& data_38C$Ploidy=="3"),]), xlab = expression("Time (log"[10]*" [minutes])"), pch=2,ylab = "Survival", col="purple",ylim = c(0,1),xlim=c(0,3.5)) 
  
  #hypoxia,2n, t 38, blue
  points(jitter(log10(data_38C$time[which(data_38C$Treatment=="hypoxia"& data_38C$Ploidy=="2" )]),a=.025),Y38[which(data_38C$Treatment=="hypoxia"& data_38C$Ploidy=="2")]/rowSums(Y38[which(data_38C$Treatment=="hypoxia"& data_38C$Ploidy=="2"),]), xlab = expression("Time (log"[10]*" [minutes])"), pch=21,ylab = "Survival", col="purple",ylim = c(0,1),xlim=c(0,3.5)) 
  
  #hyperoxia,2n, t 38, red
  points(jitter(log10(data_38C$time[which(data_38C$Treatment=="hyperoxia"& data_38C$Ploidy=="2" )]),a=.025),Y38[which(data_38C$Treatment=="hyperoxia"& data_38C$Ploidy=="2")]/rowSums(Y38[which(data_38C$Treatment=="hyperoxia"& data_38C$Ploidy=="2"),]), xlab = expression("Time (log"[10]*" [minutes])"), pch=21,ylab = "Survival", col="red",ylim = c(0,1),xlim=c(0,3.5)) 
  
  #hyperoxia,3n, t 38, red
  points(jitter(log10(data_38C$time[which(data_38C$Treatment=="hyperoxia"& data_38C$Ploidy=="3" )]),a=.025),Y38[which(data_38C$Treatment=="hyperoxia"& data_38C$Ploidy=="3")]/rowSums(Y38[which(data_38C$Treatment=="hyperoxia"& data_38C$Ploidy=="3"),]), xlab = expression("Time (log"[10]*" [minutes])"), pch=2,ylab = "Survival", col="red",ylim = c(0,1),xlim=c(0,3.5)) 
}

#make models for survival ----
Y2.a<-cbind(data$survival,data$total-data$survival) # number of alive and dead larvae of all temperatures

#models not used for making graphs ----
if(FALSE){
  m0<-glm(Y2.a~log10(time),family=binomial,data=data) #nul model 
  summary(m0)
  AIC(m0)
  Perc<-(deviance(m0)-deviance(m0))/deviance(m0) # percentage of explained deviace relative to null model
  m.0a<-glm(Y2.a~log10(time)+Ploidy,family=binomial,data=data) #log10(oxygen) is because of the non-lineair relationship between oxygen and temp 
  summary(m.0a)
  AIC(m.0a)
  Perc<-(deviance(m0)-deviance(m.0a))/deviance(m0) # percentage of explained deviace relative to null model
  m.0b<-glm(Y2.a~log10(time)+temperature,family=binomial,data=data) #log10(oxygen) is because of the non-lineair relationship between oxygen and temp 
  summary(m.0b)
  AIC(m.0b)
  Perc<-(deviance(m0)-deviance(m.0b))/deviance(m0) # percentage of explained deviace relative to null model
  m.0c<-glm(Y2.a~log10(time)+log10(oxygen),family=binomial,data=data) #log10(oxygen) is because of the non-lineair relationship between oxygen and temp 
  summary(m.0c)
  AIC(m.0c)
  Perc<-(deviance(m0)-deviance(m.0c))/deviance(m0) # percentage of explained deviace relative to null model
  
  
  m.1a<-glm(Y2.a~log10(time)+Ploidy*temperature,family=binomial,data=data) #log10(oxygen) is because of the non-lineair relationship between oxygen and temp 
  summary(m.1a)
  AIC(m.1a)
  Perc<-(deviance(m0)-deviance(m.1a))/deviance(m0) # percentage of explained deviace relative to null model
  m.1b<-glm(Y2.a~log10(time)+Ploidy*log10(oxygen),family=binomial,data=data) #log10(oxygen) is because of the non-lineair relationship between oxygen and temp 
  summary(m.1b)
  AIC(m.1b)
  Perc<-(deviance(m0)-deviance(m.1b))/deviance(m0) # percentage of explained deviace relative to null model
  m.1c<-glm(Y2.a~log10(time)+temperature*log10(oxygen),family=binomial,data=data) #log10(oxygen) is because of the non-lineair relationship between oxygen and temp 
  summary(m.1c)
  AIC(m.1c)
  Perc<-(deviance(m0)-deviance(m.1c))/deviance(m0) # percentage of explained deviace relative to null model
}

#model of interest
if(TRUE){
  m.2a<-glm(Y2.a~log10(time)+Ploidy*temperature*log10(oxygen),family=binomial,data=data) #log10(oxygen) is because of the non-lineair relationship between oxygen and temp 
  summary(m.2a)
  AIC(m.2a)
  #Perc<-(deviance(m0)-deviance(m.2a))/deviance(m0) # percentage of explained deviace relative to null model
}

#more models not used for making graphs ----
if(FALSE){
  m.2a1<-glmer(Y2.a~log10(time)+Ploidy*temperature*log10(oxygen)+(1|Cohort),family=binomial,data=data) #log10(oxygen) is because of the non-lineair relationship between oxygen and temp 
  summary(m.2a1)
  AIC(m.2a1)
  Perc<-(deviance(m0)-deviance(m.2a1))/deviance(m0) # percentage of explained deviace relative to null model
  m.3<-glm(Y2.a~log10(time)+temperature*Ploidy+log10(oxygen),family=binomial,data=data) #log10(oxygen) is because of the non-lineair relationship between oxygen and temp 
  summary(m.3)
  AIC(m.3)
  Perc<-(deviance(m0)-deviance(m.3))/deviance(m0) # percentage of explained deviace relative to null model
  m.4<-glm(Y2.a~log10(time)+log10(oxygen)*Ploidy+temperature,family=binomial,data=data) #log10(oxygen) is because of the non-lineair relationship between oxygen and temp 
  summary(m.4)
  AIC(m.4)
  Perc<-(deviance(m0)-deviance(m.4))/deviance(m0) # percentage of explained deviace relative to null model
  m.5<-glm(Y2.a~log10(time)+temperature*log10(oxygen)+Ploidy,family=binomial,data=data) #log10(oxygen) is because of the non-lineair relationship between oxygen and temp 
  summary(m.5)
  AIC(m.5)
  Perc<-(deviance(m0)-deviance(m.5))/deviance(m0) # percentage of explained deviace relative to null model
}
#####################################################################################################################
#Making prediction lines on bases of the data points and show in the graph ----
#calculating survival times and making a new dataframe

#plotting lines for 3n, 38C ----
if(TRUE){
  #plot hypoxia,3n, 38C, blue ----
  x<-seq(from=0,to=3.5,length=101) 
  
  linMod<-coefficients(m.2a)["(Intercept)"]+
    
    #linMod<-coefficients(m.2a)["(Intercept)"]+   
    coef(m.2a)["log10(time)"]*x +
    
    coef(m.2a)["Ploidy3"]+ 
    
    coef(m.2a)["temperature"]*38+ 
    
    coef(m.2a)["log10(oxygen)"]*log10(25)+ 
    
    coef(m.2a)["Ploidy3:temperature"]*38+ 
    
    coef(m.2a)["Ploidy3:log10(oxygen)"]*log10(25)+ 
    
    coef(m.2a)["temperature:log10(oxygen)"]*38*log10(25)+ 
    
    coef(m.2a)["Ploidy3:temperature:log10(oxygen)"]*38*log10(25) 
  
  yPredicted<-exp(linMod)/(1+exp(linMod)) 
  
  lines(x,yPredicted,col="purple",lwd= "3", lty="dashed")
  
  #plot 95% confidence interval without random effects (m.2a)
  new <- data.frame(time=10^(seq(from=0.0001,to=3.5,length=1001)),Ploidy=as.factor("3"),temperature=38,oxygen=25)
  res.pred <- predict(m.2a, type = "response", se.fit = TRUE,re.form=NA, newdata = new) 
  # predictions on link scale
  link.pred <- predict(m.2a, type = "link", se.fit = TRUE, re.form=NA, newdata = new) 
  inv <- family(m.2a)$linkinv # inverse of the cauchit link function
  # fit
  new$fit <- res.pred$fit # same as inv(link.pred$fit)
  # 99% Wald CI
  new$u1 <- res.pred$fit + res.pred$se.fit * qnorm(0.995)
  new$l1 <- res.pred$fit - res.pred$se.fit * qnorm(0.995)
  # 99% Wilson CI
  new$u2 <- inv(link.pred$fit + link.pred$se.fit * qnorm(0.995)) 
  new$l2 <- inv(link.pred$fit - link.pred$se.fit * qnorm(0.995))
  
  #lines(log10(new$time),new$u1, col="purple", lwd="1", lty="dashed") #plot upper 95% confidence interval
  #lines(log10(new$time),new$l1, col="purple", lwd="1", lty="dashed") #plot lower 95% confidence interval
  
  # Shade the area between the lines with a semi-transparent color
  polygon(c(log10(new$time), rev(log10(new$time))), c(new$u1, rev(new$l1)), col = adjustcolor("purple", alpha.f = 0.05), border = NA)
  
  
  #Calculate the length (time) of the heat exposure for a specific survival expectancy and save the data
  survival=0.9
  Avg<-approx(yPredicted, x, survival)$y
  up95<-approx(new$u1, log10(new$time), survival)$y
  low95<-approx(new$l1, log10(new$time), survival)$y
  group <- "hypoxia,3n, t 38"
  
  #Export data (absolute)
  Survival <- c(survival) #include survival percentage
  Average38 <- cbind(Avg,up95,low95,Survival,group)
  
  #Plotting normoxia,3n, 38C, grey ----
  x<-seq(from=0.0,to=3.5,length=101) 
  
  linMod<-coefficients(m.2a)["(Intercept)"]+ 
    
    coef(m.2a)["log10(time)"]*x+ 
    
    coef(m.2a)["Ploidy3"]+ 
    
    coef(m.2a)["temperature"]*38+ 
    
    coef(m.2a)["log10(oxygen)"]*log10(100)+ 
    
    coef(m.2a)["Ploidy3:temperature"]*38+ 
    
    coef(m.2a)["Ploidy3:log10(oxygen)"]*log10(100)+ 
    
    coef(m.2a)["temperature:log10(oxygen)"]*38*log10(100)+ 
    
    coef(m.2a)["Ploidy3:temperature:log10(oxygen)"]*38*log10(100) 
  
  yPredicted<-exp(linMod)/(1+exp(linMod)) 
  
  lines(x,yPredicted,col="grey",lwd= "3", lty="dashed") 
  
  #plot 95% confidence interval for normoxia, triploids, 38C
  new <- data.frame(time=10^(seq(from=0.0001,to=3.5,length=1001)),Ploidy=as.factor("3"),temperature=38,oxygen=100)
  res.pred <- predict(m.2a, type = "response", se.fit = TRUE, newdata = new) 
  # predictions on link scale
  link.pred <- predict(m.2a, type = "link", se.fit = TRUE, newdata = new) 
  inv <- family(m.2a)$linkinv # inverse of the cauchit link function
  # fit
  new$fit <- res.pred$fit # same as inv(link.pred$fit)
  # 99% Wald CI
  new$u1 <- res.pred$fit + res.pred$se.fit * qnorm(0.995)
  new$l1 <- res.pred$fit - res.pred$se.fit * qnorm(0.995)
  # 99% Wilson CI
  new$u2 <- inv(link.pred$fit + link.pred$se.fit * qnorm(0.995)) 
  new$l2 <- inv(link.pred$fit - link.pred$se.fit * qnorm(0.995))
  
  #lines(log10(new$time),new$u1, col="grey", lwd="1", lty="dashed") #plot upper 95% confidence interval
  #lines(log10(new$time),new$l1, col="grey", lwd="1", lty="dashed") #plot lower 95% confidence interval
  
  # Shade the area between the lines with a semi-transparent color
  polygon(c(log10(new$time), rev(log10(new$time))), c(new$u1, rev(new$l1)), col = adjustcolor("grey", alpha.f = 0.1), border = NA)
  
  #Calculate the length of the heat exposure for a specific survival expectancy, 90% survival and save the data
  Avg<-approx(yPredicted, x, survival)$y
  up95<-approx(new$u1, log10(new$time), survival)$y
  low95<-approx(new$l1, log10(new$time), survival)$y
  group <- "normoxia,3n, t 38"
  
  #Export data (absolute)
  Average38 <- rbind(Average38, list(Avg,up95,low95,Survival,group))
  
  #plot Hyperoxia,3n, 38C, red ----
  x<-seq(from=0.0,to=3.5,length=101) 
  
  linMod<-coefficients(m.2a)["(Intercept)"]+ 
    
    coef(m.2a)["log10(time)"]*x+ 
    
    coef(m.2a)["Ploidy3"]+ 
    
    coef(m.2a)["temperature"]*38+ 
    
    coef(m.2a)["log10(oxygen)"]*log10(225)+ 
    
    coef(m.2a)["Ploidy3:temperature"]*38+ 
    
    coef(m.2a)["Ploidy3:log10(oxygen)"]*log10(225)+ 
    
    coef(m.2a)["temperature:log10(oxygen)"]*38*log10(225)+ 
    
    coef(m.2a)["Ploidy3:temperature:log10(oxygen)"]*38*log10(225)                        
  
  yPredicted<-exp(linMod)/(1+exp(linMod)) 
  lines(x,yPredicted,col="red",lwd= "3", lty="dashed")
  
  #plot 95% confidence interval, tripoid, hyperoxia, 38C
  new <- data.frame(time=10^(seq(from=0.0001,to=3.5,length=1001)),Ploidy=as.factor("3"),temperature=38,oxygen=225)
  #new <- data.frame(time=seq(from=0,to=240,length=1001),Ploidy=as.factor("2"),temperature=38,oxygen=mean(data$oxygen))
  res.pred <- predict(m.2a, type = "response", se.fit = TRUE, newdata = new) 
  # predictions on link scale
  link.pred <- predict(m.2a, type = "link", se.fit = TRUE, newdata = new) 
  inv <- family(m.2a)$linkinv # inverse of the cauchit link function
  # fit
  new$fit <- res.pred$fit # same as inv(link.pred$fit)
  # 99% Wald CI
  new$u1 <- res.pred$fit + res.pred$se.fit * qnorm(0.995)
  new$l1 <- res.pred$fit - res.pred$se.fit * qnorm(0.995)
  # 99% Wilson CI
  new$u2 <- inv(link.pred$fit + link.pred$se.fit * qnorm(0.995)) 
  new$l2 <- inv(link.pred$fit - link.pred$se.fit * qnorm(0.995))
  
  #lines(log10(new$time),new$u1, col="red", lwd="1", lty="dashed") #plot upper 95% confidence interval
  #lines(log10(new$time),new$l1, col="red", lwd="1", lty="dashed") #plot lower 95% confidence interval
  
  # Shade the area between the lines with a semi-transparent color
  polygon(c(log10(new$time), rev(log10(new$time))), c(new$u1, rev(new$l1)), col = adjustcolor("red", alpha.f = 0.05), border = NA)
  
  #Calculate the length of the heat exposure for a specific survival expectancy, example 90% survival and save the data
  Avg<-approx(yPredicted, x, survival)$y
  up95<-approx(new$u1, log10(new$time), survival)$y
  low95<-approx(new$l1, log10(new$time), survival)$y
  group <- "Hyperoxia,3n, 38C"
  
  #Export data (absolute)
  Average38 <- rbind(Average38, list(Avg,up95,low95,Survival,group))
}
################################################################################################
#Plotting lines for 2n, 38C ----
if(TRUE){
  #plot hypoxia,2n, 38C, blue ----
  x<-seq(from=0,to=3.5,length=101) 
  
  linMod<-coefficients(m.2a)["(Intercept)"]+ 
    
    coef(m.2a)["log10(time)"]*x+ 
    
    #coef(m.2a)["Ploidy3"]+ 
    
    coef(m.2a)["temperature"]*38+ 
    
    coef(m.2a)["log10(oxygen)"]*log10(25)+ 
    
    #coef(m.2a)["Ploidy3:temperature"]*38+ 
    
    #coef(m.2a)["Ploidy3:log10(oxygen)"]*log10(25)+ 
    
    coef(m.2a)["temperature:log10(oxygen)"]*38*log10(25) 
  
  #coef(m.2a)["Ploidy3:temperature:log10(oxygen)"]*38*log10(25) 
  
  yPredicted<-exp(linMod)/(1+exp(linMod)) 
  
  lines(x,yPredicted,col="purple",lwd= "3") 
  
  #plot 95% confidence interval, 2n ,hypoxia, 38C
  new <- data.frame(time=10^(seq(from=0.0001,to=3.5,length=1001)),Ploidy=as.factor("2"),temperature=38,oxygen=25)
  res.pred <- predict(m.2a, type = "response", se.fit = TRUE, newdata = new) 
  # predictions on link scale
  link.pred <- predict(m.2a, type = "link", se.fit = TRUE, newdata = new) 
  inv <- family(m.2a)$linkinv # inverse of the cauchit link function
  # fit
  new$fit <- res.pred$fit # same as inv(link.pred$fit)
  # 99% Wald CI
  new$u1 <- res.pred$fit + res.pred$se.fit * qnorm(0.995)
  new$l1 <- res.pred$fit - res.pred$se.fit * qnorm(0.995)
  # 99% Wilson CI
  new$u2 <- inv(link.pred$fit + link.pred$se.fit * qnorm(0.995)) 
  new$l2 <- inv(link.pred$fit - link.pred$se.fit * qnorm(0.995))
  
  #lines(log10(new$time),new$u1, col="purple", lwd="1", lty="dashed") #plot upper 95% confidence interval
  #lines(log10(new$time),new$l1, col="purple", lwd="1", lty="dashed") #plot lower 95% confidence interval
  
  # Shade the area between the lines with a semi-transparent color
  polygon(c(log10(new$time), rev(log10(new$time))), c(new$u1, rev(new$l1)), col = adjustcolor("purple", alpha.f = 0.05), border = NA)
  
  #Calculate the length of the heat exposure for a specific survival expectancy, example 90% survival and save the data
  Avg<-approx(yPredicted, x, survival)$y
  up95<-approx(new$u1, log10(new$time), survival)$y
  low95<-approx(new$l1, log10(new$time), survival)$y
  group <- "hypoxia,2n, 38C"
  
  #Export data (absolute)
  Average38 <- rbind(Average38, list(Avg,up95,low95,Survival,group))
  
  #Plot normoxia,2n, 38C, grey ----
  x<-seq(from=0.0,to=3.5,length=101) 
  
  linMod<-coefficients(m.2a)["(Intercept)"]+ 
    
    coef(m.2a)["log10(time)"]*x+ 
    
    #coef(m.2a)["Ploidy3"]+ 
    
    coef(m.2a)["temperature"]*38+ 
    
    coef(m.2a)["log10(oxygen)"]*log10(100)+ 
    
    #coef(m.2a)["Ploidy3:temperature"]*38+ 
    
    #coef(m.2a)["Ploidy3:log10(oxygen)"]*log10(100)+ 
    
    coef(m.2a)["temperature:log10(oxygen)"]*38*log10(100)
  
  #coef(m.2a)["Ploidy3:temperature:log10(oxygen)"]*38*log10(100) 
  
  yPredicted<-exp(linMod)/(1+exp(linMod)) 
  
  lines(x,yPredicted,col="grey",lwd= "3") 
  
  #plot 95% confidence interval, 2n, normoxia, 38C
  new <- data.frame(time=10^(seq(from=0.0001,to=3.5,length=1001)),Ploidy=as.factor("2"),temperature=38,oxygen=100)
  res.pred <- predict(m.2a, type = "response", se.fit = TRUE, newdata = new) 
  # predictions on link scale
  link.pred <- predict(m.2a, type = "link", se.fit = TRUE, newdata = new) 
  inv <- family(m.2a)$linkinv # inverse of the cauchit link function
  # fit
  new$fit <- res.pred$fit # same as inv(link.pred$fit)
  # 99% Wald CI
  new$u1 <- res.pred$fit + res.pred$se.fit * qnorm(0.995)
  new$l1 <- res.pred$fit - res.pred$se.fit * qnorm(0.995)
  # 99% Wilson CI
  new$u2 <- inv(link.pred$fit + link.pred$se.fit * qnorm(0.995)) 
  new$l2 <- inv(link.pred$fit - link.pred$se.fit * qnorm(0.995))
  
  #lines(log10(new$time),new$u1, col="grey", lwd="1", lty="dashed") #plot upper 95% confidence interval
  #lines(log10(new$time),new$l1, col="grey", lwd="1", lty="dashed") #plot lower 95% confidence interval
  
  # Shade the area between the lines with a semi-transparent color
  polygon(c(log10(new$time), rev(log10(new$time))), c(new$u1, rev(new$l1)), col = adjustcolor("grey", alpha.f = 0.1), border = NA)
  
  #Calculate the length of the heat exposure for a specific survival expectancy, example 90% survival and save the data
  Avg<-approx(yPredicted, x, survival)$y
  up95<-approx(new$u1, log10(new$time), survival)$y
  low95<-approx(new$l1, log10(new$time), survival)$y
  group <- "normoxia,2n, 38C"
  
  #Export data (absolute)
  Average38 <- rbind(Average38, list(Avg,up95,low95,Survival,group))
  
  #plot Hyperoxia,2n, 38C, red ----
  x<-seq(from=0.0,to=3.5,length=101) 
  
  linMod<-coefficients(m.2a)["(Intercept)"]+ 
    
    coef(m.2a)["log10(time)"]*x+ 
    
    #coef(m.2a)["Ploidy3"]+ 
    
    coef(m.2a)["temperature"]*38+ 
    
    coef(m.2a)["log10(oxygen)"]*log10(225)+ 
    
    #coef(m.2a)["Ploidy3:temperature"]*38+ 
    
    #coef(m.2a)["Ploidy3:log10(oxygen)"]*log10(225)+ 
    
    coef(m.2a)["temperature:log10(oxygen)"]*38*log10(225) 
  
  #coef(m.2a)["Ploidy3:temperature:log10(oxygen)"]*38*log10(225)                        
  
  yPredicted<-exp(linMod)/(1+exp(linMod)) 
  lines(x,yPredicted,col="red",lwd= "3")
  
  #plot 95% confidence interval, 2n, hyperoxia, 38C
  new <- data.frame(time=10^(seq(from=0.0001,to=3.5,length=1001)),Ploidy=as.factor("2"),temperature=38,oxygen=225)
  res.pred <- predict(m.2a, type = "response", se.fit = TRUE, newdata = new) 
  # predictions on link scale
  link.pred <- predict(m.2a, type = "link", se.fit = TRUE, newdata = new) 
  inv <- family(m.2a)$linkinv # inverse of the cauchit link function
  # fit
  new$fit <- res.pred$fit # same as inv(link.pred$fit)
  # 99% Wald CI
  new$u1 <- res.pred$fit + res.pred$se.fit * qnorm(0.995)
  new$l1 <- res.pred$fit - res.pred$se.fit * qnorm(0.995)
  # 99% Wilson CI
  new$u2 <- inv(link.pred$fit + link.pred$se.fit * qnorm(0.995)) 
  new$l2 <- inv(link.pred$fit - link.pred$se.fit * qnorm(0.995))
  
  #lines(log10(new$time),new$u1, col="red", lwd="1", lty="dashed") #plot upper 95% confidence interval
  #lines(log10(new$time),new$l1, col="red", lwd="1", lty="dashed") #plot lower 95% confidence interval
  
  # Shade the area between the lines with a semi-transparent color
  polygon(c(log10(new$time), rev(log10(new$time))), c(new$u1, rev(new$l1)), col = adjustcolor("red", alpha.f = 0.05), border = NA)
  
  #Calculate the length of the heat exposure for a specific survival expectancy, example 90% survivaland save the data
  Avg<-approx(yPredicted, x, survival)$y
  up95<-approx(new$u1, log10(new$time), survival)$y
  low95<-approx(new$l1, log10(new$time), survival)$y
  group <- "Hyperoxia,2n, 38C"
  
  #Export data (absolute)
  Average38 <- rbind(Average38, list(Avg,up95,low95,Survival,group))
  
  #plot legend and title for the graph ----
  legend("topright", legend = c("Dipoid","Triploid", "Hypoxia", "Normoxia", "Hyperoxia"), pch = c(21,2,NA,NA,NA), lty= c("solid", "dashed","solid","solid","solid"), lwd = c(1,2,2,2,2),cex=1.0, col = c("black","black","purple","grey","red"))
  title("38°C")
  
  #make Average38 into a dataframe ----
  Average38_df<- data.frame(matrix(unlist(Average38), nrow=6, byrow=F),stringsAsFactors=FALSE)
  colnames(Average38_df) <- c("Avg", "up95", "low95", "Survival", "group")
  
  #change characters to numbers
  Average38_df$Avg <- as.numeric(Average38_df$Avg)
  Average38_df$up95 <- as.numeric(Average38_df$up95)
  Average38_df$low95 <- as.numeric(Average38_df$low95)
  Average38_df$Survival <- as.numeric(Average38_df$Survival)
}

# boxplot data from average 38C + 95% confidence interval ----
#ggplot(Average38_df, aes(x = group, y = Avg, fill = group)) +
#  geom_bar(stat = "identity", position = "dodge", color = "black") +
#  geom_errorbar(aes(ymin = low95, ymax = up95, color = group), 
#                position = position_dodge(width = 0.9), 
#                width = 0.25) +
#  scale_fill_manual(values = c("red","red", "purple", "purple", "grey", "grey")) +  # Specify fill colors
#  scale_color_manual(values = c("black", "black", "black","black","black","black")) +  # Specify error bar colors
#  labs(title = paste("Average exposure time for",Survival*100,"% survival at 38C (CI)"), x = "Ploidy", y = expression("Average time (log"[10]*"[minutes])")) +
#  scale_x_discrete(labels = c("2n", "3n", "2n", "3n","2n","3n")) +  # Custom x-axis labels
# scale_y_continuous(breaks = seq(0, 2.5, by = 0.5), limits = c(0, 2.5)) +  # Set specific y-axis tick positions
# theme_minimal()

####################################################################################################################################################
####################################################################################################################################################
#subsetting data for 39C and making plots ----
if(TRUE){
  data_39C<- subset(data, data$T==39) 
  Y39<-cbind(data_39C$survival,data_39C$total-data_39C$survival) # number of alive and dead larvae
  
  #normoxia,3n, t 39, groen
  plot(jitter(log10(data_39C$time[which(data_39C$Treatment=="normoxia"& data_39C$Ploidy=="3")]),a=.025),Y39[which(data_39C$Treatment=="normoxia"&data_39C$Ploidy=="3")]/rowSums(Y39[which(data_39C$Treatment=="normoxia"&data_39C$Ploidy=="3"),]), xlab = expression("Time (log"[10]*" [minutes])"), pch= 2, ylab = "Survival", col="grey",ylim = c(0,1),xlim=c(0,3.5)) 
  #normoxia,2n, t 39,groen
  points(jitter(log10(data_39C$time[which(data_39C$Treatment=="normoxia"&data_39C$Ploidy=="2")]),a=.025),Y39[which(data_39C$Treatment=="normoxia"& data_39C$Ploidy=="2")]/rowSums(Y39[which(data_39C$Treatment=="normoxia"&data_39C$Ploidy=="2"),]), xlab = expression("Time (log"[10]*" [minutes])"), pch= 21, ylab = "Survival", col="grey",ylim = c(0,1),xlim=c(0,3.5)) 
  
  #hypoxia,3n, t 39, blue
  points(jitter(log10(data_39C$time[which(data_39C$Treatment=="hypoxia"& data_39C$Ploidy=="3" )]),a=.025),Y39[which(data_39C$Treatment=="hypoxia"& data_39C$Ploidy=="3")]/rowSums(Y39[which(data_39C$Treatment=="hypoxia"& data_39C$Ploidy=="3"),]), xlab = expression("Time (log"[10]*" [minutes])"), pch=2,ylab = "Survival", col="purple",ylim = c(0,1),xlim=c(0,3.5)) 
  
  #hypoxia,2n, t 39, blue
  points(jitter(log10(data_39C$time[which(data_39C$Treatment=="hypoxia"& data_39C$Ploidy=="2" )]),a=.025),Y39[which(data_39C$Treatment=="hypoxia"& data_39C$Ploidy=="2")]/rowSums(Y39[which(data_39C$Treatment=="hypoxia"& data_39C$Ploidy=="2"),]), xlab = expression("Time (log"[10]*" [minutes])"), pch=21,ylab = "Survival", col="purple",ylim = c(0,1),xlim=c(0,3.5)) 
  
  #hyperoxia,2n, t 39, red
  points(jitter(log10(data_39C$time[which(data_39C$Treatment=="hyperoxia"& data_39C$Ploidy=="2" )]),a=.025),Y39[which(data_39C$Treatment=="hyperoxia"& data_39C$Ploidy=="2")]/rowSums(Y39[which(data_39C$Treatment=="hyperoxia"& data_39C$Ploidy=="2"),]), xlab = expression("Time (log"[10]*" [minutes])"), pch=21,ylab = "Survival", col="red",ylim = c(0,1),xlim=c(0,3.5)) 
  
  #hyperoxia,3n, t 39, red
  points(jitter(log10(data_39C$time[which(data_39C$Treatment=="hyperoxia"& data_39C$Ploidy=="3" )]),a=.025),Y39[which(data_39C$Treatment=="hyperoxia"& data_39C$Ploidy=="3")]/rowSums(Y39[which(data_39C$Treatment=="hyperoxia"& data_39C$Ploidy=="3"),]), xlab = expression("Time (log"[10]*" [minutes])"), pch=2,ylab = "Survival", col="red",ylim = c(0,1),xlim=c(0,3.5)) 
  
  #plot for 3n 39C ----
  
  #plot hypoxia,3n, 39C, blue ----
  x<-seq(from=0,to=3.5,length=101) 
  
  linMod<-coefficients(m.2a)["(Intercept)"]+ 
    
    coef(m.2a)["log10(time)"]*x+ 
    
    coef(m.2a)["Ploidy3"]+ 
    
    coef(m.2a)["temperature"]*39+ 
    
    coef(m.2a)["log10(oxygen)"]*log10(25)+ 
    
    coef(m.2a)["Ploidy3:temperature"]*39+ 
    
    coef(m.2a)["Ploidy3:log10(oxygen)"]*log10(25)+ 
    
    coef(m.2a)["temperature:log10(oxygen)"]*39*log10(25)+ 
    
    coef(m.2a)["Ploidy3:temperature:log10(oxygen)"]*39*log10(25) 
  
  yPredicted<-exp(linMod)/(1+exp(linMod)) 
  
  lines(x,yPredicted,col="purple",lwd= "3", lty="dashed") 
  
  #plot 95% confidence interval
  new <- data.frame(time=10^(seq(from=0.001,to=3.5,length=1001)),Ploidy=as.factor("3"),temperature=39,oxygen=25)
  res.pred <- predict(m.2a, type = "response", se.fit = TRUE, newdata = new) 
  # predictions on link scale
  link.pred <- predict(m.2a, type = "link", se.fit = TRUE, newdata = new) 
  inv <- family(m.2a)$linkinv # inverse of the cauchit link function
  # fit
  new$fit <- res.pred$fit # same as inv(link.pred$fit)
  # 99% Wald CI
  new$u1 <- res.pred$fit + res.pred$se.fit * qnorm(0.995)
  new$l1 <- res.pred$fit - res.pred$se.fit * qnorm(0.995)
  # 99% Wilson CI
  new$u2 <- inv(link.pred$fit + link.pred$se.fit * qnorm(0.995)) 
  new$l2 <- inv(link.pred$fit - link.pred$se.fit * qnorm(0.995))
  
  #lines(log10(new$time),new$u1, col="purple", lwd="1", lty="dashed") #plot upper 95% confidence interval
  #lines(log10(new$time),new$l1, col="purple", lwd="1", lty="dashed") #plot lower 95% confidence interval
  
  # Shade the area between the lines with a semi-transparent color
  polygon(c(log10(new$time), rev(log10(new$time))), c(new$u1, rev(new$l1)), col = adjustcolor("purple", alpha.f = 0.05), border = NA)
  
  #Calculate the length of the heat exposure for a specific survival expectancy, example 90% survival
  Avg<-approx(yPredicted, x, survival)$y
  up95<-approx(new$u1, log10(new$time), survival)$y
  low95<-approx(new$l1, log10(new$time), survival)$y
  group <- "hypoxia,3n, 39C"
  
  #Export data (absolute)
  Average39 <- cbind(Avg,up95,low95,Survival,group)
  
  #plot normoxia,3n, 39C, grey ----
  x<-seq(from=0.0,to=3.5,length=101) 
  
  linMod<-coefficients(m.2a)["(Intercept)"]+ 
    
    coef(m.2a)["log10(time)"]*x+ 
    
    coef(m.2a)["Ploidy3"]+ 
    
    coef(m.2a)["temperature"]*39+ 
    
    coef(m.2a)["log10(oxygen)"]*log10(100)+ 
    
    coef(m.2a)["Ploidy3:temperature"]*39+ 
    
    coef(m.2a)["Ploidy3:log10(oxygen)"]*log10(100)+ 
    
    coef(m.2a)["temperature:log10(oxygen)"]*39*log10(100)+ 
    
    coef(m.2a)["Ploidy3:temperature:log10(oxygen)"]*39*log10(100) 
  
  yPredicted<-exp(linMod)/(1+exp(linMod)) 
  
  lines(x,yPredicted,col="grey",lwd= "3", lty="dashed") 
  
  #plot 95% confidence interval
  new <- data.frame(time=10^(seq(from=0.0001,to=3.5,length=1001)),Ploidy=as.factor("3"),temperature=39,oxygen=100)
  res.pred <- predict(m.2a, type = "response", se.fit = TRUE, newdata = new) 
  # predictions on link scale
  link.pred <- predict(m.2a, type = "link", se.fit = TRUE, newdata = new) 
  inv <- family(m.2a)$linkinv # inverse of the cauchit link function
  # fit
  new$fit <- res.pred$fit # same as inv(link.pred$fit)
  # 99% Wald CI
  new$u1 <- res.pred$fit + res.pred$se.fit * qnorm(0.995)
  new$l1 <- res.pred$fit - res.pred$se.fit * qnorm(0.995)
  # 99% Wilson CI
  new$u2 <- inv(link.pred$fit + link.pred$se.fit * qnorm(0.995)) 
  new$l2 <- inv(link.pred$fit - link.pred$se.fit * qnorm(0.995))
  
  #lines(log10(new$time),new$u1, col="grey", lwd="1", lty="dashed") #plot upper 95% confidence interval
  #lines(log10(new$time),new$l1, col="grey", lwd="1", lty="dashed") #plot lower 95% confidence interval
  
  # Shade the area between the lines with a semi-transparent color
  polygon(c(log10(new$time), rev(log10(new$time))), c(new$u1, rev(new$l1)), col = adjustcolor("grey", alpha.f = 0.1), border = NA)
  
  #Calculate the length of the heat exposure for a specific survival expectancy, example 90% survival and save the data
  Avg<-approx(yPredicted, x, survival)$y
  up95<-approx(new$u1, log10(new$time), survival)$y
  low95<-approx(new$l1, log10(new$time), survival)$y
  group <- "normoxia,3n, t 39"
  
  #Export data (absolute)
  Average39 <- rbind(Average39, list(Avg,up95,low95,Survival,group))
  
  #Hyperoxia,3n, 39C, red ----
  x<-seq(from=0.0,to=3.5,length=101) 
  
  linMod<-coefficients(m.2a)["(Intercept)"]+ 
    
    coef(m.2a)["log10(time)"]*x+ 
    
    coef(m.2a)["Ploidy3"]+ 
    
    coef(m.2a)["temperature"]*39+ 
    
    coef(m.2a)["log10(oxygen)"]*log10(225)+ 
    
    coef(m.2a)["Ploidy3:temperature"]*39+ 
    
    coef(m.2a)["Ploidy3:log10(oxygen)"]*log10(225)+ 
    
    coef(m.2a)["temperature:log10(oxygen)"]*39*log10(225)+ 
    
    coef(m.2a)["Ploidy3:temperature:log10(oxygen)"]*39*log10(225)                        
  
  yPredicted<-exp(linMod)/(1+exp(linMod)) 
  lines(x,yPredicted,col="red",lwd= "3", lty="dashed")
  
  #plot 95% confidence interval
  new <- data.frame(time=10^(seq(from=0.0001,to=3.5,length=1001)),Ploidy=as.factor("3"),temperature=39,oxygen=225)
  res.pred <- predict(m.2a, type = "response", se.fit = TRUE, newdata = new) 
  # predictions on link scale
  link.pred <- predict(m.2a, type = "link", se.fit = TRUE, newdata = new) 
  inv <- family(m.2a)$linkinv # inverse of the cauchit link function
  # fit
  new$fit <- res.pred$fit # same as inv(link.pred$fit)
  # 99% Wald CI
  new$u1 <- res.pred$fit + res.pred$se.fit * qnorm(0.995)
  new$l1 <- res.pred$fit - res.pred$se.fit * qnorm(0.995)
  # 99% Wilson CI
  new$u2 <- inv(link.pred$fit + link.pred$se.fit * qnorm(0.995)) 
  new$l2 <- inv(link.pred$fit - link.pred$se.fit * qnorm(0.995))
  
  #lines(log10(new$time),new$u1, col="red", lwd="1", lty="dashed") #plot upper 95% confidence interval
  #lines(log10(new$time),new$l1, col="red", lwd="1", lty="dashed") #plot lower 95% confidence interval
  
  # Shade the area between the lines with a semi-transparent color
  polygon(c(log10(new$time), rev(log10(new$time))), c(new$u1, rev(new$l1)), col = adjustcolor("red", alpha.f = 0.05), border = NA)
  
  #Calculate the length of the heat exposure for a specific survival expectancy, example 90% survival and save the data
  Avg<-approx(yPredicted, x, survival)$y
  up95<-approx(new$u1, log10(new$time), survival)$y
  low95<-approx(new$l1, log10(new$time), survival)$y
  group <- "hyperoxia,3n, t 39"
  
  #Export data (absolute)
  Average39 <- rbind(Average39, list(Avg,up95,low95,Survival,group))
  ################################################################################################
  #plot 2n, 39C ----
  #plot hypoxia,2n, 39C, blue ----
  x<-seq(from=0,to=3.5,length=101) 
  
  linMod<-coefficients(m.2a)["(Intercept)"]+ 
    
    coef(m.2a)["log10(time)"]*x+ 
    
    #coef(m.2a)["Ploidy3"]+ 
    
    coef(m.2a)["temperature"]*39+ 
    
    coef(m.2a)["log10(oxygen)"]*log10(25)+ 
    
    #coef(m.2a)["Ploidy3:temperature"]*39+ 
    
    #coef(m.2a)["Ploidy3:log10(oxygen)"]*log10(25)+ 
    
    coef(m.2a)["temperature:log10(oxygen)"]*39*log10(25) 
  
  #coef(m.2a)["Ploidy3:temperature:log10(oxygen)"]*39*log10(25) 
  
  yPredicted<-exp(linMod)/(1+exp(linMod)) 
  
  lines(x,yPredicted,col="purple",lwd= "3") 
  
  #plot 95% confidence interval
  new <- data.frame(time=10^(seq(from=0.0001,to=3.5,length=1001)),Ploidy=as.factor("2"),temperature=39,oxygen=25)
  res.pred <- predict(m.2a, type = "response", se.fit = TRUE, newdata = new) 
  # predictions on link scale
  link.pred <- predict(m.2a, type = "link", se.fit = TRUE, newdata = new) 
  inv <- family(m.2a)$linkinv # inverse of the cauchit link function
  # fit
  new$fit <- res.pred$fit # same as inv(link.pred$fit)
  # 99% Wald CI
  new$u1 <- res.pred$fit + res.pred$se.fit * qnorm(0.995)
  new$l1 <- res.pred$fit - res.pred$se.fit * qnorm(0.995)
  # 99% Wilson CI
  new$u2 <- inv(link.pred$fit + link.pred$se.fit * qnorm(0.995)) 
  new$l2 <- inv(link.pred$fit - link.pred$se.fit * qnorm(0.995))
  
  #lines(log10(new$time),new$u1, col="purple", lwd="1", lty="dashed") #plot upper 95% confidence interval
  #lines(log10(new$time),new$l1, col="purple", lwd="1", lty="dashed") #plot lower 95% confidence interval
  
  # Shade the area between the lines with a semi-transparent color
  polygon(c(log10(new$time), rev(log10(new$time))), c(new$u1, rev(new$l1)), col = adjustcolor("purple", alpha.f = 0.05), border = NA)
  
  #Calculate the length of the heat exposure for a specific survival expectancy, example 90% survival and save the data
  Avg<-approx(yPredicted, x, survival)$y
  up95<-approx(new$u1, log10(new$time), survival)$y
  low95<-approx(new$l1, log10(new$time), survival)$y
  group <- "hypoxia,2n, t 39"
  
  #Export data (absolute)
  Average39 <- rbind(Average39, list(Avg,up95,low95,Survival,group))
  
  #plot normoxia,2n, 39C, grey ----
  x<-seq(from=0.0,to=3.5,length=101) 
  
  linMod<-coefficients(m.2a)["(Intercept)"]+ 
    
    coef(m.2a)["log10(time)"]*x+ 
    
    #coef(m.2a)["Ploidy3"]+ 
    
    coef(m.2a)["temperature"]*39+ 
    
    coef(m.2a)["log10(oxygen)"]*log10(100)+ 
    
    #coef(m.2a)["Ploidy3:temperature"]*39+ 
    
    #coef(m.2a)["Ploidy3:log10(oxygen)"]*log10(100)+ 
    
    coef(m.2a)["temperature:log10(oxygen)"]*39*log10(100)
  
  #coef(m.2a)["Ploidy3:temperature:log10(oxygen)"]*39*log10(100) 
  
  yPredicted<-exp(linMod)/(1+exp(linMod)) 
  
  lines(x,yPredicted,col="grey",lwd= "3") 
  
  #plot 95% confidence interval
  new <- data.frame(time=10^(seq(from=0.0001,to=3.5,length=1001)),Ploidy=as.factor("2"),temperature=39,oxygen=100)
  res.pred <- predict(m.2a, type = "response", se.fit = TRUE, newdata = new) 
  # predictions on link scale
  link.pred <- predict(m.2a, type = "link", se.fit = TRUE, newdata = new) 
  inv <- family(m.2a)$linkinv # inverse of the cauchit link function
  # fit
  new$fit <- res.pred$fit # same as inv(link.pred$fit)
  # 99% Wald CI
  new$u1 <- res.pred$fit + res.pred$se.fit * qnorm(0.995)
  new$l1 <- res.pred$fit - res.pred$se.fit * qnorm(0.995)
  # 99% Wilson CI
  new$u2 <- inv(link.pred$fit + link.pred$se.fit * qnorm(0.995)) 
  new$l2 <- inv(link.pred$fit - link.pred$se.fit * qnorm(0.995))
  
  #lines(log10(new$time),new$u1, col="grey", lwd="1", lty="dashed") #plot upper 95% confidence interval
  #lines(log10(new$time),new$l1, col="grey", lwd="1", lty="dashed") #plot lower 95% confidence interval
  
  # Shade the area between the lines with a semi-transparent color
  polygon(c(log10(new$time), rev(log10(new$time))), c(new$u1, rev(new$l1)), col = adjustcolor("grey", alpha.f = 0.1), border = NA)
  
  #Calculate the length of the heat exposure for a specific survival expectancy, example 90% survival and save the data
  Avg<-approx(yPredicted, x, survival)$y
  up95<-approx(new$u1, log10(new$time), survival)$y
  low95<-approx(new$l1, log10(new$time), survival)$y
  group <- "normoxia,2n, 39C"
  
  #Export data (absolute)
  Average39 <- rbind(Average39, list(Avg,up95,low95,Survival,group))
  
  #plot Hyperoxia,2n, 39C, red ----
  x<-seq(from=0.0,to=3.5,length=101) 
  
  linMod<-coefficients(m.2a)["(Intercept)"]+ 
    
    coef(m.2a)["log10(time)"]*x+ 
    
    #coef(m.2a)["Ploidy3"]+ 
    
    coef(m.2a)["temperature"]*39+ 
    
    coef(m.2a)["log10(oxygen)"]*log10(225)+ 
    
    #coef(m.2a)["Ploidy3:temperature"]*39+ 
    
    #coef(m.2a)["Ploidy3:log10(oxygen)"]*log10(225)+ 
    
    coef(m.2a)["temperature:log10(oxygen)"]*39*log10(225) 
  
  #coef(m.2a)["Ploidy3:temperature:log10(oxygen)"]*39*log10(225)                        
  
  yPredicted<-exp(linMod)/(1+exp(linMod)) 
  lines(x,yPredicted,col="red",lwd= "3")
  
  #plot 95% confidence interval
  new <- data.frame(time=10^(seq(from=0.0001,to=3.5,length=1001)),Ploidy=as.factor("2"),temperature=39,oxygen=225)
  res.pred <- predict(m.2a, type = "response", se.fit = TRUE, newdata = new) 
  # predictions on link scale
  link.pred <- predict(m.2a, type = "link", se.fit = TRUE, newdata = new) 
  inv <- family(m.2a)$linkinv # inverse of the cauchit link function
  # fit
  new$fit <- res.pred$fit # same as inv(link.pred$fit)
  # 99% Wald CI
  new$u1 <- res.pred$fit + res.pred$se.fit * qnorm(0.995)
  new$l1 <- res.pred$fit - res.pred$se.fit * qnorm(0.995)
  # 99% Wilson CI
  new$u2 <- inv(link.pred$fit + link.pred$se.fit * qnorm(0.995)) 
  new$l2 <- inv(link.pred$fit - link.pred$se.fit * qnorm(0.995))
  
  #lines(log10(new$time),new$u1, col="red", lwd="1", lty="dashed") #plot upper 95% confidence interval
  #lines(log10(new$time),new$l1, col="red", lwd="1", lty="dashed") #plot lower 95% confidence interval
  
  # Shade the area between the lines with a semi-transparent color
  polygon(c(log10(new$time), rev(log10(new$time))), c(new$u1, rev(new$l1)), col = adjustcolor("red", alpha.f = 0.05), border = NA)
  
  #Calculate the length of the heat exposure for a specific survival expectancy, example 90% survival and save the data
  Avg<-approx(yPredicted, x, survival)$y
  up95<-approx(new$u1, log10(new$time), survival)$y
  low95<-approx(new$l1, log10(new$time), survival)$y
  group <- "hyperoxia,2n, t 39"
  
  #Export data (absolute) ----
  Average39 <- rbind(Average39, list(Avg,up95,low95,Survival,group))
  Average39a<- data.frame(matrix(unlist(Average39), nrow=6, byrow=F),stringsAsFactors=FALSE)
  Average39_df<- data.frame(matrix(unlist(Average39), nrow=6, byrow=F),stringsAsFactors=FALSE)
  colnames(Average39_df) <- c("Avg", "up95", "low95", "Survival", "group")
  
  #plot legend and title for the graph ----
  legend("topright", legend = c("Dipoid","Triploid", "Hypoxia", "Normoxia", "Hyperoxia"), pch = c(21,2,NA,NA,NA), lty= c("solid", "dashed","solid","solid","solid"), lwd = c(1,2,2,2,2),cex=1, col = c("black","black","purple","grey","red"))
  title("39°C")
  
  #change characters to numbers
  Average39_df$Avg <- as.numeric(Average39_df$Avg)
  Average39_df$up95 <- as.numeric(Average39_df$up95)
  Average39_df$low95 <- as.numeric(Average39_df$low95)
  Average39_df$Survival <- as.numeric(Average39_df$Survival)
}

# boxplot data from average 39C + 95% confidence interval ----
#ggplot(Average39_df, aes(x = group, y = Avg, fill = group)) +
#  geom_bar(stat = "identity", position = "dodge", color = "black") +
#  geom_errorbar(aes(ymin = low95, ymax = up95, color = group), 
 #               position = position_dodge(width = 0.9), 
 #               width = 0.25) +
 # scale_fill_manual(values = c("red","red", "purple", "purple", "grey", "grey")) +  # Specify fill colors
#  scale_color_manual(values = c("black", "black", "black","black","black","black")) +  # Specify error bar colors
#  labs(title = paste("Average exposure time for",Survival*100,"% survival at 39C (CI)"), x = "Ploidy", y = expression("Average time (log"[10]*"[minutes])")) +
#  scale_x_discrete(labels = c("2n", "3n", "2n", "3n","2n","3n")) +  # Custom x-axis labels
#  scale_y_continuous(breaks = seq(0, 2.5, by = 0.5), limits = c(0, 2.5)) +  # Set specific y-axis tick positions
#  theme_minimal()

#####################################################################################################################################################################
#####################################################################################################################################################################
#subsetting data for 40C and making plots ----
if(TRUE){
  data_40C<- subset(data, data$T==40) 
  Y40<-cbind(data_40C$survival,data_40C$total-data_40C$survival) # number of alive and dead larvae
  
  #normoxia,3n, t 40, grey
  plot(jitter(log10(data_40C$time[which(data_40C$Treatment=="normoxia"& data_40C$Ploidy=="3")]),a=.025),Y40[which(data_40C$Treatment=="normoxia"&data_40C$Ploidy=="3")]/rowSums(Y40[which(data_40C$Treatment=="normoxia"&data_40C$Ploidy=="3"),]), xlab = expression("Time (log"[10]*" [minutes])"), pch= 2, ylab = "Survival", col="grey",ylim = c(0,1),xlim=c(0,3.5)) 
  #normoxia,2n, t 40,grey
  points(jitter(log10(data_40C$time[which(data_40C$Treatment=="normoxia"&data_40C$Ploidy=="2")]),a=.025),Y40[which(data_40C$Treatment=="normoxia"& data_40C$Ploidy=="2")]/rowSums(Y40[which(data_40C$Treatment=="normoxia"&data_40C$Ploidy=="2"),]), xlab = expression("Time (log"[10]*" [minutes])"), pch= 21, ylab = "Survival", col="grey",ylim = c(0,1),xlim=c(0,3.5)) 
  
  #hypoxia,3n, t 40, blue
  points(jitter(log10(data_40C$time[which(data_40C$Treatment=="hypoxia"& data_40C$Ploidy=="3" )]),a=.025),Y40[which(data_40C$Treatment=="hypoxia"& data_40C$Ploidy=="3")]/rowSums(Y40[which(data_40C$Treatment=="hypoxia"& data_40C$Ploidy=="3"),]), xlab = expression("Time (log"[10]*" [minutes])"), pch=2,ylab = "Survival", col="purple",ylim = c(0,1),xlim=c(0,3.5)) 
  
  #hypoxia,2n, t 40, blue
  points(jitter(log10(data_40C$time[which(data_40C$Treatment=="hypoxia"& data_40C$Ploidy=="2" )]),a=.025),Y40[which(data_40C$Treatment=="hypoxia"& data_40C$Ploidy=="2")]/rowSums(Y40[which(data_40C$Treatment=="hypoxia"& data_40C$Ploidy=="2"),]), xlab = expression("Time (log"[10]*" [minutes])"), pch=21,ylab = "Survival", col="purple",ylim = c(0,1),xlim=c(0,3.5)) 
  
  #hyperoxia,2n, t 40, red
  points(jitter(log10(data_40C$time[which(data_40C$Treatment=="hyperoxia"& data_40C$Ploidy=="2" )]),a=.025),Y40[which(data_40C$Treatment=="hyperoxia"& data_40C$Ploidy=="2")]/rowSums(Y40[which(data_40C$Treatment=="hyperoxia"& data_40C$Ploidy=="2"),]), xlab = expression("Time (log"[10]*" [minutes])"), pch=21,ylab = "Survival", col="red",ylim = c(0,1),xlim=c(0,3.5)) 
  
  #hyperoxia,3n, t 40, red
  points(jitter(log10(data_40C$time[which(data_40C$Treatment=="hyperoxia"& data_40C$Ploidy=="3" )]),a=.025),Y40[which(data_40C$Treatment=="hyperoxia"& data_40C$Ploidy=="3")]/rowSums(Y40[which(data_40C$Treatment=="hyperoxia"& data_40C$Ploidy=="3"),]), xlab = expression("Time (log"[10]*" [minutes])"), pch=2,ylab = "Survival", col="red",ylim = c(0,1),xlim=c(0,3.5)) 
  
  #plot 3n, 40C ----
  #plot hypoxia,3n, 40C, blue ----
  x<-seq(from=0,to=3.5,length=101) 
  
  linMod<-coefficients(m.2a)["(Intercept)"]+ 
    
    coef(m.2a)["log10(time)"]*x+ 
    
    coef(m.2a)["Ploidy3"]+ 
    
    coef(m.2a)["temperature"]*40+ 
    
    coef(m.2a)["log10(oxygen)"]*log10(25)+ 
    
    coef(m.2a)["Ploidy3:temperature"]*40+ 
    
    coef(m.2a)["Ploidy3:log10(oxygen)"]*log10(25)+ 
    
    coef(m.2a)["temperature:log10(oxygen)"]*40*log10(25)+ 
    
    coef(m.2a)["Ploidy3:temperature:log10(oxygen)"]*40*log10(25) 
  
  yPredicted<-exp(linMod)/(1+exp(linMod)) 
  
  lines(x,yPredicted,col="purple",lwd= "3", lty="dashed") 
  
  #plot 95% confidence interval
  new <- data.frame(time=10^(seq(from=0.0001,to=3.5,length=1001)),Ploidy=as.factor("3"),temperature=40,oxygen=25)
  res.pred <- predict(m.2a, type = "response", se.fit = TRUE, newdata = new) 
  # predictions on link scale
  link.pred <- predict(m.2a, type = "link", se.fit = TRUE, newdata = new) 
  inv <- family(m.2a)$linkinv # inverse of the cauchit link function
  # fit
  new$fit <- res.pred$fit # same as inv(link.pred$fit)
  # 99% Wald CI
  new$u1 <- res.pred$fit + res.pred$se.fit * qnorm(0.995)
  new$l1 <- res.pred$fit - res.pred$se.fit * qnorm(0.995)
  # 99% Wilson CI
  new$u2 <- inv(link.pred$fit + link.pred$se.fit * qnorm(0.995)) 
  new$l2 <- inv(link.pred$fit - link.pred$se.fit * qnorm(0.995))
  
  #lines(log10(new$time),new$u1, col="purple", lwd="1", lty="dashed") #plot upper 95% confidence interval
  #lines(log10(new$time),new$l1, col="purple", lwd="1", lty="dashed") #plot lower 95% confidence interval
  
  # Shade the area between the lines with a semi-transparent color
  polygon(c(log10(new$time), rev(log10(new$time))), c(new$u1, rev(new$l1)), col = adjustcolor("purple", alpha.f = 0.05), border = NA)
  
  #Calculate the length of the heat exposure for a specific survival expectancy, example 90% survival and save the data
  Avg<-approx(yPredicted, x, survival)$y
  up95<-approx(new$u1, log10(new$time), survival)$y
  low95<-approx(new$l1, log10(new$time), survival)$y
  group <- "hypoxia,3n, t 40"
  
  #Export data (absolute)
  Average40 <- cbind(Avg,up95,low95,Survival,group)
  
  #plot normoxia,3n, 40C, grey ----
  x<-seq(from=0.0,to=3.5,length=101) 
  
  linMod<-coefficients(m.2a)["(Intercept)"]+ 
    
    coef(m.2a)["log10(time)"]*x+ 
    
    coef(m.2a)["Ploidy3"]+ 
    
    coef(m.2a)["temperature"]*40+ 
    
    coef(m.2a)["log10(oxygen)"]*log10(100)+ 
    
    coef(m.2a)["Ploidy3:temperature"]*40+ 
    
    coef(m.2a)["Ploidy3:log10(oxygen)"]*log10(100)+ 
    
    coef(m.2a)["temperature:log10(oxygen)"]*40*log10(100)+ 
    
    coef(m.2a)["Ploidy3:temperature:log10(oxygen)"]*40*log10(100) 
  
  yPredicted<-exp(linMod)/(1+exp(linMod)) 
  
  lines(x,yPredicted,col="grey",lwd= "3", lty="dashed") 
  
  #plot 95% confidence interval
  new <- data.frame(time=10^(seq(from=0.0001,to=3.5,length=1001)),Ploidy=as.factor("3"),temperature=40,oxygen=100)
  res.pred <- predict(m.2a, type = "response", se.fit = TRUE, newdata = new) 
  # predictions on link scale
  link.pred <- predict(m.2a, type = "link", se.fit = TRUE, newdata = new) 
  inv <- family(m.2a)$linkinv # inverse of the cauchit link function
  # fit
  new$fit <- res.pred$fit # same as inv(link.pred$fit)
  # 99% Wald CI
  new$u1 <- res.pred$fit + res.pred$se.fit * qnorm(0.995)
  new$l1 <- res.pred$fit - res.pred$se.fit * qnorm(0.995)
  # 99% Wilson CI
  new$u2 <- inv(link.pred$fit + link.pred$se.fit * qnorm(0.995)) 
  new$l2 <- inv(link.pred$fit - link.pred$se.fit * qnorm(0.995))
  
  #lines(log10(new$time),new$u1, col="grey", lwd="1", lty="dashed") #plot upper 95% confidence interval
  #lines(log10(new$time),new$l1, col="grey", lwd="1", lty="dashed") #plot lower 95% confidence interval
  
  # Shade the area between the lines with a semi-transparent color
  polygon(c(log10(new$time), rev(log10(new$time))), c(new$u1, rev(new$l1)), col = adjustcolor("grey", alpha.f = 0.1), border = NA)
  
  #Calculate the length of the heat exposure for a specific survival expectancy, example 90% survival
  Avg<-approx(yPredicted, x, survival)$y
  up95<-approx(new$u1, log10(new$time), survival)$y
  low95<-approx(new$l1, log10(new$time), survival)$y
  group <- "normoxia,3n, 40C"
  
  #Export data (absolute)
  Average40 <- rbind(Average40, list(Avg,up95,low95,Survival,group))
  
  #plot Hyperoxia,3n, 40C, red ----
  x<-seq(from=0.0,to=3.5,length=101) 
  
  linMod<-coefficients(m.2a)["(Intercept)"]+ 
    
    coef(m.2a)["log10(time)"]*x+ 
    
    coef(m.2a)["Ploidy3"]+ 
    
    coef(m.2a)["temperature"]*40+ 
    
    coef(m.2a)["log10(oxygen)"]*log10(225)+ 
    
    coef(m.2a)["Ploidy3:temperature"]*40+ 
    
    coef(m.2a)["Ploidy3:log10(oxygen)"]*log10(225)+ 
    
    coef(m.2a)["temperature:log10(oxygen)"]*40*log10(225)+ 
    
    coef(m.2a)["Ploidy3:temperature:log10(oxygen)"]*40*log10(225)                        
  
  yPredicted<-exp(linMod)/(1+exp(linMod)) 
  lines(x,yPredicted,col="red",lwd= "3", lty="dashed")
  
  #plot 95% confidence interval
  new <- data.frame(time=10^(seq(from=0.0001,to=3.5,length=1001)),Ploidy=as.factor("3"),temperature=40,oxygen=225)
  res.pred <- predict(m.2a, type = "response", se.fit = TRUE, newdata = new) 
  # predictions on link scale
  link.pred <- predict(m.2a, type = "link", se.fit = TRUE, newdata = new) 
  inv <- family(m.2a)$linkinv # inverse of the cauchit link function
  # fit
  new$fit <- res.pred$fit # same as inv(link.pred$fit)
  # 99% Wald CI
  new$u1 <- res.pred$fit + res.pred$se.fit * qnorm(0.995)
  new$l1 <- res.pred$fit - res.pred$se.fit * qnorm(0.995)
  # 99% Wilson CI
  new$u2 <- inv(link.pred$fit + link.pred$se.fit * qnorm(0.995)) 
  new$l2 <- inv(link.pred$fit - link.pred$se.fit * qnorm(0.995))
  
  #lines(log10(new$time),new$u1, col="red", lwd="1", lty="dashed") #plot upper 95% confidence interval
  #lines(log10(new$time),new$l1, col="red", lwd="1", lty="dashed") #plot lower 95% confidence interval
  
  # Shade the area between the lines with a semi-transparent color
  polygon(c(log10(new$time), rev(log10(new$time))), c(new$u1, rev(new$l1)), col = adjustcolor("red", alpha.f = 0.05), border = NA)
  
  #Calculate the length of the heat exposure for a specific survival expectancy, example 90% survival and save the data
  Avg<-approx(yPredicted, x, survival)$y
  up95<-approx(new$u1, log10(new$time), survival)$y
  low95<-approx(new$l1, log10(new$time), survival)$y
  group <- "hyperoxia,3n, t 40"
  
  #Export data (absolute)
  Average40 <- rbind(Average40, list(Avg,up95,low95,Survival,group))
  
  ################################################################################################
  #plot for 2n 40C ----
  #plot hypoxia,2n, 40C, blue ----
  x<-seq(from=0,to=3.5,length=101) 
  
  linMod<-coefficients(m.2a)["(Intercept)"]+ 
    
    coef(m.2a)["log10(time)"]*x+ 
    
    #coef(m.2a)["Ploidy3"]+ 
    
    coef(m.2a)["temperature"]*40+ 
    
    coef(m.2a)["log10(oxygen)"]*log10(25)+ 
    
    #coef(m.2a)["Ploidy3:temperature"]*40+ 
    
    #coef(m.2a)["Ploidy3:log10(oxygen)"]*log10(25)+ 
    
    coef(m.2a)["temperature:log10(oxygen)"]*40*log10(25) 
  
  #coef(m.2a)["Ploidy3:temperature:log10(oxygen)"]*40*log10(25) 
  
  yPredicted<-exp(linMod)/(1+exp(linMod)) 
  
  lines(x,yPredicted,col="purple",lwd= "3") 
  
  #plot 95% confidence interval
  new <- data.frame(time=10^(seq(from=0.0001,to=3.5,length=1001)),Ploidy=as.factor("2"),temperature=40,oxygen=25)
  res.pred <- predict(m.2a, type = "response", se.fit = TRUE, newdata = new) 
  # predictions on link scale
  link.pred <- predict(m.2a, type = "link", se.fit = TRUE, newdata = new) 
  inv <- family(m.2a)$linkinv # inverse of the cauchit link function
  # fit
  new$fit <- res.pred$fit # same as inv(link.pred$fit)
  # 99% Wald CI
  new$u1 <- res.pred$fit + res.pred$se.fit * qnorm(0.995)
  new$l1 <- res.pred$fit - res.pred$se.fit * qnorm(0.995)
  # 99% Wilson CI
  new$u2 <- inv(link.pred$fit + link.pred$se.fit * qnorm(0.995)) 
  new$l2 <- inv(link.pred$fit - link.pred$se.fit * qnorm(0.995))
  
  #lines(log10(new$time),new$u1, col="purple", lwd="1", lty="dashed") #plot upper 95% confidence interval
  #lines(log10(new$time),new$l1, col="purple", lwd="1", lty="dashed") #plot lower 95% confidence interval
  
  # Shade the area between the lines with a semi-transparent color
  polygon(c(log10(new$time), rev(log10(new$time))), c(new$u1, rev(new$l1)), col = adjustcolor("purple", alpha.f = 0.05), border = NA)
  
  #Calculate the length of the heat exposure for a specific survival expectancy, example 90% survival and save the data
  Avg<-approx(yPredicted, x, survival)$y
  up95<-approx(new$u1, log10(new$time), survival)$y
  low95<-approx(new$l1, log10(new$time), survival)$y
  group <- "hypoxia,2n, t 40"
  
  #Export data (absolute)
  Average40 <- rbind(Average40, list(Avg,up95,low95,Survival,group))
  
  #plot normoxia,2n, 40C, grey ----
  x<-seq(from=0.0,to=3.5,length=101) 
  
  linMod<-coefficients(m.2a)["(Intercept)"]+ 
    
    coef(m.2a)["log10(time)"]*x+ 
    
    #coef(m.2a)["Ploidy3"]+ 
    
    coef(m.2a)["temperature"]*40+ 
    
    coef(m.2a)["log10(oxygen)"]*log10(100)+ 
    
    #coef(m.2a)["Ploidy3:temperature"]*40+ 
    
    #coef(m.2a)["Ploidy3:log10(oxygen)"]*log10(100)+ 
    
    coef(m.2a)["temperature:log10(oxygen)"]*40*log10(100)
  
  #coef(m.2a)["Ploidy3:temperature:log10(oxygen)"]*40*log10(100) 
  
  yPredicted<-exp(linMod)/(1+exp(linMod)) 
  
  lines(x,yPredicted,col="grey",lwd= "3") 
  
  #plot 95% confidence interval
  new <- data.frame(time=10^(seq(from=0.0001,to=3.5,length=1001)),Ploidy=as.factor("2"),temperature=40,oxygen=100)
  res.pred <- predict(m.2a, type = "response", se.fit = TRUE, newdata = new) 
  # predictions on link scale
  link.pred <- predict(m.2a, type = "link", se.fit = TRUE, newdata = new) 
  inv <- family(m.2a)$linkinv # inverse of the cauchit link function
  # fit
  new$fit <- res.pred$fit # same as inv(link.pred$fit)
  # 99% Wald CI
  new$u1 <- res.pred$fit + res.pred$se.fit * qnorm(0.995)
  new$l1 <- res.pred$fit - res.pred$se.fit * qnorm(0.995)
  # 99% Wilson CI
  new$u2 <- inv(link.pred$fit + link.pred$se.fit * qnorm(0.995)) 
  new$l2 <- inv(link.pred$fit - link.pred$se.fit * qnorm(0.995))
  
  #lines(log10(new$time),new$u1, col="grey", lwd="1", lty="dashed") #plot upper 95% confidence interval
  #lines(log10(new$time),new$l1, col="grey", lwd="1", lty="dashed") #plot lower 95% confidence interval
  
  # Shade the area between the lines with a semi-transparent color
  polygon(c(log10(new$time), rev(log10(new$time))), c(new$u1, rev(new$l1)), col = adjustcolor("grey", alpha.f = 0.1), border = NA)
  
  #Calculate the length of the heat exposure for a specific survival expectancy, example 90% survival and save the data
  Avg<-approx(yPredicted, x, survival)$y
  up95<-approx(new$u1, log10(new$time), survival)$y
  low95<-approx(new$l1, log10(new$time), survival)$y
  group <- "normoxia,2n, 40C"
  
  #Export data (absolute)
  Average40 <- rbind(Average40, list(Avg,up95,low95,Survival,group))
  
  #plot Hyperoxia,2n, 40C, red ----
  x<-seq(from=0.0,to=3.5,length=101) 
  
  linMod<-coefficients(m.2a)["(Intercept)"]+ 
    
    coef(m.2a)["log10(time)"]*x+ 
    
    #coef(m.2a)["Ploidy3"]+ 
    
    coef(m.2a)["temperature"]*40+ 
    
    coef(m.2a)["log10(oxygen)"]*log10(225)+ 
    
    #coef(m.2a)["Ploidy3:temperature"]*40+ 
    
    #coef(m.2a)["Ploidy3:log10(oxygen)"]*log10(225)+ 
    
    coef(m.2a)["temperature:log10(oxygen)"]*40*log10(225) 
  
  #coef(m.2a)["Ploidy3:temperature:log10(oxygen)"]*40*log10(225)                        
  
  yPredicted<-exp(linMod)/(1+exp(linMod)) 
  lines(x,yPredicted,col="red",lwd= "3")
  
  #plot 95% confidence interval
  new <- data.frame(time=10^(seq(from=0.0001,to=3.5,length=1001)),Ploidy=as.factor("2"),temperature=40,oxygen=225)
  res.pred <- predict(m.2a, type = "response", se.fit = TRUE, newdata = new) 
  # predictions on link scale
  link.pred <- predict(m.2a, type = "link", se.fit = TRUE, newdata = new) 
  inv <- family(m.2a)$linkinv # inverse of the cauchit link function
  # fit
  new$fit <- res.pred$fit # same as inv(link.pred$fit)
  # 99% Wald CI
  new$u1 <- res.pred$fit + res.pred$se.fit * qnorm(0.995)
  new$l1 <- res.pred$fit - res.pred$se.fit * qnorm(0.995)
  # 99% Wilson CI
  new$u2 <- inv(link.pred$fit + link.pred$se.fit * qnorm(0.995)) 
  new$l2 <- inv(link.pred$fit - link.pred$se.fit * qnorm(0.995))
  
  #lines(log10(new$time),new$u1, col="red", lwd="1", lty="dashed") #plot upper 95% confidence interval
  #lines(log10(new$time),new$l1, col="red", lwd="1", lty="dashed") #plot lower 95% confidence interval
  
  # Shade the area between the lines with a semi-transparent color
  polygon(c(log10(new$time), rev(log10(new$time))), c(new$u1, rev(new$l1)), col = adjustcolor("red", alpha.f = 0.05), border = NA)
  
  #Calculate the length of the heat exposure for a specific survival expectancy, example 90% survival and save the data
  Avg<-approx(yPredicted, x, survival)$y
  up95<-approx(new$u1, log10(new$time), survival)$y
  low95<-approx(new$l1, log10(new$time), survival)$y
  group <- "hyperoxia,2n, 40C"
  
  #Export data (absolute) ----
  Average40 <- rbind(Average40, list(Avg,up95,low95,Survival,group))
  Average40a<- data.frame(matrix(unlist(Average40), nrow=6, byrow=F),stringsAsFactors=FALSE)
  Average40_df<- data.frame(matrix(unlist(Average40), nrow=6, byrow=F),stringsAsFactors=FALSE)
  colnames(Average40_df) <- c("Avg", "up95", "low95", "Survival", "group")
  
  #plot legend and title in the graph ----
  title("40°C")
  legend("topright", legend = c("Dipoid","Triploid", "Hypoxia", "Normoxia", "Hyperoxia"), pch = c(21,2,NA,NA,NA), lty= c("solid", "dashed","solid","solid","solid"), lwd = c(1,2,2,2,2),cex=1, col = c("black","black","purple","grey","red"))
  
  #converting characters into numerics
  Average40_df$Avg <- as.numeric(Average40_df$Avg)
  Average40_df$up95 <- as.numeric(Average40_df$up95)
  Average40_df$low95 <- as.numeric(Average40_df$low95)
  Average40_df$Survival <- as.numeric(Average40_df$Survival)
} 

# boxplot data from average 40C + 95% confidence interval
#ggplot(Average40_df, aes(x = group, y = Avg, fill = group)) +
 # geom_bar(stat = "identity", position = "dodge", color = "black") +
  #geom_errorbar(aes(ymin = low95, ymax = up95, color = group), 
   #             position = position_dodge(width = 0.9), 
    #            width = 0.25) +
#  scale_fill_manual(values = c("red","red", "purple", "purple", "grey", "grey")) +  # Specify fill colors
 # scale_color_manual(values = c("black", "black", "black","black","black","black")) +  # Specify error bar colors
  #labs(title = paste("Average exposure time for",Survival*100,"% survival at 40C (CI)"), x = "Ploidy", y = expression("Average time (log"[10]*"[minutes])")) +
#  scale_x_discrete(labels = c("2n", "3n", "2n", "3n","2n","3n")) +  # Custom x-axis labels
 # scale_y_continuous(breaks = seq(0, 2.5, by = 0.5), limits = c(0, 2.5)) +  # Set specific y-axis tick positions
#theme_minimal()

#############################################################################################################################################################
#############################################################################################################################################################
# 1. Filte diploid
results <- expand.grid(Temperature = c(38, 39, 40), 
                       Ploidy = 2,           # Apenas Ploidy 3
                       Oxygen = c(25, 100, 225))

results$fit <- NA
results$upr <- NA
results$lwr <- NA
survival <- 0.5

# 2. Calculate predictions
for (i in 1:nrow(results)){
  new <- data.frame(time = 10^(seq(from = 0.0001, to = 3.5, length = 1001)),
                    Ploidy = as.factor(results$Ploidy[i]),
                    temperature = results$Temperature[i],
                    oxygen = results$Oxygen[i])
  
  res.pred <- predict(m.2a, type = "response", se.fit = TRUE, re.form = NA, newdata = new)
  
  new$fit <- res.pred$fit 
  new$u1 <- res.pred$fit + res.pred$se.fit * qnorm(0.975)
  new$l1 <- res.pred$fit - res.pred$se.fit * qnorm(0.975)
  
  results$fit[i] <- approx(new$fit, log10(new$time), survival)$y
  results$upr[i] <- approx(new$u1, log10(new$time), survival)$y
  results$lwr[i] <- approx(new$l1, log10(new$time), survival)$y
}

# 3. Define colours
cols_oxy <- c("#9932CC", "#808080", "#FF0000")

# 4. graph
plot(fit ~ Temperature, data = results, ylim = c(0.5, 2.8), col = "white",
     ylab = "Time to 50% mortality (min)", yaxt = 'n', xlab = "Temperature (°C)",
     main = "A. Diploid")

axis(side = 2, at = c(0, log10(5), 1, log10(25), 2, log10(300)), 
     labels = c(1, 5, 10, 25, 100, 300), las = 2)

# Loop 
oxy_levels <- c(25, 100, 225)

for (i in 1:length(oxy_levels)){
  temp <- results[results$Oxygen == oxy_levels[i], ]
  
  
  polygon(c(temp$Temperature, rev(temp$Temperature)), 
          c(temp$upr, rev(temp$lwr)), 
          col = adjustcolor(cols_oxy[i], alpha.f = 0.3), border = NA)
  
  lines(temp$Temperature, temp$fit, col = cols_oxy[i], lty = 1, lwd = 2)
  points(temp$Temperature, temp$fit, col = cols_oxy[i], pch = 16)
}
tabela_final_diploid <- results
tabela_final_diploid$fit_min <- 10^(results$fit)
tabela_final_diploid$upr_min <- 10^(results$upr)
tabela_final_diploid$lwr_min <- 10^(results$lwr)

print(tabela_final_diploid) 

    ### triploid
  # 1. Filter 
  results <- expand.grid(Temperature = c(38, 39, 40), 
                         Ploidy = 3,           # Apenas Ploidy 3
                         Oxygen = c(25, 100, 225))
  
  results$fit <- NA
  results$upr <- NA
  results$lwr <- NA
  survival <- 0.5
  
  # 2. Prediction
  for (i in 1:nrow(results)){
    new <- data.frame(time = 10^(seq(from = 0.0001, to = 3.5, length = 1001)),
                      Ploidy = as.factor(results$Ploidy[i]),
                      temperature = results$Temperature[i],
                      oxygen = results$Oxygen[i])
    
    res.pred <- predict(m.2a, type = "response", se.fit = TRUE, re.form = NA, newdata = new)
    
    new$fit <- res.pred$fit 
    new$u1 <- res.pred$fit + res.pred$se.fit * qnorm(0.975)
    new$l1 <- res.pred$fit - res.pred$se.fit * qnorm(0.975)
    
    results$fit[i] <- approx(new$fit, log10(new$time), survival)$y
    results$upr[i] <- approx(new$u1, log10(new$time), survival)$y
    results$lwr[i] <- approx(new$l1, log10(new$time), survival)$y
  }
  
  # 3. Colours
  cols_oxy <- c("#9932CC", "#808080", "#FF0000")
  
  # 4. Plot
  plot(fit ~ Temperature, data = results, ylim = c(0.5, 2.8), col = "white",
       ylab = "Time to 50% mortality (min)", yaxt = 'n', xlab = "Temperature (°C)",
       main = "B. Triploid")
  
  axis(side = 2, at = c(0, log10(5), 1, log10(25), 2, log10(300)), 
       labels = c(1, 5, 10, 25, 100, 300), las = 2)
  
  # Loop
  oxy_levels <- c(25, 100, 225)
  
  for (i in 1:length(oxy_levels)){
    temp <- results[results$Oxygen == oxy_levels[i], ]
    
    # Desenha o polígono de erro
    polygon(c(temp$Temperature, rev(temp$Temperature)), 
            c(temp$upr, rev(temp$lwr)), 
            col = adjustcolor(cols_oxy[i], alpha.f = 0.3), border = NA)
    
    lines(temp$Temperature, temp$fit, col = cols_oxy[i], lty = 1, lwd = 2)
    points(temp$Temperature, temp$fit, col = cols_oxy[i], pch = 16)
  }

  # Extract table
 
  
  tabela_final_triploid <- results
  tabela_final_triploid$fit_min <- 10^(results$fit)
  tabela_final_triploid$upr_min <- 10^(results$upr)
  tabela_final_triploid$lwr_min <- 10^(results$lwr)
  
  # View
  print(tabela_final_triploid)   
  
  
    
  
