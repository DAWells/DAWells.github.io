#Exploring several modelling techniques to explore a variable through time. As a demonstration of dealing with longitudinal data.

packages <- c("dplyr", "mgcv", "ggplot2", "lattice")
lapply(packages, require, character.only=T)

setwd("/Users/davidwells/Dropbox/Statistics/weight_time/data")

##########################
#                        #
#   Data preprocessing   #
#                        #
##########################

data<-read.csv("rawdata/mongoose WEIGHTS.csv")
lhdata<-read.csv("rawdata/lhdata_march_2019.csv")

#Remove the stupid names.
lhdata <- lhdata[!lhdata$indiv %in% c("UM PUP", "UM PUP1", "UM PUP2", "UM PUP3", "UM AD", "UM F", "UM SUB", "UMSUB", "UNK", "", "BP375", "BP372"),]

#Format dates.
lhdata$date <- as.Date(lhdata$date, format="%Y-%m-%d")
data$date <- as.Date(data$date, format="%d/%m/%Y")

#Sort by dates
data<-data[order(data$date),]

#Get birth days.
born <- filter(lhdata, stend == "START" & code == "BORN")

#Check duplicates.
born[duplicated(born$indiv),]

#Add age to data.
data$dob <- born$date[match(data$indiv, born$indiv)]
data$age <- data$date - data$dob
data$agen <- as.numeric(data$age)

#Get month of measurment for a seasonal effect.
data$month <- format(data$date, "%b")
data$monthn <- as.numeric(format(data$date, "%m"))

#Remove NAs.
data <- filter(data, !is.na(weight) & !is.na(age))

#To treat data as time series there can only be one record for an individual at a point in time. In this data the time point is the day so we remove any within day/individual duplications.

#Remove outright outliers, e.g. 10kg or pre-birth
data<-filter(data,age>0 & weight<4000)

#Data Single Record.
iddate <- paste(data$indiv,data$date)
dsr <- data[!duplicated(iddate),]

##########################
#                        #
#    Data exploration    #
#                        #
##########################


plot(data$age, data$weight)
xyplot(train$weight ~train$agen/365, xlab="Age in years", ylab="Weight in g")
many_records <- names(summary(data$indiv)[1:3])
mr <- filter(data, indiv %in% many_records)
ggplot(mr, aes(x=age, y=weight, colour=indiv))+geom_line()


!!!!!!!!!
#Remove outliers that are inconsistant with other measurements for that individual.
x<-filter(train,indiv %in% c("FM105","FM142"))
ggplot(x, aes(x=age/365, y=weight, colour=indiv)) + 
geom_point() +
xlab("Age in years") + ylab("Weight in g") + 
theme(panel.background=element_rect(fill="white",colour=NA),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black",size=0.5,lineend="square"),axis.ticks.length=unit(0.2,"cm"),axis.text=element_text(size=11),axis.title=element_text(size=12))+scale_colour_manual(values=c('#cc0000','#3281f7'))+ scale_fill_manual(values=c('#cc0000','#3281f7'))+guides(colour=F)




##########################
#                        #
#       Fit models       #
#                        #
##########################
weight ~ s(age)*id + seasonal_effect + id + E

#Split data sets
train<-dsr[dsr$date < as.Date("2010-01-01"),]

#Trouble with using all of the training data
d<-dsr[dsr$date < as.Date("2010-01-01") & dsr$date >= as.Date("2007-01-01"),]
train<-d[sample(nrow(d),2000),]

val<-dsr[dsr$date >= as.Date("2010-01-01") & dsr$date <= as.Date("2013-01-01"),]
test<-dsr[dsr$date > as.Date("2013-01-01"),]

#smoother for age to capture non-monotonic trend
#This relationship may differ between individuals. Therefore fit separate smoothers by id, or a random slope for id.
#Seasonal_effect based on month could be done using: factors, smoother, or a sin().
#There is obviously temporal autocorrelation in weights so this must be included in the residual variance. Fit an appropriate Covariance structure based on age alone, separate one for each individual, or the same covariance for each individual but they are only expected to covary with measurements from themself.

m6 <- gamm(weight~s(agen, k=20) + s(monthn, bs="cc"), random=list(indiv=~1), correlation = corAR1(form = ~agen|indiv), data=train[seq(1,nrow(train), length.out=2000),])

gam.check(m6$gam)
Eg<-resid(m6$gam, type="d")
plot(m6$gam$linear.predictors, Eg)
plot(Eg~dsr$age[1:1000])
plot(Eg~factor(dsr$month[1:1000]))

#gam.check shows heterogeneity, increasing variance with age.
#Therefore we allow increased residual variance with age

m7 <- gamm(weight~s(agen, k=20) + s(monthn, bs="cc"), random=list(indiv=~1), correlation = corAR1(form = ~agen|indiv), weights = varFixed(~agen), data=train)

gam.check(m7$gam)


#Comparing correlation structures

m8 <- gamm(weight~s(agen, k=20) + s(monthn, bs="cc"), random=list(indiv=~1), correlation = corARMA(p=0,q=1,form = ~agen|indiv), weights = varFixed(~agen), data=dsr[1:1000,])

m9 <- gamm(weight~s(agen, k=20) + s(monthn, bs="cc"), random=list(indiv=~1), correlation = corARMA(p=1,q=1,form = ~agen|indiv), weights = varFixed(~agen), data=dsr[1:1000,])

m10 <- gamm(weight~s(agen, k=20) + s(monthn, bs="cc"), random=list(indiv=~1), correlation = corARMA(p=1,q=2,form = ~agen|indiv), weights = varFixed(~agen), data=dsr[1:1000,])

m11 <- gamm(weight~s(agen, k=20) + s(monthn, bs="cc"), random=list(indiv=~1), correlation = corARMA(p=2,q=1,form = ~agen|indiv), weights = varFixed(~agen), data=dsr[1:1000,])

m12 <- gamm(weight~s(agen, k=20) + s(monthn, bs="cc"), random=list(indiv=~1), correlation = corARMA(p=2,q=2,form = ~agen|indiv), weights = varFixed(~agen), data=dsr[1:1000,])


plot(AIC(m7$lme,m8$lme,m9$lme,m10$lme,m11$lme))


#Visualise correlation matrix
times <- 20
H <- abs(outer(1:times, 1:times, "-"))
image(phi1 ** H)

zm<-matrix(0,nrow=times, ncol=times)
image(rbind(cbind(phi1**H,zm),cbind(zm,phi1**H)))

zm<-matrix(0,nrow=104, ncol=104)
sm<-0
times<-14
H <- abs(outer((sm+1):(sm+times), (sm+1):(sm+times), "-"))
zm[(sm+1):(sm+times),(sm+1):(sm+times)]<-phi1**H
sm<-sm+times

times<-30
H <- abs(outer((sm+1):(sm+times), (sm+1):(sm+times), "-"))
zm[(sm+1):(sm+times),(sm+1):(sm+times)]<-phi1**H
sm<-sm+times

times<-60
H <- abs(outer((sm+1):(sm+times), (sm+1):(sm+times), "-"))
zm[(sm+1):(sm+times),(sm+1):(sm+times)]<-phi1**H
sm<-sm+times

brw<-c("#000000","#440000", "#880000","#bb0000", "#dd0000", "#ffcccc", "#ffffff")
colfunc<-colorRampPalette(c("black","blue", "white"))
image(zm, col=colfunc(20), axes=F)

jpeg("correlation_matrix.jpeg")
dev.off()

#Create AR-1 covariance matrix
times <- 1:5
rho <- 0.5
sigma <- 2
###############
H <- abs(outer(times, times, "-"))
V <- sigma * rho^H
p <- nrow(V)
V[cbind(1:p, 1:p)] <- V[cbind(1:p, 1:p)] * sigma
V

################
mf <- gamm(weight~s(agen, k=20) + factor(month), random=list(indiv=~1), correlation = corAR1(form = ~agen|indiv), weights = varFixed(~agen), data=train)

m2d <- gamm(weight~s(agen, monthn), random=list(indiv=~1), correlation = corAR1(form = ~agen|indiv), weights = varFixed(~agen), data=train)

plot(m7$gam, scale=0)
plot(m2d$gam, scheme=1)
mods <- list(m7$gam,mf$gam,m2d$gam)

lapply(mods, function(mod) mean((predict(mod)-dsr$weight[1:1000])**2))


AIC(m7$lme, mf$lme, m2d$lme)

E<-resid(m7$gam)

ggplot(data=NULL,aes(x=dsr$age[1:1000], y=dsr$monthn[1:1000], size= abs(E), colour=E>0))+geom_point()

#Can AIC compare different smoothers?

##########################
#                        #
#      Predictions       #
#                        #
##########################

pred <- predict(m7$gam, type="response", newdata=val)
E <- val$weight - pred

#Get the random effect intercept
ids <- paste("1/1/", val$indiv, sep="")
id_intercept <- random.effects(m7$lme)$indiv[ids,]
unknown_id <- is.na(id_intercept)
id_intercept[unknown_id]<-0

Er <- val$weight - (pred+id_intercept)

#Predictions ignoring random effects.
xyplot(val$weight + pred ~ val$agen)
#Predictions including random intercepts.
xyplot(val$weight + pred +(pred+id_intercept)~ val$agen)

#Compare residuals with/without random intercept where known.
xyplot(E[!unknown_id] + Er[!unknown_id] ~ val$agen[!unknown_id]/365, xlab="Age in years", ylab="Residuals")

predf <- predict(mf$gam, type="response", newdata=val)
Ef <- val$weight - predf

pred2d <- predict(m2d$gam, type="response", newdata=val)
E2d <- val$weight - pred2d


val_i <- cut_interval(val$agen,n=4)

bwplot(E~val_i|val$month)

Edf<-data.frame(E=c(E,Ef,E2d), type=rep(c("E","Ef","E2d"), each=nrow(val)), val_i = rep(val_i, times=3))

ggplot(Edf, aes(x=val_i, y= E, colour=type))+geom_boxplot()

summarise(group_by(Edf, type),mean_sq = mean(E*E))

##########################
#Do predictions in gam utilise previous information correctly?
DM047

td<- filter(train, indiv=="DM047")
vd<- filter(val, indiv=="DM047")

p<-predict(m7$gam, type="response", newdata=vd)

plot(td$agen, td$weight)
xyplot(p + vd$weight ~ vd$agen)


days<-c(filter(train, indiv=="DM047")$agen, filter(val, indiv=="DM047")$agen)

plot(
ggplot(data=train
xyplot(pt + p ~ days)
xyplot(p ~ days)
xyplot(pt ~ days)
xyplot(filter(train, indiv=="DM047")$weight ~ days)

xyplot((1:5 ) + (4:0) ~1:5)

#Using autocorrelation in predictions
#The AR-1 correlation structure means the residual at a point is phi^timedif + n where n is just noise. Therefore we can use previous residuals to estimate the new ones.


##########################
#                        #
#       Structure        #
#                        #
##########################
#1. weight~age plot
#Non-monotonic relationship, heterogenity, auto-correlation

#2. Cyclical monthly variation plot
#There are both long and short term effects of time.

#3. Heterogeneity
#Fit gamm with random intercept, two smoothers, and AR-1 within indiv. Plot to show increasing variance with age.

#4. varFixed to account for increasing variance.

???

#Try other auto-correlation structures. They do not improve model fit by AIC

#Compare factor(month), s(month, bs=cc), and s(age, month)

#plot AR1 strucutre, including varFixed

#Predict test set.
#This prediction can be made using individual level random effects, should the AR-1 error structre also be incorporated? It could be used to show uncertainty growing over time.

##########################
#                        #
#         To do          #
#                        #
##########################

#Remove outliers that are inconsistent with other measurements for that individual.

#Seasonal effect could be done using sin(), factor(), or a smoother.

#Note seasonal smoother should loop back to the start.

#gam.check to ensure k is high enough

#Evaluate correlation structures

#Why is individual variance so low some times?

#Why is the month smoother 0 for some large datas?

#An issue in the model fitting is that when the data is ordered we can only get a data from a small time frame because there are so many records.