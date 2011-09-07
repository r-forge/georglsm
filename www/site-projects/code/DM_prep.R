#####################
### Begin 
#####################
rm(list=ls(all=TRUE))

source("func.R")

######################
### read original data 
######################
trans <- read.csv("beijing transition.csv")
names(trans)

#####################
#### data preparation 
#####################

tran <- trans[,2:67]
n <- nrow(tran)

### Q18: DOB converted into age (remove outliers)
age.o <- tran$Q18
age <- as.character(age.o)
age <- strsplit(age,"-")
year <- NULL
for(i in 1:n){
	temp <- age[[i]][1]
	year <- c(year,as.numeric(temp))
}
age <- 2004-year
age[age < 10 | age > 110] <- NA #age size range [10,110]
tran$Q18 <- age

### Q20: family size (remove outliers)
tran$Q20[tran$Q20 > 20] <- NA #family size range [1,20]

### check range for all columns (all in range)
apply(tran,2,function(x)range(x,na.rm=T)) 

### Q4: rating for transportation condition
tc <- tran[,4:11]
tcs <- apply(tc,1,sum)
range(tcs,na.rm=T)#[1]  8 37
sum(!is.na(tcs)) #586 valid obs
tcs <- tcs-min(tcs,na.rm=T)
tcs <- tcs/max(tcs,na.rm=T)
#transform into range [0,1]

### Q5: rating for transportation problem
tp <- tran[,12:23]
tps <- apply(tp,1,sum)
range(tps,na.rm=T)#[1] 12 72
sum(!is.na(tps)) #557 valid obs
tps <- tps-min(tps,na.rm=T)
tps <- tps/max(tps,na.rm=T)
#transform into range [0,1]

ind <- !is.na(tcs) & !is.na(tps)
cor(tcs[ind],7*12-tps[ind])^2
#cor between 2 ratings: [1] 0.3067261

### Combine results from Q4 & Q5
rate.t <- 100*(tcs-tps+1)/2
#overall rating [0,100]
sum(!is.na(rate.t))#557 valid obs
range(rate.t,na.rm=T)#[1]  0.0 92.5

### Q10&11: Driver Licence
dl <- cbind(tran$Q10,tran$Q11)

### Q15: Accident History
accid <- cbind(tran$Q15_1,tran$Q15_2,tran$Q15_3,tran$Q15_4)
accid <- (accid-1)*2.5#happen times {0,2.5,5,7.5}
acc <- accid[,1]*1000+accid[,2]*1000*30+accid[,3]*1000*60++accid[,4]*1000*120
#converted into total accident loss
range(acc,na.rm=T)#[1] 0  975000
levels(as.factor(acc))
#[1] "0"      "2500"   "5000"   "7500"   "75000"  "77500"  "80000"  "377500"
#[9] "975000"

### Q16: Driving Habits
ind.col("16_1")#37
ind.col("16_22")#58
habit <- tran[,37:58]
ind.u <- c(9,15)#2 undeterined
ind.g <- c(2,3,4,7,8,11,14,18,19,20)#10 good habits
ind.b <- c(1,5,6,10,12,13,16,17,21,22)#bad habits
hab.g <- habit[,c(ind.g)]
hab.b <- habit[,c(ind.b)]
score.dh <- apply((hab.g-3),1,sum)+apply((3-hab.b),1,sum)
#score of driving habits
range(score.dh,na.rm=T)#[1] -4 40
score.dh <- score.dh-min(score.dh,na.rm=T)
score.dh <- 100*score.dh/max(score.dh,na.rm=T)
#transform into range [0,100]
length(levels(as.factor(score.dh)))#44 levels

#### END OF PREPROCESSING 
#########################




