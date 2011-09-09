
#####################
### snowfall examples
#####################

library(snowfall)

# 1. Initialisation of snowfall.
# (if used with sfCluster, just call sfInit())
sfInit(parallel=TRUE, cpus=4, type="SOCK")

# 2. Loading data.
require(mvna)
data(sir.adm)

# 3. Wrapper, which can be parallelised.
wrapper <- function(idx) {
# Output progress in worker logfile
	#cat( "Current index: ", idx, "\n" )
index <- sample(1:nrow(sir.adm), replace=TRUE)
temp <- sir.adm[index, ]
fit <- lm(temp$time~ temp$status+ temp$pneu)
return(fit$coef)
}

# 4. Exporting needed data and loading required
# packages on workers.
sfExport("sir.adm")
sfLibrary(cmprsk)

# 5. Start network random number generator
# (as "sample" is using random numbers).
sfClusterSetupRNG()

# 6. Distribute calculation
rtime <- NULL
for(run in c(5000,1e4)){
atime <- proc.time()
result <- sfLapply(1:run, wrapper)
rtime <- rbind(rtime,c(run,(proc.time()-atime)[1:2]))
}
write.table(rtime,"rtime.txt",row.name=FALSE)

#btime <- proc.time()
#result1 <- lapply(1:5000,wrapper)
#rtime1 <- proc.time()-btime
#system.time(result1 <- lapply(1:5000,wrapper))
# Result is always in list form.
mean(unlist(result))

# 7. Stop snowfall
sfStop()


q()




