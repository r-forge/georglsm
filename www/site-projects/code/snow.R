




require(snow)

# starts up a set of R processes and assigns to the variable cl a reference to this collection of processes
cl <- makeCluster(4, "MPI")
#################################
### Intermediate Level Functions
#################################
# 1) clusterCall(cl, fun, x): asks each worker process to call fun with argument x and returns a vector of the results.
print(clusterCall(cl, function() Sys.info()[c("nodename","machine")]))
clusterCall(cl, min, rnorm(5))
clusterCall(cl, function(x) x, 1:2)
# 2.1) clusterApply(cl, vec, fun): assigns one element of the vector vec to each worker, has each worker call the function fun on its assigned element, and collects the result. If there are more elements than workers then workers are reused cyclically.
print(clusterApply(cl, 1:5, function(x) Sys.info()[c("nodename","machine")]))
clusterApply(cl, 1:10, function(x) print(1:x))
# 2.2) a load balancing version of clusterApply()
clusterApplyLB(cl, 1:10, function(x) print(1:x))
# remember to stop processes; otherwise, they keep running!
stopCluster(cl)
#################################
### High Level Functions
#################################
# parLapply(), parApply(), ...
# These higher level functions are easy to implement in terms of the intermediate level functions like clusterApply.
# parLapply <- function(cl, x, fun, ...) {
#   y <- splitList(x, length(cl)) 
#   val <- clusterApply(cl, y, lapply, fun, ...)
#   docall(c, val)
# }
# splitList(x, length(cl)): uses splitList() to split the input into P approximately equal-sized pieces for P worker processes;
splitList(1:10, 4)
#################################
### Example 
#################################
# Goal: apply func() on each row of the data
# divide data, then run each chuck on different workers
n <- 1e5
s <- matrix(rnorm(n*100),n,)
# 1) serial way
t0 <- proc.time()
v <- apply(s, 1, mean)
t1 <- proc.time() - t0
# 2) parallel way by using clusterApply() (intermediate level func)
cl <- makeCluster(4, "MPI")
t0 <- proc.time()
idx <- clusterSplit(cl, 1:dim(s)[1])
ssplt <- lapply(idx, function(i) s[i,])
v2 <- clusterApply(cl, ssplt, function(x) apply(x, 1, mean) )
t2 <- proc.time() - t0
stopCluster(cl)
# 3) parallel way by using parApply() (high level func)
# it seems slower than 2) and uses more memory!
cl <- makeCluster(4, "MPI")
t0 <- proc.time()
v3 <- parApply(cl, s, 1, mean) 
t3 <- proc.time() - t0
stopCluster(cl)
#
mean(unlist(v))
var(unlist(v))
mean(v3); var(v3)

# The power of snow lies in the ability to use the apply-style
# paradigm over a cluster of machines.
# In other words, we will be running four copies of myBigFunction() at once.
# So the snow package provides a unifying framework for
# parallelly executed apply functions.
A <- matrix(0,5000,5000)
mean(replicate(10, system.time( apply(A, 1, FUN=function(x) mean(exp(x))) ))["elapsed",])
# 3.9851
cl <- makeCluster(4, "MPI")
mean(replicate(10, system.time( parApply(cl, A, 1, FUN=function(x) mean(exp(x))) ))["elapsed",])
# 10.5032
stopCluster(cl) 


clusterSetupRNG(cl)
clusterSetupSPRNG(cl)
clusterSetupSPRNG(cl, seed=1234)



