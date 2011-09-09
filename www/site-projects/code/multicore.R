
##################################################################
# The multicore package by Simon Urbanek provides a
# convenient interface to locally running parallel computations in
# R on machines with multiple cores or CPUs. Jobs can share
# the entire initial workspace.
#
# The multicore package provides two main interfaces:
# 1) mclapply() a parallel / multicore version of lapply
# 2) the functions parallel() and collect() to launch parallel
# execution and gather results at end



require(multicore)

# 1) mclapply(): multicore version of "lapply"
  mclapply(1:30, rnorm)
  # use the same random numbers for all values
  set.seed(1)
  mclapply(1:30, rnorm, mc.preschedule=FALSE, mc.set.seed=FALSE)
  # something a bit bigger - albeit still useless :P
  unlist(mclapply(1:32, function(x) sum(rnorm(1e7))))
  
  options(cores=4) # number of cores to distribute the job to
  getOption("cores")
  system.time( unlist( mclapply(1:3, function(x) rnorm(10)) ) )
  system.time( unlist(mclapply(1:30, function(x) sum(rnorm(1e6)))) )
  #   user  system elapsed 
  # 14.450   0.490   2.577
  system.time( unlist(lapply(1:30, function(x) sum(rnorm(1e6)))) )
  #   user  system elapsed 
  # 11.060   0.000  11.067
  
  gc(verbose = getOption("verbose"), reset=TRUE)

# Actual job
system.time(do.call("cbind", 
    mclapply(subset(sampdata, select = c(a:z)), function(x) tapply(x, sampdata$groupid, sum))
))


# 2) parallel() and collect()
p <- parallel(1:10)
q <- parallel(1:20)
collect(list(p, q)) # wait for jobs to finish and collect all results
p <- parallel(1:10)
collect(p, wait=FALSE, 10) # will retrieve the result (since it's fast)
collect(p, wait=FALSE) # will signal the job as terminating
collect(p, wait=FALSE) # there is no such job
# a naive parallelized lapply can be created using parallel alone:
jobs <- lapply(1:10, function(x) parallel(rnorm(x*1e5), name=x))
collect(jobs)

# another example
wrapper <- function(x){
  A <- matrix(rnorm(1e6),1e3,); 
  det(solve(A)%*%A) }
jobs <- lapply(1:10, function(k) parallel(wrapper(k)) )
jobs <- lapply(1:10, function(k) 
            parallel({
                    A <- matrix(rnorm(100),10,)
                    det(solve(A)%*%A) }) 
               )
collect(jobs)






