
/*******************************/

proc univariate data=no2_3 normal;
var  y x1-x7;probplot y ;
run;

proc corr data=no2_3;
   var x1-x7;
   run;

proc glm data=no2_3;
model y=x1-x7 ;
run;

proc multtest data=no2_3;
model y=x1 x2 x3 x4/ bon ;
run;
 
proc iml;                                                                                                                               
use no2_3;                                                                                                                        
read  all var {x1 x2 x3 x4 x5 x6 x7}  into x;                                                                                                          
read  all var {y}  into y;                                                                                                              
/*print x y; run;*/                                                                                                                     
                                                                                                                                        
n=nrow(x);            /* Number of observations;*/
k=ncol(x);            /* Number of parameters including the intercept; */
                                                                                                                    
j=j(n,1,1); 
x10=j||x; 
          /* Display the design matrix  */                                                                                                              
cov_x=inv(x10`*x10); 
xpy=x10`*y;           /* The vector (X'X)^-1*Y  ;*/
beta=cov_x*xpy;     
PRINT beta;
          /* The estimated regression parameters;*/

/*Table 1 ANOVA for fitting regression*/

/* The fitted values, the residuals, SSE, and MSE ;*/
ssr=beta`*x10`*y;       /*SSR= sum of square of residual;*/
dfreq=k+1;             /*Degree of freedom of SSR;*/
print ssr dfreq;

msr=ssr/dfreq;         /*Mean square of residual;*/
sse=y`*y-beta`*x10`*y;  /* SSE = Sum of squares of residuals;*/
dferr=n-k-1;           /* Degrees of freedom of SSE;*/
mse=sse/dferr;        /* MSE = SSE/dferror;*/
print msr sse dferr mse;

sst=y`*y;             /* SST= sum square of total;*/
dftot=n;               /*Degrees of freedom of total;*/
fstat=msr/mse;         /* F-statistics; */ 
print sst dftot fstat; 


/*Table 2 ANOVA*/
J=j(n,n,1);
beta00=beta[2:k+1];
xbar_t=j`*x/n;
x_b=(I(n)-J/n)*x;

SSRm=beta00`*x_b`*y;
MSRm=SSRm/k;
print SSRm MSRm;

fstatRm=MSRm/MSE;
print SSE MSE fstatRm;

SSTm=y`*(I(n)-J/n)*y;
print SSTm;

/*Table 3 ANOVA showing in the term mean*/
SSM=y`*J*y/n;
MSM=SSM/1;
print SSM MSM;

fstatM=MSM/MSE;
print SSE MSE fstatM fstatRm;

SST=y`*y;
print SST;

/*Bonferroni Test of beta_j*/
var_beta=MSE*cov_x;
print var_beta;


b_1=beta[2,1]/sqrt(var_beta[2,2]);
b_2=beta[3,1]/sqrt(var_beta[3,3]);
b_3=beta[4,1]/sqrt(var_beta[4,4]);
b_4=beta[5,1]/sqrt(var_beta[5,5]);
b_5=beta[6,1]/sqrt(var_beta[6,6]);
b_6=beta[7,1]/sqrt(var_beta[7,7]);
b_7=beta[8,1]/sqrt(var_beta[8,8]);
print b_1 b_2 b_3 b_4 b_5 b_6 b_7;

p1=probt(b_1,n-k-1);
p2=2*probt(b_2,n-k-1);
p3=3*probt(b_3,n-k-1);
p4=4*probt(b_4,n-k-1);
p5=5*probt(b_5,n-k-1);
p6=6*probt(b_6,n-k-1);
p7=7*probt(b_7,n-k-1);
print p1 p2 p3 p4 p5 p6 p7;


t=tinv(.025,9);
print t;
p=probt(t,9);
print p;
p_0=probt(-4.48,9);
print p_0;


