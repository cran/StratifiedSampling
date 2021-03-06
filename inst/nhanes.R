rm(list = ls())

library(survey)
data(nhanes)

nhanes <- nhanes[-which(is.na(nhanes$HI_CHOL)),]
nhanes$SDMVSTRA <- cleanstrata(nhanes$SDMVSTRA)
nhanes$lev <- nhanes$SDMVSTRA
# nhanes$lev <- as.numeric(as.factor(as.character(paste(nhanes$race,nhanes$agecat,nha nes$RIAGENDR,sep = ""))))



nhanes <- nhanes[order(nhanes$lev),]

nhanes_split <- split(nhanes,f = nhanes$lev)
strata <- nhanes$lev


N <- nrow(nhanes)
n <- 500
Nh <- table(nhanes$lev)
sh <- as.vector(do.call(rbind,lapply(nhanes_split,FUN = function(x){ var(x$HI_CHOL) })))

den <- sum(Nh*sh)

# nh <- ceiling(n*Nh*sh/den)
nh <- n*Nh*sh/den

y <- nhanes$HI_CHOL
pik <- as.vector(inclusionprobastrata(strata,nh))
X <- as.matrix(nhanes[,c(5,7)])
s <- stratifiedcube(X,strata,pik)

# s <- fbs(X,as.matrix(strata),pik)

t(X/pik)%*%pik
t(X/pik)%*%s

varApp(X,as.matrix(strata),pik,y)



sum(y[which(s==1)]/pik[which(s == 1)])
sum(y)




  
  
N <- 1000

x1 <- rgamma(N,4,25)
x2 <- rgamma(N,4,25)
x3 <- rgamma(N,4,25)
x4 <- rgamma(N,4,25)

nstrata <- 25

strata <- as.matrix(rep((1:nstrata),each = (N/nstrata))) 
Xcat <- disjMatrix(strata)

#-------- CASE 1 pik equal and integer in each strata

pik <- inclusionprobastrata(strata,rep(2,nstrata)) #↓ one by strata
sum(pik)
X <- as.matrix(matrix(c(x1,x2),ncol = 2))

# 
# s <- fbs(X,strata,pik)

y1 <- 20*strata + rnorm(1000,120)
y2 <- 500 + 5*x1 + 5*x2 + rnorm(1000,0,270)
y3 <- 500 + 100*x1 + 100*x2 + 100*x3 + 100*x4 + rnorm(1000,0,1000)
y4 <- 500 + 200*x1 + 100*x2 + 100*x3 + 50*x4 + rnorm(1000,0,1000)


SIM = 1000
t1 <- matrix(rep(0,2*SIM),ncol = 2)
varApp(X,strata,pik,y1)



for(i in 1:SIM){
  print(i)
  
  
  # s <- stratifiedcube(cbind(pik,X),strata,pik)
  # s <- fbs(X,strata,pik)
  s <- balstrat(X,strata,pik)
  y <- y1
  y_ht <- sum(y[which(s==1)]/pik[which(s == 1)])
  t1[i,1] <- (sum(y_ht) - sum(y))^2
  t1[i,2] <- varEst(X,strata,pik,s,y)
}


colSums(t1)/SIM
varApp(X,strata,pik,y1)

