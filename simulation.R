#############################################################################
###########   Simulate Data 
#############################################################################
library(MCMCpack)
library(truncnorm)

J=60   # observation per individual
I=300  #number of individuals

### generate X
Xn=5  # number of variables
X <- matrix(NA,I*J,Xn,TRUE)
set.seed(1234)
X[,1] <- 1
X[,2] <-runif(I*J,0,9)
X[,3] <-sample(c(0,1),I*J,TRUE,c(0.95,0.05))
X[,4] <-rpois(I*J,c(1,2,3,4))  # nth session
X[,5] <-rnorm(I*J,8,5)


### generate W
Wn=4  # number of variables
W <- matrix(NA,I*J,Wn,TRUE)
set.seed(1234)
W[,1] <- 1
W[,2] <-sample(c(0,1),I*J,TRUE,c(0.7,0.3))
W[,3] <-sample(c(0,1),I*J,TRUE,c(0.6,0.4))
W[,4] <-rtruncnorm(I*J,0,1,0.02)

Z <- matrix(NA,I*J,1,TRUE)
Z[,1] <-sample(c(-1,1),I*J,TRUE,c(0.85,0.15))

# specify two states
K=2

### True parameters
sd_true <-2.5
set.seed(123)
sigma_true <- riwish(10,diag(5,5))
v_true<-0.5
mu_true<-1
delta_true <- 0.2
eta_true <- -0.6

beta_true <-matrix(c(-3,0.2,-2.5,2,2.5,-2,0.5,-3,1,0.8),Xn,K) 

thetabar_true <-c(0.1,0.5,1,-2,0.1)  ## length=Wn+1
thetai_true <- matrix(NA,Wn+1,I)
for (i in 1:I) {
  thetai_true[,i]<- mvrnorm(1,thetabar_true,sigma_true)
}

### specify initial status 
Delta0_true <- 2
A_true<-matrix(NA,I,J,TRUE)
for (i in 1:I) {
  A_true[i,1] <- rnorm(1,0,sqrt(Delta0_true))
}


### simulate state chain
Y_star <- matrix(NA,I*J,1,TRUE)
Y <- matrix(NA,I*J,1,TRUE)
U_true <- matrix(NA,I,J,TRUE)
S0_true<-matrix(NA,I,J,TRUE)

for (i in 1:I) 
{
  for (j in 1:J) 
  {
    if (j>1) {
      A_true[i,j] <- A_true[i,j-1]+S0_true[i,j-1]*delta_true*Z[J*(i-1)+j]+(1-S0_true[i,j-1])*eta_true*(Z[J*(i-1)+j])+rnorm(1,0,sqrt(v_true))
    } else {
      A_true[i,j] <- rnorm(1,0,sqrt(Delta0_true))
    }
    U_true[i,j] <- thetai_true[1:4,i]%*%W[J*(i-1)+j,]+thetai_true[5,i]*A_true[i,j]+rnorm(1)
    S0_true[i,j] <- cut(U_true[i,j],c(-Inf,0,Inf),labels = FALSE)-1
    Y_star[J*(i-1)+j] <- X[J*(i-1)+j,]%*%beta_true[,S0_true[i,j]+1] + rnorm(1,mean=0,sd=sqrt(sd_true))
    Y[J*(i-1)+j] <- max(Y_star[J*(i-1)+j],0)
  }
}
S_true <- as.vector(t(S0_true))
simulation <- cbind(Y,X,W,Z)
write.table(simulation, "simulation.csv", sep="\t")
