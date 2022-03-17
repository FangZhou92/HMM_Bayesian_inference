library(MCMCpack)
library(truncnorm)

#############################################################################
###########   0 - Read Data 
#############################################################################
d=read.csv("simulation.csv",sep = "\t")
X<- data.matrix(d[,2:6])
Xn<-5
Y<- data.matrix(d[,1])
W<- data.matrix(d[,7:12])  # variables in the past period W(i,t-1)
Wn<-6
I<-200
J<-30
K<-3
#############################################################################
###########   1 - Set Hyperparameters 
#############################################################################
alpha<-5
delta0<-10


m<-matrix(0,Xn,K,TRUE)
M<-list()
set.seed(123)
M[[1]]<-diag(x=5,Xn)
set.seed(124)
M[[2]]<-diag(x=5,Xn)
set.seed(125)
M[[3]]<-diag(x=5,Xn)

mw<-matrix(0,Wn,K,TRUE)
Mw<-list()
set.seed(123)
Mw[[1]]<-diag(x=5,Wn)
set.seed(124)
Mw[[2]]<-diag(x=5,Wn)
set.seed(125)
Mw[[3]]<-diag(x=5,Wn)

alpha0<-c(1,1,1)

# transition probability
trans_prob_f <- function (i,j) { # transition from k to l
  trans_prob <- matrix(NA,K,K)
  for (l in 1:K) 
  {
    for (k in 1:K)
    {
      if (l==1) {
        trans_prob[l,k] <- pnorm(theta4[1]-W[J*(i-1)+j,]%*%theta3[,k])
      } else if (l>1 & l<K) {
        trans_prob[l,k] <- pnorm(theta4[l]-W[J*(i-1)+j,]%*%theta3[,k])-pnorm(theta4[l-1]-W[J*(i-1)+j,]%*%theta3[,k])
      } else {
        trans_prob[l,k] <- 1- pnorm(theta4[K-1]-W[J*(i-1)+j,]%*%theta3[,k])
      }
    }
  }
  return(trans_prob)
} 

# emmision probability 
emm_prob_f <- function (i,j) {
  emm_prob <- rep(NA,K)
  delta=sqrt(theta1)
  for (k in 1:K) 
  {
    beta <- theta2[,k]
    emm_prob[k] <- (1-pnorm(X[J*(i-1)+j,]%*%beta/delta))^(ifelse(Y[J*(i-1)+j]==0,1,0))*((1/delta)*dnorm((Y[J*(i-1)+j]-X[J*(i-1)+j,]%*%beta)/delta))^(ifelse(Y[J*(i-1)+j]>0,1,0))
  }
  return(emm_prob)
}

#############################################################################
###########   2 - Gibbs sampling  
#############################################################################
# set starting value
theta1 <-4
theta2 <- matrix(c(-2,1,2,2,3,1.5,1,2,2.5,3.8,1.5,2,2,3,3),Xn,K) 
theta3 <- matrix(c(-1,1,2.5,1.8,2,3.2,-2,1.4,4.2,2,1,3.5,2,-1,1.6,1.7,1,1),Wn,K)
#theta2 <- matrix(1,Xn,K)
#theta3 <- matrix(1,Wn,K)
theta4 <-c(0,13,Inf)


S0 <- matrix(NA,I,J)
#S<-sample(c(1:K),size=I*J,replace=TRUE)  # 要 draw data 
L <- matrix(NA,I*J,1,TRUE)


# for storing estimates of interests
Iter<-5000  #number of iteration
Delta_sq <- rep(NA,Iter)  #theta1
Beta <-list()  #theta2
Ksi <- list()  #theta3
Mu <- matrix(NA,Iter,K)   #theta4

#BWD <- matrix(NA,I,K,TRUE)

for (iter in 1:Iter) 
{
  
  #step6: generate theta6={s1,s2,...sT} using forward-backward algorithm
  # initial probability
  
  start_prob <- rdirichlet(1,alpha0)
  
  for (i in 1:I)
  {
    #step 6a:forward algorithm 
    fwd<-matrix(NA,J,K)
    for (j in 1:J)
    {
      emm_prob <- emm_prob_f(i,j)
      if (j ==1) {
        prev_f_sum <- matrix(start_prob,3,1)
      } else {
        trans_prob <- trans_prob_f(i,j)
        prev_f_sum <- trans_prob%*%f_prev
      }
      f_curr <- emm_prob*prev_f_sum
      fwd[j,] <- f_curr/sum(f_curr)  # length = K
      f_prev <- f_curr/sum(f_curr)  # length = K
    }
    S0[i,J] <- sample(c(1:K),size=1,replace=TRUE,fwd[J,]/sum(fwd[J,]))
    
    # step 6b:backward algorithm
    bwd<-matrix(NA,J,K)
    for (j in (J-1):1) 
    { 
      trans_prob <- trans_prob_f(i,j+1)
      for (k in 1:K)
      {
        numerator <- trans_prob[S0[i,j+1],k] * fwd[j,k]
        denominator <- trans_prob[S0[i,j+1],]%*%fwd[j,]
        bwd[j,k] <- numerator/denominator
      }
      S0[i,j] <- sample(c(1:K),size=1,prob = bwd[j,])
    }
    #BWD[i,] <- bwd[1,]
  }
  S <- as.vector(t(S0))
  
  #step5: generate L from truncated normal distribution 
  new_mu <- c(-Inf,theta4)  #add mu_0
  S_pre <- matrix(NA,I*J,1)
  for (i in 1:I)
  {
    for (j in 1:J) {
      if (j>1) {
        S_pre[J*(i-1)+j] <- S[J*(i-1)+j-1]
      } else {
        S_pre[J*(i-1)+j] <- NA
      }
    }
  }
  
  for (i in 1:I) 
  {
    for (j in 2:J)
    {
      Left <- new_mu[S[J*(i-1)+j]]
      Right<- new_mu[S[J*(i-1)+j]+1]
      uni_mean <- W[J*(i-1)+j,]%*%theta3[,S_pre[J*(i-1)+j]]
      L[J*(i-1)+j] <- rtruncnorm(1,a=Left,b=Right,mean=uni_mean,sd=1)
    }
  }
  
  #step4: generate theta4 from uniform distribution
  for (k in 2:2)
  {
    if (k==K-1) {
      lower <- max(max(L[S==k,],na.rm = TRUE),theta4[k-1])
      upper <- min(L[S==k+1,],na.rm = TRUE)
    } else {
      lower <- max(max(L[ifelse(S==k,1,0)>0],na.rm = TRUE))
      upper <- min(min(L[ifelse(S==k+1,1,0)>0],na.rm = TRUE),theta4[k+1])
    }
    theta4[k]<-runif(1,min = lower, max = upper)
  }
  Mu[iter,] <- theta4
  
  #### step3: generate theta3 from normal distribution
  for (k in 1:K)
  {    
    a <- matrix(NA,I*J,1,TRUE)
    a <- ifelse(S_pre==k&!is.na(S_pre),1,0)
    Wnew = W[a==1,]
    Mw_star=solve(solve(Mw[[k]])+crossprod(Wnew))
    Lnew = L[a==1]
    mw_star=Mw_star%*%(solve(Mw[[k]])%*%mw[,k]+crossprod(Wnew,Lnew))
    theta3[,k]<-mvrnorm(1,mw_star,Mw_star)
  }
  Ksi[[iter]] <- theta3
  
  ########## generate y* for each observation such that y*>0 if Y>0 and y*<0 if Y=0
  y_star <- matrix(NA,I*J,TRUE)
  for (i in 1:I)
  {
    for (j in 1:J)
    {
      if (Y[J*(i-1)+j]>0) {
        y_star[J*(i-1)+j] <- Y[J*(i-1)+j]
        #y_star[J*(i-1)+j] <- rtruncnorm(1,a=0,mean=X[J*(i-1)+j,]%*%theta2[,S[J*(i-1)+j]],sd=sqrt(theta1))
      } else {
        y_star[J*(i-1)+j] <- rtruncnorm(1,b=0,mean=X[J*(i-1)+j,]%*%theta2[,S[J*(i-1)+j]],sd=sqrt(theta1))
      }
    }
  }
  
  #### step2:generate theta2 from normal distribution
  for (k in 1:K)
  {
    #b <- matrix(NA,I*J,1,TRUE)
    #b <- ifelse(S==k,1,0)
    Xnew = X[S==k,]
    M_star=solve(solve(M[[k]])+crossprod(Xnew)/theta1)
    ynew = y_star[S==k,]
    m_star=M_star%*%(crossprod(Xnew,ynew)/theta1+solve(M[[k]])%*%m[,k])
    theta2[,k]<-mvrnorm(1,m_star,M_star)
  }
  Beta[[iter]] <- theta2
  
  ######### step 1: generate theta1 from gamma distribution
  ssr = rep(0,I*J)
  for (i in 1:I)
  {
    for (j in 1:J)
    {
      ssr[J*(i-1)+j] <- y_star[J*(i-1)+j]-(X[J*(i-1)+j,]%*%theta2[,S[J*(i-1)+j]])
    }
  }
  theta1<-rinvgamma(1,I*J/2+alpha,crossprod(ssr)/2+delta0)
  Delta_sq[iter] <- theta1
  
}


