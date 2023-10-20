rm(list=ls())



minsize <- 0.01 ### make this smaller than observed value; all predicted values MUST be above this value

maxsize <- 5.75 ### make this larger than observed value; all predicted values MUST be below this value

nn <- 200 ### nn is the size of the matrix

tot.length <- 300

max.age <- 300

library(Matrix)



###################################################################



##### Compute the functions for the IPM ###########################



##### Thanks to Stephen Ellner for some of the code used here #####



####################################################################



#Survival function S(z,t)



S.fun <- function(z,intercept,z.slope,n.slope,year.eff,id.eff,N) {
  
  u<-exp(intercept+z.slope*z+n.slope*N+year.eff+id.eff)
  
  return((u/(1+u)))
  
}



# Development function G(z'|z) zz = z' for R code

G.fun <- function(z,zz,intercept.mu,z.slope.mu,n.slope.mu,year.eff.mu,id.eff.mu,intercept.va,z.slope.va,n.slope.va,year.eff.va,id.eff.va,N) {
  
  mu.z <- intercept.mu+z.slope.mu*z+n.slope.mu*N+year.eff.mu+id.eff.mu
  
  sigma.z2 <- intercept.va+z.slope.va*z+n.slope.va*N+year.eff.va+id.eff.va
  
  sigma.z2 <- ifelse(sigma.z2<0,0.0001,sigma.z2)
  
  sigma.z <- sqrt(sigma.z2)
  
  temp1 <- sqrt(2*pi)*sigma.z
  
  temp2 <- ((zz-mu.z)^2)/(2*sigma.z2)
  
  return(exp(-temp2)/temp1)
  
}



R.fun <- function(z,intercept,z.slope,n.slope,year.eff,id.eff,N){
  
  u <- exp(intercept+z.slope*z+n.slope*N+year.eff+id.eff)
  
  u <- ifelse(z<4.24,0,u)
  
  return(u)
  
}



# Inheritance function D(z'|z) zz = z' for R code. Note same structure as G in this case -- could have different probability distribution if needed

D.fun <- function(z,zz, intercept.mu,z.slope.mu,n.slope.mu,year.eff.mu,id.eff.mu,intercept.va,z.slope.va,n.slope.va,year.eff.va,id.eff.va,N) {
  
  mu.z <- intercept.mu+z.slope.mu*z+n.slope.mu*N+year.eff.mu+id.eff.mu
  
  sigma.z2 <- intercept.va+z.slope.va*z+n.slope.va*N+year.eff.va+id.eff.va
  
  sigma.z2 <- ifelse(sigma.z2<0,0.0001,sigma.z2)
  
  sigma.z <- sqrt(sigma.z2)
  
  temp1 <- sqrt(2*pi)*sigma.z
  
  temp2 <- ((zz-mu.z)^2)/(2*sigma.z2)
  
  return(exp(-temp2)/temp1)
  
}



############## The 'big matrix' M of size n x n



#### Thanks to Stephen Ellner for code -- taken from his website



bigmatrix<-function(n,s.params,r.params,g.params,d.params,N) {
  
  # boundary points b and mesh points y
  
  b <- minsize+c(0:n)*(maxsize-minsize)/n
  
  y <- 0.5*(b[1:n]+b[2:(n+1)])
  
  # create S, R, G and D matrices
  
  S <- (diag(S.fun(y,s.params[1],s.params[2],s.params[3],s.params[4], s.params[5],N)))
  
  R <- (diag(R.fun(y,r.params[1],r.params[2],r.params[3],r.params[4], r.params[5],N)))
  
  G <- (t(outer(y,y,G.fun,g.params[1],g.params[2],g.params[3], g.params[4],g.params[5],g.params[6],g.params[7],g.params[8],g.params[9], g.params[10],N)))
  
  D <- (t(outer(y,y,D.fun,d.params[1],d.params[2],d.params[3], d.params[4],d.params[5],d.params[6],d.params[7],d.params[8],d.params[9], d.params[10],N)))
  
  # scale D and G so columns sum to 1
  
  G <- G/matrix(as.vector(apply(G,2,sum)),nrow=n,ncol=n,byrow=TRUE)
  
  D <- D/matrix(as.vector(apply(D,2,sum)),nrow=n,ncol=n,byrow=TRUE)
  
  return(list(S=S,R=R,G=G,D=D,meshpts=y))
  
}



get.eigen.stuff <- function(mat){
  
  sz <- dim(mat)[1]
  
  t.now <- runif(sz)
  
  t.now <- t.now/sum(t.now)
  
  t.next <- mat%*%t.now
  
  t.next <- t.next/sum(t.next)
  
  i <- 0
  
  while (sum(abs(t.next-t.now))>0.0000001){
    
    i <- i+1
    
    print(i)
    
    t.now <- t.next
    
    t.next <- mat%*%t.now
    
    lambda <- sum(t.next)/sum(t.now)
    
    t.next <- t.next/sum(t.next)
    
  }
  
  r.now <- runif(sz)
  
  r.now <- r.now/sum(r.now)
  
  r.next <- r.now%*%mat
  
  r.next <- r.next/sum(r.next)
  
  while (sum(abs(r.next-r.now))>0.0000001){
    
    r.now <- r.next
    
    r.next <- r.now%*%mat
    
    r.next <- r.next/sum(r.next)
    
  }
  
  return(list(lambda,t.next,r.next))
  
}



b <- minsize+c(0:nn)*(maxsize-minsize)/nn

y <- 0.5*(b[1:nn]+b[2:(nn+1)])

##############################################################

# density-dependent version of code  REMEMBER LOG SCALE

s.params <- c(-0.4855082,0.3809523,0,0,0) # survival - DONE

r.params <- c(0.6931472,0,-0.001,0,0) # reproductive - DONE

#g.params <- c(2,0.654825-0.05,0,0,0,0.1,0,0,0,0) #growth - DONE

#g.params <- c(2.1373,0.5496,0,0,0,0.3,-0.1,0,0,0) #growth - DONE

g.params <- c(2.1373,0.5496,0,0,0,0.1,0,0,0,0) #growth - DONE

d.params <- c(4,0,0,0,0,0.1,0,0,0,0) # Inheritance - DONE

N.res <- array(NA,c(nn,tot.length))

N.res[,1] <- runif(nn)



for (i in 2:tot.length){
  
  N <- N.res[,i-1]
  
  N.tot <- sum(N)
  
  m <- bigmatrix(nn,s.params,r.params,g.params,d.params,N.tot)
  
  mat <- m$G%*%m$S + m$D%*%m$R
  
  N.res[,i] <- mat%*%N
  
}



pops <- apply(N.res,2,sum)