acquiextraction<-function(x,n_factors,corr,raw_data){

  if (raw_data==F)
    {
    R=x
    }
  else {

    if (corr=="Pearson"){
      R <- cor(x)
    }
    if (corr=="Polychoric"){
      R<-(psych::polychoric(x,smooth = FALSE, correct = FALSE))$rho
    }

  }

  m=size(R)[2]

  out <- EFA.MRFA::mrfa(R,n_factors,2, 0.000001, 0.000001,disp=FALSE)

  A <- out$A
  Rr<- out$Matrix

  e_variance <- sum(diag(transpose(A)%*%A)) / sum(diag(Rr))

  #centroid
  u=matrix(1,m,1)
  a=(Rr%*%u)/ matrix(((transpose(u)%*%Rr%*%u) ^(1/2)),m)

  X = solve(transpose(A)%*%A) %*% transpose(A) %*% a

  U = X %*% sqrt(solve((transpose(X)%*%X)))

  T_=U

  At<- A%*%T_


  cong=congru(a,At[,1])

  AQ_var= sum(diag(At[,1]))/sum(diag(Rr))

  acq = At[,1]
  Res = R-transpose(acq)%*%acq

  OUT<-list("Res"=Res, "acq"=transpose(acq), "AQ_var"=AQ_var)
  return(OUT)

}
