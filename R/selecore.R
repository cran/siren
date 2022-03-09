selecore<-function(Rminus,Rplus){

  m=size(Rminus)[2]
  n=size(Rplus)[2]

  umin=matrix(1,m,1)
  uplus=matrix(1,n,1)

  buff1=Rminus%*%umin
  buff2=(transpose(umin)%*% Rminus %*% umin) ^(1/2)

  vecminus = buff1 / buff2[col(buff1)]

  buff3=Rplus%*%uplus
  buff4=(transpose(uplus)%*% Rplus %*% uplus) ^(1/2)

  vecplus= buff3 / buff4[col(buff3)]

  OUT<-list("vecminus"=vecminus, "vecplus"=vecplus)
  return(OUT)

}
