extensionacq<-function(Z,LAM,PHI,nitem){

  LAM<-as.matrix(LAM)
  Z<-as.matrix(Z)

  n <- size(LAM)[1]

  corrmatrix= cor(Z)
  corrvec = corrmatrix[,nitem]

  if (nitem==1){
    corrvec=transpose(corrvec[2:n])
    TMP2 = LAM[2:n,]
  }
  else {
    if (nitem==n){
      corrvec=transpose(corrvec[1:(n-1)])
      TMP2 = LAM[2:n,]
    }
    else {
      buff1 = corrvec[1:nitem-1]

      if (nitem>2) {
        buff1 = transpose(buff1)
      }
      buff2 = transpose(corrvec[(nitem+1):n])
      corrvec = rbind(buff1,buff2)

      buff3 = LAM[1:(nitem-1),]
      buff4 = LAM[(nitem+1):n,]
      TMP2 = rbind(buff3,buff4)
    }
  }

  SRED = TMP2 %*% PHI

  estiload = (solve(transpose(SRED) %*% SRED)) %*% (transpose(SRED) %*% corrvec)


  if (nitem==1){
    LAM_nova = rbind(transpose(estiload),LAM[2:n,])
    TMP = transpose(cbind(transpose(estiload) %*% transpose(PHI) %*% transpose(LAM_nova)))
    corrvec_nova = transpose(TMP[2:n])
    }

  else {
    if (nitem==n){
      LAM_nova = rbind(transpose(estiload),LAM[1:(n-1),])
      TMP = transpose(cbind(transpose(estiload) %*% transpose(PHI) %*% transpose(LAM_nova)))
      corrvec_nova = transpose(TMP[1:(n-1)])
    }
    else {
      buff5 = LAM[1:(nitem-1),]
      buff6 = transpose(estiload)
      buff7 = LAM[(nitem+1):n,]
      LAM_nova = rbind(buff5,buff6,buff7)

      TMP = transpose((transpose(estiload) %*% transpose(PHI) %*% transpose(LAM_nova)))
      buff8 = TMP[1:(nitem-1),]
      if (nitem>2){
        buff8 = transpose(buff8)
      }
      buff9 = transpose(TMP[(nitem+1):n,])
      corrvec_nova = rbind(buff8,buff9)
    }
  }

  fit_index = sqrt(mean((corrvec - corrvec_nova)^2))

  OUT<-list("estiload"=estiload, "fit_index"=fit_index)
  return(OUT)

}
