extensionacqm<-function(Z,LAM,PHI){

  n<-size(LAM)[1]

  for (i in 1:n){
    if (i==1){
      out0<- extensionacq(Z,LAM,PHI,i)
      buff <- out0$estiload
      bufff <- out0$fit_index
    }
    else {
      out1 <- extensionacq(Z,LAM,PHI,i)
      buff1 <- out1$estiload
      buff2 <- out1$fit_index

      buff <- cbind(buff,buff1)
      bufff <- rbind(bufff,buff2)
    }


  }

  buff=transpose(buff)

  OUT<-list("estiload"=buff, "fit_index"=bufff)
  return(OUT)


}

