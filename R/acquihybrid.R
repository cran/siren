acquihybrid<-function(x,content_factors,target, corr = "Pearson", raw_data=TRUE , method="fixed", display = TRUE){


  ######################################################################
  #  x : Raw sample scores (OR corr/cov matrix)
  ######################################################################

  if (missing(x)){
    stop("The argument x is not optional, please provide a valid raw sample scores or a correlation/covariance matrix")
  }

  ######################################################################
  #  content_factors : Number of content factors to be retained
  ######################################################################

  if (missing(content_factors)){
    stop("The argument content_factors is not optional, please provide the number of content factors to be retained")
  }

  n_factors<-content_factors

  ######################################################################
  #  target : The semi-specified target if procustes rotations are selected
  ######################################################################

  if (missing(target)){
    stop("The argument target is not optional, please provide a valid target matrix. For more information, check the documentation of acquihybrid.")
  }
  else {
    t1=size(target)[1]
    t2=size(target)[2]

    #check target size
  }

  ######################################################################
  #  corr: Determine if Pearson or Polychoric matrix will be used "Pearson" or "Polychoric"
  ######################################################################



  ######################################################################
  #  raw_data: If x contains raw data
  ######################################################################

  if (is.logical(raw_data)!=TRUE){
    stop("raw_data argument should be a logical variable.")
  }


  if (raw_data==T){
    if (corr=="Pearson"){
      R <- cor(x)
    }
    if (corr=="Polychoric"){
      R<-(psych::polychoric(x,smooth = FALSE, correct = FALSE))$rho
    }
  }
  else {
    R=x

    #if raw data is false, it is a correlations/covariance matrix. Method can only be "resid"
    method="resid"
    corr = "user"
  }


  if (raw_data==T){
    n_items=size(x)[2]
    content_factors<-round(content_factors)
    if (content_factors>(n_items/4)){
      stop("The argument content_factors has to be greater than the number of items / 4")
    }
  }

  if (n_items != size(target)[1]){stop("target should have a length equal to the number of variables in x")}


  ################ BEGIN ##################


  #internally, acq is considered as a factor
  n_factors=n_factors+1

  #1: Check if the balance is complete or partial (factor by factor)

  is_balanced=matrix(1,n_factors-1)

  for (i in 1:(n_factors-1)){
    if (sum(target[,i])!=0){is_balanced[i]=0}
  }

  all_balanced=F
  if (all(is_balanced==1)){all_balanced=T}

  if (all_balanced==T){
    # All factors are fully balanced, it is not necessary to find a balanced core.
    # Applying acquiextraction for the global matrix

    cat('Computing Acquiescence extraction: Please wait                                                    \r')
    flush.console()
    OUT <- acquiextraction(x,n_factors,corr,raw_data)
    AQ_var=OUT$AQ_var
    cat('Computing Acquiescence extraction: Done!                                                          \n')
    flush.console()
    cat("\r","                                                                                                  ","\r")


  }
  else {
    #At least one factor is not fully balanced

    #Search the partially balanced factors

    ptm_one <- proc.time()

    items_neg <- matrix(nrow=n_items,ncol=n_factors)
    items_pos <- matrix(nrow=n_items,ncol=n_factors)
    for (i in 1:(n_factors-1)){

      items_neg[,i]<-target[,i]<0
      items_pos[,i]<-target[,i]>0

      which.neg=which(items_neg[,i])
      which.pos=which(items_pos[,i])

      Rminus = R[items_neg[,i],items_neg[,i]]

      #canviar diagonals pel valor mes alt de la columna

      tmp = size(Rminus)[2]

      #check
      if (tmp<2) {stop('One factor contains less than two negative items. It is not possible to obtain a balanced core.')}

      for (j in 1:tmp){
        Rminus[j,j]= sort(Rminus[,j],partial=tmp-1)[tmp-1]
      }

      Rplus = R[items_pos[,i],items_pos[,i]]

      #canviar diagonals pel valor mes alt de la columna

      tmp2= size(Rplus)[2]

      #check
      if (tmp<2) {stop('One factor contains less than two positive items. It is not possible to obtain a balanced core.')}

      for (j in 1:tmp2){
        Rplus[j,j]= sort(Rplus[,j],partial=tmp2-1)[tmp2-1]
      }

      OUT2<- selecore(Rminus,Rplus)

      vecminus=OUT2$vecminus
      vecplus=OUT2$vecplus


      if (tmp==tmp2){

        #This factor is fully balanced, extension is not required.
        # (but will be used for the global balanced core)

        core=transpose(sort(c(which.pos,which.neg)))

        if (i==1){
          core_global=core
          out_core_global=vector()
        }
        else {
          core_global=rbind(core_global,core)
        }

        # OUT <- acquiextraction(x[,core],2,corr)
        # ja no cal fer acquiextraction factor a factor, ho fem al final i
        # fem l'extension dels que quedin fora el core global

      }

      if (tmp<tmp2){
        # More positive items than negative (usual case)

        #creem un vecplus2, on anirem treient els valors que vagin sortint
        vecplus2=vecplus

        parelles = matrix(nrow=tmp,ncol=2)

        for (k in 1:tmp){
          parelles[k,1]=which.neg[k]
          parelles[k,2]=which.pos[which.min(abs(vecminus[k]-vecplus2))]

          vecplus2[which.min(abs(vecminus[k]-vecplus2))]=NA
        }

        wh=sort(which.pos[match(parelles[,2],which.pos)] )
        out_core=transpose(which.pos[is.element(which.pos,wh)==FALSE])

        if (i==1){
          out_core_global=out_core
        }
        else {
          out_core_global=rbind(out_core_global,out_core)
        }

        core.pos=sort(which.pos[match(parelles[,2],which.pos)] )
        core.neg=which.neg
        core=transpose(sort(c(core.pos,core.neg)))

        if (i==1){
          core_global=core
        }
        else {
          core_global=rbind(core_global,core)
        }

        ## in "core" we have the core items, and in "out_core" the ones outside balanced core

        # Abans feiem estimació d'acq per aquest factor i extensió per extensionacq. Ara només ens quedem
        # amb quins items formen el core, els agrupem i al final farem l'extensió global nova (més simple)


      }

      if (tmp>tmp2){
        # More negative items than positive ones

        #creem un vecminus2, on anirem treient els valors que vagin sortint
        vecminus2=vecminus

        parelles = matrix(nrow=tmp2,ncol=2)

        for (k in 1:tmp2){
          parelles[k,1]=which.pos[k]
          parelles[k,2]=which.neg[which.min(abs(vecplus[k]-vecminus2))]

          vecminus2[which.min(abs(vecplus[k]-vecminus2))]=NA
        }

        wh=sort(which.neg[match(parelles[,2],which.neg)] )
        out_core=which.neg[is.element(which.neg,wh)==FALSE]

        if (i==1){
          out_core_global=out_core
        }
        else {
          out_core_global=rbind(out_core_global,out_core)
        }

        core.neg=sort(which.neg[match(parelles[,2],which.neg)] )
        core.pos=which.pos
        core=transpose(sort(c(core.pos,core.neg)))

        if (i==1){
          core_global=core
        }
        else {
          core_global=rbind(core_global,core)
        }

      }

      if (display==TRUE){

        compT <- proc.time() - ptm_one
        compT<-compT[3]
        compT<-compT*(n_factors-i)/i

        secondsInAMinute = 60
        secondsInAnHour = 60 * secondsInAMinute
        secondsInADay = 24 * secondsInAnHour

        days <- floor(compT / secondsInADay)

        hourSeconds <- compT %% secondsInADay
        hours <- floor(hourSeconds / secondsInAnHour)

        minuteSeconds <- hourSeconds %% secondsInAnHour
        minutes <- floor(minuteSeconds / secondsInAMinute)

        remainingSeconds <- minuteSeconds %% secondsInAMinute
        seconds <- ceiling(remainingSeconds)

        if (compT > 3600){
          if (days >= 1){ #Very very rare, but just to be sure
            cat("Computing acquihybrid. Time remaining: +24 hours                                                     \r")
            flush.console()
          }
          else {
            cat("Computing acquihybrid. Time remaining: ", hours,"hours, ",minutes, "minutes and ",seconds, "seconds \r")
            flush.console()
          }
        }
        else{
          if (compT >= 60){
            cat("Computing acquihybrid. Time remaining: ", minutes, "minutes and ",seconds,"seconds \r")
            flush.console()
          }
          if (compT < 60) {
            cat("Computing acquihybrid. Time remaining ",seconds,"seconds                                                                  \r")
            flush.console()
          }
        }
      }

    }
    if (display==TRUE){
      cat("\r","                                                                                                  ","\r")

      # close(pb)
    }


    # We have the final global core (core_global) and out_core_global

    # Extract the acq from the items of the global core
    cat('Computing Acquiescence extraction: Please wait                                                    \r')
    flush.console()
    OUT <- acquiextraction(x[,core_global],3,corr,raw_data)
    cat('Computing Acquiescence extraction: Done!                                                          \n')
    flush.console()
    cat("\r","                                                                                                  ","\r")


    acqcore=OUT$acq
    AQ_var=OUT$AQ_var

    n_out_core_global=size(out_core)[1]

    for (i in 1:n_out_core_global){

      # for each item outside the core, get the correlation with the core and compute the extension

      Ri = transpose(R[out_core_global[i],core_global])

      estiload= sum(Ri) / sum(acqcore)
      if (estiload<0) {estiload=0}

      if (i==1){
        estiload_out_core=estiload
      }
      else {
        estiload_out_core=rbind(estiload_out_core,estiload)
      }

    }

    acq_loadings=vector()
    j=1
    k=1
    for (i in 1:n_items){
      if (any(i==core_global)){
        #core item, get it from acqcore
        acq_loadings[i]=acqcore[j]
        j=j+1
      }
      else {
        #out core item, get it from estiload_out_core
        acq_loadings[i]=estiload_out_core[k]
        k=k+1
      }
    }

    acq_loadings=transpose(acq_loadings)

    Res = R- acq_loadings%*%transpose(acq_loadings)

    OUT<-list("Res"=Res, "acq"=acq_loadings)

  }

  ########### CONFIRMATORY FACTOR ANALYSIS (LAVAAN)


  if (method=="resid"){
    #UTILITZAR RESIDUALS (OUT$Res)

    var_names=""
    for (i in 1:n_items){
      buff_names=paste("V",i,sep="")
      if (i==1){var_names=buff_names}
      else {var_names=cbind(var_names, buff_names)}
    }
    colnames(x)=var_names

    text=""
    for (i in 1:(n_factors-1)){
      buff=which(!target[,i]==0)
      buff2=paste("V",buff,sep="", collapse=" + ")
      text=paste(text,"factor",i," =~",buff2,"\n")
    }



    text=""
    var_names=""
    for (i in 1:n_items){
      buff_names=paste("V",i,sep="")
      if (i==1){var_names=buff_names}
      else {var_names=cbind(var_names, buff_names)}
    }
    colnames(x)=var_names

    for (i in 1:(n_factors-1)){
      buff=which(!target[,i]==0)
      f2=size(buff)[2]

      buffNA=paste("NA*V", buff[1], sep="", " +")

      buff2=paste("V",buff[2:f2],sep="", collapse=" + ")

      buff3=paste("factor",i,sep="")

      text=paste(text,buff3," =~",buffNA,buff2, "\n")
    }

    other_text="\n "
    for (i in 1:(n_factors-1)){
      other_text=paste(other_text,"factor",i,"~~1*factor",i, "\n", sep="")
    }


    text=paste(text,other_text)



    fit <- lavaan::cfa(text,sample.cov=OUT$Res,sample.nobs=1000, estimator="ULS")

    items_order=matrix("",n_items)
    j=1
    k=0
    for (i in 1:(n_factors-1)){
      buff=which(!target[,i]==0)

      f2=size(buff)[2]
      k=f2+j -1

      items_order[j:k]=paste("V",buff,sep="")
      j=j+f2

      buff_acq=which(!target[,i]==0)
      if (i==1){items_order_n=buff_acq}
      else {items_order_n=c(items_order_n,buff_acq)}
    }

    factor_names<-matrix("",n_factors)
    for (i in 1:n_factors){
      if (i==n_factors){
        factor_names[i]="ACQ"
      }
      else {
        factor_names[i]=paste("F",i,sep = "")
      }
    }

    # concat

    loadings = cbind(fit@Model@GLIST$lambda,OUT$acq[items_order_n])
    psi <- fit@Model@GLIST$psi[1:n_factors-1,1:n_factors-1]

  }
  if (method=="fixed"){

    if (corr=="Pearson"){
      acq_loadings = OUT$acq * apply(x,2,sd)
    }
    else {
      acq_loadings = OUT$acq
    }

    text=""
    var_names=""
    for (i in 1:n_items){
      buff_names=paste("V",i,sep="")
      if (i==1){var_names=buff_names}
      else {var_names=cbind(var_names, buff_names)}
    }
    colnames(x)=var_names

    for (i in 1:(n_factors-1)){
      buff=which(!target[,i]==0)
      f2=size(buff)[2]

      buffNA=paste("NA*V", buff[1], sep="", " +")

      buff2=paste("V",buff[2:f2],sep="", collapse=" + ")

      buff3=paste("factor",i,sep="")

      text=paste(text,buff3," =~",buffNA,buff2, "\n")
    }

    # reallocating items labels to the new order

    items_order=matrix("",n_items)
    j=1
    k=0
    for (i in 1:(n_factors-1)){
      buff=which(!target[,i]==0)

      f2=size(buff)[2]
      k=f2+j -1

      items_order[j:k]=paste("V",buff,sep="")
      j=j+f2

      buff_acq=which(!target[,i]==0)
      if (i==1){items_order_n=buff_acq}
      else {items_order_n=c(items_order_n,buff_acq)}
    }

    acq_text=""
    for (i in 1:(n_items)){
      if (i != n_items){
        acq_text=paste(acq_text,round(acq_loadings[items_order_n[i]],3)," * V",items_order_n[i]," +",sep="")
      }
      else {
        acq_text=paste(acq_text,round(acq_loadings[items_order_n[i]],3)," * V",items_order_n[i],sep="")
      }
    }
    acq_text=paste("ACfactor =~ ",acq_text)

    other_text="\n "
    for (i in 1:(n_factors-1)){
      other_text=paste(other_text,"factor",i,"~~1*factor",i, "\n", sep="")
    }

    other_text2="\n"
    for (i in 1:(n_factors-1)){
      other_text2=paste(other_text2,"factor",i,"~~0*ACfactor \n", sep="")
    }

    text=paste(text,acq_text,other_text,other_text2)


    if (corr=="Pearson"){
      fit <- lavaan::cfa(text,data=x,estimator="ULSMV")
      loadings = fit@Model@GLIST$lambda
      loadings[,n_factors]=OUT$acq[items_order_n]
    }
    else {
      fit <- lavaan::cfa(text,data=x, ordered=TRUE, estimator="ULSMV")
      loadings = fit@Model@GLIST$lambda
    }

    factor_names<-matrix("",n_factors)
    for (i in 1:n_factors){
      if (i==n_factors){
        factor_names[i]="ACQ"
      }
      else {
        factor_names[i]=paste("F",i,sep = "")
      }
    }

    psi <- fit@Model@GLIST$psi[1:n_factors-1,1:n_factors-1]
  }

  rownames(loadings)<-items_order
  colnames(loadings)<-factor_names
  loadings=round(loadings,3)


  log_content=matrix(FALSE,n_items,n_factors-1)
  log_AC=transpose(loadings[,n_factors]<0)
  loadings[cbind(log_content,log_AC)] <- 0.001

  # factor scores:

  pfactors=NA

  cat('Computing Factor Scores on lavaan: Please wait                                                    \r')
  flush.console()

  if (raw_data==TRUE && method=="fixed" && corr=="Pearson"){
    pfactors = lavaan::lavPredict(fit, method = "regression" )
  }
  if (raw_data==TRUE && method=="fixed" && corr=="Polychoric"){
    pfactors = lavaan::lavPredict(fit, method = "EBM")
  }
  cat('Computing Factor Scores on lavaan: Done!                                                         \n')
  flush.console()
  cat("\r","                                                                                                  ","\r")


  fit_ind=lavaan::fitmeasures(fit)

  fit_indices<-list("gfi"=fit_ind["gfi"], "srmr"=fit_ind["srmr"], "rmsea"=fit_ind["rmsea"],"cfi"=fit_ind["cfi"])

  # OUTPUT

  results<-list("loadings"=loadings, "factor_cor"=psi, "fit_indices"=fit_indices, "AQ_variance"=AQ_var, "resid_matrix"=OUT$Res , "pfactors"=pfactors)


  if (display==TRUE){return(results)}
  else {invisible(results)}

}
