## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  ?acquihybrid

## ----echo=FALSE, results='asis'-----------------------------------------------
knitr::kable(cbind(paste0("responder", 1:10),head(siren::psymas, 10)), col.names = c("Responder","I1 ","I2 ","I3 ","I4 ","I5 ","I6 ","I7 ","I8 ","I9 ","I10"))

## ----echo=FALSE, results='asis'-----------------------------------------------
knitr::kable(cbind(c("I1","I2","I3","I4","I5","I6","I7","I8","I9","I10"),c(-9,-9,0,0,0,9,0,0,9,0),c(0,0,-9,9,-9,0,9,-9,0,9)), col.names = c("Item","F1","F2"))

## ----include = FALSE----------------------------------------------------------
library(siren)
psymas_target=cbind(c(-9,-9,0,0,0,9,0,0,9,0),c(0,0,-9,9,-9,0,9,-9,0,9))
results = acquihybrid(psymas,content_factors=2,target = psymas_target, corr = "Polychoric")
psymas_loadings = results$loadings
fit_indices = results$fit_indices
psymas_pfactors = results$pfactors

## ----echo=FALSE, results='asis'-----------------------------------------------
rownames(psymas_loadings) <- c("I1", "I2", "I6", "I9", "I3", "I4", "I5", "I7", "I8", "I10")
knitr::kable(psymas_loadings, col.names = c("Item", "F1", "F2", "ACQ"))

## ----echo=FALSE, results='asis'-----------------------------------------------
knitr::kable(cbind(paste0("responder", 1:10),head(round(psymas_pfactors,3), 10)), col.names = c("Responder","F1", "F2", "ACQ"))

