## ----include=FALSE---------------------------------------------------------------------------------
knitr::opts_chunk$set(eval=TRUE, echo=TRUE, message=FALSE, warnings=FALSE)


## ---- readin, message=FALSE------------------------------------------------------------------------
library(Distance)
data("DuikerCameraTraps")


## ----smaltable, echo=FALSE-------------------------------------------------------------------------
sum(!is.na(DuikerCameraTraps$distance))
table(DuikerCameraTraps$Sample.Label)


## --------------------------------------------------------------------------------------------------
breakpoints <- c(seq(0,8,1), 10, 12, 15, 21)
hist(DuikerCameraTraps$distance, breaks=breakpoints, main="Peak activity data set",
     xlab="Radial distance (m)")


## ----fit-------------------------------------------------------------------------------------------
conversion <- convert_units("meter", NULL, "square kilometer")
trunc.list <- list(left=2, right=15)
peak.hn <- ds(DuikerCameraTraps, transect = "point", key="hn", adjustment="herm",
              cutpoints = c(seq(2,8,1), 10, 12, 15), truncation = trunc.list)
peak.unicos <- ds(DuikerCameraTraps, transect = "point", key="unif", adjustment = "cos",
                   cutpoints = c(seq(2,8,1), 10, 12, 15), truncation = trunc.list)
peak.hr <- ds(DuikerCameraTraps, transect = "point", key="hr", adjustment = "cos", 
              cutpoints = c(seq(2,8,1), 10, 12, 15), truncation = trunc.list)


## ---- sumtab---------------------------------------------------------------------------------------
knitr::kable(summarize_ds_models(peak.hn, peak.unicos, peak.hr), digits = 3, 
             caption="Model selection for three key functions fitted to duiker peak activity data")


## --------------------------------------------------------------------------------------------------
p_a <- peak.hr$ddf$fitted[1]
w <- 15
rho <- sqrt(p_a * w^2)


## --------------------------------------------------------------------------------------------------
par(mfrow=c(1,2))
plot(peak.hr, main="Peak activity", xlab="Distance (m)",
     showpoints=FALSE, lwd=3, xlim=c(0, 15))
plot(peak.hr, main="Peak activity", xlab="Distance (m)", pdf=TRUE,
     showpoints=FALSE, lwd=3, xlim=c(0, 15))


## ---- echo=FALSE-----------------------------------------------------------------------------------
par(mfrow=c(1,1))


## ---- sampfrac-------------------------------------------------------------------------------------
viewangle <- 42 # degrees
samfrac <- viewangle / 360
conversion <- convert_units("meter", NULL, "square kilometer")
peak.hr.dens <- dht2(peak.hr, flatfile=DuikerCameraTraps, strat_formula = ~1,
                     sample_fraction = samfrac, er_est = "P2", convert_units = conversion)
print(peak.hr.dens, report="density")


## ---- bootstrap, results='hide'--------------------------------------------------------------------
viewangle <- 42 # degrees
samfrac <- viewangle / 360
mysummary <- function(ests, fit){
  return(data.frame(Dhat = ests$individuals$D$Estimate))
}
duiker.boot.hr <- bootdht(model=peak.hr, flatfile=DuikerCameraTraps, resample_transects = TRUE,
                       nboot=250, summary_fun=mysummary, sample_fraction = samfrac,
                       convert.units = conversion)


## ----bootresult------------------------------------------------------------------------------------
print(summary(duiker.boot.hr))


## ---- fig.width=8, fig.cap="Distribution of density estimates from bootstrap replicates."----------
hist(duiker.boot.hr$Dhat[duiker.boot.hr$Dhat<100], breaks = 20, 
     xlab="Estimated density", main="D-hat estimates bootstraps")
abline(v=quantile(duiker.boot.hr$Dhat, probs = c(0.025,0.975)), lwd=2, lty=3)

