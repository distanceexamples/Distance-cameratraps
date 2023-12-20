## ----include=FALSE-------------------------------------------------------
knitr::opts_chunk$set(eval=TRUE, echo=TRUE, message=FALSE, warnings=FALSE)
library (kableExtra)
solution <- TRUE


## ---- readin, message=FALSE----------------------------------------------
trigger.events <- read.table(file="https://datadryad.org/stash/downloads/file_stream/73223",
                             header=TRUE)


## ----massage-------------------------------------------------------------
trigger.events$date <- paste("2014",
                       sprintf("%02i", trigger.events$month), 
                       sprintf("%02i", trigger.events$day),
                       sep="/")
trigger.events$time <- paste(sprintf("%02i", trigger.events$hour),
                       sprintf("%02i", trigger.events$minute),
                       sep=":")
trigger.events$datetime <- paste(trigger.events$date, trigger.events$time)


## ----radian, eval=solution-----------------------------------------------
## library(activity)
## trigger.events$rtime <- gettime(trigger.events$datetime,
##                                 tryFormats = "%Y/%m/%d %H:%M",
##                                 scale = "radian")


## ----activity, eval=solution---------------------------------------------
## act_result <- fitact(trigger.events$rtime, sample="data", reps=100)


## ----actplot, eval=solution, fig.cap="Fitted smooth to histogram of camera triggering times for Maxwell's duiker data."----
## plot(act_result)


## ----thenumber, eval=solution--------------------------------------------
## print(act_result@act)


## ----avmultiplier, eval=solution-----------------------------------------
## camera.operation.per.day <- 11.5
## prop.camera.time <- camera.operation.per.day / 24
## avail <- list(creation=data.frame(rate = act_result@act[1]/prop.camera.time,
##                                   SE   = act_result@act[2]/prop.camera.time))


## ----DuikerCameraTraps---------------------------------------------------
DuikerCameraTraps <- read.csv(file="https://datadryad.org/stash/downloads/file_stream/73221", 
                              header=TRUE, sep="\t")
DuikerCameraTraps$Area <- DuikerCameraTraps$Area / (1000*1000)
DuikerCameraTraps$object <- NA
DuikerCameraTraps$object[!is.na(DuikerCameraTraps$distance)] <- 1:sum(!is.na(DuikerCameraTraps$distance))


## ----smaltable-----------------------------------------------------------
sum(!is.na(DuikerCameraTraps$distance))
table(DuikerCameraTraps$Sample.Label)


## ------------------------------------------------------------------------
breakpoints <- c(seq(0,8,1), 10, 12, 15, 21)
hist(DuikerCameraTraps$distance, breaks=breakpoints, main="Peak activity data set",
     xlab="Radial distance (m)")


## ----fit-----------------------------------------------------------------
library(Distance)
trunc.list <- list(left=2, right=15)
mybreaks <- c(seq(2,8,1), 10, 12, 15)
conversion <- convert_units("meter", NULL, "square kilometer")
uni1 <- ds(DuikerCameraTraps, transect = "point", key="unif", adjustment = "cos",
           nadj=1, convert_units = conversion,
           cutpoints = mybreaks, truncation = trunc.list)
uni2 <- ds(DuikerCameraTraps, transect = "point", key="unif", adjustment = "cos",
           nadj=2, convert_units = conversion,
           cutpoints = mybreaks, truncation = trunc.list)
uni3 <- ds(DuikerCameraTraps, transect = "point", key="unif", adjustment = "cos",
           nadj=3, convert_units = conversion,
           cutpoints = mybreaks, truncation = trunc.list)

hn0 <- ds(DuikerCameraTraps, transect = "point", key="hn", adjustment = NULL,
          convert_units = conversion, cutpoints = mybreaks, truncation = trunc.list)
hn1 <- ds(DuikerCameraTraps, transect = "point", key="hn", adjustment = "cos",
          nadj=1, convert_units = conversion,
          cutpoints = mybreaks, truncation = trunc.list)
hn2 <- ds(DuikerCameraTraps, transect = "point", key="hn", adjustment = "cos",
          nadj=2, convert_units = conversion,
          cutpoints = mybreaks, truncation = trunc.list)

hr0 <- ds(DuikerCameraTraps, transect = "point", key="hr", adjustment = NULL,
          convert_units = conversion, cutpoints = mybreaks, truncation = trunc.list)
hr1 <- ds(DuikerCameraTraps, transect = "point", key="hr", adjustment = "poly",
          nadj=1, convert_units = conversion,
          cutpoints = mybreaks, truncation = trunc.list)


## ----pass1a, echo=FALSE--------------------------------------------------
knitr::kable(QAIC(uni1, uni2, uni3),
             caption="QAIC values for uniform key models.") %>%
  kable_paper(full_width = FALSE) %>%
  row_spec(3, bold=TRUE,  background = "#ff8c1a")


## ----pass1b, echo=FALSE--------------------------------------------------
knitr::kable(QAIC(hn0, hn1, hn2), 
             caption="QAIC values for half normal key models.") %>%
  kable_paper(full_width = FALSE) %>%
  row_spec(2, bold=TRUE,  background = "#ff8c1a")


## ----pass1c, echo=FALSE--------------------------------------------------
knitr::kable(QAIC(hr0, hr1),
             caption="QAIC values for hazard rate key models.") %>%
  kable_paper(full_width = FALSE) %>%
  row_spec(1, bold=TRUE,  background = "#ff8c1a")


## ----pass2---------------------------------------------------------------
chats <- chi2_select(uni3, hn1, hr0)$criteria
modnames <- unlist(lapply(list(uni3, hn1, hr0), function(x) x$ddf$name.message))
results <- data.frame(modnames, chats)
results.sort <- results[order(results$chats),]
knitr::kable(results.sort, digits=2, row.names = FALSE,
             caption="Compare with Table S5 of Howe et al. (2018)") %>%
  kable_paper(full_width = FALSE) %>%
  row_spec(1, bold=TRUE,  background = "#4da6ff")
  


## ------------------------------------------------------------------------
p_a <- uni3$ddf$fitted[1]
w <- range(mybreaks)[2] - range(mybreaks)[1]
rho <- sqrt(p_a * w^2)


## ---- fig.height=5, caption="Detection function (left) and probability density function (right) of selected (half normal + 3 cosine adjustments) detection function model."----
par(mfrow=c(1,2))
plot(uni3, main="Daytime activity", xlab="Distance (m)",
     showpoints=FALSE, lwd=3, xlim=c(0, 15))
plot(uni3, main="Daytime activity", xlab="Distance (m)", pdf=TRUE,
     showpoints=FALSE, lwd=3, xlim=c(0, 15))


## ---- echo=FALSE---------------------------------------------------------
par(mfrow=c(1,1))


## ---- sampfrac-----------------------------------------------------------
viewangle <- 42 # degrees
samfrac <- viewangle / 360
peak.uni.dens <- dht2(uni3, flatfile=DuikerCameraTraps, strat_formula = ~1,
                     sample_fraction = samfrac, er_est = "P2", multipliers = avail,
                     convert_units = conversion)
print(peak.uni.dens, report="density")


## ----mysummary-----------------------------------------------------------
mysummary <- function(ests, fit){
  return(data.frame(Label = ests$individuals$D$Label,
                    Dhat = ests$individuals$D$Estimate))
}


## ----multifunc, eval=solution--------------------------------------------
## mult <- list(availability= make_activity_fn(trigger.events$rtime, sample="data",
##                                             detector_daily_duration=camera.operation.per.day))


## ---- eval = FALSE-------------------------------------------------------
## download.file("http://distancesampling.org/R/MCDS.exe",
##               paste0(system.file(package="mrds"),"/MCDS.exe"), mode = "wb")
## #Detach and reload the Distance package to make use of it
## detach("package:Distance", unload = TRUE)
## library(Distance)


## ---- mcds.exe, eval = solution------------------------------------------
## uni3 <- ds(DuikerCameraTraps, transect = "point", key="unif", adjustment = "cos",
##            nadj=3, convert_units = conversion,
##            cutpoints = mybreaks, truncation = trunc.list,
##            optimizer = "MCDS")


## ---- bootstrap, results='hide', eval=solution---------------------------
## n.cores <- parallel::detectCores()
## daytime.boot.uni <- bootdht(model=uni3, flatfile=DuikerCameraTraps,
##                           resample_transects = TRUE, nboot=500,
##                           cores = n.cores - 1,
##                           summary_fun=mysummary, sample_fraction = samfrac,
##                           convert_units = conversion, multipliers=mult)


## ----bootresult, eval=solution-------------------------------------------
## print(summary(daytime.boot.uni))


## ---- fig.width=8, fig.cap="Distribution of density estimates from bootstrap replicates.", eval=solution----
## hist(daytime.boot.uni$Dhat, breaks = 20,
##      xlab="Estimated density", main="D-hat estimates bootstraps")
## abline(v=quantile(daytime.boot.uni$Dhat, probs = c(0.025,0.975), na.rm=TRUE), lwd=2, lty=3)


## ---- startvals, eval = solution-----------------------------------------
## uni3.with.startvals <- ds(DuikerCameraTraps, transect = "point", key="unif", adjustment = "cos",
##            nadj=3,
##            cutpoints = mybreaks, truncation = trunc.list,
##            initial_values = list(adjustment = c(as.numeric(uni2$ddf$par), 0)))


## ---- startvals2, eval = solution----------------------------------------
## print(uni2$ddf$par)


## ---- startvals3, eval = solution----------------------------------------
## uni3.with.startvals <- ds(DuikerCameraTraps, transect = "point", key="unif", adjustment = "cos",
##            nadj=3,
##            cutpoints = mybreaks, truncation = trunc.list,
##            initial_values = list(adjustment = c(0.97178178, 0.03541633, 0)))


## ---- startvals4, eval = FALSE-------------------------------------------
## print(uni3$ddf$par)
## uni3.with.startvals <- ds(DuikerCameraTraps, transect = "point", key="unif", adjustment = "cos",
##            nadj=3,
##            cutpoints = mybreaks, truncation = trunc.list,
##            optimizer = "MCDS",
##            initial_values = list(adjustment = c(0.93518220, -0.05345965, -0.08073799)))
## daytime.boot.uni <- bootdht(model=uni3.with.startvals, flatfile=DuikerCameraTraps,
##                           resample_transects = TRUE, nboot=500,
##                           cores = n.cores - 1,
##                           summary_fun=mysummary, sample_fraction = samfrac,
##                           convert_units = conversion, multipliers=mult)

