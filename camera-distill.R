## ----include=FALSE---------------------------------------------------------
knitr::opts_chunk$set(eval=TRUE, echo=TRUE, message=FALSE, warnings=FALSE)
library (kableExtra)
solution <- TRUE


## ---- readin, message=FALSE------------------------------------------------
trigger.events <- read.table(file="https://datadryad.org/stash/downloads/file_stream/73223",
                             header=TRUE)


## ----massage---------------------------------------------------------------
trigger.events$date <- paste("2014",
                       sprintf("%02i", trigger.events$month), 
                       sprintf("%02i", trigger.events$day),
                       sep="/")
trigger.events$time <- paste(sprintf("%02i", trigger.events$hour),
                       sprintf("%02i", trigger.events$minute),
                       sep=":")
trigger.events$datetime <- paste(trigger.events$date, trigger.events$time)












## ----DuikerCameraTraps-----------------------------------------------------
DuikerCameraTraps <- read.csv(file="https://datadryad.org/stash/downloads/file_stream/73221", 
                              header=TRUE, sep="\t")
DuikerCameraTraps$Area <- DuikerCameraTraps$Area / (1000*1000)
DuikerCameraTraps$object <- NA
DuikerCameraTraps$object[!is.na(DuikerCameraTraps$distance)] <- 1:sum(!is.na(DuikerCameraTraps$distance))


## ----smaltable-------------------------------------------------------------
sum(!is.na(DuikerCameraTraps$distance))
table(DuikerCameraTraps$Sample.Label)


## --------------------------------------------------------------------------
breakpoints <- c(seq(0,8,1), 10, 12, 15, 21)
hist(DuikerCameraTraps$distance, breaks=breakpoints, main="Peak activity data set",
     xlab="Radial distance (m)")


## ----fit-------------------------------------------------------------------
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


## ----pass1a, echo=FALSE----------------------------------------------------
knitr::kable(QAIC(uni1, uni2, uni3),
             caption="QAIC values for uniform key models.") %>%
  kable_paper(full_width = FALSE) %>%
  row_spec(3, bold=TRUE,  background = "#ff8c1a")


## ----pass1b, echo=FALSE----------------------------------------------------
knitr::kable(QAIC(hn0, hn1, hn2), 
             caption="QAIC values for half normal key models.") %>%
  kable_paper(full_width = FALSE) %>%
  row_spec(2, bold=TRUE,  background = "#ff8c1a")


## ----pass1c, echo=FALSE----------------------------------------------------
knitr::kable(QAIC(hr0, hr1),
             caption="QAIC values for hazard rate key models.") %>%
  kable_paper(full_width = FALSE) %>%
  row_spec(1, bold=TRUE,  background = "#ff8c1a")


## ----pass2-----------------------------------------------------------------
chats <- chi2_select(uni3, hn1, hr0)$criteria
modnames <- unlist(lapply(list(uni3, hn1, hr0), function(x) x$ddf$name.message))
results <- data.frame(modnames, chats)
results.sort <- results[order(results$chats),]
knitr::kable(results.sort, digits=2, row.names = FALSE,
             caption="Compare with Table S5 of Howe et al. (2018)") %>%
  kable_paper(full_width = FALSE) %>%
  row_spec(1, bold=TRUE,  background = "#4da6ff")
  


## --------------------------------------------------------------------------
p_a <- uni3$ddf$fitted[1]
w <- range(mybreaks)[2] - range(mybreaks)[1]
rho <- sqrt(p_a * w^2)


## ---- fig.height=5, caption="Detection function (left) and probability density function (right) of selected (half normal + 3 cosine adjustments) detection function model."----
par(mfrow=c(1,2))
plot(uni3, main="Daytime activity", xlab="Distance (m)",
     showpoints=FALSE, lwd=3, xlim=c(0, 15))
plot(uni3, main="Daytime activity", xlab="Distance (m)", pdf=TRUE,
     showpoints=FALSE, lwd=3, xlim=c(0, 15))


## ---- echo=FALSE-----------------------------------------------------------
par(mfrow=c(1,1))


## ---- sampfrac-------------------------------------------------------------
viewangle <- 42 # degrees
samfrac <- viewangle / 360
peak.uni.dens <- dht2(uni3, flatfile=DuikerCameraTraps, strat_formula = ~1,
                     sample_fraction = samfrac, er_est = "P2", multipliers = avail,
                     convert_units = conversion)
print(peak.uni.dens, report="density")


## ----mysummary-------------------------------------------------------------
mysummary <- function(ests, fit){
  return(data.frame(Label = ests$individuals$D$Label,
                    Dhat = ests$individuals$D$Estimate))
}

