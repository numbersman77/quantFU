## ---- include=FALSE, echo=FALSE-----------------------------------------------
# --------------------------------------------------------------
# generate R file with code from this file
# --------------------------------------------------------------
knitr::purl(input = "quantFU.Rmd", output = "quantFU.r")


## ---- include=TRUE, echo=TRUE-------------------------------------------------
# --------------------------------------------------------------
# packages
# --------------------------------------------------------------
packs <- c("survival", "rpact")    
for (i in 1:length(packs)){library(packs[i], character.only = TRUE)}


## ---- include=TRUE, echo=TRUE-------------------------------------------------
# --------------------------------------------------------------
# functions
# --------------------------------------------------------------

# function to compute all the different variants of follow-up quantification
quantifyFU <- function(rando, event_time, event_type, ccod){
  
  ## =================================================================
  ## compute all FU measures
  ## =================================================================

  # input arguments:
  # - rando: date of randomization
  # - event_time: time to event or time to censoring
  # - event_type: 0 = event, 1 = lost to FU, 2 = administratively censored
  # - ccod: clinical cutoff date

  require(survival)

  n <- length(rando)

  # define empirical survival function
  ecdf <- (n - 1):0 / n
  
  # objects to collect distributions
  res <- NULL
  
  # indicator for lost to follow up
  ltfu_cens <- as.numeric(event_type == 1)

  # indicator for administratively censored
  admin_cens <- as.numeric(event_type == 2)

  # usual censoring indicator: 0 = censored (for whatever reason), 1 = event
  primary_event <- as.numeric(event_type == 0)

  # indicator for censoring
  event_time_cens <- as.numeric(primary_event == 0)

  # observation time regardless of censoring
  res[[1]] <- cbind(sort(event_time), ecdf)
  m1 <- median(event_time)

  # observation time for those event-free
  d2 <- event_time[event_time_cens == 1]
  n2 <- length(d2)
  res[[2]] <- cbind(sort(d2), (n2 - 1):0 / n2)
  m2 <- median(d2)
  
  # time to censoring
  so3 <- survfit(Surv(event_time, event_time_cens) ~ 1)
  res[[3]] <- so3
  m3 <- summary(so3)$table["median"]

  # time to CCOD, potential followup
  pfu <- as.numeric((ccod - rando) / 365.25 * 12)
  res[[4]] <- cbind(sort(pfu), ecdf)
  m4 <- median(pfu)

  # known function time
  m5 <- rep(NA, n)
  m5[event_time_cens == 1] <- event_time[event_time_cens == 1]
  m5[primary_event == 1] <- pfu[primary_event == 1]
  res[[5]] <- cbind(sort(m5), ecdf)
  m5 <- median(m5)

  # Korn's potential follow-up time
  spfu <- sort(pfu)
  pt <- rep(NA, n)
  for (i in 1:n){
  
    # timepoint at which we compute potential followup (t' in Schemper et al)
    t <- spfu[i]
  
    # proportion with pfu > t
    pet <- mean(spfu > t)
  
    # time to LTFU
    so6 <- survfit(Surv(event_time, ltfu_cens) ~ 1, subset = (spfu >= t))
    pltet <- ifelse(max(so6$time > t) == 0, 0, so6$surv[so6$time > t][1])
  
    pt[i] <- pltet * pet
  }

  res[[6]] <- cbind(spfu, pt)
  
  # now take the median of the distribution, see plot(t, pt, type = "s")
  m6 <- max(spfu[pt >= 0.5])

  # Modified potential follow-up time
  m7 <- rep(NA, n)
  m7[ltfu_cens == 1] <- pmin(pfu, event_time)[ltfu_cens == 1]
  m7[ltfu_cens == 0] <- pfu[ltfu_cens == 0]
  res[[7]] <- cbind(sort(m7), ecdf)
  m7 <- median(m7)

  # Potential follow-up considering events
  m8 <- rep(NA, n)
  m8[event_time_cens == 1] <- pfu[event_time_cens == 1]
  m8[primary_event == 1] <- event_time[primary_event == 1]
  res[[8]] <- cbind(sort(m8), ecdf)
  m8 <- median(m8)

  # summarize results for medians
  dat <- matrix(NA, nrow = 8, ncol = 1)
  dat[, 1] <- c(m1, m2, m3, m4, m5, m6, m7, m8)
  rownames(dat) <- c("Observation time regardless of censoring", 
                     "Observation time for those event-free", 
                     "Time to censoring", "Time to CCOD", 
                     "Known function time", 
                     "Korn potential follow-up", 
                     "Modified potential follow-up", 
                     "Potential follow-up considering events")
  colnames(dat) <- "median"
  
  output <- list("medians" = dat, "distributions" = res)
  class(output) <- "qfu"
  return(output)
}

# function to plot distributions of various quantifiers
plot.qfu <- function(x, which = 1:8, median = TRUE){
  
  dat <- x$distributions
  
  par(las = 1)
  plot(0, 0, type = "n", xlim = c(0, 1.05 * max(dat[[1]][, 1])), ylim = c(0, 1), 
       xlab = "time-to-event endpoint", ylab = "probability", 
       main = "Distributions used for quantification of follow-up")
  if (3 %in% which){lines(dat[[3]], col = 3, conf.int = FALSE)}
  which <- which[which != 3]
  for (i in which){lines(dat[[i]], col = i)}
  
  if(isTRUE(median)){abline(h = 0.5, col = grey(0.5), lwd = 2)}
  
  legend("bottomleft", rownames(x$medians)[which], lty = 1, col = (1:8)[which], bty = "n", lwd = 2)
}


## ---- include=TRUE, echo=TRUE-------------------------------------------------
# simulate a clinical trial using rpact
# time unit is months
design <- getDesignGroupSequential(informationRates = 1,
   typeOfDesign = "asOF", sided = 1, alpha = 0.025, beta = 0.2)

simulationResult <- getSimulationSurvival(design,
    lambda2 = log(2) / 60, hazardRatio = 0.75,
    dropoutRate1 = 0.025, dropoutRate2 = 0.025, 
    dropoutTime = 12,
    accrualTime = 0:6, 
    accrualIntensity = 6 * 1:7,
    plannedEvents = 350,
    directionUpper = FALSE,
    maxNumberOfSubjects = 1000,
    maxNumberOfIterations = 1,
    maxNumberOfRawDatasetsPerStage = 1,
    seed = 2021)

# retrieve dataset
simdat <- getRawData(simulationResult)

# create variable with randomization dates
# note that addition / subtraction of date objects happens in days 
# --> multiplication by 365.25 / 12 ~= 30 below
day0 <- as.Date("2016-01-31", origin = "1899-12-30")
rando <- day0 + simdat$accrualTime * 365.25 / 12
pfs <- simdat$timeUnderObservation

# event type variable: 0 = event, 1 = lost to FU, 2 = administratively censored
event_type <- rep(NA, nrow(simdat))
event_type[simdat$event == TRUE] <- 0
event_type[simdat$event == FALSE & simdat$dropoutEvent == TRUE] <- 1
event_type[simdat$event == FALSE & simdat$dropoutEvent == FALSE] <- 2 

# PFS event
pfsevent <- as.numeric(event_type == 0)

# treatment arm
arm <- factor(simdat$treatmentGroup)

# define clinical cutoff date based on simulation result
ccod <- day0 + simdat$observationTime[1] * 365.25 / 12
  
# check
so1 <- summary(coxph(Surv(pfs, pfsevent) ~ arm))
so1


## ---- echo = TRUE, results = 'asis', message = FALSE, fig.cap = "", fig.align = "center", fig.width = 7, fig.height = 5.5----
par(las = 1)
plot(survfit(Surv(pfs, pfsevent) ~ arm), col = 2:3, mark = "'", lty = 1, xlim = c(0, 100), 
     xlab = "PFS", ylab = "probability of being event-free")


## ---- include=TRUE, echo=TRUE-------------------------------------------------
fu <- quantifyFU(rando = rando, event_time = pfs, event_type = event_type, ccod = ccod)

# medians of all these distributions:
fu$medians


## ---- echo = TRUE, results = 'asis', message = FALSE, fig.cap = "", fig.align = "center", fig.width = 7, fig.height = 5.5----
plot(fu)

