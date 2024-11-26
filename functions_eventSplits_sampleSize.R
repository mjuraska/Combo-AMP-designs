# compute total events for adequate power under a given design alternative
# (Schoenfeld, 1983, Biometrics)
((qnorm(0.975) + qnorm(0.9))^2) / ((1/4) * (log(0.7) - log(0.1))^2)

# get number of events in the numerator given
# total events 'n', hazard ratio 'hr', and probabilities 'p1' and 'p0' of
# assignment to Z=1 and Z=0, respectively;
# this is useful for obtaining event splits by arm
# for a given total number of events
n1 <- function(n, hr, p1=0.5, p0=0.5){
  p <- p1 / p0
  return(hr * p * n / (1 + p * hr))
}

# compute total target sample size 'N' for a 2-arm trial as a function of
# total events 'n', probability 'p1' of being assigned treatment Z=1,
# probability 'p0' of being assigned treatment Z=0, 
# and exponential rates (per PYR) for T|Z=1 ('rate1'),
# T|Z=0 ('rate0') and censoring time C ('rateC');
# 'tau' is the duration of follow-up (in years)
N <- function(n, p1=0.5, p0=0.5, rate1, rate0, rateC, tau=72/52){
  r1 <- rate1 / (rate1 + rateC)
  pEvent1 <- r1 - r1 * exp(-(rate1 + rateC) * tau)
  r0 <- rate0 / (rate0 + rateC)
  pEvent0 <- r0 - r0 * exp(-(rate0 + rateC) * tau)
  return(n / (p1 * pEvent1 + p0 * pEvent0))
}