out3 <- readRDS("power_expoentialModel_iterated_2_utils3.rds")

#small setting 90% two doses; both 100 and 150 gives the same target no endpoints
dftmp <- out3[[1]]

#small setting 85% two doses; 150 gives the smaller target no endpoints
dftmp <- out3[[2]]

#small setting 90% three doses
dftmp <- out3[[3]]
#small setting 85% three doses
dftmp <- out3[[4]]

#median setting 
dftmp <- out3[[5]]

#large setting 
dftmp <- out3[[6]]



#Sample size calculation based on 0.05 cutoff for one-sided tests
#PE = 90%
hr = 1 - 0.9
n = ceiling(((qnorm(0.95) + qnorm(0.9))^2) / ((1/4) * (log(0.7) - log(hr))^2))
n
n1 (n, hr, p1=0.5, p0=0.5)

#PE = 85%
hr = 1 - 0.85
n = ceiling(((qnorm(0.95) + qnorm(0.9))^2) / ((1/4) * (log(0.7) - log(hr))^2))
n
n1 (n, hr, p1=0.5, p0=0.5)

#PE = 75%
hr = 1 - 0.75
n = ceiling(((qnorm(0.95) + qnorm(0.9))^2) / ((1/4) * (log(0.7) - log(hr))^2))
n
n1 (n, hr, p1=0.5, p0=0.5)

