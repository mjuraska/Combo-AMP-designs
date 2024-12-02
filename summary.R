out1 <- readRDS("power_expoentialModel_iterated_2_utils1.rds")
out2 <- readRDS("power_expoentialModel_iterated_2_utils2.rds")
out3 <- readRDS("power_expoentialModel_iterated_2_utils3.rds")
out4 <- readRDS("power_expoentialModel_iterated_2_utils4.rds")

#small setting 90% two doses; both 100 and 150 gives the same target no endpoints
df1 <- out1[[1]]
df2 <- out3[[1]]

#small setting 85% two doses; 150 gives the smaller target no endpoints
df1 <- out1[[2]]
df2 <- out3[[2]]

#small setting 90% three doses
df1 <- out1[[3]]
df2 <- out3[[3]]
#small setting 85% three doses
df1 <- out1[[4]]
df2 <- out3[[4]]

#median setting 
df1 <- out1[[5]]
df2 <- out3[[5]]

#large setting 
df1 <- out1[[6]]
df2 <- out3[[6]]
