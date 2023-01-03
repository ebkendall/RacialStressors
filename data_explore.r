rsa_data = read.csv('Data/_FinalDataforRSASecondsStatesCovariates.csv', na.strings = "")

data_format = rsa_data[, c("ID..", "Time", "State", "RSA")]

na_patients = data_format[is.na(data_format$State), ]
na_patients = unique(na_patients$ID..)

data_format = data_format[!(data_format$ID.. %in% na_patients), ]

length(unique(data_format$ID..))

data_format$Time = data_format$Time / 100

data_format$State[data_format$State == "Baseline"] = 1
data_format$State[data_format$State == "Stress"] = 2
data_format$State[data_format$State == "Recovery"] = 3
data_format$State[data_format$State == "recovery"] = 3

data_format$State = as.numeric(data_format$State)
data_format$Time = as.numeric(data_format$Time)
data_format$RSA = as.numeric(data_format$RSA)

yes_s3 = rep(TRUE, length(unique(data_format$ID..)))
for(i in unique(data_format$ID..)) {
    sub = data_format[data_format$ID.. == i, ]
    if(sum(sub$State == 3) == 0){
        yes_s3[which(unique(data_format$ID..) == i)] = FALSE
    }
}

no_3 = unique(data_format$ID..)[!yes_s3]
data_format = data_format[!(data_format$ID.. %in% no_3), ]

save(data_format, file = 'Data/data_format.rda')

# -----------------------------------------------------------------------------
# Summary of data investigation -----------------------------------------------
# -----------------------------------------------------------------------------

means = function(df) {
    means_df = rep(0, length(unique(df$ID..)))
    for(i in 1:length(means_df)) {
        means_df[i] = mean(df[df$ID.. == unique(df$ID..)[i], 'RSA'])
    }
    return(means_df)
}
s1 = data_format[data_format$State == 1, ]
s2 = data_format[data_format$State == 2, ]
s3 = data_format[data_format$State == 3, ]


m1 = means(s1)[yes_s3]; mean(m1); var(m1)
m2 = means(s2)[yes_s3]; mean(m2); var(m2)
m3 = means(s3); mean(m3); var(m3)

cov_test = cbind(m1, cbind(m2,m3))
cov(cov_test)
colMeans(cov_test)

# hist(m2, col = "red")
# hist(m3, col = "green", add = T)
# hist(m1, add = T)
mean_RSA = matrix(0, nrow = length(unique(data_format$ID..)), ncol = 3)
sd_RSA = matrix(0, nrow = length(unique(data_format$ID..)), ncol = 3)

for(i in unique(data_format$ID..)) {
    sub_dat = data_format[data_format$ID.. == i, ]

    mean_RSA[which(unique(data_format$ID..) == i), 1] = mean(sub_dat$RSA[sub_dat$State == 1])
    mean_RSA[which(unique(data_format$ID..) == i), 2] = mean(sub_dat$RSA[sub_dat$State == 2])
    mean_RSA[which(unique(data_format$ID..) == i), 3] = mean(sub_dat$RSA[sub_dat$State == 3])
    
    sd_RSA[which(unique(data_format$ID..) == i), 1] = sd(sub_dat$RSA[sub_dat$State == 1])
    sd_RSA[which(unique(data_format$ID..) == i), 2] = sd(sub_dat$RSA[sub_dat$State == 2])
    sd_RSA[which(unique(data_format$ID..) == i), 3] = sd(sub_dat$RSA[sub_dat$State == 3])

}
save(mean_RSA, file = "Data/mean_RSA.rda")
colMeans(mean_RSA, na.rm = T)
# 
# temp = t(apply(mean_RSA, 1, diff))
# colMeans(temp, na.rm = T)
# 
# hist(mean_RSA[,2], breaks = sqrt(nrow(mean_RSA)), col = 2)
# hist(mean_RSA[,1], breaks = sqrt(nrow(mean_RSA)), col = 3, add = T)
# hist(mean_RSA[,3], breaks = sqrt(nrow(mean_RSA)), col = 4, add = T)

# -----------------------------------------------------------------------------
# Standard Deviation investigation --------------------------------------------
# -----------------------------------------------------------------------------

load('Data/data_format.rda')

s1_df = data_format[data_format$State == 1, ]
s2_df = data_format[data_format$State == 2, ]
s3_df = data_format[data_format$State == 3, ]

# Estimate of population standard deviation
sd_1 = sd(s1_df[,"RSA"])
sd_2 = sd(s2_df[,"RSA"])
sd_3 = sd(s3_df[,"RSA"])

# Individual standard deviations
sd_and_mean = matrix(ncol=9, nrow = length(unique(data_format$ID..)))
colnames(sd_and_mean) = c("m1", "m2", "m3", "sd1", "sd2", "sd3", "n1", "n2", "n3")
for(i in 1:length(unique(data_format$ID..))) {
    id = unique(data_format$ID..)[i]
    sd_and_mean[i,1] = mean(s1_df[s1_df$ID.. == id, "RSA"])
    sd_and_mean[i,2] = mean(s2_df[s2_df$ID.. == id, "RSA"])
    sd_and_mean[i,3] = mean(s3_df[s3_df$ID.. == id, "RSA"])
    
    sd_and_mean[i,4] = sd(s1_df[s1_df$ID.. == id, "RSA"])
    sd_and_mean[i,5] = sd(s2_df[s2_df$ID.. == id, "RSA"])
    sd_and_mean[i,6] = sd(s3_df[s3_df$ID.. == id, "RSA"])
    
    sd_and_mean[i,7] = sum(s1_df$ID.. == id)
    sd_and_mean[i,8] = sum(s2_df$ID.. == id)
    sd_and_mean[i,9] = sum(s3_df$ID.. == id)
}

# Sample standard deviation
sd(sd_and_mean[,"m1"]); sd_1 / sqrt(mean(sd_and_mean[,"n1"]))
sd(sd_and_mean[,"m2"]); sd_2 / sqrt(mean(sd_and_mean[,"n2"]))
sd(sd_and_mean[,"m3"]); sd_3 / sqrt(mean(sd_and_mean[,"n3"]))



