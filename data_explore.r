rsa_data = read.csv('Data/Old_data/_FinalDataforRSASecondsStatesCovariates.csv', na.strings = "")

data_format = rsa_data[, c("ID..", "Time", "State", "RSA")]

na_patients = data_format[is.na(data_format$State), ]
na_patients = unique(na_patients$ID..)

data_format = data_format[!(data_format$ID.. %in% na_patients), ]

length(unique(data_format$ID..))

data_format$Time = data_format$Time / 100

data_format$State[data_format$State == "Baseline"] = 1
data_format$State[data_format$State == "Stress"]   = 2
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

# One-way random effects ANOVA ----------------------------------------------
s1_anova = s1_df[,c("ID..", "RSA")]; s1_anova$ID.. = as.factor(s1_anova$ID..)
s2_anova = s2_df[,c("ID..", "RSA")]; s2_anova$ID.. = as.factor(s2_anova$ID..)
s3_anova = s3_df[,c("ID..", "RSA")]; s3_anova$ID.. = as.factor(s3_anova$ID..)

# https://stat.ethz.ch/~meier/teaching/anova/random-and-mixed-effects-models.html
library(lme4)
fit_1 <- lmer(RSA ~ (1 | ID..), data = s1_anova); summary(fit_1)
fit_2 <- lmer(RSA ~ (1 | ID..), data = s2_anova); summary(fit_2)
fit_3 <- lmer(RSA ~ (1 | ID..), data = s3_anova); summary(fit_3)

n_i = table(s1_df$ID..); n_i = t(n_i); colnames(n_i) = NULL
N = nrow(s1_anova)
n_0 = (1/90) * (N - (sum(n_i^2) / N))
F_hyp = qf(p = 0.95, df1 = 90, df2 = N - 91)

# STATE 1
sigma_i = 1.0365; sigma_res = 0.5693
F_stat_1 = (sigma_res + n_0 * sigma_i)/sigma_res
pf(F_stat_1, df1 = 90, df2 = N - 91, lower.tail = F)
# (F_stat_1 > F_hyp)

# STATE 2
sigma_i = 1.2908; sigma_res = 0.4845
F_stat_2 = (sigma_res + n_0 * sigma_i)/sigma_res
pf(F_stat_2, df1 = 90, df2 = N - 91, lower.tail = F)
# (F_stat_2 > F_hyp)

# STATE 3
sigma_i = 1.9359; sigma_res = 0.3202
F_stat_3 = (sigma_res + n_0 * sigma_i)/sigma_res
pf(F_stat_3, df1 = 90, df2 = N - 91, lower.tail = F)
# (F_stat_3 > F_hyp)

# Standard deviation derivations ----------------------------------------------
sd_1 = sd(s1_df[,"RSA"])
sd_2 = sd(s2_df[,"RSA"])
sd_3 = sd(s3_df[,"RSA"])

sd_and_mean = matrix(ncol=9, nrow = length(unique(data_format$ID..)))
colnames(sd_and_mean) = c("m1", "m2", "m3", 
                          "sd1", "sd2", "sd3", 
                          "n1", "n2", "n3")
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

sd(sd_and_mean[,"m1"]); sd(sd_and_mean[,"m2"]); sd(sd_and_mean[,"m3"])
mean(sd_and_mean[,"sd1"] / sqrt(sd_and_mean[,"n1"]))
mean(sd_and_mean[,"sd2"] / sqrt(sd_and_mean[,"n2"]))
mean(sd_and_mean[,"sd3"] / sqrt(sd_and_mean[,"n3"]))

m_1 = mean(s1_df[,"RSA"])
m_2 = mean(s2_df[,"RSA"])
m_3 = mean(s3_df[,"RSA"])

y_1 = dnorm(seq(0,12,by=0.1), mean = m_1, sd = sd_1)
y_2 = dnorm(seq(0,12,by=0.1), mean = m_2, sd = sd_2)
y_3 = dnorm(seq(0,12,by=0.1), mean = m_3, sd = sd_3)

# Plotting the Gaussian density
plot(x = seq(0,12,by=0.1), y = y_1, main = "Population Distribution", type = 'l', 
     xlab = "RSA", ylab = "", col = 'black', lwd = 2)
lines(x = seq(0,12,by=0.1), y = y_2, col = 'red', lwd=2)
lines(x = seq(0,12,by=0.1), y = y_3, col = 'green', lwd=2)
abline(v = m_1, col = 'black')
abline(v = m_2, col = 'red')
abline(v = m_3, col = 'green')
legend(10, 0.25, legend=c("Nominal", "Stress", "Recovery"), 
       fill = c("black","red", "green"))
    

# Sampling standard deviations
new_sd1 = sd_and_mean[,"sd1"] * sqrt(sd_and_mean[,"n1"])
new_sd2 = sd_and_mean[,"sd2"] * sqrt(sd_and_mean[,"n2"])
new_sd3 = sd_and_mean[,"sd3"] * sqrt(sd_and_mean[,"n3"])

colMeans(sd_and_mean[, 4:6])
mean(new_sd1); mean(new_sd2); mean(new_sd3)
sd_1; sd_2; sd_3

# Sample standard deviation
sd(sd_and_mean[,"m1"]); sd_1 / sqrt(mean(sd_and_mean[,"n1"]))
sd(sd_and_mean[,"m2"]); sd_2 / sqrt(mean(sd_and_mean[,"n2"]))
sd(sd_and_mean[,"m3"]); sd_3 / sqrt(mean(sd_and_mean[,"n3"]))

# Estimate of sampling distribution
sd_1 = sd(sd_and_mean[,"m1"])
sd_2 = sd(sd_and_mean[,"m2"])
sd_3 = sd(sd_and_mean[,"m3"])

m_1 = mean(sd_and_mean[,"m1"])
m_2 = mean(sd_and_mean[,"m2"])
m_3 = mean(sd_and_mean[,"m3"])

y_1 = dnorm(seq(0,12,by=0.1), mean = m_1, sd = sd_1)
y_2 = dnorm(seq(0,12,by=0.1), mean = m_2, sd = sd_2)
y_3 = dnorm(seq(0,12,by=0.1), mean = m_3, sd = sd_3)

# Plotting the Sampling distribution
plot(x = seq(0,12,by=0.1), y = y_1, main = "Sampling Distribution", type = 'l', 
     xlab = "RSA", ylab = "", col = 'black', lwd = 2)
lines(x = seq(0,12,by=0.1), y = y_2, col = 'red', lwd=2)
lines(x = seq(0,12,by=0.1), y = y_3, col = 'green', lwd=2)
abline(v = m_1, col = 'black')
abline(v = m_2, col = 'red')
abline(v = m_3, col = 'green')
legend(10, 0.25, legend=c("Nominal", "Stress", "Recovery"), 
       fill = c("black","red", "green"))


# Time series of RSA measurements with time indicators
pdf(paste0('Plots/patient_charts.pdf'))
par(mfrow=c(3, 2))

color_choice = c('blue', 'red', 'green')
for(i in 1:length(unique(data_format$ID..))) {
    sub_dat = data_format[data_format$ID.. == unique(data_format$ID..)[i], ]
    plot(sub_dat$Time, sub_dat$RSA, type = 'l', 
         xlab = "Time", ylab = "RSA", main = paste0('Patient: ', unique(data_format$ID..)[i]))
    
    # Determining state changes
    s = diff(sub_dat$State)
    abline(v = sub_dat$Time[1], col = color_choice[sub_dat$State[1]], lwd = 2)
    for(j in 1:sum(s)){
        t_switch = which(s==1)[j]+1
        abline(v = sub_dat$Time[t_switch], 
               col = color_choice[sub_dat$State[t_switch]],
               lwd = 2)
    }
}

dev.off()


load('Data/data_format.rda')

new_delta_est = matrix(nrow = length(unique(data_format$ID..)), ncol = 3)
for(i in 1:length(unique(data_format$ID..))) {
    sub_dat = data_format[data_format$ID.. == unique(data_format$ID..)[i], ]
    
    baseline = mean(sub_dat$RSA[sub_dat$State==1])
    
    s2_diff = mean(sub_dat$RSA[sub_dat$State==2] - baseline)
    s3_diff = mean(sub_dat$RSA[sub_dat$State==3] - baseline)
    new_delta_est[i,] = c(baseline, s2_diff, s3_diff)
}

colMeans(new_delta_est)

save(new_delta_est, file = "Data/new_delta_est.rda")

# Class Intro Presentation:
load('Data/data_format.rda')
sub_int = data_format[data_format$ID.. == 26746, ]

plot(sub_int$Time, sub_int$RSA, main = "Participant 26746 RSA Measurements (Before)",
     xlab = 'time', ylab = 'RSA')
abline(v = 1.69, col = "blue", lwd = 3)
abline(v = 4.69, col = "red", lwd = 3)
abline(v = 8.59, col = "lightgreen", lwd = 3)

b1 = mean(sub_int$RSA[1:10])
r1 = mean(sub_int$RSA[11:23])
g1 = mean(sub_int$RSA[24:27])

lines(sub_int$Time[1:11], rep(b1, length(sub_int$Time[1:11])), col = "blue", lwd = 3)
lines(sub_int$Time[11:24], rep(r1, length(sub_int$Time[11:24])), col = "red", lwd = 3)
lines(sub_int$Time[24:27], rep(g1, length(sub_int$Time[24:27])), col = "lightgreen", lwd = 3)
legend( 'topright', inset=c(0,-0.15), xpd=T, horiz=T, bty='n', x.intersp=.75,
        legend=c( 'Nominal', 'Stress', 'Recovery'), pch=15, pt.cex=1.5, 
        col=c( 'blue', 'red', 'lightgreen'))


plot(sub_int$Time, sub_int$RSA, main = "Participant 26746 RSA Measurements (After)",
     xlab = 'time', ylab = 'RSA')
abline(v = 1.69, col = "blue", lwd = 3)
abline(v = 4.69, col = "red", lwd = 3)
abline(v = 6.19, col = "lightgreen", lwd = 3)

b1 = mean(sub_int$RSA[1:10])
r1 = mean(sub_int$RSA[11:15])
g1 = mean(sub_int$RSA[16:27])

lines(sub_int$Time[1:11], rep(b1, length(sub_int$Time[1:11])), col = "blue", lwd = 3)
lines(sub_int$Time[11:16], rep(r1, length(sub_int$Time[11:16])), col = "red", lwd = 3)
lines(sub_int$Time[16:27], rep(g1, length(sub_int$Time[16:27])), col = "lightgreen", lwd = 3)
legend( 'topright', inset=c(0,-0.15), xpd=T, horiz=T, bty='n', x.intersp=.75,
        legend=c( 'Nominal', 'Stress', 'Recovery'), pch=15, pt.cex=1.5, 
        col=c( 'blue', 'red', 'lightgreen'))




# Sanity checking the current data with the older ones ------------------------
# First, compare the data sent to us with the old data
curr_data = read.csv('Data/_FinalDataforRSASecondsStatesCovariates.csv')
curr_id = unique(curr_data$ID..)

old_data = read.csv('Data/Study Dataset RSA Values Journal Adolescent Health Paper.csv')
old_id = unique(old_data$ID)

means_compare = matrix(data = c(curr_id, rep(NA, 125*9)), nrow = length(curr_id))
colnames(means_compare) = c('id', 'rsa_1_curr', 'rsa_1_old_no_epoch', 'rsa_1_old_30',
                            'rsa_2_curr', 'rsa_2_old_no_epoch', 'rsa_2_old_30',
                            'rsa_3_curr', 'rsa_3_old_no_epoch', 'rsa_3_old_30')
for(i in 1:length(curr_id)) {
    print(curr_id[i])
    if(curr_id[i] %in% old_id){
        sub_curr = curr_data[curr_data$ID.. == curr_id[i], ]
        sub_old = old_data[old_data$ID == curr_id[i], ]
        
        means_row = c(curr_id[i], rep(NA, 9))
        
        means_row[c(3,4,6,7,9,10)] = c(sub_old$Baseline_1_Output_no_epoch, sub_old$Baseline_1_Output_30,
                                       sub_old$Stress_1_Output_no_epoch, sub_old$Stress_1_Output_30,
                                       sub_old$Recovery_1_Output_no_epoch, sub_old$Recovery_1_Output_30)
        
        if(length(unique(sub_curr$State)) >=3) {
            curr_s1 = sub_curr$RSA[sub_curr$State == "Baseline"]
            curr_s2 = sub_curr$RSA[sub_curr$State == "Stress"]
            curr_s3 = sub_curr$RSA[sub_curr$State == "Recovery" | sub_curr$State == "recovery"]
            
            means_row[c(2,5,8)] = c(mean(curr_s1), mean(curr_s2), mean(curr_s3))
        }
        
        means_compare[i, ] = means_row
    } else {
        print('not in old dataset')
    }
}

# simplifying the means_compare
means_compare = means_compare[!is.na(means_compare[,2]), ]
means_compare = means_compare[!is.na(means_compare[,6]), ]
temp = colMeans(means_compare)
names(temp) = NULL
temp[c(2,5,8)]
temp[c(3,6,9)]
temp[c(4,7,10)]


# Second, compare the formatted data with the old data
load('Data/data_format.rda')
mean(data_format$RSA[data_format$State == 1])
mean(data_format$RSA[data_format$State == 2])
mean(data_format$RSA[data_format$State == 3])

