library(tidyverse)
library(gridExtra)

nFrames = 100

labels <- c('baseline S1 (well)   --->   S2 (mild)',
            'baseline S1 (well)   --->   S4 (dead)',
            'baseline S2 (mild)   --->   S3 (severe)',
            'baseline S2 (mild)   --->   S4 (dead)',
            'baseline S3 (severe)   --->   S4 (dead)',
            'time S1 (well)   --->   S2 (mild)',
            'time S1 (well)   --->   S4 (dead)',
            'time S2 (mild)   --->   S3 (severe)',
            'time S2 (mild)   --->   S4 (dead)',
            'time S3 (severe)   --->   S4 (dead)',
            'sex S1 (well)   --->   S2 (mild)',
            'sex S1 (well)   --->   S4 (dead)',
            'sex S2 (mild)   --->   S3 (severe)',
            'sex S2 (mild)   --->   S4 (dead)',
            'sex S3 (severe)   --->   S4 (dead)',
            'logit P( obs. state 2 | true state 1 )',
            'logit P( obs. state 1 | true state 2 )',
            'logit P( obs. state 3 | true state 2 )',
            'logit P( obs. state 2 | true state 3 )',
            'logit P( init. state 2 )','logit P( init. state 3 )')

par_index = list(beta=1:15, misclass=16:19, pi_logit=20:21)

# The true values are set as the posterior means of the thinned last 15,000 steps
# from running the MCMC routine using the numerical ODE solver (seed 10).
# The true values corresponding to the slope coefficient on time are scaled by
# 3 in order to magnify the effect of the different approaches to modelling 
# continuous time HMMs.
load('real_cav_analysis/Model_out/deSolve/mcmc_out_10.rda')
chain = mcmc_out$chain[10000:25001, ]
ind_keep = seq(1, nrow(chain), by=10)
chain = chain[ind_keep, ]
trueValues = colMeans(chain)
trueValues[6:10] = 3 * trueValues[6:10]

month_data = matrix(data=-1, nrow = nFrames, ncol = 21)
year_data = matrix(data=-1, nrow = nFrames, ncol = 21)
year_2_data = matrix(data=-1, nrow = nFrames, ncol = 21)

cred_set = vector(mode = "list", length = length(trueValues))

# For each parameter, there are 3 discretizations to consider
for(i in 1:length(cred_set)) {
    cred_set[[i]] = vector(mode = "list", length = 3)
    cred_set[[i]][[1]] = cred_set[[i]][[2]] = cred_set[[i]][[3]] =
            data.frame('lower' = c(-1), 'upper' = c(-1))
}

# -----------------------------------------------------------------------------
# Load Data and Find Credible Sets --------------------------------------------
# -----------------------------------------------------------------------------

row_ind1 = row_ind2 = row_ind3 = 1

for (i in 1:nFrames) {

  file_name = paste0("sim_cav_time_inhomog/Model_out/Month_msm/Output_msm", i, ".rda")
  if(file.exists(file_name)) {
    load(file_name)
    if (length(Output_msm$QmatricesSE) != 0) { # means a non poisitive definite matrix
        month_data[row_ind1,] = c(Output_msm$opt$par)

        # Calculating confidence intervals (Month)
        ind = 1
        for(l in 1:3) {
            for(r in c(5,13,10,14,15)){
                cred_set[[ind]][[1]][row_ind1,1] = Output_msm$Qmatrices[[l]][r] - 1.96*Output_msm$QmatricesSE[[l]][r]
                cred_set[[ind]][[1]][row_ind1,2] = Output_msm$Qmatrices[[l]][r] + 1.96*Output_msm$QmatricesSE[[l]][r]
                ind = ind + 1
            }
        }
        # Probabilities
        ind = 1
        lowerBounds = Output_msm$EmatricesL[[1]][c(5,2,10,7)]
        upperBounds = Output_msm$EmatricesU[[1]][c(5,2,10,7)]

        # Putting on the logit scale
        lowerBounds[1] = log(lowerBounds[1] / (1-lowerBounds[1]))
        lowerBounds[2:3] = log(lowerBounds[2:3] / (1-sum(lowerBounds[2:3])))
        lowerBounds[4] = log(lowerBounds[4] / (1-lowerBounds[4]))
        
        upperBounds[1] = log(upperBounds[1] / (1-upperBounds[1]))
        upperBounds[2:3] = log(upperBounds[2:3] / (1-sum(upperBounds[2:3])))
        upperBounds[4] = log(upperBounds[4] / (1-upperBounds[4]))

        for(l in par_index$misclass) {
            cred_set[[l]][[1]][row_ind1,1] = lowerBounds[ind]
            cred_set[[l]][[1]][row_ind1,2] = upperBounds[ind]
            ind = ind + 1
        }
        # Putting on the logit scale
        l_bound_init = c(Output_msm$ci[36,1], Output_msm$ci[37,1])
        u_bound_init = c(Output_msm$ci[36,2], Output_msm$ci[37,2])

        l_bound_init = log(l_bound_init / (1 - sum(l_bound_init)))
        u_bound_init = log(u_bound_init / (1 - sum(u_bound_init)))
        
        cred_set[[20]][[1]][row_ind1,] = c(l_bound_init[1], u_bound_init[1])
        cred_set[[21]][[1]][row_ind1,] = c(l_bound_init[2], u_bound_init[2])

        row_ind1 = row_ind1 + 1
    } else {
        print(paste0("Month Issue: ", i))
    }
  } else {
    print(paste0("Month File DNE: ", i))
  }

 # -----------------------------------------------------------------------------
  file_name = paste0("sim_cav_time_inhomog/Model_out/Year_msm/Output_msm", i, ".rda")
  if(file.exists(file_name)) {
    load(file_name)
    if (length(Output_msm$QmatricesSE) != 0) {
        year_data[row_ind2,] = c(Output_msm$opt$par)

        # Calculating confidence intervals (Year)
        ind = 1
        for(l in 1:3) {
            for(r in c(5,13,10,14,15)){
                cred_set[[ind]][[2]][row_ind2,1] = Output_msm$Qmatrices[[l]][r] - 1.96*Output_msm$QmatricesSE[[l]][r]
                cred_set[[ind]][[2]][row_ind2,2] = Output_msm$Qmatrices[[l]][r] + 1.96*Output_msm$QmatricesSE[[l]][r]
                ind = ind + 1
            }
        }
        # Probabilities
        ind = 1
        lowerBounds = Output_msm$EmatricesL[[1]][c(5,2,10,7)]
        upperBounds = Output_msm$EmatricesU[[1]][c(5,2,10,7)]
        
        # Putting on the logit scale
        lowerBounds[1] = log(lowerBounds[1] / (1-lowerBounds[1]))
        lowerBounds[2:3] = log(lowerBounds[2:3] / (1-sum(lowerBounds[2:3])))
        lowerBounds[4] = log(lowerBounds[4] / (1-lowerBounds[4]))

        upperBounds[1] = log(upperBounds[1] / (1-upperBounds[1]))
        upperBounds[2:3] = log(upperBounds[2:3] / (1-sum(upperBounds[2:3])))
        upperBounds[4] = log(upperBounds[4] / (1-upperBounds[4]))

        for(l in par_index$misclass) {
            cred_set[[l]][[2]][row_ind2,1] = lowerBounds[ind]
            cred_set[[l]][[2]][row_ind2,2] = upperBounds[ind]
            ind = ind + 1
        }
        
        # Putting on the logit scale
        l_bound_init = c(Output_msm$ci[36,1], Output_msm$ci[37,1])
        u_bound_init = c(Output_msm$ci[36,2], Output_msm$ci[37,2])

        l_bound_init = log(l_bound_init / (1 - sum(l_bound_init)))
        u_bound_init = log(u_bound_init / (1 - sum(u_bound_init)))
        
        cred_set[[20]][[2]][row_ind2,] = c(l_bound_init[1], u_bound_init[1])
        cred_set[[21]][[2]][row_ind2,] = c(l_bound_init[2], u_bound_init[2])

        row_ind2 = row_ind2 + 1
    } else {
        print(paste0("Year Issue: ", i))
    }
  } else {
    print(paste0("Year File DNE: ", i))
  }
 # -----------------------------------------------------------------------------

  file_name = paste0("sim_cav_time_inhomog/Model_out/YearTwo_msm/Output_msm", i, ".rda")
  if(file.exists(file_name)) {
    load(file_name)
    if (length(Output_msm$QmatricesSE) != 0) {
        year_2_data[row_ind3,] = c(Output_msm$opt$par)

        # Calculating confidence intervals (YearTwo)
        ind = 1
        for(l in 1:3) {
            for(r in c(5,13,10,14,15)){
                cred_set[[ind]][[3]][row_ind3,1] = Output_msm$Qmatrices[[l]][r] - 1.96*Output_msm$QmatricesSE[[l]][r]
                cred_set[[ind]][[3]][row_ind3,2] = Output_msm$Qmatrices[[l]][r] + 1.96*Output_msm$QmatricesSE[[l]][r]
                ind = ind + 1
            }
        }

        # Probabilities
        ind = 1
        lowerBounds = Output_msm$EmatricesL[[1]][c(5,2,10,7)]
        upperBounds = Output_msm$EmatricesU[[1]][c(5,2,10,7)]

        # Putting on the logit scale
        lowerBounds[1] = log(lowerBounds[1] / (1-lowerBounds[1]))
        lowerBounds[2:3] = log(lowerBounds[2:3] / (1-sum(lowerBounds[2:3])))
        lowerBounds[4] = log(lowerBounds[4] / (1-lowerBounds[4]))
        
        upperBounds[1] = log(upperBounds[1] / (1-upperBounds[1]))
        upperBounds[2:3] = log(upperBounds[2:3] / (1-sum(upperBounds[2:3])))
        upperBounds[4] = log(upperBounds[4] / (1-upperBounds[4]))

        for(l in par_index$misclass) {
            cred_set[[l]][[3]][row_ind3,1] = lowerBounds[ind]
            cred_set[[l]][[3]][row_ind3,2] = upperBounds[ind]
            ind = ind + 1
        }

        # Putting on the logit scale
        l_bound_init = c(Output_msm$ci[36,1], Output_msm$ci[37,1])
        u_bound_init = c(Output_msm$ci[36,2], Output_msm$ci[37,2])

        l_bound_init = log(l_bound_init / (1 - sum(l_bound_init)))
        u_bound_init = log(u_bound_init / (1 - sum(u_bound_init)))
        
        cred_set[[20]][[3]][row_ind3,] = c(l_bound_init[1], u_bound_init[1])
        cred_set[[21]][[3]][row_ind3,] = c(l_bound_init[2], u_bound_init[2])

        row_ind3 = row_ind3 + 1
    } else {
        print(paste0("YearTwo Issue: ", i))
    }
  } else {
    print(paste0("YearTwo File DNE: ", i))
  }
}

# -----------------------------------------------------------------------------
# Calculating Coverage --------------------------------------------------------
# -----------------------------------------------------------------------------

cov_df = matrix(ncol = 3, nrow = length(trueValues))
colnames(cov_df) = c("Month", "Year", "YearTwo")

for(i in 1:length(trueValues)) {
    val = trueValues[i]
    for(j in 1:ncol(cov_df)) {
        covrg = mean(cred_set[[i]][[j]]$lower <= val & val <= cred_set[[i]][[j]]$upper, na.rm=T)
        cov_df[i,j] = covrg
    }
    print(paste0("Coverage for ", round(val, 3), " is: "))
    cat('\t', '\t', '\t', "Month:   ", cov_df[i,1], '\n')
    cat('\t', '\t', '\t', "Year:    ", cov_df[i,2], '\n')
    cat('\t', '\t', '\t', "YearTwo: ", cov_df[i,3], '\n')
}

# ------------------------------------------------------------------------------

VP <- vector(mode="list", length = length(labels))
for(i in 1:length(trueValues)) {
	m = data.frame(y = month_data[,i], type = rep("Month", nrow(month_data)))
	y = data.frame(y = year_data[,i], type = rep("Year", nrow(year_data)))
	y2 = data.frame(y = year_2_data[,i], type = rep("Year_2", nrow(year_2_data)))

    plot_df = rbind(m,y,y2)
    xlabel = paste0("Month: ", round(cov_df[i,1], digits = 3),
                   ", Year: ", round(cov_df[i,2], digits = 3),
                   ", YearTwo: ", round(cov_df[i,3], digits = 3))

    VP[[i]] = ggplot(plot_df, aes(x = type, y = y)) +
      geom_violin(trim=FALSE) +
      geom_boxplot(width=0.1) +
      ggtitle(labels[i]) +
      ylab('') +
      xlab(xlabel) +
      geom_hline(yintercept=trueValues[i], linetype="dashed", color = "red") +
      theme(text = element_text(size = 7))

}

pdf("sim_cav_time_inhomog/Plots/msm.pdf", onefile = T)
grid.arrange(VP[[1]], VP[[2]], VP[[3]], VP[[4]], VP[[5]],
             VP[[6]], VP[[7]], VP[[8]], VP[[9]], ncol=3, nrow =3)
grid.arrange(VP[[10]], VP[[11]], VP[[12]], VP[[13]], VP[[14]],
             VP[[15]], VP[[16]], VP[[17]], VP[[18]], ncol=3, nrow =3)
grid.arrange(VP[[19]], VP[[20]], VP[[21]], ncol=3, nrow =3)
dev.off()

post_means = vector(mode = 'list', length = 3)
post_means[[1]] = month_data; post_means[[2]] = year_data; post_means[[3]] = year_2_data

save(post_means, file = paste0("sim_cav_time_inhomog/Plots/post_means_msm.rda"))
save(cov_df, file = paste0("sim_cav_time_inhomog/Plots/cov_df_msm.rda"))
