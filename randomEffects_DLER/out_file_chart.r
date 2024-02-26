# library(matrixStats)
library(plotrix)

# Model type -------------------------------------------------------------------
# 1: baseline only
# 2: baseline & DLER
# 3: all covariates

covariate_struct = 3
# ------------------------------------------------------------------------------

# Information defining which approach to take ----------------------------------
simulation = F
case_b = T
interm = F
it_num = 4
if(simulation) {
    trial_num = covariate_struct
    index_seeds = c(1:5)
} else {
    trial_num = covariate_struct + 6
    index_seeds = c(1:5)
}
# ------------------------------------------------------------------------------

n_post = 20000; burnin = 0; steps = 200000
# Matrix row indices for the posterior sample to use
index_post = (steps - burnin - n_post + 1):(steps - burnin)

Dir = 'Model_out/'

if(covariate_struct == 1) {
    par_index = list(zeta=1:24, misclass=38:41, delta = 25:27, tau2 = 28, 
                     sigma2 = 29:31, gamma = 32:34, delta_new = 35:37)
} else if(covariate_struct) {
    par_index = list(zeta=1:12, misclass=26:29, delta = 13:15, tau2 = 16, 
                     sigma2 = 17:19, gamma = 20:22, delta_new = 23:25)
} else {
    par_index = list(zeta=1:30, misclass=44:47, delta = 31:33, tau2 = 34, 
                     sigma2 = 35:37, gamma = 38:40, delta_new = 41:43)
}

for(seed in index_seeds) {
    file_name = NULL
    if(simulation) {
        if(case_b) {
            if(interm) {
                file_name = paste0(Dir,'mcmc_out_interm_',toString(seed),'_', trial_num, 'it', it_num, '_sim.rda')   
            } else {
                file_name = paste0(Dir,'mcmc_out_',toString(seed), '_', trial_num, '_sim_30b.rda')   
            }
        } else {
            file_name = paste0(Dir,'mcmc_out_',toString(seed), '_', trial_num, '_sim_30.rda')      
        }
    } else {
        if(case_b) {
            if(interm) {
                file_name = paste0(Dir,'mcmc_out_interm_',toString(seed),'_', trial_num, 'it', it_num, '.rda')   
            } else {
                file_name = paste0(Dir,'mcmc_out_',toString(seed), '_', trial_num, '_30b.rda')   
            }   
        } else {
            file_name = paste0(Dir,'mcmc_out_',toString(seed), '_', trial_num, '_30.rda')      
        }
    }
    
    print(file_name)
    
    load(file_name)
    b_chain_ind = 1:nrow(mcmc_out$B_chain)
    ind_B_chain = seq(1, length(b_chain_ind), by=10)
    
    if(seed == 1) {
        post_B_chain = mcmc_out$B_chain[ind_B_chain, ]
    } else {
        post_B_chain = rbind(post_B_chain, mcmc_out$B_chain[ind_B_chain, ])
    }
}

print("Dimension of posterior state sequence")
print(dim(post_B_chain))


if(simulation) {
    # Simulation
    load(paste0('Data/sim_data_', covariate_struct, '_30.rda'))
    data_format = as.data.frame(sim_data)
    
    EIDs = unique(data_format[,"ID.."])
} else {
    # Real data analysis
    load('../Data/data_format_30.rda')
    data_format = data_format_30
    
    miss_info = c(26296, 29698, 30625, 401, 423, 419, 457)
    data_format = data_format[!(data_format[,"ID.."] %in% miss_info), ]
    
    EIDs = unique(data_format[,"ID.."])
}

# New patients ---------------------------------------------------------------
if(simulation) {
    if(case_b) {
        if(interm) {
            pdf_title = paste0('Plots/chart_plot_', trial_num, '_sim_30b_it', it_num,'.pdf')
        } else {
            pdf_title = paste0('Plots/chart_plot_', trial_num, '_sim_30b', '.pdf')
        }
    } else {
        pdf_title = paste0('Plots/chart_plot_', trial_num, '_sim_30', '.pdf')
    }
} else {
    if(case_b) {
        if(interm) {
            pdf_title = paste0('Plots/chart_plot_', trial_num, '_30b_s1_it', it_num, '.pdf')
        } else {
            pdf_title = paste0('Plots/chart_plot_', trial_num, '_30b_s1', '.pdf')
        }
    } else {
        pdf_title = paste0('Plots/chart_plot_', trial_num, '_30', '.pdf')
    }
}

if(interm) {
    panels = c(4, 1)
} else {
    panels = c(3, 1)
    if(simulation) {
        load(paste0("Model_out/B_MLE_", trial_num, "_sim.rda"))
    } else {
        load(paste0("Model_out/B_MLE_", trial_num, ".rda"))   
    }
}


pdf(pdf_title)
par(mfrow=panels, mar=c(2,4,2,4))#, bg='black', fg='green')
for(i in EIDs){
	print(which(EIDs == i))
	indices_i = (data_format[,'ID..']==i)
	sub_dat = data_format[data_format[,"ID.."] == i, ]
	n_i = sum(indices_i)

    if(covariate_struct == 1) {
        plot_title = paste0('Participant: ', i)
    } else if(covariate_struct == 2) {
        
        dler_val = NULL
        for(a in EIDs) {
            dler_val = c(dler_val, unique(data_format[data_format[,"ID.."] == a, "DLER_avg"]))
        }
        mean_dler = mean(dler_val)
        
        cov_value = c(sub_dat[1,"DLER_avg"])
        cov_value = as.numeric(cov_value)
        cov_value[1] = cov_value[1] - mean_dler
        cov_value[1] = round(cov_value[1], digits = 3)

        plot_title = paste0('Participant: ', i, ', DLER: ', cov_value[1])
    } else {
        
        ages = NULL
        dler_val = NULL
        for(a in EIDs) {
            ages = c(ages, unique(data_format[data_format[,"ID.."] == a, "Age"]))
            dler_val = c(dler_val, unique(data_format[data_format[,"ID.."] == a, "DLER_avg"]))
        }
        mean_age = mean(ages)
        mean_dler = mean(dler_val)
        
        cov_value = c(sub_dat[1,c("Age", "sex1", "edu_yes", "DLER_avg")])
        cov_value = as.numeric(cov_value)
        cov_value[1] = cov_value[1] - mean_age
        cov_value[4] = cov_value[4] - mean_dler
        cov_value[1] = round(cov_value[1], digits = 3)
        cov_value[4] = round(cov_value[4], digits = 3)

        plot_title = paste0('Participant: ', i, ', sex: ', cov_value[2], 
                          ', pEdu: ', cov_value[3], ', DLER: ', cov_value[4],
                          ', age: ', cov_value[1])
    }

	t_grid = t_grid_bar = 1:n_i
	main_color = 'black'
	
	if(!simulation) {
	    x_mean_1 = c(min(which(sub_dat$State == 1)), max(which(sub_dat$State == 1))+1)
	    x_mean_2 = c(min(which(sub_dat$State == 2)), max(which(sub_dat$State == 2))+1)
	    x_mean_3 = c(min(which(sub_dat$State == 3)), max(which(sub_dat$State == 3)))
	    
	    y_mean_1 = mean(sub_dat$RSA[sub_dat$State == 1])
	    y_mean_2 = mean(sub_dat$RSA[sub_dat$State == 2])
	    y_mean_3 = mean(sub_dat$RSA[sub_dat$State == 3])   
	}
	
	b_i = as.numeric(data_format[ indices_i,"State"])
	to_s1 = c(1,(2:n_i)[diff(b_i)!=0 & b_i[-1]==1])
	to_s2 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==2]
	to_s3 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==3]

    # Plot the "ghost plot" to make everything line up
    pb = barplot(rbind( colMeans(post_B_chain[, indices_i] == 1),
			            colMeans(post_B_chain[, indices_i] == 2),
				        colMeans(post_B_chain[, indices_i] == 3)), 
            col=c( 'dodgerblue', 'firebrick1', 'yellow2'),
			xlab='time', space=0, col.main=main_color, border=NA, axes = F, plot = F) 

	plot(x=pb, y=data_format[indices_i, "RSA"], 
            xlab='time', ylab = 'RSA', col.main=main_color, 
            main = plot_title, xlim = range(pb) + c(-0.5,0.5),
	        xaxt='n', yaxt='n', col.lab = main_color)
	axis( side=1, at=pb, col.axis=main_color, labels=t_grid)
	axis( side=2, at=seq(min(data_format[indices_i, "RSA"]), 
                         max(data_format[indices_i, "RSA"])), col.axis=main_color)
	
	if(!simulation) {
	    segments(x0 = x_mean_1[1]-0.5, x1 = x_mean_1[2]-0.5, y0 = y_mean_1, y1 = y_mean_1, col = 'azure2', lwd = 3)
	    segments(x0 = x_mean_2[1]-0.5, x1 = x_mean_2[2]-0.5, y0 = y_mean_2, y1 = y_mean_2, col = 'azure2', lwd = 3)
	    segments(x0 = x_mean_3[1]-0.5, x1 = x_mean_3[2]-0.5, y0 = y_mean_3, y1 = y_mean_3, col = 'azure2', lwd = 3)   
	}
	
	if(simulation){
	    abline( v=t_grid[to_s1]-0.5, col='dodgerblue', lwd=2)
	    abline( v=t_grid[to_s2]-0.5, col='firebrick1', lwd=2)
	    abline( v=t_grid[to_s3]-0.5, col='yellow2', lwd=2)
	} else {
	    abline( v=t_grid[to_s2]-0.5, col='darkorchid4', lwd=2)
	    abline( v=t_grid[to_s3]-0.5, col='darkgrey', lwd=2)
	}

    points(x=pb, y=data_format[indices_i, "RSA"])

	barplot(rbind(  colMeans(post_B_chain[, indices_i] == 1),
			        colMeans(post_B_chain[, indices_i] == 2),
				    colMeans(post_B_chain[, indices_i] == 3)), 
            col=c( 'dodgerblue', 'firebrick1', 'yellow2'), 
			xlab='time', space=0, col.main=main_color, border=NA, 
			ylab = "Posterior probability",
            xlim=range(pb) + c(-0.5,0.5), xaxt = 'n', yaxt = 'n', col.lab = main_color) 
	legend( 'topleft', inset=c(0,-.15), xpd=T, horiz=T, bty='n', x.intersp=.75,
	        legend=c( 'Start of stressor', 'End of stressor'), pch=15, pt.cex=1.5, 
	        col=c( 'darkorchid4', 'darkgrey'))
	legend( 'topright', inset=c(0,-.15), xpd=T, horiz=T, bty='n', x.intersp=.75,
	        legend=c( 'Baseline', 'State 2', 'State 3'), pch=15, pt.cex=1.5, 
	        col=c( 'dodgerblue', 'firebrick1', 'yellow2'))
    axis( side=1, at=t_grid_bar-0.5, col.axis=main_color, labels = t_grid)
	axis( side=2, at=seq(0,1,by=0.25), col.axis=main_color)

    if(simulation){
        abline( v=t_grid[to_s1]-0.5, col='dodgerblue', lwd=2)
        abline( v=t_grid[to_s2]-0.5, col='firebrick1', lwd=2)
        abline( v=t_grid[to_s3]-0.5, col='yellow2', lwd=2)
	} else {
	    abline( v=t_grid[to_s2]-0.5, col='darkorchid4', lwd=2)
	    abline( v=t_grid[to_s3]-0.5, col='darkgrey', lwd=2)
	}
	
	if(!interm) {
	    b_i_mle = as.numeric(c(B_MLE[[which(EIDs == i)]]))
	    to_s1_mle = (2:n_i)[diff(b_i_mle)!=0 & b_i_mle[-1]==1]
	    to_s2_mle = (2:n_i)[diff(b_i_mle)!=0 & b_i_mle[-1]==2]
	    to_s3_mle = (2:n_i)[diff(b_i_mle)!=0 & b_i_mle[-1]==3]
	    
	    # Manually changing the y-axis for b_i_mle
	    b_i_mle[b_i_mle == 1] = 0
	    b_i_mle[b_i_mle == 2] = -2
	    b_i_mle[b_i_mle == 3] = -1
	    
	    plot(x=pb, y=b_i_mle, type = 's', lwd = 4, main = 'Most likely state sequence',
	         xlab='time', ylab = ' ', col.main=main_color, col.lab = main_color,
	         xlim = range(pb) + c(-0.5,0.5),
	         xaxt='n', yaxt='n', ylim = c(-2.2,0.2))
	    axis( side=1, at=pb, col.axis=main_color, labels=t_grid)
	    axis( side=2, at=-2:0, col.axis=main_color, labels = c("state 2", "state 3", "baseline"),
	          cex.axis=1)
	    
	    if(simulation){
	        abline( v=t_grid[to_s1]-0.5, col='dodgerblue', lwd=2, lty = 2)
	        abline( v=t_grid[to_s2]-0.5, col='firebrick1', lwd=2, lty = 2)
	        abline( v=t_grid[to_s3]-0.5, col='yellow2', lwd=2, lty = 2)
	    } else {
	        abline( v=t_grid[1]-0.5, col='dodgerblue', lwd=4)
	        abline( v=t_grid[to_s1_mle]-0.5, col='dodgerblue', lwd=4)
	        abline( v=t_grid[to_s2_mle]-0.5, col='firebrick1', lwd=4)
	        abline( v=t_grid[to_s3_mle]-0.5, col='yellow2', lwd=4)
	    }
	}
	
}
dev.off()

count_transitions = matrix(0, nrow=3, ncol=3)
total_trans = 0
for(i in EIDs){
    b_i_mle = as.numeric(c(B_MLE[[which(EIDs == i)]]))
	indices_i = (data_format[,'ID..']==i)
	sub_dat = data_format[data_format[,"ID.."] == i, ]
    y_1_sub = sub_dat[,"State"]

    # Check
    if(length(b_i_mle) != nrow(sub_dat)) print(paste0(i, " issue"))

    start_run = FALSE
    for(t in 1:(length(b_i_mle) - 1)) {
        if(y_1_sub[t+1] != 1 && y_1_sub[t] == 1) start_run = TRUE
        if(start_run) {
            count_transitions[b_i_mle[t], b_i_mle[t+1]] = 
                count_transitions[b_i_mle[t], b_i_mle[t+1]] + 1
            total_trans = total_trans + 1
        }
    }
}

print(paste0("Total transitions: ", total_trans))
print(nrow(data_format) - length(EIDs))

print("Transition distribution")
print(count_transitions)
print("Transition Proportions")
print(round(count_transitions / total_trans, digits = 4))

# No transition:
no_trans = NULL
for(i in EIDs) {
    if(sum(B_MLE[[which(EIDs == i)]] != 1) == 0) { no_trans = c(no_trans, i) }
}

print("No transition out of baseline: ")
print(no_trans)
print("Yes, transitioned out of baseline: ")
print(EIDs[!(EIDs %in% no_trans)])

# -----------------------------------------------------------------------------
# Blunted responses -----------------------------------------------------------
# -----------------------------------------------------------------------------
id_indiv = c(302, 33103)
pdf(paste0("Plots/blunted_resp_", trial_num, ".pdf"))
par(mfrow=c(4,1), mar=c(2,4,2,4))#, bg='black', fg='green')
for(i in id_indiv){
	print(which(EIDs == i))
	indices_i = (data_format[,'ID..']==i)
	sub_dat = data_format[data_format[,"ID.."] == i, ]
	n_i = sum(indices_i)

    if(covariate_struct == 1) {
        plot_title = paste0('Participant: ', i)
    } else if(covariate_struct == 2) {
        
        dler_val = NULL
        for(a in EIDs) {
            dler_val = c(dler_val, unique(data_format[data_format[,"ID.."] == a, "DLER_avg"]))
        }
        mean_dler = mean(dler_val)
        
        cov_value = c(sub_dat[1,"DLER_avg"])
        cov_value = as.numeric(cov_value)
        cov_value[1] = cov_value[1] - mean_dler
        cov_value[1] = round(cov_value[1], digits = 3)

        plot_title = paste0('Participant: ', i, ', DLER: ', cov_value[1])
    } else {
        
        ages = NULL
        dler_val = NULL
        for(a in EIDs) {
            ages = c(ages, unique(data_format[data_format[,"ID.."] == a, "Age"]))
            dler_val = c(dler_val, unique(data_format[data_format[,"ID.."] == a, "DLER_avg"]))
        }
        mean_age = mean(ages)
        mean_dler = mean(dler_val)
        
        cov_value = c(sub_dat[1,c("Age", "sex1", "edu_yes", "DLER_avg")])
        cov_value = as.numeric(cov_value)
        cov_value[1] = cov_value[1] - mean_age
        cov_value[4] = cov_value[4] - mean_dler
        cov_value[1] = round(cov_value[1], digits = 3)
        cov_value[4] = round(cov_value[4], digits = 3)

        plot_title = paste0('Participant: ', i, ', sex: ', cov_value[2], 
                          ', pEdu: ', cov_value[3], ', DLER: ', cov_value[4],
                          ', age: ', cov_value[1])
    }

	t_grid = t_grid_bar = 1:n_i
	main_color = 'black'
	
	if(!simulation) {
	    x_mean_1 = c(min(which(sub_dat$State == 1)), max(which(sub_dat$State == 1))+1)
	    x_mean_2 = c(min(which(sub_dat$State == 2)), max(which(sub_dat$State == 2))+1)
	    x_mean_3 = c(min(which(sub_dat$State == 3)), max(which(sub_dat$State == 3)))
	    
	    y_mean_1 = mean(sub_dat$RSA[sub_dat$State == 1])
	    y_mean_2 = mean(sub_dat$RSA[sub_dat$State == 2])
	    y_mean_3 = mean(sub_dat$RSA[sub_dat$State == 3])   
	}
	
	b_i = as.numeric(data_format[ indices_i,"State"])
	to_s1 = c(1,(2:n_i)[diff(b_i)!=0 & b_i[-1]==1])
	to_s2 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==2]
	to_s3 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==3]

    # Plot the "ghost plot" to make everything line up
    pb = barplot(rbind( colMeans(post_B_chain[, indices_i] == 1),
			            colMeans(post_B_chain[, indices_i] == 2),
				        colMeans(post_B_chain[, indices_i] == 3)), 
            col=c( 'dodgerblue', 'firebrick1', 'yellow2'),
			xlab='time', space=0, col.main=main_color, border=NA, axes = F, plot = F) 

	plot(x=pb, y=data_format[indices_i, "RSA"], 
            xlab='time', ylab = 'RSA', col.main=main_color, 
            main = plot_title, xlim = range(pb) + c(-0.5,0.5),
	        xaxt='n', yaxt='n', col.lab = main_color)
	axis( side=1, at=pb, col.axis=main_color, labels=t_grid)
	axis( side=2, at=seq(min(data_format[indices_i, "RSA"]), 
                         max(data_format[indices_i, "RSA"])), col.axis=main_color)
	
    if(!simulation) {
	    segments(x0 = x_mean_1[1]-0.5, x1 = x_mean_1[2]-0.5, y0 = y_mean_1, y1 = y_mean_1, col = 'azure2', lwd = 3)
	}

	if(simulation){
	    abline( v=t_grid[to_s1]-0.5, col='dodgerblue', lwd=2)
	    abline( v=t_grid[to_s2]-0.5, col='firebrick1', lwd=2)
	    abline( v=t_grid[to_s3]-0.5, col='yellow2', lwd=2)
	} else {
	    abline( v=t_grid[to_s2]-0.5, col='darkorchid4', lwd=2)
	    abline( v=t_grid[to_s3]-0.5, col='darkgrey', lwd=2)
	}

    points(x=pb, y=data_format[indices_i, "RSA"])

	barplot(rbind(  colMeans(post_B_chain[, indices_i] == 1),
			        colMeans(post_B_chain[, indices_i] == 2),
				    colMeans(post_B_chain[, indices_i] == 3)), 
            col=c( 'dodgerblue', 'firebrick1', 'yellow2'), 
			xlab='time', space=0, col.main=main_color, border=NA, 
			ylab = "Posterior probability",
            xlim=range(pb) + c(-0.5,0.5), xaxt = 'n', yaxt = 'n', col.lab = main_color) 
	legend( 'topleft', inset=c(0,-.15), xpd=T, horiz=T, bty='n', x.intersp=.75,
	        legend=c( 'Start of stressor', 'End of stressor'), pch=15, pt.cex=1.5, 
	        col=c( 'darkorchid4', 'darkgrey'))
	legend( 'topright', inset=c(0,-.15), xpd=T, horiz=T, bty='n', x.intersp=.75,
	        legend=c( 'Baseline', 'State 2', 'State 3'), pch=15, pt.cex=1.5, 
	        col=c( 'dodgerblue', 'firebrick1', 'yellow2'))
    axis( side=1, at=t_grid_bar-0.5, col.axis=main_color, labels = t_grid)
	axis( side=2, at=seq(0,1,by=0.25), col.axis=main_color)

    if(simulation){
        abline( v=t_grid[to_s1]-0.5, col='dodgerblue', lwd=2)
        abline( v=t_grid[to_s2]-0.5, col='firebrick1', lwd=2)
        abline( v=t_grid[to_s3]-0.5, col='yellow2', lwd=2)
	} else {
	    abline( v=t_grid[to_s2]-0.5, col='darkorchid4', lwd=2)
	    abline( v=t_grid[to_s3]-0.5, col='darkgrey', lwd=2)
	}
	
}
dev.off()

# -----------------------------------------------------------------------------
# React/recover responses -----------------------------------------------------
# -----------------------------------------------------------------------------

id_indiv = c(414, 26275)
pdf(paste0("Plots/react_recover_resp_", trial_num, ".pdf"))
par(mfrow=c(4,1), mar=c(2,4,2,4))#, bg='black', fg='green')
for(i in id_indiv){
	print(which(EIDs == i))
	indices_i = (data_format[,'ID..']==i)
	sub_dat = data_format[data_format[,"ID.."] == i, ]
	n_i = sum(indices_i)

    if(covariate_struct == 1) {
        plot_title = paste0('Participant: ', i)
    } else if(covariate_struct == 2) {
        
        dler_val = NULL
        for(a in EIDs) {
            dler_val = c(dler_val, unique(data_format[data_format[,"ID.."] == a, "DLER_avg"]))
        }
        mean_dler = mean(dler_val)
        
        cov_value = c(sub_dat[1,"DLER_avg"])
        cov_value = as.numeric(cov_value)
        cov_value[1] = cov_value[1] - mean_dler
        cov_value[1] = round(cov_value[1], digits = 3)

        plot_title = paste0('Participant: ', i, ', DLER: ', cov_value[1])
    } else {
        
        ages = NULL
        dler_val = NULL
        for(a in EIDs) {
            ages = c(ages, unique(data_format[data_format[,"ID.."] == a, "Age"]))
            dler_val = c(dler_val, unique(data_format[data_format[,"ID.."] == a, "DLER_avg"]))
        }
        mean_age = mean(ages)
        mean_dler = mean(dler_val)
        
        cov_value = c(sub_dat[1,c("Age", "sex1", "edu_yes", "DLER_avg")])
        cov_value = as.numeric(cov_value)
        cov_value[1] = cov_value[1] - mean_age
        cov_value[4] = cov_value[4] - mean_dler
        cov_value[1] = round(cov_value[1], digits = 3)
        cov_value[4] = round(cov_value[4], digits = 3)

        plot_title = paste0('Participant: ', i, ', sex: ', cov_value[2], 
                          ', pEdu: ', cov_value[3], ', DLER: ', cov_value[4],
                          ', age: ', cov_value[1])
    }

	t_grid = t_grid_bar = 1:n_i
	main_color = 'black'
	
	if(!simulation) {
	    x_mean_1 = c(min(which(sub_dat$State == 1)), max(which(sub_dat$State == 1))+1)
	    x_mean_2 = c(min(which(sub_dat$State == 2)), max(which(sub_dat$State == 2))+1)
	    x_mean_3 = c(min(which(sub_dat$State == 3)), max(which(sub_dat$State == 3)))
	    
	    y_mean_1 = mean(sub_dat$RSA[sub_dat$State == 1])
	    y_mean_2 = mean(sub_dat$RSA[sub_dat$State == 2])
	    y_mean_3 = mean(sub_dat$RSA[sub_dat$State == 3])   
	}
	
	b_i = as.numeric(data_format[ indices_i,"State"])
	to_s1 = c(1,(2:n_i)[diff(b_i)!=0 & b_i[-1]==1])
	to_s2 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==2]
	to_s3 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==3]

    # Plot the "ghost plot" to make everything line up
    pb = barplot(rbind( colMeans(post_B_chain[, indices_i] == 1),
			            colMeans(post_B_chain[, indices_i] == 2),
				        colMeans(post_B_chain[, indices_i] == 3)), 
            col=c( 'dodgerblue', 'firebrick1', 'yellow2'),
			xlab='time', space=0, col.main=main_color, border=NA, axes = F, plot = F) 

	plot(x=pb, y=data_format[indices_i, "RSA"], 
            xlab='time', ylab = 'RSA', col.main=main_color, 
            main = plot_title, xlim = range(pb) + c(-0.5,0.5),
	        xaxt='n', yaxt='n', col.lab = main_color)
	axis( side=1, at=pb, col.axis=main_color, labels=t_grid)
	axis( side=2, at=seq(min(data_format[indices_i, "RSA"]), 
                         max(data_format[indices_i, "RSA"])), col.axis=main_color)
	
    if(!simulation) {
	    segments(x0 = x_mean_1[1]-0.5, x1 = x_mean_1[2]-0.5, y0 = y_mean_1, y1 = y_mean_1, col = 'azure2', lwd = 3)
	}

	if(simulation){
	    abline( v=t_grid[to_s1]-0.5, col='dodgerblue', lwd=2)
	    abline( v=t_grid[to_s2]-0.5, col='firebrick1', lwd=2)
	    abline( v=t_grid[to_s3]-0.5, col='yellow2', lwd=2)
	} else {
	    abline( v=t_grid[to_s2]-0.5, col='darkorchid4', lwd=2)
	    abline( v=t_grid[to_s3]-0.5, col='darkgrey', lwd=2)
	}

    points(x=pb, y=data_format[indices_i, "RSA"])

	barplot(rbind(  colMeans(post_B_chain[, indices_i] == 1),
			        colMeans(post_B_chain[, indices_i] == 2),
				    colMeans(post_B_chain[, indices_i] == 3)), 
            col=c( 'dodgerblue', 'firebrick1', 'yellow2'), 
			xlab='time', space=0, col.main=main_color, border=NA, 
			ylab = "Posterior probability",
            xlim=range(pb) + c(-0.5,0.5), xaxt = 'n', yaxt = 'n', col.lab = main_color) 
	legend( 'topleft', inset=c(0,-.15), xpd=T, horiz=T, bty='n', x.intersp=.75,
	        legend=c( 'Start of stressor', 'End of stressor'), pch=15, pt.cex=1.5, 
	        col=c( 'darkorchid4', 'darkgrey'))
	legend( 'topright', inset=c(0,-.15), xpd=T, horiz=T, bty='n', x.intersp=.75,
	        legend=c( 'Baseline', 'State 2', 'State 3'), pch=15, pt.cex=1.5, 
	        col=c( 'dodgerblue', 'firebrick1', 'yellow2'))
    axis( side=1, at=t_grid_bar-0.5, col.axis=main_color, labels = t_grid)
	axis( side=2, at=seq(0,1,by=0.25), col.axis=main_color)

    if(simulation){
        abline( v=t_grid[to_s1]-0.5, col='dodgerblue', lwd=2)
        abline( v=t_grid[to_s2]-0.5, col='firebrick1', lwd=2)
        abline( v=t_grid[to_s3]-0.5, col='yellow2', lwd=2)
	} else {
	    abline( v=t_grid[to_s2]-0.5, col='darkorchid4', lwd=2)
	    abline( v=t_grid[to_s3]-0.5, col='darkgrey', lwd=2)
	}
	
}
dev.off()

# -----------------------------------------------------------------------------
# React responses -------------------------------------------------------------
# -----------------------------------------------------------------------------
id_indiv = c(29014,416) #26104,
pdf(paste0("Plots/react_resp_", trial_num, ".pdf"))
par(mfrow=c(4,1), mar=c(2,4,2,4))#, bg='black', fg='green')
for(i in id_indiv){
	print(which(EIDs == i))
	indices_i = (data_format[,'ID..']==i)
	sub_dat = data_format[data_format[,"ID.."] == i, ]
	n_i = sum(indices_i)

    if(covariate_struct == 1) {
        plot_title = paste0('Participant: ', i)
    } else if(covariate_struct == 2) {
        
        dler_val = NULL
        for(a in EIDs) {
            dler_val = c(dler_val, unique(data_format[data_format[,"ID.."] == a, "DLER_avg"]))
        }
        mean_dler = mean(dler_val)
        
        cov_value = c(sub_dat[1,"DLER_avg"])
        cov_value = as.numeric(cov_value)
        cov_value[1] = cov_value[1] - mean_dler
        cov_value[1] = round(cov_value[1], digits = 3)

        plot_title = paste0('Participant: ', i, ', DLER: ', cov_value[1])
    } else {
        
        ages = NULL
        dler_val = NULL
        for(a in EIDs) {
            ages = c(ages, unique(data_format[data_format[,"ID.."] == a, "Age"]))
            dler_val = c(dler_val, unique(data_format[data_format[,"ID.."] == a, "DLER_avg"]))
        }
        mean_age = mean(ages)
        mean_dler = mean(dler_val)
        
        cov_value = c(sub_dat[1,c("Age", "sex1", "edu_yes", "DLER_avg")])
        cov_value = as.numeric(cov_value)
        cov_value[1] = cov_value[1] - mean_age
        cov_value[4] = cov_value[4] - mean_dler
        cov_value[1] = round(cov_value[1], digits = 3)
        cov_value[4] = round(cov_value[4], digits = 3)

        plot_title = paste0('Participant: ', i, ', sex: ', cov_value[2], 
                          ', pEdu: ', cov_value[3], ', DLER: ', cov_value[4],
                          ', age: ', cov_value[1])
    }

	t_grid = t_grid_bar = 1:n_i
	main_color = 'black'
	
	if(!simulation) {
	    x_mean_1 = c(min(which(sub_dat$State == 1)), max(which(sub_dat$State == 1))+1)
	    x_mean_2 = c(min(which(sub_dat$State == 2)), max(which(sub_dat$State == 2))+1)
	    x_mean_3 = c(min(which(sub_dat$State == 3)), max(which(sub_dat$State == 3)))
	    
	    y_mean_1 = mean(sub_dat$RSA[sub_dat$State == 1])
	    y_mean_2 = mean(sub_dat$RSA[sub_dat$State == 2])
	    y_mean_3 = mean(sub_dat$RSA[sub_dat$State == 3])   
	}
	
	b_i = as.numeric(data_format[ indices_i,"State"])
	to_s1 = c(1,(2:n_i)[diff(b_i)!=0 & b_i[-1]==1])
	to_s2 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==2]
	to_s3 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==3]

    # Plot the "ghost plot" to make everything line up
    pb = barplot(rbind( colMeans(post_B_chain[, indices_i] == 1),
			            colMeans(post_B_chain[, indices_i] == 2),
				        colMeans(post_B_chain[, indices_i] == 3)), 
            col=c( 'dodgerblue', 'firebrick1', 'yellow2'),
			xlab='time', space=0, col.main=main_color, border=NA, axes = F, plot = F) 

	plot(x=pb, y=data_format[indices_i, "RSA"], 
            xlab='time', ylab = 'RSA', col.main=main_color, 
            main = plot_title, xlim = range(pb) + c(-0.5,0.5),
	        xaxt='n', yaxt='n', col.lab = main_color)
	axis( side=1, at=pb, col.axis=main_color, labels=t_grid)
	axis( side=2, at=seq(min(data_format[indices_i, "RSA"]), 
                         max(data_format[indices_i, "RSA"])), col.axis=main_color)
	
    if(!simulation) {
	    segments(x0 = x_mean_1[1]-0.5, x1 = x_mean_1[2]-0.5, y0 = y_mean_1, y1 = y_mean_1, col = 'azure2', lwd = 3)
	}

	if(simulation){
	    abline( v=t_grid[to_s1]-0.5, col='dodgerblue', lwd=2)
	    abline( v=t_grid[to_s2]-0.5, col='firebrick1', lwd=2)
	    abline( v=t_grid[to_s3]-0.5, col='yellow2', lwd=2)
	} else {
	    abline( v=t_grid[to_s2]-0.5, col='darkorchid4', lwd=2)
	    abline( v=t_grid[to_s3]-0.5, col='darkgrey', lwd=2)
	}

    points(x=pb, y=data_format[indices_i, "RSA"])

	barplot(rbind(  colMeans(post_B_chain[, indices_i] == 1),
			        colMeans(post_B_chain[, indices_i] == 2),
				    colMeans(post_B_chain[, indices_i] == 3)), 
            col=c( 'dodgerblue', 'firebrick1', 'yellow2'), 
			xlab='time', space=0, col.main=main_color, border=NA, 
			ylab = "Posterior probability",
            xlim=range(pb) + c(-0.5,0.5), xaxt = 'n', yaxt = 'n', col.lab = main_color) 
	legend( 'topleft', inset=c(0,-.15), xpd=T, horiz=T, bty='n', x.intersp=.75,
	        legend=c( 'Start of stressor', 'End of stressor'), pch=15, pt.cex=1.5, 
	        col=c( 'darkorchid4', 'darkgrey'))
	legend( 'topright', inset=c(0,-.15), xpd=T, horiz=T, bty='n', x.intersp=.75,
	        legend=c( 'Baseline', 'State 2', 'State 3'), pch=15, pt.cex=1.5, 
	        col=c( 'dodgerblue', 'firebrick1', 'yellow2'))
    axis( side=1, at=t_grid_bar-0.5, col.axis=main_color, labels = t_grid)
	axis( side=2, at=seq(0,1,by=0.25), col.axis=main_color)

    if(simulation){
        abline( v=t_grid[to_s1]-0.5, col='dodgerblue', lwd=2)
        abline( v=t_grid[to_s2]-0.5, col='firebrick1', lwd=2)
        abline( v=t_grid[to_s3]-0.5, col='yellow2', lwd=2)
	} else {
	    abline( v=t_grid[to_s2]-0.5, col='darkorchid4', lwd=2)
	    abline( v=t_grid[to_s3]-0.5, col='darkgrey', lwd=2)
	}
	
}
dev.off()