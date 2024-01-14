# library(matrixStats)
library(plotrix)

# Information defining which approach to take ----------------------------------
trial_num = 3
simulation = F
case_b = T

args = commandArgs(TRUE)
seed = as.numeric(args[1])
# ------------------------------------------------------------------------------

Dir = 'Model_out/'
load(paste0('Model_out/par_median', trial_num, '.rda'))
par_index = list(zeta=1:30, misclass=42:45, delta = 31:33, tau2 = 34, sigma2 = 35:37,
                 gamma = 38:41)

file_name = NULL
if(simulation) {
    if(case_b) {
        file_name = paste0(Dir,'mcmc_out_',toString(seed), '_', trial_num, '_sim_30b.rda')   
    } else {
        file_name = paste0(Dir,'mcmc_out_',toString(seed), '_', trial_num, '_sim_30.rda')      
    }
} else {
    if(case_b) {
        file_name = paste0(Dir,'mcmc_out_',toString(seed), '_', trial_num, '_30b.rda')   
    } else {
        file_name = paste0(Dir,'mcmc_out_',toString(seed), '_', trial_num, '_30.rda')      
    }
}

load(file_name)

if(simulation) {
    # Simulation
    load('Data/sim_data_1_30.rda')
    data_format = sim_data
    
    EIDs = unique(data_format[,"ID.."])
} else {
    # Real data analysis
    load('../Data/data_format_30.rda')
    data_format = data_format_30
    
    miss_info = c(26296, 29698, 30625, 401, 423, 419, 457)
    data_format = data_format[!(data_format[,"ID.."] %in% miss_info), ]
    
    EIDs = unique(data_format[,"ID.."])
    
    ages = NULL
    dler_val = NULL
    for(a in EIDs) {
        ages = c(ages, unique(data_format[data_format[,"ID.."] == a, "Age"]))
        dler_val = c(dler_val, unique(data_format[data_format[,"ID.."] == a, "DLER_avg"]))
    }
    mean_age = mean(ages)
    mean_dler = mean(dler_val)
}

# New patients ---------------------------------------------------------------
if(simulation) {
    if(case_b) {
        pdf_title = paste0('Plots/chart_plot_', trial_num, '_sim_30b.pdf')
    } else {
        pdf_title = paste0('Plots/chart_plot_', trial_num, '_sim_30.pdf')
    }
} else {
    if(case_b) {
        pdf_title = paste0('Plots/chart_plot_', trial_num, '_30b_s1.pdf')
    } else {
        pdf_title = paste0('Plots/chart_plot_', trial_num, '_30.pdf')
    }
}

b_chain_ind = 5000:195000

pdf(pdf_title)
panels = c(4, 1)
par(mfrow=panels, mar=c(2,4,2,4))#, bg='black', fg='green')
for(i in EIDs){
	print(which(EIDs == i))
	indices_i = (data_format[,'ID..']==i)
	sub_dat = data_format[data_format[,"ID.."] == i, ]
	n_i = sum(indices_i)

    cov_value = c(sub_dat[1,c("Age", "sex1", "edu_yes", "DLER_avg")])
    cov_value = as.numeric(cov_value)
    cov_value[1] = cov_value[1] - mean_age
    cov_value[4] = cov_value[4] - mean_dler
    baseline_mean = par_median[par_index$delta[1]] + sum(par_median[par_index$gamma] * cov_value)
    baseline_mean = round(baseline_mean, digits = 3)
    cov_value[1] = round(cov_value[1], digits = 3)
    cov_value[4] = round(cov_value[4], digits = 3)

	# t_grid = t_grid_bar = data_format[indices_i, "Time"]
	t_grid = t_grid_bar = 1:n_i
	main_color = 'black'
	
	x_mean_1 = c(min(which(sub_dat$State == 1)), max(which(sub_dat$State == 1))+1)
	x_mean_2 = c(min(which(sub_dat$State == 2)), max(which(sub_dat$State == 2))+1)
	x_mean_3 = c(min(which(sub_dat$State == 3)), max(which(sub_dat$State == 3)))
	
	y_mean_1 = mean(sub_dat$RSA[sub_dat$State == 1])
	y_mean_2 = mean(sub_dat$RSA[sub_dat$State == 2])
	y_mean_3 = mean(sub_dat$RSA[sub_dat$State == 3])
	
	if(simulation){
	    b_i = as.numeric(data_format[ indices_i,"True_state"])
	    to_s1 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==1]
	    to_s2 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==2]
	    to_s3 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==3]
    } else {
        b_i = as.numeric(data_format[ indices_i,"State"])
        to_s2 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==2]
        to_s3 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==3]
    }

    # Plot the "ghost plot" to make everything line up
    pb = barplot(rbind( colMeans(mcmc_out$B_chain[b_chain_ind, indices_i] == 1),
			            colMeans(mcmc_out$B_chain[b_chain_ind, indices_i] == 2),
				        colMeans(mcmc_out$B_chain[b_chain_ind, indices_i] == 3)), 
            col=c( 'dodgerblue', 'firebrick1', 'yellow2'),
			xlab='time', space=0, col.main=main_color, border=NA, axes = F, plot = F) 

	plot(x=pb, y=data_format[indices_i, "RSA"], 
            xlab='time', ylab = 'RSA', col.main=main_color, 
            main = paste0('Participant: ', i, ', sex: ', cov_value[2], 
                          ', pEdu: ', cov_value[3], ', DLER: ', cov_value[4],
                          ', age: ', cov_value[1], ', mean: ', baseline_mean), 
            xlim = range(pb) + c(-0.5,0.5),
	        xaxt='n', yaxt='n', col.lab = main_color)
	axis( side=1, at=pb, col.axis=main_color, labels=t_grid)
	axis( side=2, at=seq(min(data_format[indices_i, "RSA"]), 
                         max(data_format[indices_i, "RSA"])), col.axis=main_color)
	
	segments(x0 = x_mean_1[1]-0.5, x1 = x_mean_1[2]-0.5, y0 = y_mean_1, y1 = y_mean_1, col = 'azure2', lwd = 3)
	segments(x0 = x_mean_2[1]-0.5, x1 = x_mean_2[2]-0.5, y0 = y_mean_2, y1 = y_mean_2, col = 'azure2', lwd = 3)
	segments(x0 = x_mean_3[1]-0.5, x1 = x_mean_3[2]-0.5, y0 = y_mean_3, y1 = y_mean_3, col = 'azure2', lwd = 3)
	
	if(simulation){
	    abline( v=t_grid[to_s1]-0.5, col='darkolivegreen3', lwd=2)
	    abline( v=t_grid[to_s2]-0.5, col='darkorchid4', lwd=2)
	    abline( v=t_grid[to_s3]-0.5, col='darkgrey', lwd=2)
	} else {
	    abline( v=t_grid[to_s2]-0.5, col='darkorchid4', lwd=2)
	    abline( v=t_grid[to_s3]-0.5, col='darkgrey', lwd=2)
	}

    points(x=pb, y=data_format[indices_i, "RSA"])

	barplot(rbind(  colMeans(mcmc_out$B_chain[b_chain_ind, indices_i] == 1),
			        colMeans(mcmc_out$B_chain[b_chain_ind, indices_i] == 2),
				    colMeans(mcmc_out$B_chain[b_chain_ind, indices_i] == 3)), 
            col=c( 'dodgerblue', 'firebrick1', 'yellow2'), 
			xlab='time', space=0, col.main=main_color, border=NA, 
			ylab = "Posterior probability",
            xlim=range(pb) + c(-0.5,0.5), xaxt = 'n', yaxt = 'n', col.lab = main_color) 
	legend( 'topleft', inset=c(0,-.15), xpd=T, horiz=T, bty='n', x.intersp=.75,
	        legend=c( 'Start of event', 'End of event'), pch=15, pt.cex=1.5, 
	        col=c( 'darkorchid4', 'darkgrey'))
	legend( 'topright', inset=c(0,-.15), xpd=T, horiz=T, bty='n', x.intersp=.75,
	        legend=c( 'Baseline', 'State 2', 'State 3'), pch=15, pt.cex=1.5, 
	        col=c( 'dodgerblue', 'firebrick1', 'yellow2'))
    axis( side=1, at=t_grid_bar-0.5, col.axis=main_color, labels = t_grid)
	axis( side=2, at=seq(0,1,by=0.25), col.axis=main_color)

    if(simulation){
        abline( v=t_grid[to_s1]-0.5, col='darkolivegreen3', lwd=2)
        abline( v=t_grid[to_s2]-0.5, col='darkorchid4', lwd=2)
        abline( v=t_grid[to_s3]-0.5, col='darkgrey', lwd=2)
	} else {
	    abline( v=t_grid[to_s2]-0.5, col='darkorchid4', lwd=2)
	    abline( v=t_grid[to_s3]-0.5, col='darkgrey', lwd=2)
	}
	
}
dev.off()
