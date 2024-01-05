# library(matrixStats)
library(plotrix)

# Information defining which approach to take ----------------------------------
trial_num = 5
simulation = F
case_b = T
# ------------------------------------------------------------------------------

Dir = 'Model_out/'
load(paste0('Model_out/par_median', trial_num, '.rda'))
par_index = list(zeta=1:30, misclass=42:45, delta = 31:33, tau2 = 34, sigma2 = 35:37,
                 gamma = 38:41)


file_name = NULL
if(simulation) {
    if(case_b) {
        file_name = paste0("Model_out/B_chain_", trial_num, "_sim_30b.rda")
        file_name2 = paste0("Model_out/B_chain_", trial_num, "_sim_30b_MLE.rda")
    } else {
        file_name = paste0("Model_out/B_chain_", trial_num, "_sim_30.rda")   
        file_name2 = paste0("Model_out/B_chain_", trial_num, "_sim_30_MLE.rda")   
    }
} else {
    if(case_b) {
        file_name = paste0("Model_out/B_chain_", trial_num, "_30b_s1_up.rda")
        file_name2 = paste0("Model_out/B_chain_", trial_num, "_30b_s1_MLE_up.rda")
    } else {
        file_name = paste0("Model_out/B_chain_", trial_num, "_30.rda")  
        file_name2 = paste0("Model_out/B_chain_", trial_num, "_30_MLE.rda")  
    }
}

load(file_name)
load(file_name2)

if(simulation) {
    # Simulation
    load('Data/sim_data_1_30.rda')
    data_format = sim_data
} else {
    # Real data analysis
    load('Data/data_format_30.rda')
    data_format = data_format_30
    
    miss_info = c(26296, 29698, 30625, 401, 423, 419, 457)
    data_format = data_format[!(data_format[,"ID.."] %in% miss_info), ]
}

EIDs = unique(data_format[,"ID.."])

# New patients ---------------------------------------------------------------
if(simulation) {
    if(case_b) {
        pdf_title = paste0('Plots/chart_plot_', trial_num, '_sim_30b.pdf')
    } else {
        if(use_labels) {
            pdf_title = paste0('Plots/chart_plot_', trial_num, '_sim_30.pdf')
        } else {
            pdf_title = paste0('Plots/chart_plot_', trial_num, '_sim_30_nl.pdf')   
        }   
    }
} else {
    if(case_b) {
        pdf_title = paste0('Plots/chart_plot_', trial_num, '_30b_s1.pdf')
    } else {
        pdf_title = paste0('Plots/chart_plot_', trial_num, '_30.pdf')
    }
}

b_chain_ind = 45000:95000

pdf(pdf_title)
panels = c(3, 1)
par(mfrow=panels, mar=c(2,4,2,4))#, bg='black', fg='green')
for(i in EIDs){
	print(which(EIDs == i))
	indices_i = (data_format[,'ID..']==i)
	sub_dat = data_format[data_format[,"ID.."] == i, ]
	n_i = sum(indices_i)

    cov_value = c(sub_dat[1,c("Age", "sex1", "edu_yes", "DLER_avg")])
    cov_value = as.numeric(cov_value)
    cov_value[4] = round(cov_value[4], digits = 3)
    baseline_mean = par_median[par_index$delta[1]] + sum(par_median[par_index$gamma] * cov_value)
    baseline_mean = round(baseline_mean, digits = 3)

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
        
        b_i_mle = as.numeric(c(B_chain_MLE[[which(EIDs == i)]]))
	    to_s1_mle = (2:n_i)[diff(b_i_mle)!=0 & b_i_mle[-1]==1]
	    to_s2_mle = (2:n_i)[diff(b_i_mle)!=0 & b_i_mle[-1]==2]
	    to_s3_mle = (2:n_i)[diff(b_i_mle)!=0 & b_i_mle[-1]==3]
    }

    # Plot the "ghost plot" to make everything line up
    pb = barplot(rbind( colMeans(B_chain[b_chain_ind, indices_i] == 1),
			            colMeans(B_chain[b_chain_ind, indices_i] == 2),
				        colMeans(B_chain[b_chain_ind, indices_i] == 3)), 
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

	barplot(rbind(  colMeans(B_chain[b_chain_ind, indices_i] == 1),
			        colMeans(B_chain[b_chain_ind, indices_i] == 2),
				    colMeans(B_chain[b_chain_ind, indices_i] == 3)), 
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
	
	plot(x=pb, y=b_i_mle, type = 's', lwd = 4, main = 'Most likely state sequence',
	     xlab='time', ylab = 'State', col.main=main_color, col.lab = main_color,
	     xlim = range(pb) + c(-0.5,0.5),
	     xaxt='n', yaxt='n', ylim = c(0,4))
	axis( side=1, at=pb, col.axis=main_color, labels=t_grid)
	axis( side=2, at=1:3, col.axis=main_color)
	
	if(simulation){
	    abline( v=t_grid[to_s1]-0.5, col='darkolivegreen3', lwd=2)
	    abline( v=t_grid[to_s2]-0.5, col='darkorchid4', lwd=2)
	    abline( v=t_grid[to_s3]-0.5, col='darkgrey', lwd=2)
	} else {
	    abline( v=t_grid[1]-0.5, col='dodgerblue', lwd=4)
	    abline( v=t_grid[to_s1_mle]-0.5, col='dodgerblue', lwd=4)
	    abline( v=t_grid[to_s2_mle]-0.5, col='firebrick1', lwd=4)
	    abline( v=t_grid[to_s3_mle]-0.5, col='yellow2', lwd=4)
	}
	
}
dev.off()

print(pdf_title)

# Categorizing the participant encounters -------------------------------------
blunted_resp = NULL

for(i in EIDs) {
    b_i_mle = as.numeric(c(B_chain_MLE[[which(EIDs == i)]]))
    if(sum(b_i_mle == 1) == length(b_i_mle)) blunted_resp = c(blunted_resp, i)
}

print(blunted_resp)

if(case_b) {
    pdf_title2 = paste0('Plots/chart_plot_', trial_num, '_30b_s1_MLE.pdf')
} else {
    pdf_title2 = paste0('Plots/chart_plot_', trial_num, '_30_MLE.pdf')
}

pdf(pdf_title2)
panels = c(4, 1)
par(mfrow=panels, mar=c(2,4,2,4))

for(i in EIDs) {
    if(!(i %in% blunted_resp)) {
        indices_i = (data_format[,'ID..']==i)
        sub_dat = data_format[data_format[,"ID.."] == i, ]
        n_i = sum(indices_i)

        t_grid = t_grid_bar = 1:n_i
        main_color = 'black'
        
        if(simulation){
            b_i = as.numeric(data_format[ indices_i,"True_state"])
            to_s1 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==1]
            to_s2 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==2]
            to_s3 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==3]
        } else {
            b_i = as.numeric(data_format[ indices_i,"State"])
            to_s2 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==2]
            to_s3 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==3]
            
            b_i_mle = as.numeric(c(B_chain_MLE[[which(EIDs == i)]]))
            to_s1_mle = (2:n_i)[diff(b_i_mle)!=0 & b_i_mle[-1]==1]
            to_s2_mle = (2:n_i)[diff(b_i_mle)!=0 & b_i_mle[-1]==2]
            to_s3_mle = (2:n_i)[diff(b_i_mle)!=0 & b_i_mle[-1]==3]
        }

        # Plot the "ghost plot" to make everything line up
        pb = barplot(rbind( colMeans(B_chain[b_chain_ind, indices_i] == 1),
                            colMeans(B_chain[b_chain_ind, indices_i] == 2),
                            colMeans(B_chain[b_chain_ind, indices_i] == 3)), 
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
        
        if(simulation){
            abline( v=t_grid[to_s1]-0.5, col='darkolivegreen3', lwd=2)
            abline( v=t_grid[to_s2]-0.5, col='darkorchid4', lwd=2)
            abline( v=t_grid[to_s3]-0.5, col='darkgrey', lwd=2)
        } else {
            abline( v=t_grid[to_s2]-0.5, col='darkorchid4', lwd=2)
            abline( v=t_grid[to_s3]-0.5, col='darkgrey', lwd=2)
        }

        plot(x=pb, y=b_i_mle, type = 's', lwd = 4, main = 'Most likely state sequence',
            xlab='time', ylab = 'State', col.main=main_color, col.lab = main_color,
            xlim = range(pb) + c(-0.5,0.5),
            xaxt='n', yaxt='n', ylim = c(0,4))
        legend( 'topleft', inset=c(0,-.15), xpd=T, horiz=T, bty='n', x.intersp=.75,
                legend=c( 'Start of event', 'End of event'), pch=15, pt.cex=1.5, 
                col=c( 'darkorchid4', 'darkgrey'))
        legend( 'topright', inset=c(0,-.15), xpd=T, horiz=T, bty='n', x.intersp=.75,
                legend=c( 'Baseline', 'State 2', 'State 3'), pch=15, pt.cex=1.5, 
                col=c( 'dodgerblue', 'firebrick1', 'yellow2'))
        axis( side=1, at=pb, col.axis=main_color, labels=t_grid)
        axis( side=2, at=1:3, col.axis=main_color)
        
        if(simulation){
            abline( v=t_grid[to_s1]-0.5, col='darkolivegreen3', lwd=2)
            abline( v=t_grid[to_s2]-0.5, col='darkorchid4', lwd=2)
            abline( v=t_grid[to_s3]-0.5, col='darkgrey', lwd=2)
        } else {
            abline( v=t_grid[1]-0.5, col='dodgerblue', lwd=4)
            abline( v=t_grid[to_s1_mle]-0.5, col='dodgerblue', lwd=4)
            abline( v=t_grid[to_s2_mle]-0.5, col='firebrick1', lwd=4)
            abline( v=t_grid[to_s3_mle]-0.5, col='yellow2', lwd=4)
        }
    }
}
dev.off()

# Blunted response:
# c(25897,  26104,  26116,  26194,  26242,  26275,  26371,  26452,  26527,  26701,
#   26734,  26746,  26872,  26956,  27070,  27076,  27205,  27289,  27685,  27832,
#   28651,  28792,  29275,  29293,  29341,  29419,  29617,  29686,  29761,  29791,
#     302,  30550,  30559,  31588,  31750,  32260,  32305,  32419,  33478,  33718,
#   34105,    403,    404,    406,    407,    408,    410,    412,    413,    414,
#     415,    416,    417,    424,    425,    426,    427, 433001,    435,    436,
#     441,    445,    448,    450,    451,    452,    456,    459,    500,    501,
#     506,    508,    510,    516,    517)

# Saved ones for the paper -----------------------------------------------------
main_id = c(26116, 26269, 29014, 31039, 31471)
for(i in main_id) {
    pdf(paste0("Plots/", i, ".pdf"), paper="a4r")
    par(mfrow=c(3,1), mar=c(2,5,2,5))
    indices_i = (data_format[,'ID..']==i)
    sub_dat = data_format[data_format[,"ID.."] == i, ]
    n_i = sum(indices_i)

    t_grid = t_grid_bar = 1:n_i
    main_color = 'black'
    
    if(simulation){
        b_i = as.numeric(data_format[ indices_i,"True_state"])
        to_s1 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==1]
        to_s2 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==2]
        to_s3 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==3]
    } else {
        b_i = as.numeric(data_format[ indices_i,"State"])
        to_s2 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==2]
        to_s3 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==3]
        
        b_i_mle = as.numeric(c(B_chain_MLE[[which(EIDs == i)]]))
        to_s1_mle = (2:n_i)[diff(b_i_mle)!=0 & b_i_mle[-1]==1]
        to_s2_mle = (2:n_i)[diff(b_i_mle)!=0 & b_i_mle[-1]==2]
        to_s3_mle = (2:n_i)[diff(b_i_mle)!=0 & b_i_mle[-1]==3]
    }

    # Plot the "ghost plot" to make everything line up
    pb = barplot(rbind( colMeans(B_chain[b_chain_ind, indices_i] == 1),
                        colMeans(B_chain[b_chain_ind, indices_i] == 2),
                        colMeans(B_chain[b_chain_ind, indices_i] == 3)), 
            col=c( 'dodgerblue', 'firebrick1', 'yellow2'),
            xlab='time', space=0, col.main=main_color, border=NA, axes = F, plot = F) 

    plot(x=pb, y=data_format[indices_i, "RSA"], 
            xlab='time', ylab = 'RSA', col.main=main_color, 
            main = paste0('Participant: ', i, ', sex: ', cov_value[2], 
                          ', pEdu: ', cov_value[3], ', DLER: ', cov_value[4],
                          ', age: ', cov_value[1], ', mean: ', baseline_mean), 
            xlim = range(pb) + c(-0.5,0.5),
            xaxt='n', yaxt='n', col.lab = main_color, 
            cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    axis( side=1, at=pb, col.axis=main_color, labels=t_grid)
    axis( side=2, at=seq(min(data_format[indices_i, "RSA"]), 
                        max(data_format[indices_i, "RSA"])), 
                        col.axis=main_color, 
                        cex.lab=1.5, cex.axis=1.5)
    
    if(simulation){
        abline( v=t_grid[to_s1]-0.5, col='darkolivegreen3', lwd=2)
        abline( v=t_grid[to_s2]-0.5, col='darkorchid4', lwd=2)
        abline( v=t_grid[to_s3]-0.5, col='darkgrey', lwd=2)
    } else {
        abline( v=t_grid[to_s2]-0.5, col='darkorchid4', lwd=2)
        abline( v=t_grid[to_s3]-0.5, col='darkgrey', lwd=2)
    }

    barplot(rbind(  colMeans(B_chain[b_chain_ind, indices_i] == 1),
                    colMeans(B_chain[b_chain_ind, indices_i] == 2),
                    colMeans(B_chain[b_chain_ind, indices_i] == 3)), 
            col=c( 'dodgerblue', 'firebrick1', 'yellow2'), 
            xlab='time', space=0, col.main=main_color, border=NA, 
            ylab = "Posterior probability",
            xlim=range(pb) + c(-0.5,0.5), xaxt = 'n', yaxt = 'n', 
            col.lab = main_color, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5) 
    axis( side=1, at=t_grid_bar-0.5, col.axis=main_color, labels = t_grid)
    axis( side=2, at=c(0,1), col.axis=main_color,
            cex.lab=1.5, cex.axis=1.5)

    if(simulation){
        abline( v=t_grid[to_s1]-0.5, col='darkolivegreen3', lwd=2)
        abline( v=t_grid[to_s2]-0.5, col='darkorchid4', lwd=2)
        abline( v=t_grid[to_s3]-0.5, col='darkgrey', lwd=2)
    } else {
        abline( v=t_grid[to_s2]-0.5, col='darkorchid4', lwd=2)
        abline( v=t_grid[to_s3]-0.5, col='darkgrey', lwd=2)
    }
    legend( 'topleft', inset=c(0,-0.16), xpd=T, horiz=T, bty='n', x.intersp=.75,
            legend=c( 'Baseline', 'State 2', 'State 3'), pch=15, pt.cex=1.5, 
            col=c( 'dodgerblue', 'firebrick1', 'yellow2'), cex = 1.5)
    legend( 'topright', inset=c(0,-0.16), xpd=T, horiz=T, bty='n', x.intersp=.75,
            legend=c( 'Start of event', 'End of event'), pch=15, pt.cex=1.5, 
            col=c( 'darkorchid4', 'darkgrey'), cex = 1.5)
    
    plot(x=pb, y=b_i_mle, type = 's', lwd = 4,
        xlab='time', ylab = 'State', col.main=main_color, col.lab = main_color,
        xlim = range(pb) + c(-0.5,0.5),
        xaxt='n', yaxt='n', ylim = c(0,4), cex.lab=1.5, cex.axis=1.5, 
        cex.main=1.5, cex.sub=1.5)
    axis( side=1, at=pb, col.axis=main_color, labels=t_grid)
    axis( side=2, at=1:3, col.axis=main_color, cex.lab=1.5, cex.axis=1.5)
    
    if(simulation){
        abline( v=t_grid[to_s1]-0.5, col='darkolivegreen3', lwd=2)
        abline( v=t_grid[to_s2]-0.5, col='darkorchid4', lwd=2)
        abline( v=t_grid[to_s3]-0.5, col='darkgrey', lwd=2)
    } else {
        abline( v=t_grid[1]-0.5, col='dodgerblue', lwd=4)
        abline( v=t_grid[to_s1_mle]-0.5, col='dodgerblue', lwd=4)
        abline( v=t_grid[to_s2_mle]-0.5, col='firebrick1', lwd=4)
        abline( v=t_grid[to_s3_mle]-0.5, col='yellow2', lwd=4)
    }
    dev.off()
}