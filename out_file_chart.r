# library(matrixStats)
library(plotrix)

# Information defining which approach to take ----------------------------------
trial_num = 3
simulation = F
case_b = T
# ------------------------------------------------------------------------------

Dir = 'Model_out/'

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
        file_name = paste0("Model_out/B_chain_", trial_num, "_30b_s1.rda")
        file_name2 = paste0("Model_out/B_chain_", trial_num, "_30b_s1_MLE.rda")
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

pdf(pdf_title)
panels = c(4, 1)
par(mfrow=panels, mar=c(2,4,2,4), bg='black', fg='green')
for(i in EIDs){
	print(which(EIDs == i))
	indices_i = (data_format[,'ID..']==i)
	sub_dat = data_format[data_format[,"ID.."] == i, ]
	n_i = sum(indices_i)

	# t_grid = t_grid_bar = data_format[indices_i, "Time"]
	t_grid = t_grid_bar = 1:n_i
	
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
	
    b_chain_ind = 50000:95000

    # Plot the "ghost plot" to make everything line up
    pb = barplot(rbind( colMeans(B_chain[b_chain_ind, indices_i] == 1),
			            colMeans(B_chain[b_chain_ind, indices_i] == 2),
				        colMeans(B_chain[b_chain_ind, indices_i] == 3)), 
            col=c( 'dodgerblue', 'firebrick1', 'yellow2'), 
			xlab='time', space=0, col.main='green', border=NA, axes = F, plot = F) 

	plot(x=pb, y=data_format[indices_i, "RSA"], 
            xlab='time', ylab = 'RSA', col.main='green', 
            main = paste0('Participant: ', i), xlim = range(pb) + c(-0.5,0.5),
	        xaxt='n')
	# grid( nx=20, NULL, col='white')
	axis( side=1, at=pb, col.axis='green', labels=t_grid)
	axis( side=2, at=seq(min(data_format[indices_i, "RSA"]), 
                         max(data_format[indices_i, "RSA"])), col.axis='green')
	
	if(simulation){
	    abline( v=t_grid[to_s1]-0.5, col='darkolivegreen3', lwd=4)
	    abline( v=t_grid[to_s2]-0.5, col='darkorchid4', lwd=4)
	    abline( v=t_grid[to_s3]-0.5, col='darkgrey', lwd=4)
	} else {
	    abline( v=t_grid[to_s2]-0.5, col='darkorchid4', lwd=4)
	    abline( v=t_grid[to_s3]-0.5, col='darkgrey', lwd=4)
	    
	    abline( v=t_grid[1]-0.5, col='dodgerblue', lwd=1)
        abline( v=t_grid[to_s1_mle]-0.5, col='dodgerblue', lwd=1)
	    abline( v=t_grid[to_s2_mle]-0.5, col='firebrick1', lwd=1)
	    abline( v=t_grid[to_s3_mle]-0.5, col='yellow2', lwd=1)
	}

	barplot(rbind(  colMeans(B_chain[b_chain_ind, indices_i] == 1),
			        colMeans(B_chain[b_chain_ind, indices_i] == 2),
				    colMeans(B_chain[b_chain_ind, indices_i] == 3)), 
            col=c( 'dodgerblue', 'firebrick1', 'yellow2'), 
			xlab='time', space=0, col.main='green', border=NA, 
            xlim=range(pb) + c(-0.5,0.5)) 
    # grid( nx=20, NULL, col='white')
	legend( 'topright', inset=c(0,-.22), xpd=T, horiz=T, bty='n', x.intersp=.75,
	        legend=c( 'Baseline', 'Stress', 'Recovery'), pch=15, pt.cex=1.5, 
	        col=c( 'dodgerblue', 'firebrick1', 'yellow2'))
	legend( 'topleft', inset=c(0,-.22), xpd=T, horiz=T, bty='n', x.intersp=.75,
	        legend=c( 'Start of stress period', 'End of stress period'), pch=15, pt.cex=1.5, 
	        col=c( 'darkorchid4', 'darkgrey'))
    axis( side=1, at=t_grid_bar-0.5, col.axis='green', labels = t_grid)
	axis( side=2, at=0:1, col.axis='green')

    if(simulation){
        abline( v=t_grid[to_s1]-0.5, col='darkolivegreen3', lwd=4)
        abline( v=t_grid[to_s2]-0.5, col='darkorchid4', lwd=4)
        abline( v=t_grid[to_s3]-0.5, col='darkgrey', lwd=4)
	} else {
	    abline( v=t_grid[to_s2]-0.5, col='darkorchid4', lwd=4)
	    abline( v=t_grid[to_s3]-0.5, col='darkgrey', lwd=4)
	    
	    abline( v=t_grid[1]-0.5, col='dodgerblue', lwd=1)
	    abline( v=t_grid[to_s1_mle]-0.5, col='dodgerblue', lwd=1)
	    abline( v=t_grid[to_s2_mle]-0.5, col='firebrick1', lwd=1)
	    abline( v=t_grid[to_s3_mle]-0.5, col='yellow2', lwd=1)
	}
}
dev.off()

print(pdf_title)
