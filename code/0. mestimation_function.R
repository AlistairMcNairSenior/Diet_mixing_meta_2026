
# Clean up
rm(list=ls())

# Function for meta-estimation plots
gg_mestimation<-function(data, group, stat, n, control_mu, mu, ci_l=NA, ci_u=NA, se=NA, tau, n_breaks=waiver()){
		
		# Source packages
		require(ggplot2)
		require(ggbeeswarm)
		
		# Set up data
		data$group<-data[,group]
		data$stat<-data[,stat]
		data$n<-data[,n]
		
		# name axes etc
		if(stat == "lnX"){
			es<-"lnRR"
			ylab<-"Log Sample Mean"
		}else{
			es<-"lnCVR"
			ylab<-"Log Sample CV"
		}
		
		# Calculate CIs and PIs
		if((is.na(ci_l) + is.na(ci_u)) > 0){
			CI<-c(mu - se*1.96, mu + se*1.96)
		} else {
			CI<-c(ci_l, ci_u)
		}
		PI2<-c(mu - tau*1.96, mu + tau*1.96)
		
		plot<-ggplot(data, aes(x=group, y=stat, group=group, size=n)) +
			annotate("rect", xmin=0.15, xmax=2.85, 
			ymin=control_mu+PI2[1], 
			ymax=control_mu+PI2[2], 
			fill='cornflowerblue', alpha=0.3) +
			annotate("rect", xmin=0.15, xmax=2.85, 
			ymin=control_mu+CI[1], 
			ymax=control_mu+CI[2], 
			fill='red', alpha=0.3) +
			geom_hline(yintercept=control_mu+mu, col="red", size=0.5) +
			geom_hline(yintercept=control_mu, col="black", size=0.5) +
			geom_beeswarm(alpha=0.5) +
			theme_bw() +
			ylab(ylab) + xlab("Treatment Condition") +
			labs(size="Sample Size") +
			scale_x_discrete(labels=c("Mix Food", "Single Food")) +
			scale_y_continuous(sec.axis = sec_axis(~ . - control_mu, breaks=c(-2, -1, 0, 1, 2), name=" ")) + 
			scale_size_continuous(range=c(0.01, 3.5), breaks=n_breaks) +
			annotate("text", label=es, x=3.05, y=control_mu, angle=90) +
			coord_cartesian(xlim=c(0.75, 2.25), clip="off") +
			theme(legend.position="inside", legend.position.inside=c(0,0), legend.justification=c(0,0), legend.background=element_blank())

		return(plot)

}

