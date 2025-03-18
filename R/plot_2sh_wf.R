# Code to analyze output data files from SLiM simulations and generate plots
# 2 serial homlogs with conserved identities and evolving states (z1 and z2), WF populations

library(ggplot2)

setwd("your_working_directory/2sh_wf")

# Generate a data matrix containing parameter combinations to be examined
n_share=c(0,5,10,15,20) # Number of shared targets of 2 IDGs (all value to consider)
cor_selection=c(0,0.9,-0.9) # Coefficient of correlational selection (all value to consider)
par.comb.all=matrix(0,nrow=1,ncol=2)
for(i in 1:length(n_share)){
  for(j in 1:length(cor_selection)){
    row=c(n_share[i],cor_selection[j])
    par.comb.all=rbind(par.comb.all,row)
  }
}
par.comb.all=par.comb.all[2:nrow(par.comb.all),]
colnames(par.comb.all)=c("n_share","cor_selection")
rownames(par.comb.all)=NULL

# Plot results for directional selection (one type of scenario each time)
# Go to the corresponding directory and set optimum of z2
setwd("./dir_for_selection_on_1_trait");opt2=0
setwd("./dir_for_divergent selection");opt2=-10
setwd("./dir_for_concordant_selection");opt2=10

# Open files and extract information
prefix=paste("sim_out_end_10_",opt2,"_",sep="") # Common part shared by all file names
out=matrix(0,nrow=nrow(par.comb.all),ncol=4)
for(c in 1:nrow(par.comb.all)){
  fn=paste(prefix,par.comb.all[c,1],"_",par.comb.all[c,2],".txt",sep="")
  d<-read.table(fn,sep="\t")
  out[c,1]=mean(d[,2]-40) # Mean of z1 (converted to difference relative to ancestral state)
  out[c,2]=mean(d[,3]-40) # Mean of z2 (converted to difference relative to ancestral state
  out[c,3]=sd(d[,2]-40)/(sqrt(nrow(d))) # SE of z1
  out[c,4]=sd(d[,3]-40)/(sqrt(nrow(d))) # SE of z2
}
colnames(out)=c("z1","z2","se1","se2") # Name the columns

d=data.frame(par.comb.all,out) # Data frame used for plotting
setwd("..")
fn=paste("sum_2sh_dir_",opt2,".txt",sep="")
write.table(d,file=fn,sep="\t") # Write summary data file

# Plot 1 curve in each panel (each figure will be 3x3; axis labels are empty, preparing for merging panels afterwards)
# Subsets based on correlational selection
d1=d[which(d$cor_selection==0),]
d2=d[which(d$cor_selection==0.9),]
d3=d[which(d$cor_selection==-0.9),]

dsub=d1 # Decide which subset to plot
# Plot z1 divergence
g<-ggplot(dsub,aes(x=n_share,y=z1))
g=g+geom_point(size=3)+geom_line(lwd=2)+geom_errorbar(aes(ymin=z1-se1,ymax=z1+se1,width=0.2))
g=g+theme_classic()+theme(axis.text=element_text(size=12))
g=g+xlab("")+ylab("")+ylim(c(0,11))
# Plot divergence between SHs
h<-ggplot(dsub,aes(x=n_share,y=abs(z1-z2)))
lb=abs(dsub$z1-dsub$z2)-(dsub$se1+dsub$se2);lb[which(lb<0)]=0 # If lower bound for error bar is negative, set it to be 0
h=h+geom_point(size=3)+geom_line(lwd=2)+geom_errorbar(aes(ymin=lb, ymax=abs(z1-z2)+(se1+se2)),width=0.2)
h=h+theme_classic()+theme(axis.text=element_text(size=12))
h=h+xlab("")+ylab("")+ylim(c(0,20))
# Write file
fn1=paste("slim_out_z1_",opt2,"_",dsub[1,2],".pdf",sep="")
fn2=paste("slim_out_div_",opt2,"_",dsub[1,2],".pdf",sep="")
ggsave(plot=g,file=fn1,width=7,height=5)
ggsave(plot=h,file=fn2,width=7,height=5)


# Plot for stabilizing selection
setwd("./dir_for_stabilizing_selection")
prefix="sim_out_end_0_0_" # Common part shared by all file names
out=matrix(0,nrow=nrow(par.comb.all),ncol=11)
for(c in 1:nrow(par.comb.all)){
  fn=paste(prefix,par.comb.all[c,1],"_",par.comb.all[c,2],".txt",sep="")
  d<-read.table(fn,sep="\t")
  out[c,1]=mean(d[,2]-40) # Mean of z1 (converted to difference relative to ancestral state)
  out[c,2]=mean(d[,3]-40) # Mean of z2 (converted to difference relative to ancestral state)
  out[c,3]=var(d[,2]) # Variance of z1
  out[c,4]=var(d[,3]) # Variance of z2
  out[c,5]=sd(d[,2]-40)/(sqrt(nrow(d))) # SE of z1
  out[c,6]=sd(d[,3]-40)/(sqrt(nrow(d))) # SE of z2
  out[c,7]=sqrt(2*var(d[,2])^2/(nrow(d)-1)) # SE of var1
  out[c,8]=sqrt(2*var(d[,3])^2/(nrow(d)-1)) # SE of var2
  out[c,9]=sqrt(out[c,7]+out[c,8]) # SE of variance of divergence between SHs
  out[c,10]=cor(d[,2],d[,3]) # Correlation between z1 and z2
  out[c,11]=(1-out[c,10]^2)/sqrt(nrow(d)-3)
}
colnames(out)=c("z1","z2","var1","var2","se_z1","se_z2","se_v1","se_v2","se_vd","cor","se_cor") # Name the columns

d=data.frame(par.comb.all,out) # Data frame used for plotting
setwd("..")
write.table(d,file="sum_2sh_stab.txt",sep="\t") # Write summary data file

d1=d[which(d$cor_selection==0),]
d2=d[which(d$cor_selection==0.9),]
d3=d[which(d$cor_selection==-0.9),]

dsub=d3 # Decide which subset to plot
# Plot mean divergence between SHs
h<-ggplot(dsub,aes(x=n_share,y=abs(z1-z2)))
lb=abs(dsub$z1-dsub$z2)-(dsub$se_z1+dsub$se_z2);lb[which(lb<0)]=0 # If lower bound for error bar is negative, set it to be 0
h=h+geom_point(size=3)+geom_line(lwd=2)+geom_errorbar(aes(ymin=lb, ymax=abs(z1-z2)+(se_z1+se_z2)),width=0.2)
h=h+theme_classic()+theme(axis.text=element_text(size=12))
h=h+xlab("")+ylab("")+ylim(c(0,5))
fn=paste("slim_stab_out_div_",dsub[1,2],".pdf",sep="")
ggsave(plot=h,file=fn,width=7,height=5)

# Plot variance of divergence between SHs
h<-ggplot(dsub,aes(x=n_share,y=var1+var2))
h=h+geom_point(size=3)+geom_line(lwd=2)+geom_errorbar(aes(ymin=var1+var2-se_vd, ymax=var1+var2+se_vd),width=0.2)
h=h+theme_classic()+theme(axis.text=element_text(size=12))
h=h+xlab("")+ylab("")+ylim(c(0,12.5))
fn=paste("slim_stab_out_vd_",dsub[1,2],".pdf",sep="")
ggsave(plot=h,file=fn,width=7,height=5)

# Plot correlation between SH
h<-ggplot(dsub,aes(x=n_share,y=cor))
h=h+geom_point(size=3)+geom_line(lwd=2)+geom_errorbar(aes(ymin=cor-se_cor, ymax=cor+se_cor),width=0.2)
h=h+theme_classic()+theme(axis.text=element_text(size=12))
h=h+xlab("")+ylab("")+ylim(c(-1,1))
fn=paste("slim_stab_out_cor_",dsub[1,2],".pdf",sep="")
ggsave(plot=h,file=fn,width=7,height=5)

