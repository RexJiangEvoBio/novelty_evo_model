# Code to analyze output data files from SLiM simulations and generate plots
# 2 serial homlogs with conserved identities and evolving states (z1 and z2), non-WF populations

library(ggplot2)

setwd("your_working_directory/2sh_nwf")

# Generate a data matrix containing parameter combinations to be examined
n_share=c(0,5,10,15,20)
cor_selection=c(0,0.9,-0.9)
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
prefix=paste("sim_nwf_out_end_10_",opt2,"_",sep="") # Common part shared by all file names
out=matrix(0,nrow=nrow(par.comb.all),ncol=6)
for(c in 1:nrow(par.comb.all)){
  fn=paste(prefix,par.comb.all[c,1],"_",par.comb.all[c,2],".txt",sep="")
  d<-read.table(fn,sep="\t")
  out[c,1]=mean(d[,3]-40) # Mean of z1 (converted to difference relative to ancestral state)
  out[c,2]=mean(d[,4]-40) # Mean of z2 (converted to difference relative to ancestral state)
  out[c,3]=sd(d[,3]-40)/(sqrt(nrow(d))) # SE of z1
  out[c,4]=sd(d[,4]-40)/(sqrt(nrow(d))) # SE of z2
  out[c,5]=mean(d[,2]) # Population size
  out[c,6]=sd(d[,2])/(sqrt(nrow(d))) # SE of population size
}
colnames(out)=c("z1","z2","se1","se2","N","se_N") # Name the columns

setwd("..") # Change working directory back before plotting
d=data.frame(par.comb.all,out) # Data frame used for plotting

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
g=g+xlab("")+ylab("")+ylim(c(0,10))
# Plot divergence between SHs
h<-ggplot(dsub,aes(x=n_share,y=abs(z1-z2)))
lb=abs(dsub$z1-dsub$z2)-(dsub$se1+dsub$se2);lb[which(lb<0)]=0 # If lower bound for error bar is negative, set it to be 0
h=h+geom_point(size=3)+geom_line(lwd=2)+geom_errorbar(aes(ymin=lb, ymax=abs(z1-z2)+(se1+se2)),width=0.2)
h=h+theme_classic()+theme(axis.text=element_text(size=12))
h=h+xlab("")+ylab("")+ylim(c(0,20))
# Plot population size
k<-ggplot(dsub,aes(x=n_share,y=N))
k=k+geom_point(size=3)+geom_line(lwd=2)+geom_errorbar(aes(ymin=N-se_N, ymax=N+se_N),width=0.2)
k=k+theme_classic()+theme(axis.text=element_text(size=12))
k=k+xlab("")+ylab("")+ylim(c(0,2100))
# Write file
fn1=paste("slim_nwf_out_z1_",opt2,"_",dsub[1,2],".pdf",sep="")
fn2=paste("slim_nwf_out_div_",opt2,"_",dsub[1,2],".pdf",sep="")
fn3=paste("slim_nwf_out_N_",opt2,"_",dsub[1,2],".pdf",sep="")
ggsave(plot=g,file=fn1,width=7,height=5)
ggsave(plot=h,file=fn2,width=7,height=5)
ggsave(plot=k,file=fn3,width=7,height=5)


