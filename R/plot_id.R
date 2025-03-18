# Code to analyze output data files from SLiM simulations and generate plots
# Character with shifting identity

library(ggplot2)

# Generate a data matrix containing parameter combinations to be examined
opt.all=c(-10,-5,0,5,10) # Optima
l.C=c(1) # Number of loci affecting
par.comb.all=matrix(0,nrow=1,ncol=2)
for(i in 1:length(opt.all)){
  for(j in 1:length(l.C)){
    row=c(opt.all[i],l.C[j])
    par.comb.all=rbind(par.comb.all,row)
  }
}
par.comb.all=par.comb.all[2:nrow(par.comb.all),]
colnames(par.comb.all)=c("opt","l_C")
rownames(par.comb.all)=NULL

# Open files and extract information
# Plot results for non-WF populations
setwd("your_working_directory/id_wf")
prefix=paste("sim_id_out_end_",sep="") # Common part shared by all file names
out=matrix(0,nrow=nrow(par.comb.all),ncol=6)
for(c in 1:nrow(par.comb.all)){
  fn=paste(prefix,par.comb.all[c,1],"_",par.comb.all[c,2],".txt",sep="")
  d<-read.table(fn,sep="\t")
  out[c,1]=length(which(d[,5]==1))/nrow(d) # Proportion of populations where identity 2 is fixed
  out[c,2]=mean(d[,5]) # Mean C_2
  out[c,3]=mean(d[,2]) # Mean phenotype
  out[c,4]=sqrt(out[c,1]*(1-out[c,1])/nrow(d)) # SE of fixation proportion
  out[c,5]=sd(d[,5])/sqrt(nrow(d)) # SE of C_2
  out[c,6]=sd(d[,2])/sqrt(nrow(d)) # SE of z
}
colnames(out)=c("frac","c2","z","se_frac","se_c2","se_z") # Name the columns

d=data.frame(par.comb.all,out) # Data frame used for plotting
#d$l_C=factor(d$l_C,levels=sort(unique(d$l_C)),ordered=TRUE) # Factorize l_C for plotting (run only if there are >1 levels)
# Plot fixation fraction
#g<-ggplot(d,aes(x=opt,y=frac,group=l_C,col=l_C))
g<-ggplot(d,aes(x=opt,y=frac)) # l_C has only 1 level
g=g+geom_bar(stat="identity",color="grey",fill="grey")
g=g+geom_errorbar(aes(ymin=frac-se_frac,ymax=frac+se_frac,width=0.2))
g=g+theme_classic()+theme(axis.text=element_text(size=12))
g=g+xlab("")+ylab("")+ylim(c(0,1))
ggsave(plot=g,file="plot_frac.pdf",width=7,height=5)
# Plot C_2
#g<-ggplot(d,aes(x=opt,y=c2,group=l_C,col=l_C))
g<-ggplot(d,aes(x=opt,y=c2)) # l_C has only 1 level
g=g+geom_bar(stat="identity",color="grey",fill="grey")
ub=d$c2+d$se_c2;lb=d$c2-d$se_c2;ub[which(ub>1)]=1;lb[which(lb<0)]=0
g=g+geom_errorbar(aes(ymin=lb,ymax=ub,width=0.2))
g=g+theme_classic()+theme(axis.text=element_text(size=12))
g=g+xlab("")+ylab("")+ylim(c(0,1))
ggsave(plot=g,file="plot_c2.pdf",width=7,height=5)

# Plot results for non-WF populations
setwd("your_working_directory/id_nwf")
prefix=paste("sim_nwf_id_out_end_",sep="") # Common part shared by all file names
out=matrix(0,nrow=nrow(par.comb.all),ncol=8)
for(c in 1:nrow(par.comb.all)){
  fn=paste(prefix,par.comb.all[c,1],"_",par.comb.all[c,2],".txt",sep="")
  d<-read.table(fn,sep="\t")
  out[c,1]=length(which(d[,6]==1))/nrow(d) # Proportion of populations where identity 2 is fixed
  out[c,2]=mean(d[,6]) # Mean C_2
  out[c,3]=mean(d[,4]) # Mean phenotype
  out[c,4]=mean(d[,3]) # Mean population size
  out[c,5]=sqrt(out[c,1]*(1-out[c,1])/nrow(d)) # SE of fixation proportion
  out[c,6]=sd(d[,6])/sqrt(nrow(d)) # SE of C_2
  out[c,7]=sd(d[,4])/sqrt(nrow(d)) # SE of z
  out[c,8]=sd(d[,3])/sqrt(nrow(d)) # # SE of population size
}
colnames(out)=c("frac","c2","z","N","se_frac","se_c2","se_z","se_N") # Name the columns

d=data.frame(par.comb.all,out) # Data frame used for plotting

#d$l_C=factor(d$l_C,levels=sort(unique(d$l_C)),ordered=TRUE) # Factorize l_C for plotting (run only if there are >1 levels)
# Plot fixation fraction
#g<-ggplot(d,aes(x=opt,y=frac,group=l_C,col=l_C))
g<-ggplot(d,aes(x=opt,y=frac)) # l_C has only 1 level
g=g+geom_bar(stat="identity",color="grey",fill="grey")
ub=d$frac+d$se_frac;lb=d$frac-d$se_frac;ub[which(ub>1)]=1;lb[which(lb<0)]=0
g=g+geom_errorbar(aes(ymin=lb,ymax=ub,width=0.2))
g=g+theme_classic()+theme(axis.text=element_text(size=12))
g=g+xlab("")+ylab("")+ylim(c(0,1))
ggsave(plot=g,file="plot_frac_nwf.pdf",width=7,height=5)
# Plot C_2
#g<-ggplot(d,aes(x=opt,y=c2,group=l_C,col=l_C))
g<-ggplot(d,aes(x=opt,y=c2)) # l_C has only 1 level
g=g+geom_bar(stat="identity",color="grey",fill="grey")
ub=d$c2+d$se_c2;lb=d$c2-d$se_c2;ub[which(ub>1)]=1;lb[which(lb<0)]=0
g=g+geom_errorbar(aes(ymin=lb,ymax=ub,width=0.2))
g=g+theme_classic()+theme(axis.text=element_text(size=12))
g=g+xlab("")+ylab("")+ylim(c(0,1))
ggsave(plot=g,file="plot_c2_nwf.pdf",width=7,height=5)
# Plot N
h<-ggplot(d,aes(x=opt,y=N))
#h=h+geom_point(size=3)+geom_line(lwd=2)
h=h+geom_bar(stat="identity",fill="grey")
lb=d$N-d$se_N;lb[which(lb<0)]=0
h=h+geom_errorbar(aes(ymin=lb,ymax=N+se_N,width=0.4))
h=h+theme_classic()+theme(axis.text=element_text(size=12))
h=h+xlab("")+ylab("")
ggsave(plot=h,file="plot_N_nwf.pdf",width=7,height=5)

# Generate a bar plot for DFE of mutations affecting morphogen gradient
d<-read.table("dfe_id.txt",header=TRUE,sep="\t")
d$opt=factor(d$opt,levels=sort(unique(d$opt)),ordered=TRUE)
d$effect=factor(d$effect,levels=c("beneficial","neutral","deleterious"),ordered=TRUE)
g<-ggplot(d,aes(x=opt,y=frac,fill=effect))
g=g+geom_bar(stat="identity")+scale_fill_manual(values=c("orange","darkgrey","purple"))
g=g+theme_classic()+theme(axis.text=element_text(size=12),legend.position="none")
g=g+xlab("")+ylab("")
ggsave(plot=g,file="plot_dfe_id.pdf",width=7,height=2)