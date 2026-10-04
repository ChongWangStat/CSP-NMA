suppressPackageStartupMessages(library(netmeta))
covid_data <- data.frame(
  study=c("RECOVERY","RECOVERY","RECOVERY","RECOVERY","SOLIDARITY","SOLIDARITY","SOLIDARITY","SOLIDARITY","ACTT-1","ACTT-1","CoDEX","CoDEX","REMAP-CAP","REMAP-CAP"),
  t=c("A","B","C","D","A","C","D","E","A","E","A","B","A","B"),
  r=c(1110,482,353,396,440,104,148,301,75,59,91,85,33,78),
  n=c(4321,2104,1596,1561,4088,947,1399,2743,521,541,148,151,58,136)
)
pw <- pairwise(treat=t,event=r,n=n,studlab=study,data=covid_data,sm="OR",incr=0.5,method.incr="all",allstudies=TRUE)
fit <- netmeta(pw,common=TRUE,random=TRUE,reference.group="A")
cat("R_VERSION\t",R.version.string,"\n",sep="")
cat("NETMETA_VERSION\t",as.character(packageVersion("netmeta")),"\n",sep="")
cat("TAU\t",fit$tau,"\n",sep="")
cat("TAU2\t",fit$tau^2,"\n",sep="")
cat("Q\t",fit$Q,"\n",sep="")
cat("DF_Q\t",fit$df.Q,"\n",sep="")
cat("P_Q\t",fit$pval.Q,"\n",sep="")
write.csv(data.frame(tau=fit$tau,tau2=fit$tau^2,Q=fit$Q,df=fit$df.Q,p=fit$pval.Q),"covid_tau_probe.csv",row.names=FALSE)
