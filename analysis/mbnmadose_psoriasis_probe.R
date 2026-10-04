options(width=220)
suppressPackageStartupMessages(library(MBNMAdose))
for (nm in c("psoriasis75","psoriasis90","psoriasis100")) {
  data(list=nm, package="MBNMAdose")
  d <- get(nm)
  cat("\n===",nm,"===\n")
  cat("rows",nrow(d)," studies_raw",length(unique(d$studyID))," agents",length(unique(d$agent)),"\n")
  cat("agents:",paste(sort(unique(d$agent)),collapse=", "), "\n")
  cc <- d[complete.cases(d[,c("n","r","agent","studyID")]),]
  armn <- aggregate(agent~studyID,data=cc,FUN=length)
  cat("complete rows",nrow(cc)," studies",length(unique(cc$studyID)),
      " full_contrasts",sum(choose(armn$agent,2)),
      " rank_contrasts",sum(armn$agent-1),"\n")
  print(table(armn$agent))
  print(cc)
  write.csv(cc,paste0(nm,"_complete.csv"),row.names=FALSE)
}
sessionInfo()
