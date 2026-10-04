options(width=220)
suppressPackageStartupMessages({
  library(multinma)
  library(dplyr)
})

data("plaque_psoriasis_agd", package="multinma")
data("plaque_psoriasis_ipd", package="multinma")

ipd_arms <- plaque_psoriasis_ipd |>
  group_by(studyc, trtc) |>
  summarise(pasi75_r=sum(pasi75), pasi75_n=dplyr::n(), .groups="drop")

agd <- plaque_psoriasis_agd |>
  select(studyc, trtc, pasi75_r, pasi75_n)

dat <- bind_rows(agd, ipd_arms) |>
  arrange(studyc, trtc)

trial_arms <- dat |>
  count(studyc, name="n_arms")
Tn <- n_distinct(dat$trtc)
Mn <- n_distinct(dat$studyc)
full_contrasts <- sum(choose(trial_arms$n_arms,2))
rank_contrasts <- sum(trial_arms$n_arms-1)

cat("R_VERSION\t",R.version.string,"\n",sep="")
cat("MULTINMA_VERSION\t",as.character(packageVersion("multinma")),"\n",sep="")
cat("TRIALS\t",Mn,"\n",sep="")
cat("TREATMENTS\t",Tn,"\n",sep="")
cat("FULL_PAIRWISE_CONTRASTS\t",full_contrasts,"\n",sep="")
cat("RANK_CONTRASTS\t",rank_contrasts,"\n",sep="")
cat("\n=== ARM COUNT DISTRIBUTION ===\n")
print(table(trial_arms$n_arms))
cat("\n=== TREATMENTS ===\n")
print(sort(unique(dat$trtc)))
cat("\n=== TRIALS AND ARMS ===\n")
print(trial_arms,n=Inf)
cat("\n=== DATA ===\n")
print(dat,n=Inf)

write.csv(dat,"multinma_plaque_psoriasis_pasi75.csv",row.names=FALSE)
write.csv(trial_arms,"multinma_plaque_psoriasis_trial_arms.csv",row.names=FALSE)
sessionInfo()
