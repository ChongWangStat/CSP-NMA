options(width = 220)
suppressPackageStartupMessages(library(netmeta))

covid_data <- data.frame(
  study = c("RECOVERY","RECOVERY","RECOVERY","RECOVERY",
            "SOLIDARITY","SOLIDARITY","SOLIDARITY","SOLIDARITY",
            "ACTT-1","ACTT-1","CoDEX","CoDEX","REMAP-CAP","REMAP-CAP"),
  t = c("A","B","C","D","A","C","D","E","A","E","A","B","A","B"),
  r = c(1110,482,353,396,440,104,148,301,75,59,91,85,33,78),
  n = c(4321,2104,1596,1561,4088,947,1399,2743,521,541,148,151,58,136),
  stringsAsFactors = FALSE
)

cat("R_VERSION\t", R.version.string, "\n", sep="")
cat("NETMETA_VERSION\t", as.character(packageVersion("netmeta")), "\n", sep="")

pw <- pairwise(treat=t, event=r, n=n, studlab=study, data=covid_data,
               sm="OR", incr=0.5, method.incr="all", allstudies=TRUE)
nma <- netmeta(pw, common=TRUE, random=FALSE, reference.group="A")

cat("\n=== PAIRWISE_ROWS ===\n")
pair_df <- data.frame(idx=seq_along(pw$TE), studlab=pw$studlab, treat1=pw$treat1,
                      treat2=pw$treat2, TE=pw$TE, seTE=pw$seTE)
print(pair_df, row.names=FALSE, digits=10)

cat("\n=== TE_COMMON ===\n")
print(nma$TE.common, digits=10)

cat("\n=== NETMETA_COMPARISON_VECTORS ===\n")
for (nm in c("comparison","TE.direct.common","seTE.direct.common","TE.common","seTE.common")) {
  if (!is.null(nma[[nm]])) {
    cat("--", nm, "--\n")
    print(nma[[nm]], digits=10)
  }
}

cat("\n=== HAT_MATRIX_FULL_SHORTESTPATH ===\n")
hsp <- hatmatrix(nma, method="shortestpath", type="full")
print(names(hsp))
if (!is.null(hsp$common)) print(hsp$common, digits=10) else print(hsp, digits=10)

cat("\n=== HAT_MATRIX_FULL_RANDOMWALK ===\n")
hrw <- hatmatrix(nma, method="randomwalk", type="full")
print(names(hrw))
if (!is.null(hrw$common)) print(hrw$common, digits=10) else print(hrw, digits=10)

cat("\n=== NETCONTRIB_SHORTESTPATH ===\n")
ncs <- netcontrib(nma, method="shortestpath", common=TRUE, random=FALSE,
                  study=TRUE, path=TRUE)
str(ncs, max.level=2)
cat("\n-- common --\n")
print(ncs$common, digits=10)
cat("\n-- study.common --\n")
print(ncs$study.common, digits=10)
cat("\n-- path.common --\n")
print(ncs$path.common, digits=10)

cat("\n=== NETCONTRIB_RANDOMWALK ===\n")
ncr <- netcontrib(nma, method="randomwalk", common=TRUE, random=FALSE)
str(ncr, max.level=2)
cat("\n-- common --\n")
print(ncr$common, digits=10)

cat("\n=== DIMNAMES ===\n")
cat("shortest common\n"); print(dimnames(ncs$common))
cat("shortest study.common\n"); print(dimnames(ncs$study.common))
cat("shortest path.common\n"); print(dimnames(ncs$path.common))
cat("randomwalk common\n"); print(dimnames(ncr$common))

# Attempt a transparent reconstruction audit when netcontrib columns align with
# aggregate direct comparison estimates. We only compute where exact matching
# names are available; no silent matching.
cat("\n=== RECONSTRUCTION_AUDIT_CANDIDATES ===\n")
cat("The following prints are diagnostic. Final manuscript values must use an explicitly documented mapping.\n")

write.csv(pair_df, "covid_pairwise_rows.csv", row.names=FALSE)
write.csv(as.data.frame(ncs$common), "covid_netcontrib_shortest_common.csv")
write.csv(as.data.frame(ncr$common), "covid_netcontrib_randomwalk_common.csv")
saveRDS(list(pw=pw,nma=nma,ncs=ncs,ncr=ncr,hsp=hsp,hrw=hrw),
        "covid_netmeta_audit_v2.rds")
sessionInfo()
