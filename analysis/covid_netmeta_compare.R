options(width = 200)
suppressPackageStartupMessages({
  library(netmeta)
})
source("csp_functions.R")

covid_data <- data.frame(
  study = c(
    "RECOVERY","RECOVERY","RECOVERY","RECOVERY",
    "SOLIDARITY","SOLIDARITY","SOLIDARITY","SOLIDARITY",
    "ACTT-1","ACTT-1",
    "CoDEX","CoDEX",
    "REMAP-CAP","REMAP-CAP"
  ),
  id = c(1,1,1,1,2,2,2,2,3,3,4,4,5,5),
  t = c("A","B","C","D","A","C","D","E","A","E","A","B","A","B"),
  r = c(1110,482,353,396,440,104,148,301,75,59,91,85,33,78),
  n = c(4321,2104,1596,1561,4088,947,1399,2743,521,541,148,151,58,136),
  stringsAsFactors = FALSE
)

cat("R VERSION:", R.version.string, "\n")
cat("NETMETA VERSION:", as.character(packageVersion("netmeta")), "\n")

# Exact CSP analysis used in the preprint example.
csp <- fit_csp_nma(covid_data, model = "fixed")

# Use the same arm data and fixed-effect log-odds-ratio analysis in netmeta.
# method.incr='all' matches the CSP 0.5 continuity-correction convention
# (no arm in this dataset has zero events, so it does not alter these data).
pw <- pairwise(
  treat = t, event = r, n = n, studlab = study,
  data = covid_data, sm = "OR",
  incr = 0.5, method.incr = "all",
  allstudies = TRUE
)
nma <- netmeta(pw, common = TRUE, random = FALSE, reference.group = "A")

cat("\n=== PAIRWISE DATA ===\n")
print(data.frame(studlab=pw$studlab, treat1=pw$treat1, treat2=pw$treat2,
                 TE=pw$TE, seTE=pw$seTE), row.names=FALSE)

cat("\n=== CSP THETA ===\n")
print(csp$theta_hat)

cat("\n=== NETMETA TE.COMMON MATRIX ===\n")
print(nma$TE.common)

cat("\n=== NETMETA COMPARISONS / DIRECT ESTIMATES ===\n")
print(data.frame(comparison=nma$comparison, TE.direct.common=nma$TE.direct.common,
                 seTE.direct.common=nma$seTE.direct.common), row.names=FALSE)

cat("\n=== DAVIES FULL H COMMON ===\n")
hm <- hatmatrix(nma, method = "davies", type = "full")
print(names(hm))
if (!is.null(hm$common)) print(hm$common) else print(hm)

cat("\n=== NETCONTRIB SHORTESTPATH ===\n")
nc_s <- netcontrib(nma, method = "shortestpath", common = TRUE, random = FALSE,
                   study = TRUE, path = TRUE)
print(names(nc_s))
print(nc_s$common)
cat("\n--- study.common ---\n")
print(nc_s$study.common)
cat("\n--- path.common ---\n")
print(nc_s$path.common)

cat("\n=== NETCONTRIB RANDOMWALK ===\n")
nc_r <- netcontrib(nma, method = "randomwalk", common = TRUE, random = FALSE)
print(names(nc_r))
print(nc_r$common)

cat("\n=== NETMETA OBJECT RELEVANT NAMES ===\n")
print(grep("direct|TE|H|hat|comparison", names(nma), value=TRUE, ignore.case=TRUE))

cat("\n=== NETCONTRIB DIMNAMES ===\n")
print(dimnames(nc_s$common))
print(dimnames(nc_r$common))

# Save key objects for a second-stage audit if needed.
saveRDS(list(csp=csp, pw=pw, nma=nma, hm=hm, nc_s=nc_s, nc_r=nc_r),
        file="covid_netmeta_audit_objects.rds")
