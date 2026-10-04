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

sep <- nma$sep.trts
split_comp <- function(x) strsplit(x, sep, fixed=TRUE)[[1]]

get_mat_value <- function(M, comp) {
  uv <- split_comp(comp)
  if (length(uv) != 2L) stop("Cannot parse comparison: ", comp)
  if (!(uv[1] %in% rownames(M)) || !(uv[2] %in% colnames(M)))
    stop("Comparison treatments not found: ", comp)
  as.numeric(M[uv[1], uv[2]])
}

# Fitted NMA effects on the exact orientation used by netmeta rows.
all_comps <- combn(nma$trts, 2, FUN=function(z) paste(z, collapse=sep))
fit_vec <- setNames(vapply(all_comps, function(z) get_mat_value(nma$TE.common, z), numeric(1)), all_comps)

# Aggregate two-step hat matrix: signed coefficients on aggregate direct effects.
Hlong_obj <- hatmatrix(nma, method="Davies", type="long", common=TRUE, random=FALSE)
Hlong <- Hlong_obj$common
direct_comps <- colnames(Hlong)
direct_TE <- setNames(vapply(direct_comps, function(z) get_mat_value(nma$TE.direct.common, z), numeric(1)), direct_comps)
fit_from_H <- as.numeric(Hlong %*% direct_TE)
names(fit_from_H) <- rownames(Hlong)

# Nonnegative contribution matrices shown to users.
ncs <- netcontrib(nma, method="shortestpath", common=TRUE, random=FALSE,
                  study=TRUE, path=TRUE)
ncr <- netcontrib(nma, method="randomwalk", common=TRUE, random=FALSE)

reconstruct_from_contrib <- function(W, label) {
  if (!identical(colnames(W), direct_comps)) {
    # align by names, but refuse silent omissions
    if (!all(colnames(W) %in% names(direct_TE)))
      stop(label, ": contribution columns do not map to direct comparisons")
    d <- direct_TE[colnames(W)]
  } else {
    d <- direct_TE
  }
  vals <- as.numeric(W %*% d)
  names(vals) <- rownames(W)
  vals
}

sp_recon <- reconstruct_from_contrib(ncs$common, "shortestpath")
rw_recon <- reconstruct_from_contrib(ncr$common, "randomwalk")

# Align fitted estimates to output row orientations.
fit_for_rows <- function(rows) setNames(vapply(rows, function(z) get_mat_value(nma$TE.common, z), numeric(1)), rows)
fit_sp <- fit_for_rows(rownames(ncs$common))
fit_rw <- fit_for_rows(rownames(ncr$common))
fit_h <- fit_for_rows(rownames(Hlong))

audit_sp <- data.frame(
  comparison=rownames(ncs$common),
  fitted_NMA=as.numeric(fit_sp),
  contribution_weighted_direct=as.numeric(sp_recon),
  reconstruction_error=as.numeric(sp_recon-fit_sp),
  row_sum=rowSums(ncs$common),
  method="shortestpath",
  stringsAsFactors=FALSE
)
audit_rw <- data.frame(
  comparison=rownames(ncr$common),
  fitted_NMA=as.numeric(fit_rw),
  contribution_weighted_direct=as.numeric(rw_recon),
  reconstruction_error=as.numeric(rw_recon-fit_rw),
  row_sum=rowSums(ncr$common),
  method="randomwalk",
  stringsAsFactors=FALSE
)
audit_h <- data.frame(
  comparison=rownames(Hlong),
  fitted_NMA=as.numeric(fit_h),
  signed_hat_reconstruction=as.numeric(fit_from_H),
  reconstruction_error=as.numeric(fit_from_H-fit_h),
  stringsAsFactors=FALSE
)

cat("\n=== SIGNED_AGGREGATE_HAT_RECONSTRUCTION ===\n")
print(audit_h, row.names=FALSE, digits=12)
cat("\nmax_abs_signed_hat_error\t", format(max(abs(audit_h$reconstruction_error)), digits=16), "\n", sep="")

cat("\n=== SHORTESTPATH_REPORTED_CONTRIBUTION_RECONSTRUCTION ===\n")
print(audit_sp, row.names=FALSE, digits=12)
cat("\nmax_abs_shortestpath_error\t", format(max(abs(audit_sp$reconstruction_error)), digits=16), "\n", sep="")
cat("median_abs_shortestpath_error\t", format(median(abs(audit_sp$reconstruction_error)), digits=16), "\n", sep="")

cat("\n=== RANDOMWALK_REPORTED_CONTRIBUTION_RECONSTRUCTION ===\n")
print(audit_rw, row.names=FALSE, digits=12)
cat("\nmax_abs_randomwalk_error\t", format(max(abs(audit_rw$reconstruction_error)), digits=16), "\n", sep="")
cat("median_abs_randomwalk_error\t", format(median(abs(audit_rw$reconstruction_error)), digits=16), "\n", sep="")

cat("\n=== A:E AND B:E DETAIL ===\n")
sel <- c(paste("A","E",sep=sep), paste("B","E",sep=sep))
print(rbind(audit_sp[audit_sp$comparison %in% sel,],
            audit_rw[audit_rw$comparison %in% sel,]), row.names=FALSE, digits=12)

cat("\n=== CONTRIBUTION MATRICES ===\n")
cat("-- shortestpath common --\n"); print(ncs$common, digits=8)
cat("-- randomwalk common --\n"); print(ncr$common, digits=8)

cat("\n=== DIRECT EFFECT VECTOR USED IN RECONSTRUCTION ===\n")
print(direct_TE, digits=12)

cat("\n=== SHORTESTPATH STUDY CONTRIBUTION OBJECT ===\n")
print(ncs$study.common, digits=8)
cat("\n=== SHORTESTPATH PATH CONTRIBUTION OBJECT ===\n")
print(ncs$path.common, digits=8)

write.csv(audit_h, "covid_signed_hat_reconstruction.csv", row.names=FALSE)
write.csv(audit_sp, "covid_shortestpath_reconstruction.csv", row.names=FALSE)
write.csv(audit_rw, "covid_randomwalk_reconstruction.csv", row.names=FALSE)
write.csv(data.frame(comparison=names(direct_TE), direct_TE=as.numeric(direct_TE)),
          "covid_aggregate_direct_effects.csv", row.names=FALSE)
write.csv(as.data.frame(ncs$common), "covid_netcontrib_shortest_common.csv")
write.csv(as.data.frame(ncr$common), "covid_netcontrib_randomwalk_common.csv")
if (!is.null(ncs$study.common)) write.csv(ncs$study.common, "covid_netcontrib_shortest_study.csv", row.names=FALSE)
if (!is.null(ncs$path.common)) write.csv(ncs$path.common, "covid_netcontrib_shortest_path.csv", row.names=FALSE)
saveRDS(list(pw=pw,nma=nma,Hlong=Hlong_obj,ncs=ncs,ncr=ncr),
        "covid_netmeta_audit_v3.rds")

cat("\n=== SESSION INFO ===\n")
sessionInfo()
