options(width = 240)
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

# The aggregate two-step H matrix acts on direct-comparison estimates formed
# from the SAME multi-arm-adjusted pairwise weights used to construct H.
aggregate_direct <- function(x, model="common") {
  comps <- x$comparisons
  se_adj <- if (model=="common") x$seTE.adj.common else x$seTE.adj.random
  out <- setNames(numeric(length(comps)), comps)
  out_se <- setNames(numeric(length(comps)), comps)
  for (comp in comps) {
    uv <- split_comp(comp)
    forward <- x$treat1 == uv[1] & x$treat2 == uv[2]
    reverse <- x$treat1 == uv[2] & x$treat2 == uv[1]
    idx <- which(forward | reverse)
    if (!length(idx)) stop("No pairwise rows for ", comp)
    sign <- ifelse(forward[idx], 1, -1)
    w <- 1 / se_adj[idx]^2
    out[comp] <- sum(w * sign * x$TE[idx]) / sum(w)
    out_se[comp] <- sqrt(1/sum(w))
  }
  list(TE=out, se=out_se)
}
agg <- aggregate_direct(nma, "common")
direct_TE <- agg$TE

get_mat_value <- function(M, comp) {
  uv <- split_comp(comp)
  as.numeric(M[uv[1], uv[2]])
}
fit_for_rows <- function(rows)
  setNames(vapply(rows, function(z) get_mat_value(nma$TE.common, z), numeric(1)), rows)

Hlong_obj <- hatmatrix(nma, method="Davies", type="long", common=TRUE, random=FALSE)
Hlong <- Hlong_obj$common
if (!identical(colnames(Hlong), names(direct_TE))) {
  if (!all(colnames(Hlong) %in% names(direct_TE))) stop("H columns cannot be aligned")
  direct_for_H <- direct_TE[colnames(Hlong)]
} else direct_for_H <- direct_TE

fit_H <- as.numeric(Hlong %*% direct_for_H)
names(fit_H) <- rownames(Hlong)
fit_H_truth <- fit_for_rows(rownames(Hlong))
audit_H <- data.frame(
  comparison=rownames(Hlong),
  fitted_NMA=as.numeric(fit_H_truth),
  signed_hat_reconstruction=fit_H,
  reconstruction_error=fit_H-as.numeric(fit_H_truth),
  stringsAsFactors=FALSE
)

# Four current netcontrib methods.
methods <- c("shortestpath","randomwalk","cccp","pseudoinverse")
nc <- lapply(methods, function(m)
  netcontrib(nma, method=m, common=TRUE, random=FALSE,
             study=(m=="shortestpath"), path=(m=="shortestpath")))
names(nc) <- methods

reconstruct_weights <- function(W, method) {
  if (!all(colnames(W) %in% names(direct_TE)))
    stop(method, ": contribution columns cannot be aligned to aggregate direct effects")
  d <- direct_TE[colnames(W)]
  recon <- as.numeric(W %*% d)
  fit <- fit_for_rows(rownames(W))
  data.frame(
    comparison=rownames(W),
    fitted_NMA=as.numeric(fit),
    contribution_weighted_direct=recon,
    reconstruction_error=recon-as.numeric(fit),
    method=method,
    stringsAsFactors=FALSE
  )
}
audits <- lapply(methods, function(m) reconstruct_weights(nc[[m]]$common,m))
audit_all <- do.call(rbind,audits)
rownames(audit_all) <- NULL

# Inspect the exact path-factor identity underlying the L2 and L1 methods.
Hfull <- hatmatrix(nma, method="Davies", type="full", common=TRUE, random=FALSE)$common
cm_pi <- netmeta:::contribution.matrix.ruecker.pseudoinv(nma, "common")
cm_l1 <- netmeta:::contribution.matrix.ruecker.cccp(nma, "common")

phi_resid <- function(cm,label) {
  vals <- numeric(length(cm$phi))
  estimator_resid <- numeric(length(cm$phi))
  # full direct vector with zero entries for missing edges
  yfull <- setNames(rep(0,ncol(Hfull)),colnames(Hfull))
  yfull[names(direct_TE)] <- direct_TE
  for (i in seq_along(cm$phi)) {
    ph <- cm$phi[[i]]
    Z <- cm$zlist[[i]]
    vals[i] <- max(abs(as.numeric(ph %*% Z) - Hfull[i,]))
    estimator_resid[i] <- as.numeric(ph %*% (Z %*% yfull) - Hfull[i,] %*% yfull)
  }
  data.frame(comparison=rownames(Hfull),method=label,
             max_coefficient_factorization_error=vals,
             estimator_factorization_error=estimator_resid)
}
phi_audit <- rbind(phi_resid(cm_l1,"cccp_L1_path_factor"),
                   phi_resid(cm_pi,"pseudoinverse_L2_path_factor"))

cat("\n=== AGGREGATE_DIRECT_EFFECTS_USED_BY_H ===\n")
print(data.frame(comparison=names(direct_TE),TE=as.numeric(direct_TE),
                 se=as.numeric(agg$se)),row.names=FALSE,digits=12)

cat("\n=== SIGNED_H_RECONSTRUCTION ===\n")
print(audit_H,row.names=FALSE,digits=12)
cat("max_abs_signed_H_error\t",format(max(abs(audit_H$reconstruction_error)),digits=16),"\n",sep="")

cat("\n=== FINAL_NONNEGATIVE_CONTRIBUTION_RECONSTRUCTION ===\n")
print(audit_all,row.names=FALSE,digits=12)
for (m in methods) {
  z <- audit_all[audit_all$method==m,]
  cat(m,"_max_abs_error\t",format(max(abs(z$reconstruction_error)),digits=16),"\n",sep="")
  cat(m,"_median_abs_error\t",format(median(abs(z$reconstruction_error)),digits=16),"\n",sep="")
}

cat("\n=== EXACT_PATH_FACTOR_CHECKS_FOR_L1_L2 ===\n")
print(phi_audit,row.names=FALSE,digits=12)
cat("max_phiZ_minus_H\t",format(max(phi_audit$max_coefficient_factorization_error),digits=16),"\n",sep="")
cat("max_phiZy_minus_Hy\t",format(max(abs(phi_audit$estimator_factorization_error)),digits=16),"\n",sep="")

cat("\n=== A:E_AND_B:E ===\n")
sel <- c(paste("A","E",sep=sep),paste("B","E",sep=sep))
print(audit_all[audit_all$comparison %in% sel,],row.names=FALSE,digits=12)

write.csv(data.frame(comparison=names(direct_TE),direct_TE=as.numeric(direct_TE),direct_se=as.numeric(agg$se)),
          "covid_aggregate_direct_effects_v4.csv",row.names=FALSE)
write.csv(audit_H,"covid_signed_hat_reconstruction_v4.csv",row.names=FALSE)
write.csv(audit_all,"covid_netcontrib_reconstruction_all_methods_v4.csv",row.names=FALSE)
write.csv(phi_audit,"covid_exact_path_factor_audit_v4.csv",row.names=FALSE)
for (m in methods)
  write.csv(as.data.frame(nc[[m]]$common),
            paste0("covid_netcontrib_",m,"_common_v4.csv"))
if (!is.null(nc$shortestpath$study.common))
  write.csv(nc$shortestpath$study.common,"covid_shortestpath_study_v4.csv",row.names=FALSE)
if (!is.null(nc$shortestpath$path.common))
  write.csv(nc$shortestpath$path.common,"covid_shortestpath_path_v4.csv",row.names=FALSE)

cat("\n=== SESSION INFO ===\n")
sessionInfo()
