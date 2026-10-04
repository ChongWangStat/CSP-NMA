options(width=220)
suppressPackageStartupMessages({
  library(NMAforest)
  library(netmeta)
  library(igraph)
})

cat("R_VERSION\t", R.version.string, "\n", sep="")
cat("NMAFOREST_VERSION\t", as.character(packageVersion("NMAforest")), "\n", sep="")
cat("NETMETA_VERSION\t", as.character(packageVersion("netmeta")), "\n", sep="")
cat("IGRAPH_VERSION\t", as.character(packageVersion("igraph")), "\n", sep="")

data(example_data, package="NMAforest")

res <- NMAforest(
  data=example_data,
  sm="OR",
  reference="x",
  model="random",
  comparison=c("x","y"),
  study="study",
  treat="t",
  event="r",
  N="n",
  study_id="id",
  study_path=FALSE
)

out <- res$output
cat("\n=== NMAFOREST_BINARY_OUTPUT ===\n")
print(out, row.names=FALSE)

getnum <- function(label, col) {
  x <- out[out$Label == label, col]
  if (length(x) != 1L) stop("Expected one row for ", label)
  as.numeric(x)
}

theta_D <- getnum("Overall Direct Effect","EffectSize")
p_D <- getnum("Overall Direct Effect","proportion")
theta_I <- getnum("Overall Indirect Effect","EffectSize")
p_I <- getnum("Overall Indirect Effect","proportion")
theta_N <- getnum("Overall NMA Effect","EffectSize")

path_rows <- grepl("^Path ", out$Label)
path_theta <- as.numeric(out$EffectSize[path_rows])
path_p <- as.numeric(out$proportion[path_rows])

recon_DI <- p_D*theta_D + p_I*theta_I
recon_paths <- p_D*theta_D + sum(path_p*path_theta)
err_DI <- recon_DI - theta_N
err_paths <- recon_paths - theta_N

# Conservative upper bound on discrepancy caused solely by 0.001 rounding
# of each displayed effect and proportion, plus 0.001 rounding of the NMA row.
round_bound <- function(theta, p) {
  sum(abs(p)*0.0005 + abs(theta)*0.0005 + 0.0005*0.0005) + 0.0005
}
rb_DI <- round_bound(c(theta_D,theta_I), c(p_D,p_I))
rb_paths <- round_bound(c(theta_D,path_theta), c(p_D,path_p))

summary <- data.frame(
  reconstruction=c("direct_plus_indirect","direct_plus_reported_paths"),
  reported_NMA=c(theta_N,theta_N),
  reconstructed=c(recon_DI,recon_paths),
  error=c(err_DI,err_paths),
  conservative_rounding_bound=c(rb_DI,rb_paths),
  abs_error_over_rounding_bound=c(abs(err_DI)/rb_DI,abs(err_paths)/rb_paths)
)

cat("\n=== RECONSTRUCTION_AUDIT ===\n")
print(summary, row.names=FALSE, digits=12)

write.csv(out,"nmaforest_binary_output.csv",row.names=FALSE)
write.csv(summary,"nmaforest_binary_reconstruction.csv",row.names=FALSE)
write.csv(example_data,"nmaforest_binary_example_data.csv",row.names=FALSE)

cat("\n=== SESSION INFO ===\n")
sessionInfo()
