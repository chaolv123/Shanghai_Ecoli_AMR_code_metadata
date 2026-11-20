repos <- c(CRAN = "https://cloud.r-project.org")
required <- c("glmmTMB", "sandwich")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = repos)
stopifnot(vapply(required, requireNamespace, logical(1), quietly = TRUE))
cat("Required R packages are available:\n")
for (pkg in required) cat(pkg, as.character(packageVersion(pkg)), "\n")
