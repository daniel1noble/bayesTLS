# Build the committed phylogenetic test fixture:
#   tests/testthat/fixtures/phylo_truth.rds
#
# Contents: an 8-species tree-structured CORRELATION matrix A (unit diagonal)
# plus the fixed "true" species CTmax deviations drawn once from
# MVN(0, sigma^2 * A). Committing the draw makes the gated recovery tests
# (tests/testthat/test-phylo-fit.R) deterministic without an `ape` dependency:
# A is constructed by hand as shared-branch-length correlations on a balanced
# ultrametric tree (root edge 0.2 -> clade split 0.6 -> sister split 0.85 ->
# tips at 1), which is exactly what ape::vcv(phy, corr = TRUE) would return
# for that tree. Regenerate only if the fixture design changes.
#
# Run from the package root:  Rscript data-raw/make_phylo_test_fixture.R

sp <- paste0("sp", 1:8)

# Shared heights: sisters 0.85, within-clade cousins 0.6, across clades 0.2.
A <- matrix(0.2, 8, 8, dimnames = list(sp, sp))
for (block in list(1:4, 5:8)) A[block, block] <- 0.6
for (pair in list(1:2, 3:4, 5:6, 7:8)) A[pair, pair] <- 0.85
diag(A) <- 1
stopifnot(isSymmetric(A), all(eigen(A, only.values = TRUE)$values > 0))

# One fixed draw of species CTmax deviations (degrees C) with real
# phylogenetic signal.
sigma_phylo <- 1.5
set.seed(20260928)
dev <- drop(chol(A) %*% rnorm(8)) * sigma_phylo
names(dev) <- sp

truth <- list(
  A            = A,
  species      = sp,
  ctmax_global = 41,     # degrees C at the 1 h reference
  z            = 4,      # common across species in the fixture
  sigma_phylo  = sigma_phylo,
  ctmax_dev    = dev,
  k            = 6,      # 4PL steepness (log10-time scale)
  low          = 0.01,
  up           = 0.99
)

dir.create("tests/testthat/fixtures", showWarnings = FALSE, recursive = TRUE)
saveRDS(truth, "tests/testthat/fixtures/phylo_truth.rds")
cat("Wrote tests/testthat/fixtures/phylo_truth.rds\n")
print(round(dev, 3))
