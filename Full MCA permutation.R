library(FactoMineR)

# --- SETUP ----
xdata <- read.delim("xdata.txt", header=TRUE, stringsAsFactors = TRUE)
xdata <- xdata[, -1]  # remove recording name column for clarity

id_orig <- xdata$Combination    # original call-type labels

# Compute weights as in original analysis
tab <- table(id_orig)
weights <- 1 / (length(tab) * tab[id_orig])
weights <- as.numeric(weights)

# --- 1. OBSERVED MCA AND INERTIA ----

library("FactoMineR")
res_obs <- MCA(xdata, 
               quali.sup = 1:4, 
               quanti.sup = c(5,6),
               row.w = weights,
               ncp = 5)

coords_obs <- res_obs$ind$coord[, 1:5]
w_obs <- res_obs$call$row.w

# compute observed barycenters
unique_ids <- unique(id_orig)
n_groups <- length(unique_ids)
k <- ncol(coords_obs)

global_center_obs <- colSums(coords_obs * w_obs) / sum(w_obs)

bary_obs <- matrix(NA, n_groups, k)
group_weights <- numeric(n_groups)

for (g in seq_len(n_groups)) {
  grp <- unique_ids[g]
  sel <- which(id_orig == grp)
  wg <- w_obs[sel]
  group_weights[g] <- sum(wg)
  bary_obs[g, ] <- colSums(coords_obs[sel, , drop = FALSE] * wg) / sum(wg)
}

I_between_obs <- sum(group_weights * rowSums((bary_obs - 
                                                matrix(global_center_obs, n_groups, k, byrow = TRUE))^2))

I_total_obs <- sum(w_obs * rowSums((coords_obs - 
                                      matrix(global_center_obs, nrow(coords_obs), k, byrow = TRUE))^2))

R_obs <- I_between_obs / I_total_obs
cat("Observed R:", R_obs, "\n")

# --- 2. FULL PERMUTATION TEST ----

nperm <- 10000
R_perm <- numeric(nperm)

set.seed(123)

for (p in 1:nperm) {
  
  # permute labels
  id_perm <- sample(id_orig)
  
  # rebuild data with permuted Combination labels
  xperm <- xdata
  xperm$Combination <- id_perm
  
  # recompute MCA
  res_perm <- MCA(xperm, 
                  quali.sup = 1:4, 
                  quanti.sup = c(5,6),
                  row.w = weights,
                  ncp = 5)
  
  coords <- res_perm$ind$coord[, 1:5]
  w <- res_perm$call$row.w
  
  # recompute global center
  global_center <- colSums(coords * w) / sum(w)
  
  # recompute barycenters
  bary <- matrix(NA, n_groups, k)
  group_weights <- numeric(n_groups)
  
  for (g in seq_len(n_groups)) {
    grp <- unique_ids[g]
    sel <- which(id_perm == grp)
    wg <- w[sel]
    group_weights[g] <- sum(wg)
    bary[g, ] <- colSums(coords[sel, , drop = FALSE] * wg) / sum(wg)
  }
  
  I_between <- sum(group_weights * rowSums((bary - 
                                              matrix(global_center, n_groups, k, byrow = TRUE))^2))
  
  I_total <- sum(w * rowSums((coords - 
                                matrix(global_center, nrow(coords), k, byrow = TRUE))^2))
  
  R_perm[p] <- I_between / I_total
}

# --- 3. P-VALUE ----
pval <- mean(R_perm >= R_obs)
cat("Full MCA permutation p-value:", pval, "\n")
