##########################
# 期待値 W/2、分散 W^2/12 の一様分布からの i.i.d. 標本 が n個
# その平均を n_trial 回計算し、W/2を引く
##########################
# Load library
library(purrr)

# Set seed
set.seed(11111)

#############################################
# Simulation setting
#############################################
n = 100 # The count for a k-mer
W = 250
n_trial = 10000

#############################################
# Empirical null distribution of AUC scores
#############################################
# Sampling n times from a uniform distribution U(0,W) for `n_trial` times
rep(n, n_trial) %>% 
  map(function(x) mean(runif(x, min=0, max=W))) %>% 
  unlist() -> sampled_means

# Calculate AUC scores
sampled_aucs <- sampled_means - W/2

# Visualize
plot(density(sampled_aucs))
rug(jitter(sampled_aucs))

#############################################
# Theoretical null distribution of AUC scores
#############################################
# According to the central limit theorem, 
# when n is sufficiently large, the sample mean of U(0,W) - W/2
# approximately follows the normal distribution N(W/2-W/2,W^2/12/n)
x <- seq(min(sampled_aucs), max(sampled_aucs), by = .1)
y <- dnorm(x, mean = 0, sd = sqrt(W^2/12/n))
lines(x, y, col=3)

#############################################
# Compare the two null distributions
#############################################
# Difference of S.D. of AUC score
ratio <- sd(sampled_aucs)/sqrt(W^2/12/n)

pdf("pval_moccs_auc_density_plot.pdf")
plot(density(sampled_aucs), main = sprintf("ratio of S.D. of empirical over theoretical\n=%.4f", ratio))
rug(jitter(sampled_aucs))
x <- seq(min(sampled_aucs), max(sampled_aucs), by = .1)
y <- dnorm(x, mean = 0, sd = sqrt(W^2/12/n))
lines(x, y, col=3)
dev.off()

sink(file = "pval_moccs_auc_ratio.txt")
sprintf("ratio of S.D. of empirical over theoretical: %.4f", ratio)
sink()


