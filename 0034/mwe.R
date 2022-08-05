library(iterpc)
library(brunnermunzel)
library(testthat)

#normal Wilcoxon test statistic
calcW <- function(x, y) {
  m <- length(x)
  pooledSample <- c(x, y)
  ranks <- rank(pooledSample)
  ranksX <- ranks[1:m]
  return(sum(ranksX))
}

#Brunner-Munzel test statistic
BrunnerMunzelStat <- function(x, y) {
  ns <- c(length(x), length(y))
  N <- sum(ns)
  pooledSample <- c(x, y)
  ranks <- rank(pooledSample)
  avRanks <- rep(NA, 2)
  avRanks[1] <- mean(ranks[1:ns[1]])
  avRanks[2] <- mean(ranks[(ns[1] + 1):N])
  
  pHat <- (avRanks[2] - avRanks[1])
  
  
  # sigmas
  sigmas <- rep(NA, 2)
  for (i in 1:2) {
    constant <- (ns[i] + 1) / 2
    if (i == 1) {
      cRanks <- rank(x)
      offset <- 0
    } else {
      cRanks <- rank(y)
      offset <- ns[1]
    }
    sumsOfSquares <- 0
    
    for (k in 1:ns[i]) {
      sumsOfSquares <- sumsOfSquares + (ranks[offset + k] - cRanks[k] - avRanks[i] + constant)^2
    }
    sigmas[i] <- 1 / (ns[i] - 1) * sumsOfSquares
  }
  
  Vn <- sqrt(N) * sqrt(1 / ns[2] * sigmas[1] + 1 / ns[1] * sigmas[2])
  
  
  Tn <- pHat / Vn * sqrt(prod(ns) / N)
  return(Tn)
}




perm_test <- function(x, y, testStat) {
  W <- testStat(x, y)
  allValues <- c(x, y)
  groupMemberShips <- c(rep(TRUE, length(x)), rep(FALSE, length(y)))
  
  possiblePerms <- getall(iterpc(table(groupMemberShips), order = TRUE))
  n1 <- length(x)
  n2 <- length(y)
  permutations <- choose(n1 + n2, n1)
  hypothetical <- rep(NA, permutations)
  permIter <- function(i) c(TRUE, FALSE)[possiblePerms[i, ]]
  for (i in 1:permutations) {
    curGroupMembersShips <- permIter(i)
    curGroup1 <- allValues[curGroupMembersShips]
    curGroup2 <- allValues[!curGroupMembersShips]
    hypothetical[i] <- testStat(curGroup1, curGroup2)
  }
  pValSmaller <- sum(W<=hypothetical)/permutations
  pValGreater <- sum(W>=hypothetical)/permutations
  pVal <- min(pValSmaller,pValGreater)*2
  return(list(W = W, pVal = pVal))
}

#test statistics are the same
reps <- 10^3
for (i in 1:reps){
  x <- rnorm(5)
  y <- rnorm(5)
  res_my <- BrunnerMunzelStat(x,y)
  res_bm <- brunnermunzel.test(x,y)
  testthat::expect_equivalent(res_bm$statistic,res_my)
}

# my permutation test seems to work. It leads to the same results
# as the wilcoxon test if used with the W test statistic
reps <- 100
p_vals_R <- rep(NA,reps)
p_vals_me <- rep(NA,reps)
for (i in 1:reps){
  x <- rnorm(5)
  y <- rnorm(5)
  p_vals_R[i] <- wilcox.test(x,y, exact = TRUE)$p.value
  p_vals_me[i] <- perm_test(x,y, calcW)$pVal
}
testthat::expect_equivalent(p_vals_R, p_vals_me)

## but using the Brunner Munzel test statistic does not always lead to the same results
# as the permutation test from the package. 20 out of 100 results doe not agree
set.seed(1290)
reps <- 100
p_val_jul <- p_val_bm <- p_val_lawstat <- p_val_bm_perm <- numeric(reps)
for (i in seq_len(reps)){
  x <- rnorm(5)
  y <- rnorm(5)
  
  res_jul <- perm_test(x,y,BrunnerMunzelStat)
  res_bm_perm <- brunnermunzel.permutation.test(x,y)
  
  p_val_jul[i] <- res_jul$pVal
  p_val_bm_perm[i] <- res_bm_perm$p.value
}
testthat::expect_equivalent(p_val_jul, p_val_bm_perm)



