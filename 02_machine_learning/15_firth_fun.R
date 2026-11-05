source("config.R")  ## 统一路径配置，见 config.R`r`n## Firth 惩罚 logistic 的最小实现（供 15 / 17 共用）
## 单变量：X = [1, x]
loglik_star <- function(beta, X, y) {
  eta <- drop(X %*% beta)
  ll <- sum(y * eta - log1p(exp(eta)))
  p <- pmin(pmax(plogis(eta), 1e-10), 1 - 1e-10)
  W <- p * (1 - p)
  I <- crossprod(X * sqrt(W))
  detI <- determinant(I, logarithm = TRUE)
  if (!is.finite(as.numeric(detI$modulus))) return(-Inf)
  ll + 0.5 * as.numeric(detI$modulus)
}

firth_fit <- function(X, y, maxit = 300, tol = 1e-9) {
  beta <- rep(0, ncol(X)); ll_cur <- loglik_star(beta, X, y)
  for (it in seq_len(maxit)) {
    eta <- drop(X %*% beta); p <- plogis(eta); W <- p * (1 - p)
    I <- crossprod(X * sqrt(W)); Iinv <- solve(I)
    XW <- X * sqrt(W); h <- rowSums((XW %*% Iinv) * XW)
    U <- crossprod(X, y - p + h * (0.5 - p))
    delta <- drop(Iinv %*% U)
    k <- 1; ok <- FALSE
    for (s in 1:25) {
      bn <- beta + k * delta; ll_new <- loglik_star(bn, X, y)
      if (is.finite(ll_new) && ll_new >= ll_cur) { ok <- TRUE; break }
      k <- k / 2
    }
    if (!ok) break
    beta <- beta + k * delta
    if (max(abs(k * delta)) < tol) { ll_cur <- loglik_star(beta, X, y); break }
    ll_cur <- loglik_star(beta, X, y)
  }
  eta <- drop(X %*% beta); p <- plogis(eta); W <- p * (1 - p)
  V <- solve(crossprod(X * sqrt(W)))
  list(beta = beta, vcov = V, ll = ll_cur, fitted = p, se = sqrt(diag(V)))
}

firth_lrt <- function(fit, X, y) {
  Xr <- X[, 1, drop = FALSE]
  b0 <- 0; ll0 <- loglik_star(matrix(b0), Xr, y)
  for (it in 1:200) {
    eta <- drop(Xr %*% b0); p <- plogis(eta); W <- p * (1 - p)
    I <- crossprod(Xr * sqrt(W)); Iinv <- solve(I)
    XW <- Xr * sqrt(W); h <- rowSums((XW %*% Iinv) * XW)
    U <- crossprod(Xr, y - p + h * (0.5 - p))
    d <- drop(Iinv %*% U); b0 <- b0 + d
    if (max(abs(d)) < 1e-9) break
  }
  ll0 <- loglik_star(matrix(b0), Xr, y)
  stat <- 2 * (fit$ll - ll0)
  list(stat = stat, p = pchisq(stat, ncol(X) - 1, lower.tail = FALSE))
}

firth_profile_ci <- function(fit, X, y, level = 0.95) {
  crit <- qchisq(level, 1); target <- fit$ll - crit / 2
  prof <- function(b) {
    obj <- function(b0) -loglik_star(c(b0, b), X, y)
    v <- tryCatch(-optimize(obj, interval = c(-60, 60))$objective, error = function(e) NA_real_)
    v
  }
  se <- fit$se[2]; bhat <- fit$beta[2]
  ## 逐步向外扩展，直到剖面值跌破阈值
  grid_lo <- seq(bhat, bhat - 12 * se, length.out = 120)
  grid_hi <- seq(bhat, bhat + 12 * se, length.out = 120)
  lo <- NA_real_; hi <- NA_real_
  find_edge <- function(grid) {
    v <- sapply(grid, prof) - target
    v[!is.finite(v)] <- -Inf
    bad <- which(v < 0)
    if (!length(bad)) return(NA_real_)
    i <- min(bad)
    if (i == 1) return(grid[1])
    (grid[i - 1] + grid[i]) / 2
  }
  c(lo = find_edge(grid_lo), hi = find_edge(grid_hi))
}

firth_report <- function(case, ctrl, label) {
  x <- c(case, ctrl); y <- c(rep(1, length(case)), rep(0, length(ctrl)))
  X <- cbind(1, x)
  g <- suppressWarnings(glm(y ~ x, family = binomial))
  f <- firth_fit(X, y)
  lrt <- firth_lrt(f, X, y)
  ci <- firth_profile_ci(f, X, y)
  data.frame(cohort = label, n = length(y), events = sum(y),
             glm_OR = exp(coef(g)[2]),
             firth_OR = exp(f$beta[2]),
             lo_wald = exp(f$beta[2] - 1.96 * f$se[2]),
             hi_wald = exp(f$beta[2] + 1.96 * f$se[2]),
             lo_prof = exp(ci["lo"]), hi_prof = exp(ci["hi"]),
             LRT_p = lrt$p, stringsAsFactors = FALSE)
}
