#' Predict p-Value from Test Statistic
#'
#' @param bh.test The value of the test statistic of the Bayer Hanck Test
#' @param trendtype Type of deterministic component
#' @param test.type Selection of tests to choose from
#' @param k degrees of freedom - 1
#'
#' @noRd
#'
get_p_value <- function(bh.test, trendtype, test.type, k, ...){
  get_lambda <- function(data, case_w, art, test_w){
    data %>%
      dplyr::filter(case == case_w,
                    test.type == test_w) %>%
      dplyr::select(dplyr::all_of(art)) %>%
      dplyr::pull()
  }

  get_model <- function(trendtype, test){
    models %>%
      dplyr::filter(case == trendtype,
                    test.type == test) %>%
      dplyr::select(model) %>%
      dplyr::pull() %>%
      purrr::pluck(1)
  }

  get_critical_val <- function(trendtype, k_s, test){
    models %>%
      dplyr::filter(case == trendtype,
                    test.type == test) %>%
      tidyr::unnest(critical) %>%
      dplyr::filter(k == k_s) %>%
      dplyr::pull(crit_val)
  }

  get_p_trans <- function(model){
    model %>%
      purrr::pluck('terms') %>%
      purrr::pluck(2)
  }

  # getting critical val
  crit_val <- get_critical_val(trendtype, k, test.type)

  if (crit_val <= bh.test) {
    p.value <- 1e-12
  }else{
    # getting lambda for p_values
    lambda_p <- get_lambda(models, trendtype, 'p', 'all')
    # getting lambda for stat
    lambda_stat <- get_lambda(models, trendtype, 'stat', test.type)
    # saving model
    model <- get_model(trendtype, test.type)
    # dependent var
    dep_var <- get_p_trans(model) %>% as.character()

    # generating data set
    new_data <- tibble::tibble(dep = 1L,
                       stat_Fisher_all_bc = ((bh.test^lambda_stat)-1)/lambda_stat,
                       stat_Fisher_E_J_bc = ((bh.test^lambda_stat)-1)/lambda_stat,
                       k_dummy = as.factor(k))
    colnames(new_data)[1] <- dep_var
    # approximation of the model
    p.value_raw <- as.vector(stats::model.matrix(object = model, data = new_data) %*% coef(model))

    p.value_trans <- if(stringr::str_detect(dep_var, '_bc')){
      Re((lambda_p*as.complex(p.value_raw) + 1)^(1/lambda_p))
    } else {p.value_raw}

    p.value <- ifelse(p.value_trans >= 1, 9.9999e-1, ifelse(p.value_trans <= 0, 1e-12, p.value_trans))
  }

  p.value <- ifelse(p.value <= 1e-12, paste(c("<", 1e-12), collapse = " "),
                    ifelse(p.value <= 1e-4, format(p.value, scientific = TRUE, digits = 4),
                           format(round(p.value, 4), scientific = FALSE)))

  return(p.value)
}

#' Johansen ML estimation of a VECM
#'
#' Reproduces VECM(x, lag = p, include = ..., estim = "ML") from tsDyn package.
#'
#' @param x       numeric matrix or data.frame of levels, one column per series.
#' @param lag     number of lagged differences p (VAR order is p + 1).
#' @param include deterministic terms entered *unrestricted* in the short-run
#'   part. "both" means constant and trend. A trend without a constant is not
#'   a standard Johansen case and is not offered.
#' @param LRinclude deterministic terms *restricted* to the cointegration space.
#' @param r       cointegration rank used for the coefficient estimates.
#'
#' @return list with eigenvalues, trace / max-eigenvalue statistics, beta,
#'   alpha, Pi, the full coefficient matrix, residuals, the log-likelihood and
#'   \code{npar}, the parameter count needed by \code{rank_test()}.
#' @noRd
vecm_ml <- function(x,
                    lag = 1L,
                    include = c("const", "none", "both"),
                    LRinclude = c("none", "const", "trend", "both"),
                    r = 1L) {

  include   <- match.arg(include)
  LRinclude <- match.arg(LRinclude)

  x <- as.matrix(x)
  N <- nrow(x)
  K <- ncol(x)
  if (is.null(colnames(x))) colnames(x) <- paste0("y", seq_len(K))
  p <- as.integer(lag)
  if (p < 0L) stop("'lag' must be non-negative.")
  if (N <= K + p + 1L) stop("Not enough observations.")
  if (r < 0L || r > K) stop("'r' must lie between 0 and ncol(x).")

  ## ---- design matrices ------------------------------------------------
  ## effective sample: t = p + 2, ..., N   =>   Tn = N - p - 1
  dx <- diff(x)
  E  <- stats::embed(dx, p + 1L)
  Z0 <- E[, seq_len(K), drop = FALSE]                 # Delta x_t
  Z2 <- if (p > 0L) E[, -seq_len(K), drop = FALSE] else NULL
  Z1 <- x[(p + 1L):(N - 1L), , drop = FALSE]          # x_{t-1}
  Tn <- nrow(Z0)

  colnames(Z0) <- colnames(x)
  colnames(Z1) <- colnames(x)
  if (p > 0L)
    colnames(Z2) <- paste0(colnames(x), " -", rep(seq_len(p), each = K))

  det_mat <- function(kind) {
    D <- matrix(numeric(0), nrow = Tn, ncol = 0L)
    if (kind %in% c("const", "both")) D <- cbind(D, Intercept = rep(1, Tn))
    if (kind %in% c("trend", "both")) D <- cbind(D, Trend = seq_len(Tn))
    D
  }

  ## Deterministics come BEFORE the lagged differences, matching tsDyn's
  ## coefficient layout (ECT, Intercept, Trend, then lags). The column order of
  ## W is irrelevant to the concentration step: qr.resid projects onto the
  ## column span, which is unchanged by permutation. So lambda, beta and the
  ## test statistics are untouched by this ordering.
  W  <- cbind(det_mat(include), Z2)       # unrestricted regressors
  Z1 <- cbind(Z1, det_mat(LRinclude))     # levels + restricted deterministics
  K1 <- ncol(Z1)

  ## ---- concentrate out the short-run part ------------------------------
  qrW <- if (ncol(W) > 0L) qr(W) else NULL
  conc <- function(Y) if (is.null(qrW)) Y else qr.resid(qrW, Y)
  R0 <- conc(Z0)
  R1 <- conc(Z1)

  S00 <- crossprod(R0)     / Tn
  S01 <- crossprod(R0, R1) / Tn
  S11 <- crossprod(R1)     / Tn

  ## ---- reduced rank regression -----------------------------------------
  ## solve |lambda S11 - S10 S00^{-1} S01| = 0 via a symmetric eigenproblem
  C  <- chol(S11)                      # S11 = t(C) %*% C
  Ci <- backsolve(C, diag(K1))         # C^{-1}
  M  <- crossprod(Ci, t(S01)) %*% solve(S00, S01) %*% Ci
  M  <- (M + t(M)) / 2                 # kill numerical asymmetry
  ev <- eigen(M, symmetric = TRUE)

  lambda <- ev$values
  V <- Ci %*% ev$vectors               # normalised so that t(V) S11 V = I
  rownames(V) <- colnames(Z1)
  colnames(V) <- paste0("ev", seq_len(K1))

  ## only K eigenvalues can be non-zero, even with restricted deterministics
  lam <- pmax(pmin(lambda[seq_len(K)], 1 - .Machine$double.eps), 0)

  ## ---- test statistics --------------------------------------------------
  eigen_stat <- -Tn * log(1 - lam)             # H0: rank = r vs rank = r + 1
  trace_stat <- rev(cumsum(rev(eigen_stat)))   # H0: rank <= r
  names(trace_stat) <- names(eigen_stat) <- paste0("r = ", seq_len(K) - 1L)

  ## ---- coefficients for the chosen rank ---------------------------------
  beta <- alpha <- Pi <- NULL
  if (r > 0L) {
    beta <- V[, seq_len(r), drop = FALSE]
    top  <- beta[seq_len(r), , drop = FALSE]
    if (abs(det(top)) < .Machine$double.eps^0.75)
      warning("Leading r x r block of beta is near singular; ",
              "Phillips normalisation may be unstable.")
    beta <- beta %*% solve(top)                # beta = [I_r ; beta_2]
    ect  <- Z1 %*% beta
    ## tsDyn drops the suffix for a single cointegrating vector
    colnames(ect) <- if (r == 1L) "ECT" else paste0("ECT", seq_len(r))
  } else {
    ect <- matrix(numeric(0), nrow = Tn, ncol = 0L)
  }

  RHS <- cbind(ect, W)
  if (ncol(RHS) > 0L) {
    fit    <- stats::lm.fit(RHS, Z0)
    coefs  <- t(fit$coefficients)
    resid  <- Z0 - RHS %*% t(coefs)
  } else {
    coefs  <- matrix(numeric(0), nrow = K, ncol = 0L)
    resid  <- Z0
  }
  rownames(coefs) <- colnames(x)
  npar <- K * ncol(RHS)          # total parameters; npar / K = per equation

  if (r > 0L) {
    alpha <- coefs[, seq_len(r), drop = FALSE]
    Pi    <- alpha %*% t(beta)
  }

  Sigma  <- crossprod(resid) / Tn
  loglik <- -Tn * K / 2 * (log(2 * pi) + 1) -
    Tn / 2 * log(det(S00)) -
    Tn / 2 * sum(log(1 - lam[seq_len(r)]))

  structure(list(
    lambda      = lambda,
    eigenvector = V,
    trace       = trace_stat,
    lambda_max  = eigen_stat,
    beta        = beta,
    alpha       = alpha,
    Pi          = Pi,
    coefficients = coefs,
    residuals   = resid,
    Sigma       = Sigma,
    logLik      = loglik,
    T           = Tn,
    K           = K,
    npar        = npar,
    lag         = p,
    r           = r,
    include     = include,
    LRinclude   = LRinclude,
    S00 = S00, S01 = S01, S11 = S11
  ), class = "vecm_ml")
}


#' Test of the cointegrating rank
#'
#' Maximum likelihood trace and maximum eigenvalue tests for the cointegration
#' rank, with p-values from the Doornik (1998, 1999) gamma approximation.
#'
#'
#' If \code{r_null} is supplied, the p-value for that specific null is
#' returned. Otherwise the rank of the first non-rejected test is selected.
#'
#' @param vecm   object returned by \code{vecm_ml()}.
#' @param type   test driving the automatic rank selection.
#' @param r_null optional specific rank to test.
#' @param cval   level for the automatic selection.
#'
#' @return object of class \code{"rank_test"} with \code{print} and
#'   \code{summary} methods.
#' @noRd
rank_test <- function(vecm, type = c("eigen", "trace"), r_null, cval = 0.05) {

  type <- match.arg(type)
  if (!is.list(vecm) || is.null(vecm$lambda))
    stop("'vecm' must be the result of vecm_ml().")
  if (!is.numeric(cval) || length(cval) != 1L || cval <= 0 || cval >= 1)
    stop("'cval' must be a single number strictly between 0 and 1.")

  Tn <- vecm$T
  K  <- vecm$K

  ## only the first K eigenvalues are non-zero when deterministics are
  ## restricted to the cointegration space
  lambda <- vecm$lambda[seq_len(K)]

  case <- doornik_case(vecm$include, vecm$LRinclude)

  ## ---- statistics -------------------------------------------------------
  ll    <- log(1 - lambda)
  trace <- -Tn * rev(cumsum(rev(ll)))   # H0: rank <= r
  eigen <- -Tn * ll                     # H0: rank = r vs rank = r + 1

  ## ---- p-values ---------------------------------------------------------
  ## n = K - r, so n runs K, K-1, ..., 1 as r runs 0, ..., K-1
  n_seq <- K:1L

  ## adjusted sample size: effective T minus parameters per equation.
  ## The small-sample correction adjusts the moments of the gamma, NOT the
  ## statistic, so the uncorrected 'trace' is passed in either case.
  if (is.null(vecm$npar))
    stop("'vecm$npar' missing. Add it to vecm_ml(): npar <- K * ncol(RHS).")
  adjT <- Tn - floor(vecm$npar / K)

  pv <- function(stat, ty, small = FALSE, Tuse = NULL)
    vapply(seq_along(stat), function(i)
      doornik_p(stat[i], n = n_seq[i], case = case, type = ty,
                T = Tuse, smallSamp = small),
      numeric(1))

  trace_pval   <- pv(trace, "trace")
  trace_pval_T <- if (adjT > 0) pv(trace, "trace", TRUE, adjT) else rep(NA_real_, K)
  eigen_pval   <- pv(eigen, "eigen")   # no published small-sample correction

  ## ---- rank selection ---------------------------------------------------
  pvals <- switch(type, trace = trace_pval, eigen = eigen_pval)

  if (all(is.na(pvals))) {
    ## Doornik tables not populated, or no p-value available for this case.
    ## Fail loudly rather than letting NA reach the if() below.
    #warning("No p-values available; rank selection skipped.")
    pval <- NA_real_
    rank <- NA_integer_
    spec <- "unavailable"
  } else if (missing(r_null)) {
    w    <- if (all(pvals < cval)) length(pvals) else which(pvals > cval)[1L]
    pval <- pvals[w]
    rank <- if (pval > cval) w - 1L else w
    spec <- "unspec"
  } else {
    if (!is.numeric(r_null) || length(r_null) != 1L ||
        r_null < 0L || r_null > K - 1L)
      stop("'r_null' must be a single integer between 0 and ", K - 1L, ".")
    pval <- pvals[r_null + 1L]
    rank <- as.integer(r_null)
    spec <- "specified"
  }

  res_df <- data.frame(
    r            = 0:(K - 1L),
    trace        = trace,
    trace_pval   = trace_pval,
    trace_pval_T = trace_pval_T,
    eigen        = eigen,
    eigen_pval   = eigen_pval,
    row.names    = NULL
  )

  structure(list(
    res_df = res_df,
    r      = rank,
    pval   = pval,
    case   = case,
    T      = Tn,
    T_adj  = adjT,
    K      = K,
    type   = type,
    cval   = cval,
    r_null = spec,
    call   = match.call()
  ), class = "rank_test")
}


#' Gamma-approximation p-value for one statistic
#'
#' approach and tables borrowed from tsDyn package
#'
#' @param x     the test statistic (UNCORRECTED, even when smallSamp = TRUE).
#' @param n     n = K - r, the number of common trends under the null.
#' @param case  one of Hnames.
#' @param type  "trace" or "eigen".
#' @param T     adjusted sample size, only used when smallSamp = TRUE.
#' @param smallSamp apply the Table 9 correction. Trace only.
#' @noRd
doornik_p <- function(x, n, case, type = c("trace", "eigen"),
                      T = NULL, smallSamp = FALSE) {
  type <- match.arg(type)
  case <- match.arg(case, Hnames)
  if (is.na(x)) return(NA_real_)
  if (type == "trace") {
    v_mean <- c(n^2, n, 1, as.numeric(n == 1), as.numeric(n == 2), sqrt(n))
    v_var  <- c(n^2, n, 1, as.numeric(n == 1), as.numeric(n == 2))
    mu <- drop(v_mean %*% tab7_mean[, case])
    v  <- drop(v_var  %*% tab7_var[, case])
    if (smallSamp) {
      if (is.null(T)) stop("'T' required when smallSamp = TRUE.")
      ## NOTE: the correction multiplies the MOMENTS, it does not rescale the
      ## statistic. The statistic passed in stays the uncorrected one.
      vs <- c(sqrt(n)/T, n/T, n^2/T^2, as.numeric(n == 1)/T,
              as.numeric(n == 1), as.numeric(n == 2), as.numeric(n == 3))
      mu <- exp(log(mu) + drop(vs %*% tab9_mean[, case]))
      v  <- exp(log(v)  + drop(vs %*% tab9_var[, case]))
    }
  } else {
    ## no small-sample correction is published for the eigen statistic
    v_eig <- c(n, 1, as.numeric(n == 1), as.numeric(n == 2), sqrt(n))
    mu <- drop(v_eig %*% tab8_mean[, case])
    v  <- drop(v_eig %*% tab8_var[, case])
  }
  shape <- mu^2 / v
  rate  <- mu   / v
  stats::pgamma(x, shape = shape, rate = rate, lower.tail = FALSE)
}

#' Map a vecm_ml() specification onto Doornik's case labels
#' @noRd
doornik_case <- function(include, LRinclude) {
  key <- paste(include, LRinclude, sep = "/")
  switch(key,
         "none/none"   = "H_z",
         "none/const"  = "H_c",
         "const/none"  = "H_lc",
         "const/trend" = "H_l",
         "both/none"   = "H_ql",
         stop("No Doornik case for include = '", include,
              "', LRinclude = '", LRinclude, "'."))
}
