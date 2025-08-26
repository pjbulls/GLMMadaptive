mixed_model <- function (fixed, random, data, family, weights = NULL,
                         na.action = na.exclude, n_phis = NULL, initial_values = NULL, 
                         control = list(), ...) {
    call <- match.call()
    # set family
    if (is.character(family))
        family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
        family <- family()
    if (is.null(family$family)) {
        stop("'family' not recognized.\n")
    }
    if (family$family == "gaussian")
        stop("use function lme() from package 'nlme' or function lmer() from ",
             "package 'lme4'.\n")
    if (length(grep("Negative Binomial", family$family))) {
        stop("Because the namespace of the MASS package seems also to be loaded\n",
             "  use 'family = GLMMadaptive::negative.binomial()'.")
    }
    if (family$family == "Gamma" && is.null(family$log_dens)) {
        if (family$link != "log") {
            warning("with the Gamma family currently only the log link ",
                    "function works. It has been reset automatically.")
        }
        family <- Gamma.fam()
    }
    known_families <- c("binomial", "poisson", "negative binomial", "Gamma")
    if (inherits(data, "tbl_df") || inherits(data, "tbl"))
        data <- as.data.frame(data) # in case 'data' is a tibble
    orig_data <- data
    groups <- unique(all.vars(getID_Formula(random)))
    data[groups] <- lapply(data[groups], 
                           function (x) if (!is.factor(x)) factor(x, levels = unique(x)) else x)
    data <- data[order(data[[groups[1L]]]), ]
    # drop unused levels of factors
    factors <- sapply(data, is.factor)
    data[factors] <- lapply(data[factors], function (x) x[, drop = TRUE])
    # construct model frames
    # fixed effects
    mfX <- model.frame(terms(fixed, data = data), data = data, na.action = na.pass)
    termsX <- terms(mfX)
    # random effects
    form_random <- constructor_form_random(random, data)
    mfZ <- lapply(form_random, function (form, data) 
        model.frame(terms(form, data = data), data = data, na.action = na.pass), 
        data = data)
    termsZ <- lapply(mfZ, terms)
    # fixed effects ZI part
    mfX_zi <- termsX_zi <- NULL
    # random effects ZI part
    mfZ_zi <- termsZ_zi <- NULL
    # delete missing data if na.action is na.fail or na.omit
    chr_na_action <- as.character(body(na.action))[2]
    if (chr_na_action == "na.exclude" || chr_na_action == "na.omit") {
        complete_cases <- cbind(complete.cases(mfX), sapply(mfZ, complete.cases))
        keep <- apply(complete_cases, 1, all)
        mfX <- mfX[keep, , drop = FALSE]
        mfZ[] <- lapply(mfZ, function (mf) mf[keep, , drop = FALSE])
    }
    # id variable
    id_nam <- all.vars(getID_Formula(random))
    id_orig <- model.frame(terms(getID_Formula(random), data = data), 
                           data = data, na.action = na.pass)
    if (chr_na_action == "na.exclude" || chr_na_action == "na.omit")
        id_orig <- id_orig[keep, , drop = FALSE]
    id <- match(id_orig[[1L]], unique(id_orig[[1L]]))
    # weights
    if (!is.null(weights)) {
        if (!is.numeric(weights))
            stop("'weights' must be a numeric vector.")
        if (length(weights) != length(unique(id)))
            stop("the length of 'weights' does not match with the number of groups in 'data'.")
    }
    # model response
    y <- model.response(mfX)
    if (is.factor(y)) {
        if (family$family == "binomial")
            y <- as.numeric(y != levels(y)[1L])
        else
            stop("the response variable should not be a factor.\n")
    }
    # construct model matrices
    # fixed effects
    offset <- model.offset(mfX)
    X <- model.matrix(termsX, mfX)
    # random effects
    Z <- mapply(constructor_Z, termsZ, mfZ, MoreArgs = list(id = id), SIMPLIFY = FALSE)
    Z <- do.call("cbind", Z)
    # fixed effects ZI part
    offset_zi <- X_zi <- NULL
    # random effects ZI part
    Z_zi <- NULL
    nRE <- ncol(Z)
    ###########################
    # control settings
    con <- list(iter_EM = 30, iter_qN_outer = 15, iter_qN = 10, iter_qN_incr = 10,
                optimizer = "optim", optim_method = "BFGS", parscale_betas = 0.1, 
                parscale_D = 0.01, parscale_phis = 0.01, parscale_gammas = 0.01, 
                tol1 = 1e-04, tol2 = 1e-05, tol3 = 1e-08, numeric_deriv = "fd", 
                nAGQ = if (nRE < 3) 11 else 7, update_GH_every = 10, max_coef_value = 12, 
                max_phis_value = exp(10), verbose = FALSE)
    control <- c(control, list(...))
    namC <- names(con)
    con[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% namC]) > 0)
        warning("unknown names in control: ", paste(noNms, collapse = ", "))
    ###########################
    # initial values
    diag_D <- (random[[length(random)]])[[1]] == as.name("||")
    inits <- if (family$family %in% known_families || (is.list(initial_values) &&
                                                       inherits(initial_values$betas, 'family'))) {
        betas <- if (family$family %in% known_families) {
            if (family$family == "negative binomial")
                glm.fit(X, y, family = poisson(), offset = offset)$coefficients
            else if (family$family == "Gamma")
                glm.fit(X, y, family = Gamma(), offset = offset)$coefficients
            else if (family$family == "censored normal")
                glm.fit(X, y, family = gaussian(), offset = offset)$coefficients
            else if (family$family == "zero-inflated binomial")
                glm.fit(X, y, family = binomial(), offset = offset)$coefficients
            else
                glm.fit(X, y, family = family, offset = offset)$coefficients
        } else {
            glm.fit(X, y, family = initial_values$betas, offset = offset)$coefficients
        }
        list(betas = betas * sqrt(1.346), D = if (diag_D) rep(1, nRE) else diag(nRE))
    } else {
        list(betas = rep(0, ncol(X)), D = if (diag_D) rep(1, nRE) else diag(nRE))
    }
    ##########################
    # penalized
    ##########################
    # Functions
    Funs <- list(
        mu_fun = family$linkinv,
        var_fun = family$variance,
        mu.eta_fun = family$mu.eta
    )
    if (family$family %in% known_families && is.null(family$log_dens)) {
        Funs$log_dens <- switch(family$family,
                                'binomial' = binomial_log_dens,
                                'poisson' = poisson_log_dens,
                                'negative binomial' = negative.binomial_log_dens)
    } else if (family$family %in% known_families && !is.null(family$log_dens)) {
        Funs$log_dens <- family$log_dens
    } else if (!family$family %in% known_families && !is.null(family$log_dens)) {
        Funs$log_dens <- family$log_dens
    } else {
        stop("'log_dens' component of the 'family' argument is NULL with no default.\n")
    }
    if (!is.function(Funs$log_dens)) {
        stop("'log_dens' component of the 'family' argument must be a function.\n")
    }
    if (!is.function(Funs$mu_fun)) {
        stop("'linkinv' component of the 'family' argument must be a function.\n")
    }
    if (!is.function(Funs$mu_fun)) {
        stop("'linkinv' component of the 'family' argument must be a function.\n")
    }
    if (!is.null(family$score_eta_fun) && is.function(family$score_eta_fun)) {
        Funs$score_eta_fun <- family$score_eta_fun
    }
    if (!is.null(family$score_phis_fun) && is.function(family$score_phis_fun)) {
        Funs$score_phis_fun <- family$score_phis_fun
    }
    has_phis <- inherits(try(Funs$log_dens(y, rep(0, length(y)), Funs$mu_fun, 
                                           phis = NULL), TRUE),
                         "try-error")
    if (has_phis) {
        if (family$family %in% c("negative binomial", "zero-inflated negative binomial",
                                 "hurdle negative binomial", "hurdle log-normal",
                                 "beta", "hurdle beta", "Conway Maxwell Poisson", 
                                 "Gamma", "censored normal", "beta binomial",
                                 "Student's-t")) {
            n_phis <- 1
        } else if (is.null(n_phis)) {
            stop("argument 'n_phis' needs to be specified.\n")
        }
        inits$phis <- rep(0.0, n_phis)
    }
    if (!is.null(initial_values) && is.list(initial_values) &&
        !inherits(initial_values$betas, 'family')) {
        lngths <- lapply(inits[(nams.initial_values <- names(initial_values))], length)
        if (!isTRUE(all.equal(lngths, lapply(initial_values, length)))) {
            warning("'initial_values' is not a list with elements of appropriate ",
                    "length; default initial_values are used instead.\n")
        } else {
            inits[nams.initial_values] <- initial_values
        }
    }
    ###############
    # Fit the model
    out <- mixed_fit(y, X, Z, id, offset, family, inits, Funs, 
                     con, weights)
    # check whether Hessian is positive definite at convergence
    H <- out$Hessian
    if (any(is.na(H) | !is.finite(H))) {
        warning("infinite or missing values in Hessian at convergence.\n")
    } else {
        ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
        if (!all(ev >= -1e-06 * abs(ev[1L]))) 
            warning("Hessian matrix at convergence is not positive definite; ", 
                    "unstable solution.\n")
    }
    # fix names
    names(out$coefficients) <- colnames(X)
    RE_nams <- colnames(Z)
    dimnames(out$D) <- list(RE_nams, RE_nams)
    if (!is.null(out$phis))
        names(out$phis) <- paste0("phi_", seq_along(out$phis))
    if (!is.null(out$gammas))
        names(out$gammas) <- colnames(X_zi)
    all_nams <- if (diag_D) {
        nams_D <- paste0("D_", seq_len(nRE), seq_len(nRE))
        c(names(out$coefficients), nams_D, names(out$phis), 
          if (!is.null(out$gammas)) paste0("zi_", names(out$gammas)))
    } else {
        nams_D <- paste0("D_", apply(which(upper.tri(out$D, TRUE), arr.ind = TRUE), 1, 
                                     paste0, collapse = ""))
        c(names(out$coefficients), nams_D, names(out$phis), 
          if (!is.null(out$gammas)) paste0("zi_", names(out$gammas)))
    }
    dimnames(out$Hessian) <- list(all_nams, all_nams)
    out$data <- orig_data
    out$id <- id_orig
    out$id_name <- id_nam 
    out$offset <- offset
    dimnames(out$post_modes) <- list(unique(id_orig[[1L]]), RE_nams)
    names(out$post_vars) <- unique(id_orig[[1L]])
    out$post_vars[] <- lapply(out$post_vars, function (v) {
        dimnames(v) <- list(RE_nams, RE_nams)
        v
    })
    out$Terms <- list(termsX = termsX, termsZ = termsZ)
    out$model_frames <- list(mfX = mfX, mfZ = mfZ)
    out$control <- con
    out$Funs <- Funs
    out$family <- family
    out$weights <- weights
    out$na.action <- na.action
    out$contrasts <- attr(X, "contrasts")
    out$call <- call
    class(out) <- "MixMod"
    out
}
