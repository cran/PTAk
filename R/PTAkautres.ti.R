"CANDPARA" <-
function (X, dim = 3, test = 1e-08, Maxiter = 1000, smoothing = FALSE, 
    smoo = list(NA), verbose = getOption("verbose"), file = NULL, 
    modesnam = NULL, addedcomment = "") 
{
    datanam <- substitute(X)
    sym <- NULL
    if (is.list(X)) {
        if (is.list(X$met)) 
            metrics <- TRUE
        else stop(paste("------with metrics X must be a list with $data and $met----"))
    }
    else metrics <- FALSE
    if (metrics) {
        nam <- dimnames(X$data)
        diX <- length(dim(X$data))
        for (d in 1:diX) {
            if (length(X$met[[d]]) > 1) {
                if (length(X$met[[d]]) == dim(X$data)[d]^2) {
                  tempp <- d
                  t12 <- CONTRACTION(X$data, Powmat(X$met[[d]], 
                    1/2), Xwiz = d, zwiX = 1)
                  d <- tempp
                  lacola <- (1:diX)[-d]
                  laperm <- c(lacola, d)
                }
                else {
                  lacola <- (1:diX)[-d]
                  laperm <- c(d, lacola)
                  lacol <- (dim(X$data))[lacola]
                  pt12 <- matrix(aperm(X$data, laperm), ncol = prod(lacol))
                  t12 <- sqrt(X$met[[d]]) * pt12
                }
                t12 <- array(t12, (dim(X$data))[laperm])
                X$data <- aperm(t12, match(1:diX, laperm))
            }
            else X$data <- X$data * X$met[[d]]
        }
        met <- X$met
        X <- X$data
        dimnames(X) <- nam
    }
    debtime <- proc.time()
    pass. <- function(a, r) {
        pasta <- a
        for (i in 2:r) pasta <- paste(pasta, a, sep = "")
        return(pasta)
    }
    if (verbose) {
        cat("\n", "       ----------+++++++++++------------", 
            "\n", ifelse(smoothing, paste("Smoothed ", "\n"), 
                ""), "              PARAFAC/CANDECOMP ", "\n", 
            file = ifelse(is.null(file), "", file), append = TRUE)
        cat("       ----------+++++++++++------------", "\n", 
            file = ifelse(is.null(file), "", file), append = TRUE)
        cat(" Data is ... ", deparse(datanam), "...", file = ifelse(is.null(file), 
            "", file), append = TRUE)
        cat("  .... Tensor of order ", length(dim(X)), file = ifelse(is.null(file), 
            "", file), append = TRUE)
        cat("  ....  with dimensions: ", dim(X), "\n", file = ifelse(is.null(file), 
            "", file), append = TRUE)
        if (!is.null(modesnam)) 
            cat("modes are ", modesnam, "\n", file = ifelse(is.null(file), 
                "", file), append = TRUE)
        if (metrics) 
            cat("---Analysis with non-Identity metrics  ------", 
                "\n", file = ifelse(is.null(file), "", file), 
                append = TRUE)
        if (!addedcomment == "") 
            cat("\n", addedcomment, "\n", file = ifelse(is.null(file), 
                "", file), append = TRUE)
    }
    if (!is.array(X)) {
        stop(paste("--- X must be an array  ! ---"))
    }
    ord <- length(dim(X))
    if (is.null(modesnam)) {
        modesnam <- paste(rep("mo", ord), 1:ord)
    }
    if (smoothing) {
        if (length(smoo) < ord) 
            smoo <- rep(list(smoo[[1]]), ord)
    }
    else smoo <- list(NULL)
    sval0 <- INITIA(X, modesnam = modesnam, method = "svd", dim = dim)
    test0 <- 1
    atest <- 0
    sval <- sval0
    iter <- 0
    if (smoothing) {
        for (j in 1:ord) if (!is.list(smoo[[j]])) 
            smoo[[j]] <- list(smoo[[j]])
        for (a in 2:dim) {
            for (j in 1:ord) if (length(smoo[[j]]) < a) 
                smoo[[j]][[a]] <- smoo[[j]][[a - 1]]
        }
    }
    while (test0 > test) {
        iter <- iter + 1
        if (verbose & iter%%100 == 1) 
            cat("\n", " ----------- iteration-", iter, file = ifelse(is.null(file), 
                "", file), append = TRUE)
        for (i in 1:ord) {
            if (iter == 1) {
                if (verbose) 
                  cat("\n", i, "^", sval0[[i]]$d, file = ifelse(is.null(file), 
                    "", file), append = TRUE)
                sval[[i]]$d <- NULL
            }
        }
        if (i == 1) {
            tzz <- 1
            Z <- 1
        }
        else {
            ifelse(dim == 1, tzz <- 1, tzz <- sval[[1]]$v %*% 
                t(sval[[1]]$v))
            Z <- t(sval[[1]]$v)
        }
        for (j in 2:ord) {
            if (!j == i) {
                if (dim > 1) 
                  tzz <- tzz * (sval[[j]]$v %*% t(sval[[j]]$v))
                Z <- RaoProd(t(sval[[j]]$v), Z)
            }
        }
        sval[[i]]$v <- t(matrix(aperm(X, c(i, (1:length(dim(X)))[-i])), 
            nrow = dim(X)[i]) %*% Z %*% Ginv(tzz))
        if (smoothing) {
            for (a in 1:dim) if (is.function(smoo[[i]][[a]])) 
                sval[[i]]$v[a, ] <- smoo[[i]][[a]](sval[[i]]$v[a, 
                  ])
        }
        sval[[i]]$d <- sqrt(diag(sval[[i]]$v %*% t(sval[[i]]$v)))
        sval[[i]]$v <- sval[[i]]$v/sval[[i]]$d
        atest <- atest + sum((sval[[i]]$v - sval0[[i]]$v)^2)
        if (!is.null(sym)) {
            for (i in ord:1) {
                if (!i == sym[i]) 
                  sval[[sym[i]]] <- sval[[i]]
            }
        }
        sval0 <- sval
        test0 <- sqrt(atest)
        atest <- 0
        if (verbose & (iter%%100) == 1) 
            cat("\n", "----------- test =         ", test0, "\n", 
                file = ifelse(is.null(file), "", file), append = TRUE)
        if (iter > (Maxiter - 1) & (iter - Maxiter)%%100 == 0) {
            cat("\n \n \n \n \n ", " WARNING ****** Iteration already =  ", 
                iter, "\n")
            cat(" ** type anything to STOP ** just RETURN to carry on **", 
                "\n")
            conti <- scan("", what = "", n = 1, quiet = TRUE, 
                flush = TRUE, )
            if (length(conti) > 0) 
                stop(paste(" ---- Aborted by request ---- "))
        }
    }
    pourRR2 <- function() {
        tens <- t(sval[[1]]$v) %*% diag(sval[[ord]]$d)
        for (r in 2:ord) {
            tens <- RaoProd(t(sval[[r]]$v), tens)
        }
        return(summary(lm(as.vector(X) ~ tens - 1)))
    }
    pass. <- function(a, r) {
        pasta <- a
        for (i in 2:r) pasta <- paste(pasta, a, sep = "")
        return(pasta)
    }
    PCnam <- paste("v", pass.(1:dim, ord), sep = "")
    ssX <- sum(X^2)
    sstens <- (sval[[i]]$d^2)
    PCT <- 100 * sstens/ssX
    sval[[i]]$lm <- pourRR2()
    if (verbose) {
        cat(" --------optimisation  done ", "\n", file = ifelse(is.null(file), 
            "", file), append = TRUE)
        cat(" --------Final iteration----", iter, "\n", file = ifelse(is.null(file), 
            "", file), append = TRUE)
        cat(" ----------- test =         ", test0, "\n", file = ifelse(is.null(file), 
            "", file), append = TRUE)
        cat("\n", " --Norms-- ", sval[[i]]$d, "\n", " --Percent-- ", 
            PCT)
        cat("\n", " ---Total R2 ", sval[[i]]$lm$r.squared * 100, 
            "%", "\n")
    }
    cat("-----Execution Time-----", (proc.time() - debtime)[3], 
        "\n")
    if (metrics) {
        for (d in 1:length(sval)) {
            if (length(met[[d]]) > 1) {
                if (length(met[[d]]) == dim(X)[d]^2) {
                  sval[[d]]$v <- sval[[d]]$v %*% Powmat(met[[d]], 
                    -1/2)
                }
                else {
                  sval[[d]]$v <- t(1/sqrt(met[[d]]) * t(sval[[d]]$v))
                }
            }
            else sval[[d]]$v <- sval[[d]]$v * 1/sqrt(met[[d]])
        }
    }
    sval[[i]]$pct <- as.vector(PCT)
    sval[[i]]$ssX <- as.vector(ssX)
    sval[[i]]$vsnam <- PCnam
    sval[[i]]$datanam <- datanam
    sval[[i]]$method <- match.call()
    class(sval) <- c("solutions.CANDPARA", "solutions.PTAk")
    invisible(return(sval))
}
"PCAn" <-
function (X, dim = c(2, 2, 2, 3), test = 1e-12, Maxiter = 400, 
    smoothing = FALSE, smoo = list(NA), verbose = getOption("verbose"), 
    file = NULL, modesnam = NULL, addedcomment = "") 
{
    datanam <- substitute(X)
    sym <- NULL
    if (is.list(X)) {
        if (is.list(X$met)) 
            metrics <- TRUE
        else stop(paste("------with metrics X must be a list with $data and $met----"))
    }
    else metrics <- FALSE
    if (metrics) {
        nam <- dimnames(X$data)
        diX <- length(dim(X$data))
        for (d in 1:diX) {
            if (length(X$met[[d]]) > 1) {
                if (length(X$met[[d]]) == dim(X$data)[d]^2) {
                  tempp <- d
                  t12 <- CONTRACTION(X$data, Powmat(X$met[[d]], 
                    1/2), Xwiz = d, zwiX = 1)
                  d <- tempp
                  lacola <- (1:diX)[-d]
                  laperm <- c(lacola, d)
                }
                else {
                  lacola <- (1:diX)[-d]
                  laperm <- c(d, lacola)
                  lacol <- (dim(X$data))[lacola]
                  pt12 <- matrix(aperm(X$data, laperm), ncol = prod(lacol))
                  t12 <- sqrt(X$met[[d]]) * pt12
                }
                t12 <- array(t12, (dim(X$data))[laperm])
                X$data <- aperm(t12, match(1:diX, laperm))
            }
            else X$data <- X$data * X$met[[d]]
        }
        met <- X$met
        X <- X$data
        dimnames(X) <- nam
    }
    debtime <- proc.time()
    pass. <- function(a, r) {
        pasta <- a
        for (i in 2:r) pasta <- paste(pasta, a, sep = "")
        return(pasta)
    }
    if (verbose) {
        cat("----------+++++++++++------------", "\n", ifelse(smoothing, 
            paste("Smoothed ", "\n"), ""), " PCA-n modes  ", 
            "\n", file = ifelse(is.null(file), "", file), append = TRUE)
        cat(" Data is ... ", deparse(datanam), "...", "\n", file = ifelse(is.null(file), 
            "", file), append = TRUE)
        cat("  .... Tensor of order ", length(dim(X)), file = ifelse(is.null(file), 
            "", file), append = TRUE)
        cat("  ....  with dimensions: ", dim(X), "\n", file = ifelse(is.null(file), 
            "", file), append = TRUE)
        if (!is.null(modesnam)) 
            cat("modes are ", modesnam, "\n", file = ifelse(is.null(file), 
                "", file), append = TRUE)
        if (metrics) 
            cat("---Analysis with non-Identity metrics  ------", 
                "\n", file = ifelse(is.null(file), "", file), 
                append = TRUE)
        if (!addedcomment == "") 
            cat("\n", addedcomment, "\n", file = ifelse(is.null(file), 
                "", file), append = TRUE)
    }
    if (!is.array(X)) {
        stop(paste("--- X must be an array  ! ---"))
    }
    solutions <- NULL
    ord <- length(dim(X))
    if (smoothing) {
        if (length(smoo) < ord) 
            smoo <- rep(list(smoo[[1]]), ord)
    }
    else smoo <- list(NULL)
    if (is.null(modesnam)) {
        modesnam <- paste(rep("mo", ord), 1:ord)
    }
    if (!length(dim) == ord) 
        stop(" Wrong length for dim argument (= Rank-spaces !)")
    for (j in 1:ord) if (dim[j] > dim(X)[j]) 
        stop(" (dim argument) some Rank-spaces are too big!")
    sval0 <- INITIA(X, modesnam = modesnam, method = "svd", dim = dim)
    test0 <- 1
    atest <- 0
    sval <- sval0
    iter <- 0
    if (smoothing) {
        for (j in 1:ord) {
            if (!is.list(smoo[[j]])) 
                smoo[[j]] <- list(smoo[[j]])
            for (a in 2:dim[j]) if (length(smoo[[j]]) < a) 
                smoo[[j]][[a]] <- smoo[[j]][[a - 1]]
        }
    }
    while (test0 > test) {
        iter <- iter + 1
        if (verbose & iter%%100 == 1) {
            cat(" ----------- iteration-", iter, "\n", file = ifelse(is.null(file), 
                "", file), append = TRUE)
        }
        for (i in 1:ord) {
            if (iter == 1) {
                if (verbose) {
                  cat(" ", i, "^", sval0[[i]]$d, file = ifelse(is.null(file), 
                    "", file), append = TRUE)
                }
            }
            sval[[i]]$d <- NULL
        }
        corematv <- X
        for (j in 1:ord) {
            if (j < i) 
                corematv <- CONTRACTION(corematv, sval[[j]]$v, 
                  Xwiz = 1, zwiX = 2)
            if (j > i) 
                corematv <- CONTRACTION(corematv, sval[[j]]$v, 
                  Xwiz = 2, zwiX = 2)
        }
        corematv <- matrix(corematv, nrow = dim(X)[i])
        if (smoothing) 
            svdcormatv <- svdsmooth(corematv, nomb = dim[i], 
                smooth = list(NA, smoo[[i]]))
        else svdcormatv <- svd(corematv)
        sval[[i]]$dopt <- svdcormatv$d[1:dim[i]]
        sval[[i]]$v <- t(svdcormatv$u[, 1:dim[i]])
        if (all(svdcormatv$u[, 1] < 0)) 
            sval[[i]]$v <- -sval[[i]]$v
        coremat <- array(t(corematv) %*% t(sval[[i]]$v), c(dim[-i], 
            dim[i]))
        atest <- atest + sum((sval[[i]]$v - sval0[[i]]$v)^2)
        if (!is.null(sym)) {
            for (i in ord:1) {
                if (!i == sym[i]) 
                  sval[[sym[i]]] <- sval[[i]]
            }
        }
        sval0 <- sval
        if (verbose & iter%%100 == 0) {
            cat(" --", coremat, file = ifelse(is.null(file), 
                "", file), append = TRUE)
        }
        test0 <- sqrt(atest)
        atest <- 0
        if (verbose & (iter%%100) == 1) {
            cat("\n", "----------- test =         ", test0, "\n", 
                file = ifelse(is.null(file), "", file), append = TRUE)
        }
        if (iter > (Maxiter - 1) & (iter - Maxiter)%%100 == 0) {
            cat("\n \n \n \n \n ", " WARNING ****** Iteration already =  ", 
                iter, "\n")
            cat(" ** type anything to STOP ** just RETURN to carry on **", 
                "\n")
            conti <- scan("", what = "", n = 1, quiet = TRUE, 
                flush = TRUE, )
            if (length(conti) > 0) 
                stop(paste(" ---- Aborted by request ---- "))
        }
    }
    PCnam <- outer(outer(1:dim[1], 1:dim[2], FUN = "paste", sep = ""), 
        1:dim[3], FUN = "paste", sep = "")
    if (ord > 3) {
        for (t in 4:ord) {
            PCnam <- outer(PCnam, 1:dim[q], FUN = "paste", sep = "")
        }
    }
    PCnam <- paste("v", PCnam, sep = "")
    ssX <- sum(X^2)
    sstens <- sum(coremat^2)
    totPCT <- 100 * sstens/ssX
    if (verbose) {
        cat(" --------optimisation  done ", "\n", file = ifelse(is.null(file), 
            "", file), append = TRUE)
        cat(" --------Final iteration----", iter, "\n", file = ifelse(is.null(file), 
            "", file), append = TRUE)
        cat(" ----------- test =         ", test0, "\n", file = ifelse(is.null(file), 
            "", file), append = TRUE)
        cat("\n", " --Core Matrix-- ", coremat, "\n", " --  Percent -- ", 
            totPCT, "%", "\n")
    }
    if (metrics) {
        for (d in 1:length(sval)) {
            if (length(met[[d]]) > 1) {
                if (length(met[[d]]) == dim(X)[d]^2) {
                  sval[[d]]$v <- sval[[d]]$v %*% Powmat(met[[d]], 
                    -1/2)
                }
                else {
                  sval[[d]]$v <- t(1/sqrt(met[[d]]) * t(sval[[d]]$v))
                }
            }
            else sval[[d]]$v <- sval[[d]]$v * 1/sqrt(met[[d]])
        }
    }
    cat("-----Execution Time-----", (proc.time() - debtime)[3], 
        "\n")
    sval[[i]]$d <- as.vector(coremat)
    sval[[i]]$coremat <- coremat
    sval[[i]]$pct <- as.vector(totPCT)
    sval[[i]]$ssX <- as.vector(ssX)
    sval[[i]]$vsnam <- PCnam
    sval[[i]]$datanam <- datanam
    sval[[i]]$method <- match.call()
    class(sval) <- c("solutions.PCAn", "solutions.PTAk")
    invisible(return(sval))
}
"REBUILDPCAn" <-
function (solu) 
{
    lo <- length(solu)
    recon <- t(solu[[1]]$v)
    for (k in 2:lo) {
        recon <- recon %o% t(essreconf[[k]]$v)
    }
    reconf <- CONTRACTION(recon, solu[[lo]]$coremat, Xwiz = (1:(lo - 
        1)) * 2, zwiX = 1:(lo - 1))
    cat("\n", "--- RMSerror---", sqrt(mean((eval(solu[[lo]]$datanam) - 
        reconf)^2)))
    invisible(return(reconf))
}
