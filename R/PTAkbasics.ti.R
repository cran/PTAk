"CONTRACTION" <-
function (X, z, Xwiz = NULL, zwiX = NULL, rezwiX = FALSE, usetensor = TRUE) 
{
    if (usetensor) {
        if (is.null(zwiX)) {
            if (is.vector(z)) 
                zwiX <- 1
            else zwiX <- 1:length(dim(z))
        }
        if (is.null(Xwiz)) {
            if (is.vector(X)) 
                Xwiz <- 1
            else Xwiz <- 1:length(dim(X))
        }
        return(tensor(X, z, Xwiz, zwiX))
    }
    else {
        non <- function(A, awi) {
            (1:length(dim(A)))[!(1:length(dim(A))) %in% awi]
        }
        zbX <- FALSE
        if (length(dim(as.array(X))) < length(dim(as.array(z)))) {
            zbX <- TRUE
            temp <- X
            X <- z
            z <- temp
            temp <- Xwiz
            Xwiz <- zwiX
            zwiX <- temp
        }
        if (is.vector(z)) {
            zwiX <- 1
            zz <- z
            lacolz <- NULL
            lacolaz <- NULL
            if (is.null(Xwiz)) 
                Xwiz <- (1:length(dim(X)))[dim(X) %in% length(z)]
        }
        else {
            if (is.null(zwiX)) {
                zwiX <- 1:length(dim(z))
            }
            if (is.null(Xwiz)) 
                Xwiz <- match(dim(z)[zwiX], dim(X))
            Xwiz <- Xwiz[!is.na(Xwiz)]
            if (rezwiX) 
                zwiX <- match(dim(X)[Xwiz], dim(z))
            zwiX <- zwiX[!is.na(zwiX)]
            if (!all(dim(X)[Xwiz] == dim(z)[zwiX])) 
                stop(paste(" @@@@@ WRONG matching for contraction!@@@@@"))
            if (all(dim(z) %in% dim(X)[Xwiz]) & length(dim(z)) < 
                length(dim(X)[Xwiz])) {
                zz <- as.vector(aperm(z, zwiX))
                lacolz <- NULL
                lacolaz <- NULL
            }
            else {
                czwiX <- non(z, zwiX)
                lacolaz <- (1:length(dim(z)))[czwiX]
                lacolz <- (dim(z))[lacolaz]
                zz <- matrix(aperm(z, c(zwiX, lacolaz)), ncol = prod(lacolz))
            }
        }
        lacola <- (1:length(dim(X)))[non(X, Xwiz)]
        laperm <- c(Xwiz, lacola)
        lacol <- (dim(X))[lacola]
        toconz <- matrix(aperm(X, laperm), ncol = prod(lacol))
        Xz <- t(toconz) %*% zz
        dinam <- function(A) {
            namA <- rep(paste(1:length(dim(A))), dim(A))
            dinamA <- list(namA[1:dim(A)[1]])
            for (e in 2:length(dim(A))) {
                dinamA <- c(dinamA, list(namA[(sum(dim(A)[1:(e - 
                  1)]) + 1):sum(dim(A)[1:e])]))
            }
            return(dinamA)
        }
        if (!is.null(dimnames(X))) 
            dimnamX <- dimnames(X)
        else dimnamX <- dinam(X)
        if (!is.vector(z)) {
            if (!is.null(dimnames(z))) 
                dimnamz <- dimnames(z)
            else dimnamz <- dinam(z)
        }
        if (is.null(lacolaz)) 
            ladim <- dimnamX[lacola]
        else ladim <- c(dimnamX[lacola], dimnamz[lacolaz])
        Xz <- array(Xz, c(lacol, lacolz), dimnames = ladim)
        if (zbX) 
            Xz <- aperm(Xz, c((1:length(dim(Xz)))[-lacola], lacola))
        return(Xz)
    }
}
"CONTRACTION.list" <-
function (X, zlist, moins = 1, zwiX = NULL, usetensor = TRUE, 
    withapply = FALSE) 
{
    mplu <- 1
    lz <- length(zlist)
    for (tu in 1:lz) {
        if (!tu == moins) {
            if (withapply) 
                X <- apply(X, (1:length(dim(X)))[-mplu], FUN = function(x) as.vector(x) %*% 
                  zlist[[tu]]$v)
            else X <- CONTRACTION(X, zlist[[tu]]$v, Xwiz = mplu, 
                zwiX = zwiX[tu], usetensor = usetensor)
        }
        else mplu <- mplu + 1
    }
    return(X)
}
"INITIA" <-
function (X, modesnam = NULL, method = "Presvd", dim = 1, ...) 
{
    if (!is.array(X)) {
        stop(paste("--- X must be an array  ! ---"))
    }
    VV <- list(NULL)
    if (is.null(modesnam)) 
        modesnam <- paste(rep("m", length(dim(X))), 1:length(dim(X)))
    if (!is.function(method) && method == "Presvd") 
        dim <- 1
    if (length(dim) == 1) 
        dim <- rep(dim, length(dim(X)))
    for (i in 1:length(dim(X))) {
        cci <- (1:length(dim(X)))[-i]
        if (is.function(method)) 
            VV[[i]] <- method(matrix(aperm(X, c(cci, i)), ncol = dim(X)[i]), 
                ...)
        else {
            if (method == "Presvd") 
                VV[[i]] <- PPMA(matrix(aperm(X, c(cci, i)), ncol = dim(X)[i]), 
                  pena = list(NULL, NULL))
            if (method == "svd") 
                VV[[i]] <- svd(matrix(aperm(X, c(cci, i)), ncol = dim(X)[i]))
        }
        if (dim[i] > dim(X)[i]) 
            dimi <- dim(X)[i]
        else dimi <- dim[i]
        VV[[i]]$d <- VV[[i]]$d[1:dimi]
        VV[[i]]$modesnam <- modesnam[[i]]
        VV[[i]]$n <- dimnames(X)[[i]]
        VV[[i]]$v <- t(VV[[i]]$v[, 1:dimi])
        if (dimi == 1) 
            VV[[i]]$v <- as.vector(VV[[i]]$v)
        VV[[i]]$u <- NULL
    }
    return(VV)
}
"PROJOT" <-
function (X, solu, numo = 1, bortho = TRUE, Ortho = TRUE, metrics = NULL) 
{
    txDy <- function(x, D, y) {
        if (!is.null(D)) {
            if (is.vector(D)) 
                y <- D * y
            if (is.matrix(D)) 
                y <- D %*% y
        }
        if (is.vector(x) & is.vector(y)) 
            return(sum(x * y))
        else return(t(x) %*% y)
    }
    projmat <- function(Y, x, D = NULL, bortho = TRUE) {
        if (is.vector(x)) 
            Y <- x %*% (1/txDy(x, D, x) * txDy(x, D, Y))
        if (is.matrix(x)) {
            if (!bortho) 
                Y <- x %*% Powmat(txDy(x, D, x), -1) %*% txDy(x, 
                  D, Y)
            else Y <- x %*% ((1/diag(txDy(x, D, x))) * txDy(x, 
                D, Y))
        }
        return(Y)
    }
    ldx <- length(dim(X))
    if (!is.list(numo)) 
        numo <- rep(list(numo), ldx)
    if (!is.list(bortho)) 
        bortho <- rep(list(bortho), ldx)
    if (!is.list(Ortho)) 
        Ortho <- rep(list(Ortho), ldx)
    if (!is.list(metrics)) 
        metrics <- rep(list(metrics), ldx)
    for (i in 1:ldx) {
        if (!is.null(numo[[i]])) {
            z <- solu[[i]]$v
            if (is.matrix(z)) {
                if (!dim(X)[i] == dim(z)[2]) 
                  stop("----WRONG DIMENSIONS----")
                else z <- z[numo[[i]], ]
            }
            else {
                if (!dim(X)[i] == length(z)) 
                  stop("----WRONG DIMENSIONS----")
            }
            lacola <- (1:length(dim(X)))[-i]
            laperm <- c(i, lacola)
            lacol <- (dim(X))[lacola]
            toconz <- matrix(aperm(X, laperm), ncol = prod(lacol))
            if (!is.vector(z)) 
                z <- t(z)
            PXz <- projmat(toconz, z, D = metrics[[i]], bortho = bortho[[i]])
            PXz <- array(PXz, c(dim(X)[i], lacol))
            if (i == ldx) 
                PXz <- aperm(PXz, c(2:ldx, 1))
            if (!i == 1 & !i == ldx) 
                PXz <- aperm(PXz, c(2:(i), 1, (i + 1):ldx))
            dimnames(PXz) <- dimnames(X)
            if (Ortho[[i]]) 
                X <- X - PXz
            else X <- PXz
        }
    }
    return(X)
}
"REBUILD" <-
function (solutions, nTens = 1:2, testvar = 1, redundancy = FALSE) 
{
    if (!is.list(solutions)) {
        stop(" should be a solutions object see PTA3")
    }
    ord <- length(solutions)
    if (as.character(solutions[[ord]]$method)[1] == "PCA") 
        REBUILDPCAn(solutions)
    else {
        tensfin <- 0
        deja <- NULL
        dejaTP <- NULL
        testpass <- length(nTens)
        for (cp in nTens) {
            if (100 * (solutions[[ord]]$d[cp]^2/solutions[[ord]]$ssX[1]) > 
                testvar) {
                if (!solutions[[ord]]$d[cp] %in% deja || (redundancy & 
                  substr(solutions[[ord]]$vsnam[cp], 2, 1) == 
                    "t")) {
                  if (!substr(solutions[[ord]]$vsnam[cp], 1, 
                    1) == "*") {
                    deja <- c(deja, solutions[[ord]]$d[cp])
                    dejaTP <- c(dejaTP, cp)
                  }
                  tens <- solutions[[1]]$v[cp, ] * solutions[[ord]]$d[cp]
                  names(tens) <- solutions[[1]]$n
                  for (d in 2:ord) {
                    atens <- solutions[[d]]$v[cp, ]
                    names(atens) <- solutions[[d]]$n
                    tens <- tens %o% atens
                  }
                  tensfin <- tensfin + tens
                }
            }
        }
        pcre <- 100 * sum(deja^2)/solutions[[ord]]$ssX[1]
        cat("-- Variance Percent rebuilt", solutions[[ord]]$datanam, 
            " at ", pcre, "% ", "\n")
        cat("-- MSE ", mean((eval(solutions[[ord]]$datanam) - 
            tensfin)^2), "\n")
        cat("-- with ", length(deja), " Principal Tensors out of ", 
            length(nTens), " given", "\n")
        if (pcre > 100) {
            cat("\n", "--WARNING !--- redundancy in choice of solutions to rebuild !!!", 
                "\n")
            print(pcre, digits = 20)
        }
        comp <- 100 - 100 * (sum(dim(eval(solutions[[ord]]$datanam))) + 
            1) * length(deja)/prod(dim(eval(solutions[[ord]]$datanam)))
        cat("-- compression    ", comp, " %", "\n")
        if (comp < 0) {
            cat("******no compression ....", "\n")
            dejadedans <- cbind(dejaTP, deja)
            rownames(dejadedans) <- solutions[[ord]]$vsnam[dejaTP]
            print(dejadedans)
        }
        return(tensfin)
    }
}
"RESUM" <-
function (solb, sola = NULL, numass = NULL, verbose = getOption("verbose"), 
    file = NULL, summary = FALSE, testvar = 0.1, not = NULL) 
{
    if (!is.null(sola)) {
        numlast <- length(sola[[length(sola)]]$d)
        if (is.null(numass)) 
            num <- numlast
        if (!is.null(numass)) 
            num <- numass
        for (i in 1:length(solb)) {
            for (j in 1:length(sola)) {
                if (as.vector(solb[[i]]$modesnam) == as.vector(sola[[j]]$modesnam)) {
                  sola[[j]]$v <- rbind(sola[[j]]$v, solb[[i]]$v)
                  if ("iter" %in% names(solb[[i]])) 
                    sola[[j]]$iter <- c(sola[[j]]$iter, solb[[i]]$iter)
                  if ("test" %in% names(solb[[i]])) 
                    sola[[j]]$test <- c(sola[[j]]$test, solb[[i]]$test)
                }
            }
        }
        for (k in 1:length(sola)) {
            if (is.matrix(sola[[k]]$v)) {
                if ((dim(sola[[k]]$v)[[1]]) == numlast) {
                  sola[[k]]$v <- rbind(sola[[k]]$v, rep(1, length(solb[[length(solb)]]$d)) %x% 
                    t(sola[[k]]$v[num, ]))
                }
            }
            if (!is.matrix(sola[[k]]$v)) {
                sola[[k]]$v <- rbind(sola[[k]]$v, rep(1, length(solb[[length(solb)]]$d)) %x% 
                  t(sola[[k]]$v))
            }
        }
        sola[[k]]$d <- c(sola[[k]]$d, solb[[i]]$d)
        sola[[k]]$pct <- c(sola[[k]]$pct, solb[[i]]$pct)
        sola[[k]]$ssX <- c(sola[[k]]$ssX, solb[[i]]$ssX)
        if ("smoocheck" %in% names(solb[[i]])) 
            sola[[k]]$smoocheck <- cbind(sola[[k]]$smoocheck, 
                solb[[i]]$smoocheck)
        for (m in 1:length(solb[[i]]$vsnam)) {
            for (n in 1:length(sola[[k]]$vsnam)) {
                if ((round(sola[[k]]$d[n], digits = 10) == round(solb[[i]]$d[m], 
                  digits = 10)) & (!substr(solb[[i]]$vsnam[m], 
                  1, 1) == "*")) {
                  solb[[i]]$vsnam[m] <- paste("*t", solb[[i]]$vsnam[m], 
                    sep = "")
                }
            }
        }
        sola[[k]]$vsnam <- c(sola[[k]]$vsnam, solb[[i]]$vsnam)
    }
    else {
        sola <- solb
        k <- length(sola)
    }
    pctota <- (100 * (sola[[k]]$d)^2)/sola[[k]]$ssX[1]
    if (verbose & !summary) {
        cat("                ------Percent Rebuilt from Selected ----", 
            sum(pctota[!substr(sola[[k]]$vsnam, 1, 1) == "*" & 
                pctota > testvar]), "%", "\n")
    }
    if (!is.null(file)) {
        cat("               ------Percent Rebuilt from Slected ----", 
            sum(pctota[!substr(sola[[k]]$vsnam, 1, 1) == "*" & 
                pctota > testvar]), "%", "\n", file = file, append = TRUE)
        if (verbose) {
            sink(file = file, append = TRUE)
            summ <- as.matrix(cbind(1:length(sola[[k]]$d), sola[[k]]$d, 
                sola[[k]]$ssX, sola[[k]]$pct, pctota))
            dimnames(summ) <- list(sola[[k]]$vsnam, c("-no-", 
                "--Sing Val--", "--ssX--", "--local Pct--", "--Global Pct--"))
            summ <- summ[pctota > testvar, ]
            print(summ, digits = 5)
            sink()
        }
    }
    if (summary) {
        if (is.null(not)) 
            not <- TRUE
        cat("\n", "++++ PTA- ", k, "modes ++++ ", "\n")
        di <- NULL
        for (r in 1:length(sola)) di <- c(di, length(sola[[r]]$v[1, 
            ]))
        nostar <- !substr(sola[[k]]$vsnam, 1, 1) == "*"
        cat("               data= ", deparse(sola[[k]]$datanam), 
            " ", di, "\n")
        cat("   ", sola[[k]]$addedcomment, "\n")
        cat("                ------Percent Rebuilt----", sum(pctota[nostar]), 
            "%", "\n")
        summ <- matrix(cbind(1:length(sola[[k]]$d), sola[[k]]$d, 
            sola[[k]]$ssX, sola[[k]]$pct, pctota), ncol = 5)
        summ <- summ[(pctota > testvar) & not, ]
        summ <- matrix(summ, ncol = 5)
        dimnames(summ) <- list(sola[[k]]$vsnam[pctota > testvar & 
            not], c("-no-", "--Sing Val--", "--ssX--", "--local Pct--", 
            "--Global Pct--"))
        sumex <- sum(summ[!substr(dimnames(summ)[[1]], 1, 1) == 
            "*", 5])
        cat("                ------Percent Rebuilt from Selected ----", 
            sumex, "%", "\n")
        print(summ, digits = 5)
        cat("\n", "++++               ++++", "\n")
        if (!is.null(testvar) & !testvar == 0) 
            cat(" Shown are selected over ", length(sola[[k]]$vsnam[nostar]), 
                " PT  with var>", testvar, "% total", "\n")
        else cat(" over ", length(sola[[k]]$vsnam[nostar]), " PT ", 
            "\n")
    }
    else invisible(sola)
}
"TENSELE" <-
function (T, moins = NULL, asarray = TRUE, order = NULL, id = NULL) 
{
    dim <- NULL
    if (is.null(order)) 
        order <- length(T):1
    if (is.list(id)) 
        asarray <- TRUE
    vu <- 0
    for (i in order) {
        if (!(i %in% moins)) {
            if (is.null(id[[i]])) 
                Tv <- T[[i]]$v
            else Tv <- T[[i]]$v[id[[i]], ]
            if (vu == 0) 
                tensel <- Tv
            if (!vu == 0) {
                if (asarray) {
                  tensel <- Tv %o% tensel
                }
                if (!asarray) {
                  tensel <- as.vector(tensel %x% Tv)
                  dim <- c(length(Tv), dim)
                }
            }
            vu <- 1
        }
    }
    return(tensel)
}
"is.solutions.PTAk" <-
function (solut) 
{
    if (isol <- is.list(solut)) {
        kor <- length(solut)
        if (is.vector(solut[[kor]]$d)) {
            nbTens <- length(solut[[kor]]$d)
            for (i in 1:kor) {
                if (is.list(solut[[i]])) {
                  if ((nbTens == 1 & is.vector(solut[[i]]$v)) || 
                    (length(solut[[kor]]$d) > 1 & ifelse(is.matrix(solut[[i]]$v), 
                      dim(solut[[i]]$v)[1] == nbTens, FALSE))) {
                    if (!is.character(solut[[i]]$modesnam)) 
                      isol <- FALSE
                  }
                }
            }
            if (isol) {
                if (all(c(is.vector(solut[[kor]]$pct), is.vector(solut[[kor]]$ssX), 
                  is.vector(solut[[kor]]$vsnam), is.character(solut[[kor]]$vsnam), 
                  is.array(eval(solut[[kor]]$datanam)), is.call(solut[[kor]]$method)))) {
                  if (!length(solut[[kor]]$pct) == nbTens || 
                    !length(solut[[kor]]$ssX) == nbTens || !length(solut[[kor]]$vsnam) == 
                    nbTens) 
                    isol <- FALSE
                }
                else isol <- FALSE
            }
        }
    }
    return(isol)
}
"summary.solutions.FCAk" <-
function (sola, testvar = 0.5, dontshow = "*") 
{
    nostar <- (!substr(sola[[length(sola)]]$vsnam, 1, 1) == "*")
    if (dontshow == "*") 
        show <- (!substr(sola[[length(sola)]]$vsnam, 1, 1) == 
            "*")
    else if (!is.null(dontshow)) 
        show <- dontshow & nostar
    else show <- TRUE
    k <- length(sola)
    pctotafc <- (100 * (sola[[k]]$d)^2)/(sola[[k]]$ssX[1] - 1)
    pctota <- (100 * (sola[[k]]$d)^2)/sola[[k]]$ssX[1]
    cat("\n", "++++ FCA- ", k, "modes++++ ", "\n")
    di <- NULL
    for (r in 1:length(sola)) di <- c(di, length(sola[[r]]$v[1, 
        ]))
    cat("     ++ Contingency Table ", deparse(sola[[k]]$datanam), 
        " ", di, " ++", "\n")
    cat("   ", sola[[k]]$addedcomment, "\n")
    cat("                -----Total Percent Rebuilt----", sum(pctota[nostar]), 
        "%", "\n")
    cat("     ++ Percent of lack of complete independence rebuilt  ++ ", 
        sum(pctotafc[show][-1]), "%", "\n")
    cat("                                    selected pctoafc > ", 
        testvar, "%  total= ", sum(pctotafc[show & pctotafc > 
            testvar][-1]), "\n")
    summ <- matrix(cbind(1:length(sola[[k]]$d), sola[[k]]$d, 
        sola[[k]]$ssX, pctota, pctotafc), ncol = 5)
    summ <- summ[pctotafc > testvar & show, ]
    summ <- matrix(summ, ncol = 5)
    summ[1, 5] <- NA
    dimnames(summ) <- list(sola[[k]]$vsnam[pctotafc > testvar & 
        show], c("-no-", "--Sing Val--", "--ssX--", "--Global Pct--", 
        "--FCA--"))
    print(summ, digits = 5)
    cat("\n", "++++               ++++", "\n")
    if (!is.null(testvar) & !testvar == 0) 
        cat(" Shown are selected  over ", length(sola[[k]]$vsnam[nostar]) - 
            1, " PT  with pct AFC >", testvar, "% ", "\n")
    else cat(" over ", length(sola[[k]]$vsnam), " PT (with*)", 
        "\n")
}
"summary.solutions.PTAk" <-
function (solution, testvar = 1, dontshow = "*") 
{
    if (!is.null(dontshow)) 
        if (dontshow == "*") 
            dontshow <- (!substr(solution[[length(solution)]]$vsnam, 
                1, 1) == "*")
    RESUM(solution, summary = TRUE, testvar = testvar, not = dontshow)
}
