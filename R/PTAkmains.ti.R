"APSOLU3" <-
function (X, solu, pt3 = NULL, nbPT2 = 1, smoothing = FALSE,
    smoo = list(NA), verbose = getOption("verbose"), file = NULL)
{
    if (is.list(X)) {
        if (is.list(X$met))
            metrics <- TRUE
        else stop(paste("------with metrics X must be a list with $data and $met----"))
    }
    else metrics <- FALSE
    if (metrics) {
        nam <- dimnames(X$data)
        for (d in 1:3) {
            if (length(X$met[[d]]) > 1) {
                if (length(X$met[[d]]) == dim(X$data)[d]^2) {
                  tempp <- d
                  t12 <- CONTRACTION(X$data, Powmat(X$met[[d]],
                    1/2), Xwiz = d, zwiX = 1)
                  d <- tempp
                  lacola <- (1:3)[-d]
                  laperm <- c(lacola, d)
                }
                else {
                  lacola <- (1:3)[-d]
                  laperm <- c(d, lacola)
                  lacol <- (dim(X$data))[lacola]
                  pt12 <- matrix(aperm(X$data, laperm), ncol = prod(lacol))
                  t12 <- sqrt(X$met[[d]]) * pt12
                }
                t12 <- array(t12, (dim(X$data))[laperm])
                X$data <- aperm(t12, match(1:3, laperm))
            }
            else X$data <- X$data * X$met[[d]]
        }
        met <- X$met
        X <- X$data
        dimnames(X) <- nam
        for (d in 1:length(solu)) {
            if (length(met[[d]]) > 1) {
                if (length(met[[d]]) == dim(X)[d]^2) {
                  solu[[d]]$v <- solu[[d]]$v %*% Powmat(met[[d]],
                    1/2)
                }
                else {
                  solu[[d]]$v <- t(sqrt(met[[d]]) * t(solu[[d]]$v))
                }
            }
            else solu[[d]]$v <- solu[[d]]$v * sqrt(met[[d]])
        }
    }
    if (!is.null(pt3))
        numsol <- pt3
    else numsol <- length(solu[[length(solu)]]$d)
    Zsol <- list(NULL, NULL)
    if (smoothing & !length(smoo) == 3)
        stop(paste("--- Smoothing list must be of length 3! ---"))
    for (i in 1:3) {
        if (verbose) {
            cat(" ----------APSOLU3------------", file = ifelse(is.null(file),
                "", file), append = TRUE)
            cat(" ---- Associated solution to entry ---", i,
                file = ifelse(is.null(file), "", file), append = TRUE)
            cat("  ....  of dimension: ", dim(X)[i], "\n", file = ifelse(is.null(file),
                "", file), append = TRUE)
        }
        tracei <- i
        Z <- CONTRACTION(X, if (is.matrix(solu[[i]]$v))
            solu[[i]]$v[numsol, ]
        else solu[[i]]$v, Xwiz = i, zwiX = ifelse(length(numsol) ==
            1, 1, 2))
        i <- tracei
        if (nbPT2 == 1)
            nomb <- min(dim(Z))
        else nomb <- min(dim(Z), nbPT2)
        if (smoothing == TRUE)
            solq <- svdsmooth(Z, nomb = nomb, smooth = smoo[-i])
        else solq <- svd(Z)
        Zsol[[1]]$modesnam <- solu[[((1:3)[-i])[1]]]$modesnam
        nomb <- min(nomb, length(solq$d))
        Zsol[[1]]$v <- t(solq$u[, 1:nomb])
        Zsol[[2]]$modesnam <- solu[[((1:3)[-i])[2]]]$modesnam
        Zsol[[2]]$v <- t(solq$v[, 1:nomb])
        if (smoothing == TRUE) {
            Zsol[[2]]$smoocheck <- array(NA, c(3, nomb))
            Zsol[[2]]$smoocheck[(1:3)[-i], ] <- solq$smoocheck
        }
        if (!dim(Z)[2] > length(solq$d))
            ssX <- sum(solq$d^2)
        else ssX <- sum(as.vector(Z)^2)
        Zsol[[2]]$d <- solq$d[1:nomb]
        Zsol[[2]]$pct <- (100 * (solq$d^2)/ssX)[1:nomb]
        Zsol[[2]]$ssX <- rep(ssX, nomb)
        Zsol[[2]]$vsnam <- c(paste("*", dim(X)[i], solu[[length(solu)]]$vsnam[numsol],
            dim(X)[-i][1], dim(X)[-i][2], sep = ""), rep(paste(dim(X)[i],
            solu[[length(solu)]]$vsnam[numsol], dim(X)[-i][1],
            dim(X)[-i][2]), nomb - 1))
        solu <- RESUM(Zsol, solu, numass = numsol, verbose = verbose,
            file = file)
    }
    if (metrics) {
        for (d in 1:length(solu)) {
            if (length(met[[d]]) > 1) {
                if (length(met[[d]]) == dim(X)[d]^2) {
                  solu[[d]]$v <- solu[[d]]$v %*% Powmat(met[[d]],
                    -1/2)
                }
                else {
                  solu[[d]]$v <- t(1/sqrt(met[[d]]) * t(solu[[d]]$v))
                }
            }
            else solu[[d]]$v <- solu[[d]]$v * 1/sqrt(met[[d]])
        }
    }
    return(solu)
}
"APSOLUk" <-
function (X, solu, nbPT, nbPT2 = 1, smoothing = FALSE, smoo = list(NA),
    minpct = 0.1, ptk = NULL, verbose = getOption("verbose"),
    file = NULL, modesnam = NULL)
{
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
        for (d in 1:length(solu)) {
            if (length(met[[d]]) > 1) {
                if (length(met[[d]]) == dim(X)[d]^2) {
                  solu[[d]]$v <- solu[[d]]$v %*% Powmat(met[[d]],
                    1/2)
                }
                else {
                  solu[[d]]$v <- t(sqrt(met[[d]]) * t(solu[[d]]$v))
                }
            }
            else solu[[d]]$v <- solu[[d]]$v * sqrt(met[[d]])
        }
    }
    numsol <- length(solu[[length(solu)]]$d)
    if (!is.null(ptk))
        numsol <- ptk
    Zsol <- list(NULL, NULL)
    kor <- length(dim(X))
    if (smoothing & !length(smoo) == kor)
        stop(paste("--- Smoothing list must be of length ", kor,
            "! ---"))
    for (i in 1:kor) {
        if (verbose) {
            cat("\n", "\n", "            ++++++++++++++++ --APSOLUk-- ",
                file = ifelse(is.null(file), "", file), append = TRUE)
            cat(solu[[kor]]$vsnam[numsol], " Associated solution to entry ---",
                i, file = ifelse(is.null(file), "", file), append = TRUE)
            cat("  ....  of dimension: ", dim(X)[i], "\n", file = ifelse(is.null(file),
                "", file), append = TRUE)
        }
        tracei <- i
        Z <- CONTRACTION(X, matrix(solu[[i]]$v, ncol = dim(X)[i])[numsol,
            ], Xwiz = i)
        i <- tracei
        if (length(dim(Z)) == 3) {
            solZ <- PTA3(Z, nbPT = nbPT[1], nbPT2 = nbPT2, smoothing = smoothing,
                smoo = smoo[-i], minpct = minpct, verbose = verbose,
                file = file, modesnam = modesnam[-i])
        }
        if (length(dim(Z)) > 3) {
            solZ <- PTAk(Z, nbPT = nbPT, nbPT2 = nbPT2, smoothing = smoothing,
                smoo = smoo[-i], minpct = minpct, verbose = verbose,
                file = file, modesnam = modesnam[-i])
        }
        nno <- length(solZ[[length(solZ)]]$vsnam)
        for (n in 1:nno) {
            if (!substr(solZ[[length(solZ)]]$vsnam[n], 1, 1) ==
                "*") {
                solZ[[length(solZ)]]$vsnam[n] <- paste(dim(X)[i],
                  solZ[[length(solZ)]]$vsnam[n], sep = "-")
            }
            else {
                solZ[[length(solZ)]]$vsnam[n] <- paste("*", paste(dim(X)[i],
                  substr(solZ[[length(solZ)]]$vsnam[n], 2, 100),
                  sep = "-"), sep = "")
            }
        }
        solZ[[length(solZ)]]$vsnam[1] <- paste("*", solZ[[length(solZ)]]$vsnam[1],
            sep = "")
        if (smoothing == TRUE) {
            smooche <- solZ[[length(solZ)]]$smoocheck
            solZ[[length(solZ)]]$smoocheck <- array(NA, c(kor,
                nno))
            solZ[[length(solZ)]]$smoocheck[(1:kor)[-i], ] <- smooche
        }
        solu <- RESUM(solZ, solu, numass = numsol, verbose = verbose,
            file = file)
    }
    if (metrics) {
        for (d in 1:length(solu)) {
            if (length(met[[d]]) > 1) {
                if (length(met[[d]]) == dim(X)[d]^2) {
                  solu[[d]]$v <- solu[[d]]$v %*% Powmat(met[[d]],
                    -1/2)
                }
                else {
                  solu[[d]]$v <- t(1/sqrt(met[[d]]) * t(solu[[d]]$v))
                }
            }
            else solu[[d]]$v <- solu[[d]]$v * 1/sqrt(met[[d]])
        }
    }
    return(solu)
}
"FCAk" <-
function (X, nbPT = 3, nbPT2 = 1, minpct = 0.01, smoothing = FALSE,
    smoo = rep(list(function(u) ksmooth(1:length(u), u, kernel = "normal",
        bandwidth = 3, x.points = (1:length(u)))$y), length(dim(X))),
    verbose = getOption("verbose"), file = NULL, modesnam = NULL,
    addedcomment = "", chi2 = TRUE, E = NULL)
{
    ldx <- length(dim(X))
    if (verbose) {
        cat("\n", "       ----------+++++++++++------------",
            "\n", ifelse(smoothing, paste("Penalised ", "\n"),
                ""), "       Correspondence Analysis on ", ldx,
            " modes ", "\n", file = ifelse(is.null(file), "",
                file), append = TRUE)
        if (!is.null(modesnam))
            cat("modes are ", modesnam, "\n", file = ifelse(is.null(file),
                "", file), append = TRUE)
        cat("       ----------+++++++++++------------", "\n",
            fil = ifelse(is.null(file), "", file), append = TRUE)
        cat("Data   = complete independence    + lack of independence ...",
            "\n", file = ifelse(is.null(file), "", file), append = TRUE)
        cat(" lack of independence = partial independence + lack of independence ... etc ...",
            "\n", file = ifelse(is.null(file), "", file), append = TRUE)
    }
    Y <- FCAmet(X, chi2 = chi2, E = E)
    solutions <- PTAk(Y, nbPT = nbPT, nbPT2 = nbPT2, smoothing = smoothing,
        smoo = smoo, minpct = minpct, verbose = verbose, file = file,
        modesnam = modesnam, addedcomment = addedcomment)
    solutions[[ldx]]$datanam <- substitute(X)
    solutions[[ldx]]$method <- match.call()
    solutions[[ldx]]$addedcomment <- addedcomment
    class(solutions) <- c("solutions.FCAk", "solutions.PTAk")
    invisible(solutions)
}
"PPMA" <-
function (X, test = 1e-10, pena = list(function(u) ksmooth(1:length(u),
    u, kernel = "normal", bandwidth = 3, x.points = (1:length(u)))$y,
    NA), ini = mean, vsmin = 1e-20, Maxiter = 2000)
{
    v0 <- apply(X, 2, FUN = ini)
    if (all(v0 < 1e-04)) {
        v0 <- (X[sample(1:dim(X)[1], 1), ])
        if (max(abs(X)) < test * 1e-08) {
            cat(" Sum of squares veryyyy smallll  .......", "\n")
            return(list(u = rep(0, dim(X)[1]), v = rep(0, dim(X)[2]),
                d = 0, iter = 0, test = NA))
        }
    }
    test0 <- 1
    ite <- 1
    while (test0 > test) {
        u <- as.vector(X %*% v0)
        if (is.function(pena[[1]]))
            u <- pena[[1]](u)
        d <- sqrt(u %*% u)
        if (!d == 0)
            u <- u/d
        v <- as.vector(u %*% X)
        if (is.function(pena[[2]]))
            v <- pena[[2]](v)
        d <- sqrt(v %*% v)
        if (!d == 0)
            v <- v/d
        if (test0 == 1)
            u0 <- u
        if (!d < vsmin)
            test0 <- sum((u - u0)^2) + sum((v - v0)^2)
        else test0 <- 0
        v0 <- v
        u0 <- u
        ite <- ite + 1
        if (ite > (Maxiter - 1) && (ite - Maxiter)%%200 == 0) {
            cat("\n \n \n \n \n ", " WARNING ****** Iteration already =  ",
                ite, "test= ", test0, "\n")
            cat(" ** type  999  to STOP ** just RETURN to carry on **",
                "\n")
            cat(" or type a new test value initial was", test,
                "\n")
            conti <- scan("", n = 1, quiet = TRUE, flush = TRUE)
            if (length(conti) > 0) {
                if (conti == 999)
                  stop(paste(" ---- Aborted by request ---- "))
                if (is.numeric(conti))
                  test <- conti
            }
        }
    }
    return(list(u = as.matrix(u), v = as.matrix(v), d = as.vector(d),
        iter = ite, test = test0))
}
"PTA3" <-
function (X, nbPT = 2, nbPT2 = 1, smoothing = FALSE, smoo = list(function(u) ksmooth(1:length(u),
    u, kernel = "normal", bandwidth = 4, x.points = (1:length(u)))$y,
    function(u) smooth.spline(u, df = 3)$y, NA), minpct = 0.1,
    verbose = getOption("verbose"), file = NULL, modesnam = NULL,
    addedcomment = "")
{
    datanam <- substitute(X)
    if (is.list(X)) {
        if (is.list(X$met))
            metrics <- TRUE
        else stop(paste("------with metrics X must be a list with $data and $met----"))
    }
    else metrics <- FALSE
    if (metrics) {
        nam <- dimnames(X$data)
        for (d in 1:3) {
            if (length(X$met[[d]]) > 1) {
                if (length(X$met[[d]]) == dim(X$data)[d]^2) {
                  tempp <- d
                  t12 <- CONTRACTION(X$data, Powmat(X$met[[d]],
                    1/2), Xwiz = d, zwiX = 1)
                  d <- tempp
                  lacola <- (1:3)[-d]
                  laperm <- c(lacola, d)
                }
                else {
                  lacola <- (1:3)[-d]
                  laperm <- c(d, lacola)
                  lacol <- (dim(X$data))[lacola]
                  pt12 <- matrix(aperm(X$data, laperm), ncol = prod(lacol))
                  t12 <- sqrt(X$met[[d]]) * pt12
                }
                t12 <- array(t12, (dim(X$data))[laperm])
                X$data <- aperm(t12, match(1:3, laperm))
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
                ""), "              PTA 3modes ", "\n", file = ifelse(is.null(file),
                "", file), append = TRUE)
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
    solutions <- NULL
    if (smoothing) {
        if (smoothing & !length(smoo) == 3)
            stop(paste("--- Smoothing list must be of length 3! ---"))
        for (j in 1:3) if (!is.list(smoo[[j]]))
            smoo[[j]] <- list(smoo[[j]])
    }
    for (t in 1:nbPT) {
        if (verbose) {
            cat("----- Principal Tensor ---- ", paste("vs", pass.(t,
                3), sep = ""), file = ifelse(is.null(file), "",
                file), append = TRUE)
        }
        if (smoothing) {
            if (t > 1) {
                for (j in 1:3) if (length(smoo[[j]]) == t - 1)
                  smoo[[j]][[t]] <- smoo[[j]][[t - 1]]
            }
            tosmoo <- list(smoo[[1]][[t]], smoo[[2]][[t]], smoo[[3]][[t]])
        }
        else tosmoo <- list(NA)
        solut <- SINGVA(X, verbose = verbose, file = file,
            PTnam = paste("vs", pass.(t, 3), sep = ""),
            smoothing = smoothing, smoo = tosmoo, modesnam = modesnam)
        if (is.null(solutions) & verbose)
            cat(" --- GLobal Percent --- ", (100 * solut[[3]]$d^2)/solut[[3]]$ssX[1],
                "%", "\n", file = ifelse(is.null(file), "", file),
                append = TRUE)
        if (verbose & !is.null(solutions)) {
            cat("                 -- GLobal Percent -- ", (100 *
                solut[[3]]$d^2)/solutions[[3]]$ssX[1], "%", "\n",
                file = ifelse(is.null(file), "", file), append = TRUE)
        }
        if (!is.null(solutions)) {
            if (100 * solut[[length(solut)]]$d^2/solutions[[length(solutions)]]$ssX[1] <
                minpct) {
                cat("\n", "\n", " ++ Last 3-modes vs < ", minpct,
                  "% stopping this level and under ++", "\n")
                solutions <- RESUM(solut, solutions, verbose = verbose,
                  file = file)
                break
            }
        }
        if (nbPT2 >= 1)
            solut <- APSOLU3(X, solut, pt3 = NULL, nbPT2 = nbPT2,
                smoothing = smoothing, smoo = tosmoo, verbose = verbose,
                file = file)
        if (verbose)
            cat("\n", "+++ PTA 3modes  ------After ---", paste("vs",
                pass.(t, 3), sep = ""), file = ifelse(is.null(file),
                "", file), append = TRUE)
        solutions <- RESUM(solut, solutions, verbose = verbose,
            file = file)
        if (t < nbPT)
            X <- PROJOT(X, solut)
    }
    if (metrics) {
        for (d in 1:3) {
            if (length(met[[d]]) > 1) {
                if (length(met[[d]]) == dim(X)[d]^2) {
                  solutions[[d]]$v <- solutions[[d]]$v %*% Powmat(met[[d]],
                    -1/2)
                }
                else {
                  solutions[[d]]$v <- t(1/sqrt(met[[d]]) * t(solutions[[d]]$v))
                }
            }
            else solutions[[d]]$v <- solutions[[d]]$v * 1/sqrt(met[[d]])
        }
    }
    solutions[[3]]$method <- match.call()
    solutions[[3]]$addedcomment <- addedcomment
    solutions[[length(solutions)]]$datanam <- datanam
    cat("\n", "-----Execution Time-----", (proc.time() - debtime)[3],
        "\n")
    class(solutions) <- c("solutions.PTAk")
    invisible(solutions)
}
"PTAk" <-
function (X, nbPT = 2, nbPT2 = 1, minpct = 0.1, smoothing = FALSE,
    smoo = list(, NA), verbose = getOption("verbose"), file = NULL,
    modesnam = NULL, addedcomment = "")
{
    datanam <- substitute(X)
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
            paste("Penalised ", "\n"), ""), " Principal Tensor Analysis on k modes ",
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
    kor <- length(dim(X))
    if (length(nbPT) < kor - 2) {
        nbPT <- rep(nbPT[1], kor - 2)
    }
    if (is.null(modesnam)) {
        modesnam <- paste(rep("mo", kor), 1:kor)
    }
    solutions <- NULL
    if (smoothing) {
        if (smoothing & !length(smoo) == kor)
            stop(paste("--- Smoothing list must be of length ",
                kor, "! ---"))
        for (j in 1:kor) if (!is.list(smoo[[j]]))
            smoo[[j]] <- list(smoo[[j]])
    }
    for (t in 1:(nbPT[kor - 2])) {
        if (verbose)
            cat("\n", "\n", "                ++++++  k-modes Solutions  ---- k=",
                kor, paste(", vs", pass.(t, kor), sep = ""),
                "++++++", "\n", "\n", file = ifelse(is.null(file),
                  "", file), append = TRUE)
        tosmoo <- list(NA)
        if (smoothing) {
            if (t > 1) {
                for (j in 1:kor) if (length(smoo[[j]]) == t -
                  1)
                  smoo[[j]][[t]] <- smoo[[j]][[t - 1]]
            }
            if (kor == 3)
                tosmoo <- list(smoo[[1]][[t]], smoo[[2]][[t]],
                  smoo[[3]][[t]])
            if (kor == 4)
                tosmoo <- list(smoo[[1]][[t]], smoo[[2]][[t]],
                  smoo[[3]][[t]], smoo[[4]][[t]])
            if (kor == 5)
                tosmoo <- list(smoo[[1]][[t]], smoo[[2]][[t]],
                  smoo[[3]][[t]], smoo[[4]][[t]], smoo[[5]][[t]])
            if (kor == 6)
                tosmoo <- list(smoo[[1]][[t]], smoo[[2]][[t]],
                  smoo[[3]][[t]], smoo[[4]][[t]], smoo[[5]][[t]],
                  smoo[[6]][[t]])
            if (kor == 7)
                tosmoo <- list(smoo[[1]][[t]], smoo[[2]][[t]],
                  smoo[[3]][[t]], smoo[[4]][[t]], smoo[[5]][[t]],
                  smoo[[6]][[t]], smoo[[7]][[t]])
            if (kor == 8)
                tosmoo <- list(smoo[[1]][[t]], smoo[[2]][[t]],
                  smoo[[3]][[t]], smoo[[4]][[t]], smoo[[5]][[t]],
                  smoo[[6]][[t]], smoo[[7]][[t]], smoo[[8]][[t]])
        }

        solut <- SINGVA(X, verbose = verbose, file = file,
            PTnam = paste("vs", pass.(t, kor), sep = ""), 
            smoothing = smoothing, smoo = tosmoo, modesnam = modesnam)
        if (is.null(solutions) & verbose)
            cat("                 -- GLobal Percent -- ", solut[[kor]]$pct,
                "%", "\n", file = ifelse(is.null(file), "", file),
                append = TRUE)
        if (verbose & !is.null(solutions)) {
            cat("                 -- GLobal Percent -- ", (100 *
                solut[[length(solut)]]$d^2)/solutions[[length(solutions)]]$ssX[1],
                "%", "\n", file = ifelse(is.null(file), "", file),
                append = TRUE)
        }
        if (!is.null(solutions)) {
            if (100 * solut[[length(solut)]]$d^2/solutions[[length(solutions)]]$ssX[1] <
                minpct) {
                cat("\n", "\n", " ++ Last ", kor, "-modes vs < ",
                  minpct, "% stopping this level and under ++",
                  "\n", file = ifelse(is.null(file), "", file),
                  append = TRUE)
                solutions <- RESUM(solut, solutions, verbose = verbose,
                  file = file)
                break
            }
        }
        if (kor - 3 > 0) {
            if (!nbPT[kor - 3] == 0) {
                solut <- APSOLUk(X, solut, nbPT = nbPT, nbPT2 = nbPT2,
                  smoothing = smoothing, smoo = tosmoo, minpct = minpct,
                  ptk = NULL, verbose = verbose, file = file,
                  modesnam = modesnam)
            }
        }
        if (kor == 3 & nbPT2 >= 1) {
            ptk <- NULL
            solut <- APSOLU3(X, solut, pt3 = ptk, nbPT2 = nbPT2,
                smoothing = smoothing, smoo = tosmoo, verbose = verbose,
                file = file)
        }
        solutions <- RESUM(solut, solutions, verbose = verbose,
            file = file)
        if (is.null(solutions[[length(solutions)]]$datanam))
            solutions[[length(solutions)]]$datanam <- datanam
        if (t < nbPT[kor - 2])
            X <- PROJOT(X, solut)
    }
    if (metrics) {
        for (d in 1:length(solutions)) {
            if (length(met[[d]]) > 1) {
                if (length(met[[d]]) == dim(X)[d]^2) {
                  solutions[[d]]$v <- solutions[[d]]$v %*% Powmat(met[[d]],
                    -1/2)
                }
                else {
                  solutions[[d]]$v <- t(1/sqrt(met[[d]]) * t(solutions[[d]]$v))
                }
            }
            else solutions[[d]]$v <- solutions[[d]]$v * 1/sqrt(met[[d]])
        }
    }
    solutions[[kor]]$method <- match.call()
    solutions[[kor]]$addedcomment <- addedcomment
    cat("-----Execution Time-----", (proc.time() - debtime)[3],
        "\n")
    class(solutions) <- c("solutions.PTAk")
    invisible(solutions)
}
"SINGVA" <-
function (X, test = 1e-12, PTnam = "vs111", Maxiter = 2000, verbose = getOption("verbose"),
    file = NULL, smoothing = FALSE, smoo = list(NA), modesnam = NULL,
    Ini = "Presvd", sym = NULL)
{
    datanam <- substitute(X)
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
    if (!is.array(X)) {
        stop(paste("--- X must be an array  ! ---"))
    }
    ord <- length(dim(X))
    if (verbose) {
        cat("\n", "       ----------+++++++++++------------ RPVSCC algorithm ",
            "\n", file = ifelse(is.null(file), "", file), append = TRUE)
        cat("                             ------------ Singular Value  ",
            PTnam, "\n", file = ifelse(is.null(file), "", file),
            append = TRUE)
        cat("                                       ----  dimensions:  ",
            dim(X), "\n", file = ifelse(is.null(file), "", file),
            append = TRUE)
    }
    sval0 <- INITIA(X, modesnam = modesnam, method = Ini)
    if (!is.null(sym)) {
        if (!ord == length(sym))
            stop(paste("--- Wrong length for parameter sym ! ---"))
        for (i in 1:ord) {
            if (!i == sym[i])
                sval0[[i]] <- sval0[[sym[i]]]
        }
    }
    if (verbose) {
        cat(" ----------------------", "\n", "Initialisation  done",
            "\n", file = ifelse(is.null(file), "", file), append = TRUE)
    }
    sval <- sval0
    if (smoothing) {
        sval[[ord]]$smoocheck <- array(FALSE, c(ord, 1))
        if (!length(smoo) == ord)
            stop(paste("--- Smoothing list must be of length ",
                ord, "! ---"))
        for (i in 1:ord) smoo[[i]] <- toplist(smoo[[i]])
    }
    test0 <- 1
    atest <- 0
    iter <- 0
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
                sval[[i]]$d <- NULL
            }
            tracei <- i
            v <- CONTRACTION.list(X, sval0, moins = i)
            i <- tracei
            if (smoothing) {
                if (is.function(smoo[[i]])) {
                  v <- smoo[[i]](as.vector(v))
                  sval[[ord]]$smoocheck[i, 1] <- TRUE
                }
            }
            sval[[i]]$v <- as.vector(v)
            sval0[[i]]$v <- as.vector(sval0[[i]]$v)
            sd <- sqrt(sval[[i]]$v %*% sval[[i]]$v)
            if (verbose & iter%%100 == 0) {
                cat(" --", sd, file = ifelse(is.null(file), "",
                  file), append = TRUE)
            }
            if (!sd == 0)
                sval[[i]]$v <- (sval[[i]]$v)/sd
            atest <- atest + (sval[[i]]$v - sval0[[i]]$v) %*%
                (sval[[i]]$v - sval0[[i]]$v)
            if (!is.null(sym)) {
                for (i in ord:1) {
                  if (!i == sym[i])
                    sval[[sym[i]]] <- sval[[i]]
                }
            }
            sval0 <- sval
        }
        test0 <- sqrt(atest)
        atest <- 0
        if (verbose & (iter%%100) == 1) {
            cat("\n", "----------- test =         ", test0, "\n",
                file = ifelse(is.null(file), "", file), append = TRUE)
        }
        if (iter > (Maxiter - 1) && (iter - Maxiter)%%200 == 0) {
            cat("\n \n \n \n \n ", " WARNING ****** Iteration already =  ",
                iter, "test= ", test0, "\n")
            cat(" ** type  999  to STOP ** just RETURN to carry on **",
                "\n")
            cat(" or type a new test value initial was", test,
                "\n")
            conti <- scan("", n = 1, quiet = TRUE, flush = TRUE)
            if (length(conti) > 0) {
                if (conti == 999)
                  stop(paste(" ---- Aborted by request ---- "))
                if (is.numeric(conti))
                  test <- conti
            }
        }
        sval0 <- sval
    }
    ssX <- sum(X^2)
    sstens <- sd^2
    totPCT <- 100 * sstens/ssX
    if (!verbose) {
        cat(" ---Final iteration--- ", iter, "\n")
        cat(" --Singular Value-- ", sd, " -- Local Percent -- ",
            totPCT, "%", "\n")
    }
    else {
        cat(" --------Final iteration----", iter, "\n", file = ifelse(is.null(file),
            "", file), append = TRUE)
        cat(" ----------- test =         ", test0, "\n", file = ifelse(is.null(file),
            "", file), append = TRUE)
        cat("\n", " --Singular Value-- ", sd, " -- Local Percent -- ",
            totPCT, "%", "\n", file = ifelse(is.null(file), "",
                file), append = TRUE)
    }
    sval[[i]]$iter <- iter
    sval[[i]]$test <- test
    sval[[i]]$d <- as.vector(sd)
    sval[[i]]$pct <- as.vector(totPCT)
    sval[[i]]$ssX <- as.vector(ssX)
    sval[[i]]$vsnam <- PTnam
    if (metrics) {
        for (d in 1:length(sval)) {
            if (length(met[[d]]) > 1) {
                if (length(met[[d]]) == dim(X)[d]^2) {
                  sval[[d]]$v <- as.vector(sval[[d]]$v %*% Powmat(met[[d]],
                    -1/2))
                }
                else {
                  sval[[d]]$v <- 1/sqrt(met[[d]]) * sval[[d]]$v
                }
            }
            else sval[[d]]$v <- sval[[d]]$v * 1/sqrt(met[[d]])
        }
    }
    class(sval) <- c("solutions.PTAk")
    return(sval)
}
"SVDgen" <-
function (Y, D2 = 1, D1 = 1, smoothing = FALSE, nomb = min(dim(Y)),
    smoo = list(function(u) ksmooth(1:length(u), u, kernel = "normal",
        bandwidth = 3, x.points = (1:length(u)))$y))
{
    nomb <- min(nomb, dim(Y))
    datanam <- substitute(Y)
    dinam <- dimnames(Y)
    if (length(D1) == dim(Y)[1]^2) {
        Y <- Powmat(D1, 1/2) %*% Y
    }
    else if (length(D1) == dim(Y)[1]) {
        Y <- sqrt(D1) * Y
    }
    else if (length(D1) == 1) {
        Y <- sqrt(D1) * Y
    }
    else stop(paste(" ----- Wrong DIMENSION for Metric of the First Entry ------ !!!!@@@@"))
    if (length(D2) == dim(Y)[2]^2) {
        Y <- Y %*% Powmat(D2, 1/2)
    }
    else if (length(D2) == dim(Y)[2]) {
        Y <- t(sqrt(D2) * t(Y))
    }
    else if (length(D2) == 1) {
        Y <- sqrt(D2) * Y
    }
    else stop(paste(" ----- Wrong DIMENSION for Metric of the First Entry ------ !!!!@@@@"))
    dimnames(Y) <- dinam
    if (smoothing == TRUE)
        result <- svdsmooth(Y, nomb = nomb, smooth = smoo)
    else result <- svd(Y)
    if (length(D1) == dim(Y)[1]^2) {
        result$u <- Powmat(D1, -1/2) %*% result$u
    }
    else if (length(D1) == dim(Y)[1]) {
        result$u <- 1/sqrt(D1) * result$u
    }
    else if (length(D1) == 1) {
        result$u <- (1/sqrt(D1)) * result$u
    }
    if (length(D2) == dim(Y)[2]^2) {
        result$v <- Powmat(D2, -1/2) %*% result$v
    }
    else if (length(D2) == dim(Y)[2]) {
        result$v <- 1/sqrt(D2) * result$v
    }
    else if (length(D2) == 1) {
        result$v <- 1/sqrt(D2) * result$v
    }
    if (all(result$u[, 1] < 0)) {
        result$u <- result$u * (-1)
        result$v <- result$v * (-1)
    }
    solutions <- list(NULL, NULL)
    solutions[[1]]$v <- t(result$u)[1:nomb, ]
    solutions[[1]]$modesnam <- "lignes"
    solutions[[1]]$n <- dimnames(Y)[[1]]
    solutions[[2]]$v <- t(result$v)[1:nomb, ]
    solutions[[2]]$modesnam <- "colonnes"
    solutions[[2]]$n <- dimnames(Y)[[2]]
    solutions[[2]]$d <- result$d[1:nomb]
    if (smoothing | nomb < min(dim(Y)))
        ssX <- sum(as.vector(Y)^2)
    else ssX <- sum(result$d^2)
    solutions[[2]]$pct <- (100 * result$d^2/ssX)[1:nomb]
    solutions[[2]]$ssX <- rep(ssX, nomb)
    solutions[[2]]$vsnam <- paste("vs", 1:nomb, sep = "")
    solutions[[2]]$datanam <- datanam
    solutions[[2]]$method <- match.call()
    if (smoothing)
        solutions[[2]]$smoocheck <- result$smoocheck
    class(solutions) <- c("solutions.PTAk")
    return(solutions)
}
"svdsmooth" <-
function (X, nomb = min(dim(X)), smooth = list(function(u) ksmooth(1:length(u),
    u, kernel = "normal", bandwidth = 3, x.points = (1:length(u)))$y),
    vsmin = 1e-16)
{
    if (!is.list(smooth[[1]]))
        smooth[[1]] <- list(smooth[[1]])
    if (length(smooth) < 2)
        smooth <- rep(list(smooth[[1]]), 2)
    if (!is.list(smooth[[2]]))
        smooth[[2]] <- list(smooth[[2]])
    solu <- list(NULL, NULL)
    solu[[1]]$v <- array(0, c(nomb, dim(X)[1]))
    solu[[2]]$v <- array(0, c(nomb, dim(X)[2]))
    solu[[2]]$d <- rep(0, nomb)
    solu[[2]]$smoocheck <- array(NA, c(2, nomb))
    solu[[2]]$smoocheck[1:2, 1] <- c(is.function(smooth[[1]][[1]]),
        is.function(smooth[[2]][[1]]))
    fi <- PPMA(X, pena = list(smooth[[1]][[1]], smooth[[2]][[1]]))
    solu[[1]]$v[1, ] <- fi$u
    solu[[2]]$v[1, ] <- fi$v
    solu[[2]]$d[1] <- fi$d
    for (qi in 2:nomb) {
        X <- PROJOT(X, solu, numo = (qi - 1))
        if (length(smooth[[1]]) == qi - 1)
            smooth[[1]][[qi]] <- smooth[[1]][[qi - 1]]
        if (length(smooth[[2]]) == qi - 1)
            smooth[[2]][[qi]] <- smooth[[2]][[qi - 1]]
        tempi <- list(toplist(smooth[[1]][[qi]]), toplist(smooth[[2]][[qi]]))
        fi <- PPMA(X, pena = tempi)
        solu[[1]]$v[qi, ] <- fi$u
        solu[[2]]$v[qi, ] <- fi$v
        solu[[2]]$d[qi] <- fi$d
        solu[[2]]$smoocheck[1:2, qi] <- c(is.function(smooth[[1]][[qi]]),
            is.function(smooth[[2]][[qi]]))
        if (fi$d < vsmin)
            break
    }
    return(list(u = t(solu[[1]]$v), d = solu[[2]]$d, v = t(solu[[2]]$v),
        smoocheck = solu[[2]]$smoocheck))
}
"toplist" <-
function (li)
{
    while (is.list(li)) {
        li <- li[[1]]
    }
    return(li)
}
