"CauRuimet" <-
function (Z, ker = 1, m0 = 1, withingroup = TRUE, loc = substitute(apply(Z,
    2, mean, trim = 0.1)), matrixmethod = TRUE)
{
    debtime <- proc.time()
    if (m0 == "tridiag") {
        m0 <- array(as.integer(0), c(dim(Z)[1], dim(Z)[1]))
        m0[1:2, 1] <- c(1, 1)
        m0[(dim(Z)[1] - 1):dim(Z)[1], dim(Z)[1]] <- c(1, 1)
        for (j in 2:(dim(Z)[1] - 1)) {
            m0[j - 1, j] <- 1
            m0[j, j] <- 1
            m0[j + 1, j] <- 1
        }
    }
    mz <- eval(loc)
    Sz <- sweep(Z, 2, mz)
    Sz <- t(Sz) %*% Sz/(dim(Z)[1] - 1)
    norm2S <- function(u, S = Powmat(Sz, (-1))) {
        return(t(u) %*% S %*% u)
    }
    if (is.numeric(ker)) {
        g <- ker
        ker <- function(t) {
            return(exp(-(g * t)))
        }
    }
    if (withingroup) {
        if (matrixmethod) {
            distZiZj <- norm2S(t(Z))
            diadis <- diag(distZiZj)/2
            distZiZj <- 2 * sweep(sweep(-distZiZj, 2, -diadis),
                1, -diadis)
            M <- m0 * ker(distZiZj)
            sumM <- (sum(as.vector(M)) - dim(Z)[1])/2
            M <- diag(apply(M, 2, sum)) - M
            W <- norm2S(Z, M)/sumM
        }
        else {
            W <- matrix(0, nrow = dim(Z)[2], ncol = dim(Z)[2])
            totad <- 0
            for (i in 1:(dim(Z)[1] - 1)) for (j in (i + 1):dim(Z)[1]) {
                ad <- as.double(ker(norm2S(Z[i, ] - Z[j, ])))
                if (is.matrix(m0))
                  ad <- ad * m0[i, j]
                W <- W + ad * ((Z[i, ] - Z[j, ]) %o% (Z[i, ] -
                  Z[j, ]))
                totad <- totad + ad
            }
            totad <- totad
            W <- W/totad
        }
    }
    else {
        W <- matrix(0, nrow = dim(Z)[2], ncol = dim(Z)[2])
        totad <- 0
        for (i in 1:(dim(Z)[1])) {
            ad <- as.double(ker(norm2S(Z[i, ] - mz)))
            W <- W + ad * ((Z[i, ] - mz) %o% (Z[i, ] - mz))
            if (is.matrix(m0))
                ad <- ad * m0[i, j]
            totad <- totad + ad
        }
        totad <- totad * dim(Z)[1]^2
        W <- W/totad
    }
    cat("-----Execution Time-----", (proc.time() - debtime)[3],
        "\n")
    return(W)
}
"Detren" <-
function (dat, Mm = c(1, 3), rsd = TRUE, tren = function(x) smooth.spline(as.vector(x),
    df = 5)$y)
{
    tre <- apply(dat, Mm, FUN = tren)
    dimi <- c(dim(dat)[-Mm], dim(dat)[Mm])
    tre <- aperm(array(tre, dimi), match(dimi, dim(dat)))
    if (rsd)
        return(dat - tre)
    else return(tre)
}
"FCAmet" <-
function (X, chi2 = FALSE, E = NULL)
{
    if (!is.array(X)) {
        stop(paste("--- X must be an array  ! ---"))
    }
    datanam <- substitute(X)
    ord <- length(dim(X))
    N <- sum(X)
    metafc <- rep(list(NULL), ord)
    Indep <- metafc[[1]] <- apply(X, 1, sum)/N
    for (t in 2:ord) {
        metafc[[t]] <- apply(X, t, sum)/N
        Indep <- Indep %o% metafc[[t]]
    }
    if (chi2) {
        Indep <- array(Indep, dim(X))
        Chi2 <- N * sum((X/N - Indep)^2/Indep)
        cat("\n", " --")
        cat("\n", "++++ Data is            ", deparse(datanam),
            "        +++++++")
        cat("\n", "-------------- Multiple contingency Table of dimensions ",
            dim(X), "  ----", "\n")
        cat("\n", "-------------- Chi2 = ", Chi2, " with ddl = ",
            prod(dim(X) - 1))
        cat("\n", " ------------- p(>Chi2)= ", pchisq(Chi2, df = prod(dim(X) -
            1), lower.tail = FALSE), "\n")
        cat("\n", " --", "\n")
    }
    if (!is.null(E))
        invisible(list(data = (X/N - E)/Indep, met = metafc,
            count = N))
    else invisible(list(data = (X/N)/Indep, met = metafc, count = N))
}
"Ginv" <-
function (A)
{
    Powmat(A, -1)
}
"IterMV" <-
function (n = 10, dat = J12.dat, Mm = c(1, 3), Vm = c(2, 3),
    fFUN = mean, usetren = FALSE, tren = function(x) smooth.spline(as.vector(x),
        df = 5)$y, rsd = TRUE)
{
    sdev <- function(x) {
        sd(as.vector(x))
    }
    for (i in 1:n) {
        if (usetren) {
            dat <- Detren(dat = dat, Mm = Mm, tren = tren, rsd = rsd)
        }
        else {
            mean.dat <- apply(dat, Mm, FUN = fFUN)
            dat <- sweep(dat, Mm, mean.dat)
        }
        sd.dat <- apply(dat, Vm, sdev)
        if (sd.dat == 1)
            warning("zero variances were replaced by 1")
        sd.dat <- ifelse(sd.dat == 0, 1, sd.dat)
        dat <- sweep(dat, Vm, sd.dat, FUN = "/")
    }
    cat("\n", "---Max of the means: ", max(apply(dat, Mm, mean),
        nam = TRUE), "\n")
    return(dat)
}
"Multcent" <-
function (dat = J12.dat, bi = c(1, 2), by = 3, centre = mean,
    centrebyBA = c(TRUE, FALSE), scalebyBA = c(TRUE, FALSE))
{
    if (centrebyBA[1]) {
        me <- apply(dat, by, FUN = centre)
        dat <- sweep(dat, by, me)
    }
    sdev <- function(x) {
        sd(as.numeric(x))
    }
    if (scalebyBA[1]) {
        sca <- apply(dat, by, sdev)
        if (sca == 1)
            warning("zero variances were replaced by 1")
        sca <- ifelse(sca == 0, 1, sca)
        dat <- sweep(dat, by, sca, FUN = "/")
    }
    if (!is.null(bi)) {
        for (g in 1:length(bi)) {
            me <- apply(dat, c(bi[g], by), FUN = centre)
            dat <- sweep(dat, c(bi[g], by), me)
        }
    }
    if (centrebyBA[2]) {
        me <- apply(dat, by, FUN = centre)
        dat <- sweep(dat, by, me)
    }
    if (scalebyBA[2]) {
        sca <- apply(dat, by, sdev)
        if (sca == 1)
            warning("zero variances were replaced by 1")
        sca <- ifelse(sca == 0, 1, sca)
        dat <- sweep(dat, by, sca, FUN = "/")
    }
    return(dat)
}
"Powmat" <-
function (A, pw, eltw = FALSE)
{
    A <- as.matrix(A)
    if (eltw) {
        dimA <- dim(A)
        A <- as.vector(A)
        RR <- A^pw
        RR[abs(RR) == Inf] <- A[abs(RR) == Inf]
        if (dimA[2] > 1)
            RR <- matrix(RR, ncol = dimA[2])
    }
    else {
        valsin <- svd(A)
        diago <- valsin$d[valsin$d > 1e-06]
        diago <- diago^pw
        if (length(diago) == 0) {
            RR <- matrix(0, ncol(A), nrow(A))
            return(RR)
        }
        if (length(diago) == 1)
            RR <- t(as.matrix(valsin$v[, 1]) %*% t(as.matrix(valsin$u[,
                1]))) * diago
        else RR <- valsin$u[, 1:length(diago)] %*% diag(diago) %*%
            t(valsin$v[, 1:length(diago)])
        RR <- as.matrix(RR)
        if (pw < 0 & (!min(dim(RR)) == 1))
            RR <- t(RR)
        if (length(RR) == 1)
            RR <- as.numeric(RR)
        else if (dim(RR)[1] == 1)
            RR <- as.vector(RR)
    }
    return(RR)
}
"RaoProd" <-
function (A, B)
{
    A <- as.matrix(A)
    B <- as.matrix(B)
    if (min(dim(A)) == 1 & min(dim(B)) == 1)
        return(as.vector(A) %x% as.vector(B))
    else {
        if (length(A) == 1 || length(B) == 1) {
            ifelse(length(B) == 1, return(A * as.vector(B)),
                return(as.vector(A) * B))
        }
        else {
            if (!dim(A)[2] == dim(B)[2])
                stop("Wrong number of columns")
            f <- dim(A)[2]
            re <- array(0, c(dim(A)[1] * dim(B)[1], f))
            for (w in 1:f) {
                re[, w] <- A[, w] %x% B[, w]
            }
            return(re)
        }
    }
}
"RiskJack.plot" <-
function (solution, nbvs = 1:20, mod = NULL, max = NULL, rescaled = TRUE,
    ...)
{
    qchoix <- nbvs
    ord <- length(solution)
    if (is.null(max))
        max <- length(solution[[ord]]$d)
    if (is.null(mod))
        mod <- 1:length(solution)
    val <- solution[[ord]]$d[!substr(solution[[ord]]$vsnam, 1,
        1) == "*"]
    if (ord > 2)
        iden <- order(val)[length(val):1]
    else iden <- 1:max
    covalid <- function() {
        mindiff <- min((val[iden][-length(solution[[ord]]$d)]^2 -
            val[iden][-1]^2))
        for (mode in mod) {
            if (length(solution[[mode]]$v[1, ]) < solution[[ord]]$ssX[1]^2/mindiff/length(solution[[mode]]$v[1,
                ]^2)) {
                cat(" WARNING ..mode ", mode, " ..n= ", length(solution[[mode]]$v[1,
                  ]), " validity condition  >", solution[[ord]]$ssX[1]^2/mindiff/length(solution[[mode]]$v[1,
                  ]^2), "\n")
            }
        }
    }
    RJack <- matrix(rep(0, max(mod) * length(qchoix)), c(max(mod),
        length(qchoix)))
    for (m in mod) {
        for (q in qchoix) {
            tl <- 0
            if (q > (max - 1))
                q <- max - 1
            for (k in 1:q) {
                for (j in (q + 1):max(iden)) {
                  l1 <- solution[[ord]]$d[iden[j]]^2
                  l2 <- solution[[ord]]$d[iden[k]]^2
                  tjk <- mean(solution[[m]]$v[iden[j], ]^2 *
                    solution[[m]]$v[iden[k], ]^2) * l1 * l2
                  diff <- (l1 - l2)^2
                  tl <- tl + tjk/diff
                }
            }
            RJack[m, match(q, qchoix)] <- tl * 1/(length(solution[[m]]$v[j,
                ]) - 1)
            if (q == (max - 1))
                q <- max(qchoix)
        }
    }
    for (u in mod) {
        if (rescaled)
            RJack[u, ] <- (RJack[u, ] - min(RJack[u, ]))/(max(RJack[u,
                ]) - min(RJack[u, ]))
        plot(qchoix, RJack[u, ], xlab = "Nb of dimensions", ylab = "Risk's approx",
            lty = u, col = u, type = "b", ...)
        par(new = TRUE)
    }
    legend(max(qchoix) - 1.5, max(RJack)/2, paste("Risk-mode",
        mod), col = mod, lty = mod, bty = "n", cex = 0.7)
    invisible(par(new = FALSE))
}
"Susan1D" <-
function (y, x = NULL, sigmak = NULL, sigmat = NULL, ker = list(function(u) return(exp(-0.5 *
    u^2))))
{
    if (is.null(x))
        x <- 1:length(y)
    else {
        if (!length(x) == length(y))
            stop("Wrong length for x")
        y <- y[order(x)]
        x <- sort(x)
    }
    if (is.null(sigmat))
        sigmat <- 8 * (length(y)^(-1/5))
    if (is.null(sigmak))
        sigmak <- 1/2 * (range(y)[2] - range(y)[1])
    if (length(ker) < 2)
        ker <- list(t = ker[[1]], k = ker[[1]])
    knei <- max(1, round(2 * sigmat))
    resul <- y
    for (t in 1:length(y)) {
        xt <- 0
        wjt <- 0
        for (j in max(1, t - knei):min(length(y), t + knei)) {
            wj <- ker$t((x[j] - x[t])/sigmat) * ker$k((y[j] -
                y[t])/sigmak)
            xt <- xt + wj * y[j]
            wjt <- wjt + wj
        }
        resul[t] <- xt/wjt
    }
    return(resul)
}
"plot.solutions.PTAk" <-
function (solution, labels = TRUE, mod = 1, nb1 = 1, nb2 = NULL,
    coefi = list(NULL, NULL), xylab = TRUE, ppch = (1:length(solution)),
    lengthlabels = 2, ylimit = NULL, scree = FALSE, ordered = TRUE,
    nbvs = 40, RiskJack = NULL, method = "", ...)
{
    if (is.null(coefi[[1]]))
        coefi[[1]] <- rep(1, length(solution))
    if (is.null(coefi[[2]]))
        coefi[[2]] <- rep(1, length(solution))
    if (is.null(lengthlabels))
        lengthlabels <- rep(10, length(solution))
    if (length(lengthlabels) == 1)
        lengthlabels <- rep(lengthlabels, length(solution))
    ord <- length(solution)
    if (as.character(solution[[ord]]$method)[1] == "FCA" | method ==
        "FCA") {
        divv <- solution[[ord]]$ssX[1] - 1
        perclab <- "% FCA"
        if (length(nbvs) == 1)
            nbvs <- 2:nbvs
        else if (1 %in% nbvs)
            nbvs <- nbvs[-match(1, nbvs)]
    }
    else {
        divv <- solution[[ord]]$ssX[1]
        perclab <- "% global"
    }
    di <- NULL
    for (r in 1:length(solution)) di <- c(di, length(solution[[r]]$v[1,
        ]))
    if (!scree) {
        ylim <- ylimit
        xlab <- ""
        ylab <- ""
        if (is.null(nb2)) {
            xlim <- c(1, max(di[mod]) + 1)
            if (xylab)
                ylab <- paste(solution[[ord]]$vsnam[nb1], " local",
                  round(solution[[ord]]$pct[nb1], 2), "% ", round((100 *
                    (solution[[ord]]$d[nb1])^2)/divv, 2), perclab)
        }
        else {
            if (xylab)
                xlab <- paste(solution[[ord]]$vsnam[nb1], " local",
                  round(solution[[ord]]$pct[nb1], 2), "% ", round((100 *
                    (solution[[ord]]$d[nb1])^2)/divv, 2), perclab)
            if (xylab)
                ylab <- paste(solution[[ord]]$vsnam[nb2], " local",
                  round(solution[[ord]]$pct[nb2], 2), "% ", round((100 *
                    (solution[[ord]]$d[nb2])^2)/divv, 2), perclab)
        }
        for (u in mod) {
            if (!is.null(nb2)) {
                xy <- t(solution[[u]]$v[c(nb1, nb2), ]) %*% diag(c(coefi[[1]][u],
                  coefi[[2]][u]))
                xyn <- t(solution[[u]]$v[c(nb1, nb2), ]) %*%
                  diag(c(coefi[[1]][u], coefi[[2]][u]))
                ylim <- c(min(ylim, xyn), max(ylim, xyn))
            }
            else {
                xy <- solution[[u]]$v[nb1, ] * coefi[[1]][u]
                xyn <- solution[[u]]$v[nb1, ] * coefi[[1]][u]
                ylim <- c(min(ylim, xyn), max(ylim, xyn))
            }
        }
        for (u in mod) {
            if (!is.null(nb2)) {
                xy <- t(solution[[u]]$v[c(nb1, nb2), ]) %*% diag(c(coefi[[1]][u],
                  coefi[[2]][u]))
                xyn <- t(solution[[u]]$v[c(nb1, nb2), ]) %*%
                  diag(c(coefi[[1]][u], coefi[[2]][u]))
                ylim <- c(min(ylim, xyn), max(ylim, xyn))
                xlim <- ylim
                xaxt <- "s"
            }
            else {
                xy <- solution[[u]]$v[nb1, ] * coefi[[1]][u]
                xyn <- solution[[u]]$v[nb1, ] * coefi[[1]][u]
                ylim <- c(min(ylim, xyn), max(ylim, xyn))
                if (!"xaxt" %in% names(list(...)))
                  xaxt <- "n"
            }
            if (labels) {
                if ("xlab" %in% names(list(...))) {
                  if ("ylab" %in% names(list(...)))
                    plot(xy, xlim = xlim, ylim = ylim, pch = ppch[u],
                      xaxt = xaxt, ...)
                  else plot(xy, xlim = xlim, ylim = ylim, ylab = ylab,
                    pch = ppch[u], xaxt = xaxt, ...)
                  if (is.null(nb2))
                    axis(1, 1:length(xy))
                }
                else {
                  if ("ylab" %in% names(list(...)))
                    plot(xy, xlim = xlim, ylim = ylim, pch = ppch[u],
                      xaxt = xaxt, xlab = xlab, ...)
                  else plot(xy, xlim = xlim, ylim = ylim, ylab = ylab,
                    pch = ppch[u], xaxt = xaxt, xlab = xlab,
                    ...)
                }
                if (!is.null(solution[[u]]$n)) {
                  if (is.factor(solution[[u]]$n)) {
                    if ("cex" %in% names(list(...)))
                      cex <- list(...)$cex
                    else cex <- par("cex")
                    text(xy, labels = substr(levels(solution[[u]]$n)[as.numeric(solution[[u]]$n)],
                      1, lengthlabels[u]), col = as.numeric(solution[[u]]$n),
                      pos = 4, cex = cex)
                    if (is.null(nb2)) {
                      par(new = TRUE)
                      plot(xy ~ solution[[u]]$n, xlab = "", ylab = "",
                        ylim = ylim, cex = cex)
                      par(new = FALSE)
                    }
                  }
                  else text(xy, labels = substr(solution[[u]]$n,
                    1, lengthlabels[u]), col = u, pos = 4)
                }
            }
            else if ("xlab" %in% names(list(...))) {
                if ("ylab" %in% names(list(...)))
                  plot(xy, xlim = xlim, ylim = ylim, pch = ppch[u],
                    col = u, xaxt = xaxt, ...)
                else plot(xy, xlim = xlim, ylim = ylim, ylab = ylab,
                  pch = ppch[u], col = u, xaxt = xaxt, ...)
                if (is.null(nb2))
                  axis(1, 1:length(xy))
            }
            else {
                if ("ylab" %in% names(list(...)))
                  plot(xy, xlim = xlim, ylim = ylim, pch = ppch[u],
                    col = u, xaxt = xaxt, xlab = xlab, ...)
                else plot(xy, xlim = xlim, ylim = ylim, ylab = ylab,
                  pch = ppch[u], col = u, xaxt = xaxt, xlab = xlab,
                  ...)
            }
            abline(h = 0, col = "green", lty = 2)
            abline(v = 0, col = "green", lty = 2)
            par(new = TRUE)
        }
        invisible(par(new = FALSE))
    }
    else {
        if (!is.null(ordered)) {
            if (ordered == TRUE) {
                ld <- length(solution[[ord]]$d[!substr(solution[[ord]]$vsnam,
                  1, 1) == "*"])
                if (length(nbvs) == 1) {
                  nbvs <- min(max(5, nbvs), ld)
                  nbvs <- 1:nbvs
                }
                scre <- 100 * ((solution[[ord]]$d[!substr(solution[[ord]]$vsnam,
                  1, 1) == "*"])^2)/divv
                scre <- (sort(scre[nbvs]))
                scre <- scre[length(scre):1]
                nbvs <- nbvs[1:length(scre)]
                plot(nbvs, scre, xlab = "Ordered ", ylab = "Squared Singular Values (%)",
                  xaxt = "n", ...)
                axis(1, at = nbvs)
                par(new = TRUE)
                plot(nbvs, ylim = c(0, 100), cumsum(scre), axes = FALSE,
                  lwd = 2, lty = 1, type = "b", pch = "c", col = 3,
                  xlab = "", ylab = "")
                axis(4, at = atpc <- seq(0, 100, 10), labels = formatC(atpc,
                  format = "fg"), col.axis = 3)
                par(new = TRUE)
                if (!is.null(RiskJack))
                  RiskJack.plot(solution, nbvs = nbvs, mod = mod,
                    max = min(RiskJack + length(nbvs), ld), rescaled = TRUE,
                    axes = FALSE, ann = FALSE, pch = "r")
                par(new = FALSE)
            }
            if (ordered == FALSE) {
                ld <- length(solution[[ord]]$d[!substr(solution[[ord]]$vsnam,
                  1, 1) == "*"])
                if (length(nbvs) == 1) {
                  nbvs <- min(max(5, nbvs), ld)
                  nbvs <- 1:nbvs
                }
                scre <- ((solution[[ord]]$d)^2)[nbvs]
                plot(nbvs, scre, xlab = "Unordered with redundancy",
                  ylab = "Squared Singular Values", ...)
            }
        }
    }
    invisible(par(new = FALSE))
}
