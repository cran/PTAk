"demo.FCAk" <-
function (data) 
{
    cat(" If you've got a nice little dataset say with 4 or 5 entries", 
        "\n", " send it to me  and I'll put the demo in the next release ", 
        "\n")
    cat(" didier@fmrib.ox.ac.uk ", "\n")
}
"demo.PTA3" <-
function (bootn = 1, show = 5, toshow = c(1, 3, 9, 12), yli = list(c(-1, 
    +1)), openX11s = FALSE) 
{
    data(iris3)
    cat(" A little fun using iris3 and matching randomly 15 for each iris sample!", 
        "\n")
    cat("   then performing a PTA-3modes.  If many draws are done, plots")
    cat("   show the stability of the first and third Principal Tensors.", 
        "\n")
    cat("iris3 is centered and reduced beforehand for each original variables.", 
        "\n")
    iris3 <- (iris3 - rep(1, 50) %o% apply(iris3, 2, mean) %o% 
        rep(1, 3))
    iris3 <- iris3 * rep(1, 50) %o% ((apply(matrix(aperm(iris3, 
        c(1, 3, 2)), ncol = 4), 2, var))^(-1/2)) %o% rep(1, 3)
    vs <- rep(list(NULL), 12)
    dimnames(iris3) <- list(NULL, c("SL", "SW", "PL", "PW"), 
        c("Set", "Ver", "Vir"))
    if (length(yli) < length(toshow)) 
        yli <- rep(list(yli[[1]]), length(toshow))
    for (combien in 1:bootn) {
        iridem <- iris3[1:15, , ]
        iridem[, , 1] <- iris3[sample(1:50)[1:15], , 1]
        iridem[, , 2] <- iris3[sample(1:50)[1:15], , 2]
        iridem[, , 3] <- iris3[sample(1:50)[1:15], , 3]
        result1 <- PTA3(iridem, verbose = !(combien == 1))
        if (combien%%show == 1) {
            RESUM(result1, summary = TRUE)
            if (openX11s) 
                X11(width = 4, height = 4, pointsize = 12)
            plot(result1, scree = TRUE, lty = 3, type = "b", 
                main = paste("Bootstrap ", combien))
            if (openX11s) 
                X11(width = 15, height = 5, pointsize = 12)
            par(mfrow = c(1, 3))
            plot(result1, labels = TRUE, mod = c(1, 2, 3), nb1 = toshow[2], 
                nb2 = NULL, ylimit = yli[[match(toshow[2], toshow)]])
            plot(result1, labels = TRUE, main = paste("Bootstrap ", 
                combien), mod = c(1, 2, 3), nb1 = toshow[3], 
                nb2 = NULL, ylimit = yli[[match(toshow[3], toshow)]])
            plot(result1, labels = TRUE, mod = c(1, 2, 3), nb1 = toshow[1], 
                nb2 = toshow[4], ylimi = yli[[match(toshow[1], 
                  toshow)]])
        }
        for (s in toshow) {
            for (u in 1:3) {
                if (combien == 1) 
                  vs[[s]][[u]] <- list(d = NULL, v = NULL, n = result1[[u]]$n)
                vs[[s]][[u]]$d <- c(vs[[s]][[u]]$d, result1[[u]]$d[s])
                vs[[s]][[u]]$ssX <- c(vs[[s]][[u]]$ssX, result1[[u]]$ssX[s])
                vs[[s]][[u]]$v <- rbind(vs[[s]][[u]]$v, result1[[u]]$v[s, 
                  ])
                vs[[s]][[u]]$pct <- c(vs[[s]][[u]]$pct, result1[[u]]$pct[s])
            }
        }
    }
    if (combien > 2) {
        if (openX11s) 
            X11(width = 8, height = 8, pointsize = 12)
        par(mfrow = c(2, 2))
        for (s in toshow[1:4]) {
            for (i in 1:combien) {
                plot.solutions.PTAk(vs[[s]], cex = 0.7, , ylimit = yli[[match(s, 
                  toshow)]], ppch = rep(i%%25, 3), labels = TRUE, 
                  mod = c(1, 2, 3), nb1 = i, nb2 = NULL, axes = FALSE, 
                  main = paste(combien, " PT", s), xylab = FALSE)
                par(new = TRUE)
            }
            par(new = FALSE)
        }
        if (openX11s) 
            X11(width = 5, height = 5, pointsize = 12)
        par(new = FALSE, mfrow = c(1, 1))
        varia <- as.ts(cbind(vs[[toshow[1]]][[3]]$d, vs[[toshow[2]]][[3]]$d, 
            vs[[toshow[3]]][[3]]$d, vs[[toshow[4]]][[3]]$d))
        varia <- varia - rep(1, combien) %o% apply(varia, 2, 
            mean)
        varia <- varia * rep(1, combien) %o% (apply(varia, 2, 
            var)^(-1/2))
        plot.ts(varia, main = paste("Std Error from mean Vsing "))
    }
    cat("***what was done is: ******", "\n")
    
    print(demo.PTA3)
    cat("**************", "\n")
    invisible(list(vs13912 = vs, lastresult = result1))
}
"demo.PTAk" <-
function (data) 
{
    cat(" If you've got a nice little dataset say with 4 or 5 entries", 
        "\n", " send it to me  and I'll put the demo in the next release ", 
        "\n")
    cat(" didier@fmrib.ox.ac.uk ", "\n")
}
"demo.SVD.CauRui" <-
function (ker = 0.84, withingroup = TRUE, openX11s = FALSE) 
{
    data(iris)
    iris2 <- as.matrix(iris[, 1:4])
    dimnames(iris2)[[1]] <- iris[, 5]
    D2 <- CauRuimet(iris2, ker = ker, withingroup = withingroup)
    D2 <- Powmat(D2, (-1))
    iris2 <- sweep(iris2, 2, apply(iris2, 2, mean))
    res <- SVDgen(iris2, D2 = D2, D1 = 1)
    res2 <- SVDgen(iris2)
    if (openX11s) 
        X11(width = 6, height = 4)
    par(mfrow = c(1, 2))
    if (withingroup) 
        plot(res, nb1 = 1, nb2 = 2, cex = 0.5, main = expression(paste("Unknown ", 
            Wg^{
                -1
            })), type = "n")
    else plot(res, nb1 = 1, nb2 = 2, cex = 0.5, main = expression(paste(Wo^{
        -1
    })), type = "n")
    mtext(paste("ker= ", ker), side = 3, cex = 1, line = 0)
    plot(res2, nb1 = 1, nb2 = 2, cex = 0.5, main = "Id", type = "n")
    summary(res, testvar = 0)
    cat("\n", "SVD with metrics in variables space", "\n")
    browser()
    print(demo.SVD.CauRui)
    invisible(res)
}
"demo.svdsmoo1" <-
function (snr = 3, openX11s = FALSE) 
{
    library(modreg)
    data(timage12)
    resx <- 91
    resy <- 109
    if (openX11s) 
        X11(width = 4, height = 4)
    filled.contour(1:resx, 1:resy, timage12, main = "ori")
    bsr <- range(as.vector(timage12))
    bs <- bsr[2] - bsr[1]
    cat("----", "\n")
    cat(" range of the image = ", bs, "\n", "noise added with sd =  range /snr = ", 
        bs/snr, "\n")
    bruitim <- array(rnorm(resx * resy, mean = 0, sd = bs/snr), 
        c(resx, resy))
    timage12 <- timage12 + bruitim
    d.svd <- SVDgen(timage12, nomb = 25)
    debtime <- proc.time()
    smoofun <- rep(list(Susan1D), 7)
    smoofun[[8]] <- NA
    d.svdo <- SVDgen(timage12, nomb = 25, smoothing = TRUE, smoo = list(smoofun))
    duree <- (proc.time() - debtime)[3]
    cat("duree ", duree, "s", "\n")
    summary(d.svd, testvar = 0.5)
    summary(d.svdo, testvar = 0.5)
    if (openX11s) 
        X11(width = 6, height = 4)
    par(mfrow = c(1, 2))
    plot(d.svd, scree = TRUE, RiskJack = 0, nbvs = 1:12, main = "svd")
    plot(d.svdo, scree = TRUE, nbvs = 1:12, RiskJack = 0, main = "svdo")
    d.svd.re <- REBUILD(d.svd, nTens = 1:7, testvar = 0)
    d.svdo.re <- REBUILD(d.svdo, nTens = 1:7, testvar = 0)
    if (openX11s) 
        X11(width = 9, height = 4)
    par(mfrow = c(1, 3))
    zlimi <- bsr
    image(1:resx, 1:resy, timage12, col = cm.colors(20), zlim = zlimi)
    image(1:resx, 1:resy, d.svd.re, col = cm.colors(20), zlim = zlimi)
    image(1:resx, 1:resy, d.svdo.re, col = cm.colors(20), zlim = zlimi)
    par(mfrow=c(1,1))
    if (openX11s) 
        X11(width = 4, height = 4)
    filled.contour(1:resx, 1:resy, timage12, main = paste("ori+bruit (snr=", 
        snr, ")"))
    if (openX11s) 
        X11(width = 4, height = 4)
    filled.contour(1:resx, 1:resy, d.svd.re, main = "svd 7")
    if (openX11s) 
        X11(width = 4, height = 4)
    filled.contour(1:resx, 1:resy, d.svdo.re, main = "svdsmooth 7")
}
