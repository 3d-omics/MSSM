################
# Extended RLQ #
################

rlqESLTP <- function(dudiE, dudiS, dudiL, dudiT, dudiP, ...){
	tabE <- dudiE$li/sqrt(dudiE$eig[1])
	tabS <- dudiS$li/sqrt(dudiS$eig[1])
	tabP <- dudiP$li/sqrt(dudiP$eig[1])
	tabT <- dudiT$li/sqrt(dudiT$eig[1])
	tabES <- cbind.data.frame(tabE, tabS)
	names(tabES) <- c(paste("E", 1:ncol(tabE), sep = ""),
		paste("S", 1:ncol(tabS), sep = ""))
	tabTP <- cbind.data.frame(tabT, tabP)
	names(tabTP) <- c(paste("T", 1:ncol(tabT), sep = ""),
		paste("P", 1:ncol(tabP), sep = ""))
	pcaES <- dudi.pca(tabES, scale = F, row.w = dudiL$lw, scan = FALSE,
		nf = (length(dudiE$eig) + length(dudiS$eig)))
	pcaTP <- dudi.pca(tabTP, scale = F, row.w = dudiL$cw, scan = FALSE,
		nf = (length(dudiT$eig) + length(dudiP$eig)))

	X <- rlq(pcaES, dudiL, pcaTP, ...)

	U <- as.matrix(X$l1) * unlist(X$lw)
	U <- data.frame(as.matrix(pcaES$tab[, 1:ncol(tabE)]) %*% U[1:ncol(tabE), 1:X$nf])
	row.names(U) <- row.names(pcaES$tab)
	names(U) <- names(X$lR)
	X$lR_givenE <- U
	

	U <- as.matrix(X$l1) * unlist(X$lw)
	U <- data.frame(as.matrix(pcaES$tab[, -(1:ncol(tabE))]) %*% U[-(1:ncol(tabE)), 1:X$nf])
	row.names(U) <- row.names(pcaES$tab)
	names(U) <- names(X$lR)
	X$lR_givenS <- U


	U <- as.matrix(X$c1) * unlist(X$cw)
	U <- data.frame(as.matrix(pcaTP$tab[, 1:ncol(tabT)]) %*% U[1:ncol(tabT), 1:X$nf])
	row.names(U) <- row.names(pcaTP$tab)
	names(U) <- names(X$lQ)
	X$lQ_givenT <- U

	U <- as.matrix(X$c1) * unlist(X$cw)
	U <- data.frame(as.matrix(pcaTP$tab[, -(1:ncol(tabT))]) %*% U[-(1:ncol(tabT)), 1:X$nf])
	row.names(U) <- row.names(pcaTP$tab)
	names(U) <- names(X$lQ)
	X$lQ_givenP <- U

	X$row.w <- dudiL$lw

	X$col.w <- dudiL$cw
	
	class(X) <- c("rlqESLTP", "rlq", "dudi")

	return(X)
}


plot.rlqESLTP <- function(X, which = "E", phy = NULL, xy = NULL, traits = NULL,
	env = NULL, type = NULL, ax = 1){

	
	if(which == "S"){
	if(is.null(xy)) stop("xy required")
	par(mfrow = c(1, 3))
	s.value(xy, X$lR_givenE[, ax], zmax = max(X$lR[, ax]), addaxes = F, clegend = 2,
		sub = "environment-based", csub=2)
	s.value(xy, X$lR_givenS[, ax], zmax = max(X$lR[, ax]), addaxes = F, clegend = 2,
		sub = "space-based", csub = 2)
	s.value(xy, X$lR[,ax], zmax = max(X$lR[, ax]), addaxes = F, clegend = 2,
		sub = "global", csub=2)
	}

	if(which == "P"){
	if(is.null(phy)) stop("phy required")
	par(mfrow = c(1, 3))
	dotchart.phylog(phy, X$lQ_givenT[names(phy$leaves), ax], cleav=0, cdot = 1,
		scaling = F, yjoi = 0, cex.axis = 1.5, sub="trait-based", csub=0)
	dotchart.phylog(phy, X$lQ_givenP[names(phy$leaves), ax], cleav = 0, cdot = 1,
		scaling = F, yjoi = 0, cex.axis = 1.5, sub = "phylogeny-based", csub = 0)
	dotchart.phylog(phy, X$lQ[names(phy$leaves), ax], cleav = 0, cdot = 1,
		scaling = F, yjoi = 0, cex.axis = 1.5, sub = "global", csub = 0)
	}

	if(which == "T" | which == "E"){
		
        mfrow = n2mfrow(length(type))
        par(mfrow = mfrow)

		if(which == "T"){
			ltab <- traits
			w <- X$col.w
            sco1 <- X$lQ
            sco2 <- X$mQ
		}
		else{
			ltab <- env
			w <- X$row.w
            sco1 <- X$lR
            sco2 <- X$mR
		}
		if (is.data.frame(ltab)) ltab <- list(ltab)
		for(i in 1:length(ltab)){
			if(type[i] == "Q"){
				thetab <- ltab[[i]]
                if(!any(is.na(thetab))){
                    thetabS <- scalewt(thetab, w)
                    corS <- (t(thetabS)%*%diag(w)%*%sco2[, ax])[, 1]
                }
                else{
                    funcorS <- function(j){
                        x <- thetab[, j]
                        xsna <- x[!is.na(x)]
                        sco2sna <- sco2[!is.na(x), ax]
                        wsna <- w[!is.na(x)]
                        thetabSsna <- scalewt(xsna, wsna)
                        corSsna <- t(thetabSsna)%*%diag(wsna)%*%sco2sna
                        return(corSsna)
                    }
                    corS <- sapply(1:ncol(thetab), funcorS)
                    names(corS) <- names(thetab)
                }
				dotchart(sort(corS), lab = rownames(corS)[order(corS)],
					main = "Pearson correlation")
				abline(v = 0)

			}
			if(type[i] == "O"){
				thetab <- ltab[[i]]
				thetab <- as.data.frame(apply(thetab, 2, rank))

                if(!any(is.na(thetab))){
                    thetabS <- scalewt(thetab, w)
                    corS <- t(thetabS)%*%diag(w)%*%scalewt(rank(sco2[, ax]), w)
                }
                else{
                    funcorS <- function(j){
                        x <- thetab[, j]
                        xsna <- x[!is.na(x)]
                        wsna <- w[!is.na(x)]
                        sco2sna <- scalewt(rank(sco2[!is.na(x), ax]), wsna)
                        thetabSsna <- scalewt(xsna, wsna)
                        corSsna <- t(thetabSsna)%*%diag(wsna)%*%sco2sna
                        return(corSsna)
                    }
                    corS <- sapply(1:ncol(thetab), funcorS)
                    names(corS) <- names(thetab)
                }
				dotchart(sort(corS), lab = rownames(corS)[order(corS)],
					main = "Spearman correlation")
				abline(v = 0)
				
			}
			if(type[i] == "N"){
                thetab <- ltab[[i]]
                funmod <- function(unx){

                    if(!any(is.na(unx))){
                    mod <- model.matrix(~-1+factor(unx))
                    colnames(mod) <- levels(factor(unx))
                    rownames(mod) <- rownames(thetab)
                    return(as.data.frame(mod))
                    }
                    else{
                        mod <- model.matrix(~-1+factor(unx))
                        correctedtab <- matrix(NA, nrow(thetab), ncol(mod))
                        correctedtab[as.numeric(rownames(mod)), ] <- mod
                        colnames(correctedtab) <- levels(factor(unx))
                        rownames(correctedtab) <- rownames(thetab)
                        return(as.data.frame(correctedtab))
                    }
                    }
                    res <- cbind.data.frame(apply(thetab, 2, funmod))
                    sco.distrina(sco1[, ax], res)
			}
			if(type[i] == "F" | type[i] == "B" | type[i] == "D"){
                thetab <- ltab[[i]]
				sco.distrina(sco1[, ax], thetab)
			}
            if(type[i] == "C"){
                thetab <- ltab[[i]]
                if(!any(is.na(thetab))){
                    alphat <- t(t(thetab * 2 * pi)/attributes(thetab)$max)
                    alphatcos <- scalewt(cos(alphat), w)
                    alphatsin <- scalewt(sin(alphat), w)
                    rxc <- t(alphatcos)%*%diag(w)%*%sco2[, ax]
                    rxs <- t(alphatsin)%*%diag(w)%*%sco2[, ax]
                    rcs <- diag(t(alphatsin)%*%diag(w)%*%alphatcos)
                    corC <- (sqrt((rxc^2 + rxs^2 - 2*rxc*rxs*rcs)/(1 -
                    rcs^2)))[, 1]
                }
                else{
                    funcorC <- function(j){
                        x <- thetab[, j]
                        xsna <- x[!is.na(x)]
                        sco2sna <- sco2[!is.na(x), ax]
                        wsna <- w[!is.na(x)]
                        alphat <- xsna * 2 * pi/attributes(thetab)$max[j]
                        alphatcos <- scalewt(cos(alphat), wsna)
                        alphatsin <- scalewt(sin(alphat), wsna)
                        rxc <- t(alphatcos)%*%diag(wsna)%*%sco2sna
                        rxs <- t(alphatsin)%*%diag(wsna)%*%sco2sna
                        rcs <- diag(t(alphatsin)%*%diag(wsna)%*%alphatcos)
                        corCsna <- sqrt((rxc^2 + rxs^2 - 2*rxc*rxs*rcs)/(1 - rcs^2))
                        return(corCsna)
                    }
                    corC <- sapply(1:ncol(thetab), funcorC)
                    names(corC) <- names(thetab)
                }
           		dotchart(sort(corC), lab = rownames(corC)[order(corC)],
					main = "Circular correlation")
				abline(v = 0)

            }
		}
	}

    par(mfrow=c(1, 1))

}

sco.distrina <- function (score, df, y.rank = TRUE, csize = 1, labels =
names(df),
    clabel = 1, xlim = NULL, grid = TRUE, cgrid = 0.75, include.origin = TRUE,
    origin = 0, sub = NULL, csub = 1)
{
    if (!is.vector(score))
        stop("vector expected for score")
    if (!is.numeric(score))
        stop("numeric expected for score")
    if (!is.data.frame(df))
        stop("data.frame expected for df")
    n <- length(score)
    if ((nrow(df) != n))
        stop("Non convenient match")
    n <- length(score)
    nvar <- ncol(df)
    opar <- par(mar = par("mar"))
    on.exit(par(opar))
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    ymin <- scoreutil.base(y = score, xlim = xlim, grid = grid,
        cgrid = cgrid, include.origin = include.origin, origin = origin,
        sub = sub, csub = csub)
    ymax <- par("usr")[4]
    ylabel <- strheight("A", cex = par("cex") * max(1, clabel)) *
        1.4
    xmin <- par("usr")[1]
    xmax <- par("usr")[2]
    xaxp <- par("xaxp")
    nline <- xaxp[3] + 1
    v0 <- seq(xaxp[1], xaxp[2], le = nline)
    if (grid) {
        segments(v0, rep(ymin, nline), v0, rep(ymax, nline),
            col = gray(0.5), lty = 1)
    }
    rect(xmin, ymin, xmax, ymax)
    sum.col <- apply(df, 2, sum, na.rm = TRUE)
    labels <- labels[sum.col > 0]
    df <- df[, sum.col > 0]
    nvar <- ncol(df)
    sum.col <- apply(df, 2, sum, na.rm = TRUE)
    df <- sweep(df, 2, sum.col, "/")
    y.distri <- (nvar:1)
    if (y.rank) {
        y.distri <- drop(apply(df, 2, function(x) sum(x[!is.na(x)] * score[!is.na(x)])))
        y.distri <- rank(y.distri, ties.method = "first")
    }
    ylabel <- strheight("A", cex = par("cex") * max(1, clabel)) *
        1.4
    y.distri <- (y.distri - min(y.distri))/(max(y.distri) - min(y.distri))
    y.distri <- ymin + ylabel + (ymax - ymin - 2 * ylabel) *
        y.distri
    res <- matrix(0, nvar, 2)
    for (i in 1:nvar) {
        w <- df[, i]
	wna <- w[!is.na(w)]
	scorena <- score[!is.na(w)]
        y0 <- y.distri[i]
        x.moy <- sum(wna * scorena)
        x.et <- sqrt(sum(wna * (scorena - x.moy)^2))
        res[i, 1] <- x.moy
        res[i, 2] <- x.et * x.et
        x1 <- x.moy - x.et * csize
        x2 <- x.moy + x.et * csize
        etiagauche <- TRUE
        if ((x1 - xmin) < (xmax - x2))
            etiagauche <- FALSE
        segments(x1, y0, x2, y0)
        if (clabel > 0) {
            cha <- labels[i]
            cex0 <- par("cex") * clabel
            xh <- strwidth(cha, cex = cex0)
            xh <- xh + strwidth("x", cex = cex0)
            if (etiagauche)
                x0 <- x1 - xh/2
            else x0 <- x2 + xh/2
            text(x0, y0, cha, cex = cex0)
        }
        points(x.moy, y0, pch = 20, cex = par("cex") * 2)
    }
    res <- as.data.frame(res)
    names(res) <- c("mean", "var")
    rownames(res) <- names(df)
    invisible(res)
}


###################################################################
# Tests for phylogenetic and traits clustering vs overdispersion  #
###################################################################

TPQE <- function(df, dis, nrep = 999, popw = NULL, alter = "two-sided"){

    if(is.null(popw))
        popw <- apply(df, 2, sum) / sum(df)

    fun1 <- function(df1, dis){
    stati <- divc(df1, dis)
    statimean <- sum(popw * stati[, 1])
    vtot <- apply(sweep(sweep(df1, 2, apply(df1, 2, sum), "/"), 2, popw, "*"), 1, sum)
    statpop <- divc(cbind.data.frame(vtot), dis)[, 1]

    targstat <- (statpop - statimean) / statpop

    return(targstat)

    }

    valobs <- fun1(df, dis)
    valsim <- sapply(1:nrep, function(i) fun1(df[sample(1:nrow(df)), ], dis))

    test1 <- as.randtest(valsim, valobs, alter = alter)
    test1$call <- "TPQE"
    return(test1)

}

################################
# Test for phylogenetic signal #
################################

# See Pavoine, S., Baguette, M., & Bonsall, M.B. (2010) Decomposition of trait diversity among the nodes of a phylogenetic tree. Ecological Monographs, In press.
# For further use of this function
# Last update: June 14 2010

rtest.decdiv <- function(phy, freq, dis = NULL, nrep = 99, vranking = "complexity", ties.method = "average", option = 1:3, optiontest = NULL, tol = 1e-08)
{

    #*******************************************************************************#
    #                         Checking of the parameters                            #
    #*******************************************************************************#
    
    if(!is.vector(freq)) stop("freq must be a unique vector")
    if (!is.numeric(nrep) | nrep <= 1) 
        stop("Non convenient nrep")
    if(sum(freq) < tol) stop("empty sample")
    if(any(freq < -tol)) stop("negative values in df")
    
    #*******************************************************************************#
    #                               Basic notations                                 #
    #*******************************************************************************#
    
    freq[freq < tol] <- 0
    freq <- freq / sum(freq)
    
    nsp <- length(phy$leaves)
  	nnodes <- length(phy$nodes)
    if(is.null(dis))
    dis <- as.dist(sqrt(2*matrix(1, nsp, nsp) - diag(rep(1, nsp))))
    
    #*******************************************************************************#
    #                               Node ranking                                    #
    #*******************************************************************************#

    complexity <- function(phy){   

	    matno <- as.data.frame(matrix(0, nnodes, nnodes))
	    rownames(matno) <- names(phy$nodes)
	    names(matno) <- names(phy$nodes)
        pathnodes <- phy$path[-(1:nsp)]
	    for(i in 1:nnodes){
	        matno[pathnodes[[i]], i] <- 1
	    }
        listno <- lapply(1:nnodes, function(i) names(matno)[matno[i, ] > 0])
        names(listno) <- names(phy$nodes)
        nbdes <- cbind.data.frame(lapply(phy$parts, function(x) prod(1:length(x))))
        compl <- lapply(listno, function(x) prod(nbdes[x]))
        compltab <- cbind.data.frame(compl)
        compltab <- cbind.data.frame(t(compltab))
        names(compltab) <- "complexity"
        return(compltab)
        
    }

    droot <- function(phy){
        roottab <- cbind.data.frame(phy$droot[-(1:nsp)])
        names(roottab) <- "droot"
        return(roottab)
    }

    if(is.numeric(vranking)){
        vrank <- as.data.frame(rank(vranking, ties.method = ties.method))
        names(vrank) <- "free"
    }
    else
    vrank <- sapply(vranking, function(x) rank(get(x)(phy), ties.method = ties.method))
   
    if(!any(option == 3))
        r1 <- length(option)
    else
        r1 <- length(option) + length(vranking) - 1

    #*******************************************************************************#
    #                       Field observations                                      #
    #*******************************************************************************#
    
    vobs <- decdiv(phy, freq, dis, tol = 1e-08)

    #*******************************************************************************#
    #                       Statistics for the four tests                           #
    #*******************************************************************************#    
    
    namnodes <- rownames(vobs)
    stat1 <- function(v){
        v <- v/sum(v)
        return(max(v))
    }
    stat2 <- function(v){
        v <- v/sum(v)
        fun1 <- function(m){
            return(abs(sum(sort(v)[1:m]) - m/nnodes))
        }
        return(max(unlist(lapply(1:nnodes, fun1))))
    }
    stat3 <- function(v){
        # this statistics has been sightly changed because the consideration of ties in Ollier et al.
        # was not explained, althought ties always happen with such a methodology. 
        funstat3 <- function(vrank1){
            v <- v/sum(v)
            return(sum(rank(vrank1, ties.method = ties.method)*v)/nnodes)
        }
        return(apply(vrank, 2, funstat3))
    }

    methods <- c("stat1", "stat2", "stat3")[option]
    
    #*******************************************************************************#
    #                      Statistics on field observations                         #
    #*******************************************************************************#
    
    statobs <- unlist(sapply(methods, function(x) get(x)(vobs[, 1])))

    #*******************************************************************************#
    #                              Permutation scheme                               #
    #*******************************************************************************#     
    
    funperm <- function(i){
        e <- sample(1:nsp)
        vtheo <- decdiv(phy, freq[e], as.dist(as.matrix(dis)[e, e]), tol = tol)
        stattheo <- unlist(sapply(methods, function(x) get(x)(vtheo[, 1])))
        return(stattheo)
    }
    
    tabsimu <- as.data.frame(t(cbind.data.frame(lapply(1:nrep, funperm))))
    rownames(tabsimu) <- paste("t", 1:nrep, sep="")
    if(r1 == 2 & methods[1] == "stat3")
        names(tabsimu) <- paste("stat3", names(tabsimu), sep=".")
    
    #*******************************************************************************#
    #                                     End                                       #
    #*******************************************************************************# 
    
    optiondefault <- c("greater", "greater", "two-sided", "two-sided", "two-sided")
    names(optiondefault) <- c("stat1", "stat2", "stat3.complexity", "stat3.droot", "stat3.free")
    
    if(r1 == 1)
    {
    if(!is.null(optiontest))
    return(as.randtest(obs = statobs, sim = tabsimu[, 1], alter = optiontest, call = "rtest.decdiv"))
    else
    return(as.randtest(tabsimu[, 1], statobs, alter = optiondefault[names(tabsimu)], call = "rtest.decdiv"))
    }
    
    if(!is.null(optiontest))
    return(as.krandtest(obs = statobs, sim = tabsimu, alter = optiontest, call = "rtest.decdiv"))
    else
    return(as.krandtest(obs = statobs, sim = tabsimu, alter = optiondefault[names(tabsimu)], 
        call = "rtest.decdiv"))

}

decdiv <- function(phy, df, dis = NULL, tol = 1e-08){
	
	if(is.vector(df)){
        df <- cbind.data.frame(df)
    }
    if(!is.data.frame(df)) stop("df should be a data frame")
    if(any(apply(df, 2, sum)<tol)) stop("null column in df")
    if(any(df < -tol)) stop("negative values in df")
    df[df < tol] <- 0
    df <- as.data.frame(apply(df, 2, function(x) x/sum(x)))
    
    disc2 <- function(samples, dis = NULL, structures = NULL, tol = 1e-08) 
    {
        if (!inherits(samples, "data.frame")) 
            stop("Non convenient samples")
        if (any(samples < 0)) 
            stop("Negative value in samples")
        if (any(apply(samples, 2, sum) < 1e-16)) 
            stop("Empty samples")
        if (!is.null(dis)) {
            if (!inherits(dis, "dist")) 
                stop("Object of class 'dist' expected for distance")
            if (!is.euclid(dis)) 
                warning("Euclidean property is expected for distance")
            dis <- as.matrix(dis)
            if (nrow(samples) != nrow(dis)) 
                stop("Non convenient samples")
        }
        if (is.null(dis)) 
            dis <- (matrix(1, nrow(samples), nrow(samples)) - diag(rep(1, 
                nrow(samples)))) * sqrt(2)
        if (!is.null(structures)) {
            if (!inherits(structures, "data.frame")) 
                stop("Non convenient structures")
            m <- match(apply(structures, 2, function(x) length(x)), 
                ncol(samples), 0)
            if (length(m[m == 1]) != ncol(structures)) 
                stop("Non convenient structures")
            m <- match(tapply(1:ncol(structures), as.factor(1:ncol(structures)), 
                function(x) is.factor(structures[, x])), TRUE, 0)
            if (length(m[m == 1]) != ncol(structures)) 
                stop("Non convenient structures")
        }
        Structutil <- function(dp2, Np, unit) {
            if (!is.null(unit)) {
                modunit <- model.matrix(~-1 + unit)
                sumcol <- apply(Np, 2, sum)
                Ng <- modunit * sumcol
                lesnoms <- levels(unit)
            }
            else {
                Ng <- as.matrix(Np)
                lesnoms <- colnames(Np)
            }
            sumcol <- apply(Ng, 2, sum)
            Lg <- t(t(Ng)/sumcol)
            colnames(Lg) <- lesnoms
            Pg <- as.matrix(apply(Ng, 2, sum)/nbhaplotypes)
            rownames(Pg) <- lesnoms
            deltag <- as.matrix(apply(Lg, 2, function(x) t(x) %*% 
                dp2 %*% x))
            ug <- matrix(1, ncol(Lg), 1)
            dg2 <- t(Lg) %*% dp2 %*% Lg - 1/2 * (deltag %*% t(ug) + 
                ug %*% t(deltag))
            colnames(dg2) <- lesnoms
            rownames(dg2) <- lesnoms
            return(list(dg2 = dg2, Ng = Ng, Pg = Pg))
        }
        Diss <- function(dis, nbhaplotypes, samples, structures) {
            structutil <- list(0)
            structutil[[1]] <- Structutil(dp2 = dis, Np = samples, 
                NULL)

            ###
            diss <- list(as.dist(structutil[[1]]$dg2))
            fun1 <- function(x){
                y <- x
                y[y<tol] <- 0
                return(y)
            }
            diss <- lapply(diss, fun1)
            diss <- lapply(diss, function(x) sqrt(2*x))
            ###

            if (!is.null(structures)) {
                for (i in 1:length(structures)) {
                    structutil[[i + 1]] <- Structutil(structutil[[1]]$dg2, 
                    structutil[[1]]$Ng, structures[, i])
                }
                ###
                diss <- c(diss, tapply(1:length(structures), factor(1:length(structures)),
                    function(x) (as.dist(structutil[[x + 1]]$dg2))))
                diss <- lapply(diss, fun1)
                diss <- lapply(diss, function(x) sqrt(2*x))
                ###
            }
            return(diss)
        }
        nbhaplotypes <- sum(samples)
        diss <- Diss(dis^2, nbhaplotypes, samples, structures)
        names(diss) <- c("samples", names(structures))
        if (!is.null(structures)) {
            return(diss)
        }
        return(diss$samples)
    }
   
    decdivV <- function(freq){
       nsp <- length(phy$leaves)
	   nnodes <- length(phy$nodes)
	   matno <- as.data.frame(matrix(0, nnodes, nsp))
	   rownames(matno) <- names(phy$nodes)
	   names(matno) <- names(phy$leaves)
	   for(i in 1:nsp){
		  matno[phy$path[[i]][-length(phy$path[[i]])], i] <- 1
	   }
	   matfr <- as.matrix(matno) %*% diag(freq)
	   matfr2 <- as.data.frame(t(matfr))
    	divno <- divc(matfr2, dis)
        matfr3 <- cbind.data.frame(matfr2, diag(freq))
        names(matfr3) <- c(names(matfr2), names(phy$leaves))
        matfr4 <- matfr3[, apply(matfr3, 2, sum)!=0]
        if(ncol(matfr4)==0) stop("only one species considered")
        discno <- disc2(matfr4, dis, tol = tol)
        lambdano <- apply(matfr4, 2, sum)
        prdist <- diag(lambdano)%*%as.matrix(discno^2/2)%*%diag(lambdano)
        colnames(prdist) <- rownames(prdist) <- names(matfr4)
        fun1 <- function(x){
            x <- x[apply(matfr3[x], 2, sum)!=0]
            if(length(x) == 1) return(0)
            else return(sum(prdist[x, x])/2)
        }
        res <- unlist(lapply(phy$parts, fun1))
        lambdano <- apply(matfr3, 2, sum)
        lambdano[lambdano < tol] <- 1
        res <- res * 1/as.vector(lambdano)[1:nnodes]
	   return(res)
    }
    return(apply(df, 2, decdivV))
}
