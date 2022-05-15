#  ----------------------------------------------------------------------
#  Contains simple helper functions that are used in other scripts.
#  ----------------------------------------------------------------------

library(ggplot2)
library(gridExtra)
library(grid)
library(scico)

minmax <- function(x) {
    if (max(x) != min(x)) {
        return((x - min(x)) / (max(x) - min(x)))
    } else {
        return(x * 0)
    }
}

weightMatrix.nD <- function(x, span = 0.5) {
    ncells <- nrow(x)
    coords <- as.matrix(x)
    d <- as.matrix(dist(coords))

    # extract a weights vector per cell

    W.raw <- sapply(seq_len(ncells), function(cell) {
        dvec <- d[cell, ]
        vals <- rep(0, ncells)
        vals[order(dvec)[1:ceiling(span * ncells)]] <- seq(1, 0, length.out = ceiling(span * ncells))
        # vals[order(dvec)[1:ceiling(span*ncells)]] <- 1
        return(vals)
    }, simplify = FALSE)

    W <- do.call(rbind, W.raw)

    return(W)
}

weightMatrix.gaussian <- function(coords, l = 20, cell = 0) {
    ncells <- nrow(coords)
    coords <- as.matrix(coords)
    d <- as.matrix(dist(coords))

    W.raw <- sapply(seq_len(ncells), function(cell) {
        dvec <- d[cell, ]
        vals <- rep(0, ncells)
        for (i in 1:length(vals)) {
            # vals[i] = exp(-d[cell,i]^2 / l^2)
            vals[i] <- (1 / (l * sqrt(2 * pi))) * exp(-(1 / 2) * d[cell, i]^2 / l^2)
        }
        return(vals)
    }, simplify = FALSE)

    W <- do.call(rbind, W.raw)

    return(W)
}


corTaylor <- function(x, w = 1) {
    if (!inherits(x, "matrix")) {
        stop("Input must be inherit ’matrix’ class.")
    }

    if (length(w) == 1) {
        w <- rep(1, ncol(x))
    }

    x <- apply(x, 1, function(i) w * i)
    d <- ncol(x)

    return((1 / sqrt(d)) * sd(eigen(cor(x), only.values = TRUE)$values))
}


tmaxEigenVal <- function(x, w = 1) {
    if (length(w) == 1) {
        w <- rep(1, ncol(x))
    }

    x <- apply(x, 1, function(i) w * i)

    covmat <- cov(x)

    return((sum(colSums(covmat))) / ncol(x))
}

maxEigenVal <- function(x, w = 1) {
    if (length(w) == 1) {
        w <- rep(1, ncol(x))
    }

    x <- apply(x, 1, function(i) w * i)

    eigvals <- eigen(t(x) %*% x, only.values = TRUE)$values
    # eigvals <- eigen(cov(x),only.values=TRUE)$values

    return(max(eigvals))
}

maxSingVal <- function(x, w = 1) {
    if (!inherits(x, "matrix")) {
        stop("Input must be inherit ’matrix’ class.")
    }

    if (length(w) == 1) {
        w <- rep(1, ncol(x))
    }

    x <- apply(x, 1, function(i) w * i)

    return(max(svd(t(x))$d))
}
minSingVal <- function(x, w = 1) {
    if (!inherits(x, "matrix")) {
        stop("Input must be inherit ’matrix’ class.")
    }

    if (length(w) == 1) {
        w <- rep(1, ncol(x))
    }

    x <- apply(x, 1, function(i) w * i)

    return(min(svd(t(x))$d))
}

maxByMinSingVal <- function(x, w = 1) {
    if (!inherits(x, "matrix")) {
        stop("Input must be inherit ’matrix’ class.")
    }

    if (length(w) == 1) {
        w <- rep(1, ncol(x))
    }

    x <- apply(x, 1, function(i) w * i)
    maxVal <- max(svd(t(x))$d)
    minVal <- min(svd(t(x))$d)
    if (minVal == 0) {
        minVal <- 1^(-120)
    }
    maxByMin <- maxVal / minVal
    if (maxByMin < 0) {
        message("negative value ")
    }
    # message("min val is", min(svd(t(x))$d))
    return(log10(maxVal / minVal))

    # return(max(svd(t(x))$d) / min(svd(t(x))$d)))
}

plotmeigs <- function(df, meigs) {
    n <- length(meigs)
    i <- 1
    plot.vals <- lapply(1:n, function(i) {
        ggplot(df, aes(x = x, y = y)) +
            geom_point(aes(colour = meigs[[i]]), size = 3) +
            theme_minimal() +
            theme(panel.grid = element_blank()) +
            theme(axis.text = element_blank()) +
            xlab("") +
            ylab("") +
            labs(colour = "") +
            theme(legend.position = "bottom") +
            theme(plot.title = element_text(hjust = 0.5, face = "italic")) +
            # scale_color_gradient(low = "#ffffd2", high = "#7c1824") +
            # scale_color_viridis_c(option = "lajolla", na.value="red", direction=-1) +
            scale_color_scico(palette = "lajolla") +
            coord_fixed() +
            guides(colour = guide_colourbar(title.position = "top", title.hjust = 0.5)) +
            theme(legend.key.width = unit(0.5, "inches")) +
            theme(legend.text = element_text(size = 12)) +
            theme(plot.margin = margin(-10, 0, -10, 0)) +
            theme(plot.title = element_text(size = 20)) +
            theme(axis.title = element_text(size = 15)) +
            theme(legend.title = element_text(size = 20)) +
            labs(colour = colnames(meigs)[i]) +
            NULL
    })

    grid.arrange(grobs = plot.vals, nrow = 1, ncol = n)

}

plotvals <- function(n, df, title, vals, labels, size, nrow, ncol) {
    plot.vals <- lapply(1:n, function(i) {
        ggplot(df, aes(x = x, y = y)) +
            geom_point(aes(colour = vals[[i]]), size = size) +
            theme_minimal() +
            theme(panel.grid = element_blank()) +
            theme(axis.text = element_blank()) +
            xlab("") +
            ylab("") +
            labs(colour = "") +
            theme(legend.position = "bottom") +
            theme(plot.title = element_text(hjust = 0.5, face = "italic")) +
            scale_color_viridis_c(option = "plasma", na.value="red", direction=-1) +
            # scale_color_binned() +
            # scale_color_binned(n.breaks=4, nice.breaks=F, labels=function(x) {
            #                        if (x >= 10) round(x)
            #                        else sprintf("%.2f", x)
            #     }) +
            coord_fixed() +
            guides(colour = guide_colourbar(title.position = "top", title.hjust = 0.5)) +
            theme(legend.key.width = unit(0.5, "inches")) +
            theme(legend.text = element_text(size = 12)) +
            theme(plot.margin = margin(-10, 0, 10, 0)) +
            theme(plot.title = element_text(size = 20)) +
            theme(axis.title = element_text(size = 15)) +
            theme(legend.title = element_text(size = 20)) +
            labs(colour = labels[i]) +
            NULL
    })


    title <- textGrob(title, gp = gpar(fontsize = 20))

    # Add a zeroGrob of height 2cm on top of the title
    title <- arrangeGrob(zeroGrob(), title,
        widths = unit(1, "npc"),
        heights = unit(c(2, 1), c("cm", "npc")),
        as.table = FALSE
    )

    grid.arrange(grobs = plot.vals, nrow = nrow, ncol = ncol, top = title)
}

scatterPlot <- function(scatterdf) {
    plotvar <- ggplot(scatterdf, aes(x = x, y = y)) +
        geom_point(size = 2, shape = 23) +
        labs(x = "Cell Count ", y = "Singular value")

    plotbox <- ggplot(scatterdf, aes(x = x, y = y, group = x)) +
        geom_boxplot() +
        labs(x = "Cell Count ", y = "Singular value")

    do.call("grid.arrange", c(list(plotvar, plotbox), ncol = 2))
}

scatterPlotSpotlight <- function(scatterdf) {
    plotvar <- ggplot(scatterdf, aes(x = x, y = y)) +
        geom_point(size = 2, shape = 23) +
        scale_y_continuous(limits = c(0, 1)) +
        labs(x = "Cell Count ", y = "Fraction")

    plotbox <- ggplot(scatterdf, aes(x = x, y = y, group = x)) +
        geom_boxplot() +
        scale_y_continuous(limits = c(0, 1)) +
        labs(x = "Cell Count ", y = "Fraction")

    do.call("grid.arrange", c(list(plotvar, plotbox), ncol = 2))
}
