minmax <- function(x) {
    (x - min(x)) / (max(x) - min(x))
}

weightMatrix.nD <- function(x, span = 0.5)
{

    ncells <- nrow(x)
    coords <- as.matrix(x)
    d <- as.matrix(dist(coords))

    # extract a weights vector per cell

    W.raw <- sapply(seq_len(ncells), function(cell) {
                       dvec <- d[cell,]
                       vals <- rep(0, ncells)
                       vals[order(dvec)[1:ceiling(span*ncells)]] <- seq(1, 0, length.out = ceiling(span*ncells))
                       # vals[order(dvec)[1:ceiling(span*ncells)]] <- 1
                       return(vals)
            }, simplify = FALSE)

    W <- do.call(rbind, W.raw)

    return(W)
}

weightMatrix.gaussian  <- function(coords, l=20, cell=0)
{

    ncells <- nrow(coords)
    coords <- as.matrix(coords)
    d <- as.matrix(dist(coords))

    W.raw <- sapply(seq_len(ncells), function(cell) {
        dvec <- d[cell,]
        vals <- rep(0, ncells)
        for(i in 1:length(vals)) {
            vals[i] <- (1 / (l * sqrt(2*pi))) * exp(-(1/2) * d[cell,i]^2 / l^2)
        }
        return(vals)
    }, simplify = FALSE)

    W <- do.call(rbind, W.raw)

    return(W)
}


corTaylor <- function(x, w = 1)
{
    if(!inherits(x, "matrix")) {
        stop("Input must be inherit ’matrix’ class.")
    }

    if (length(w) == 1) {
        w <- rep(1, ncol(x))
    }

    x <- apply(x, 1, function(i) w*i)
    d <- ncol(x)

    return ((1/sqrt(d)) * sd(eigen(cor(x),only.values=TRUE)$values))
}


maxEigenVal <- function(x, w=1)
{
    if(!inherits(x,"matrix")) {
        stop("Input must be inherit ’matrix’ class.")
    }

    if (length(w) == 1) {
        w <- rep(1, ncol(x))
    }

    d <- nrow(x)
    x <- apply(x, 1, function(i) w*i)

    return (max(eigen(t(x) %*% x,only.values=TRUE)$values))
}

maxSingVal <- function(x, w=1) {
    if(!inherits(x,"matrix")) {
        stop("Input must be inherit ’matrix’ class.")
    }

    if (length(w) == 1) {
        w <- rep(1, ncol(x))
    }

    d <- nrow(x)
    x <- apply(x, 1, function(i) w*i)

    return (max(svd(t(x))$d))
}

zScore <- function(x, index, meanOfZenes, sdOfZenes)
{
    zscore_for_single_gene <- matrix(nrow = nrow(x), ncol = 1);

    sm <- 0
    for(i in 1:(nrow(x))) {
        zscore_for_single_gene[i, 1] <- (x[i, index] - meanOfZenes[i]) / sdOfZenes[i];
        sm <- sm + zscore_for_single_gene[i, 1]
    }
    aggrZscore <- sm / sqrt(nrow(x));

    return (aggrZscore)
}


plotvals <- function(n, df, vals, labels, size) {
    plot.vals <- lapply(1:n, function(i) {
        ggplot(df, aes(x = -x, y = -y)) +
            geom_point(aes(colour = vals[[i]]), size = size) +
            theme_minimal() +
            theme(panel.grid = element_blank()) +
            theme(axis.text = element_blank()) +
            xlab("") +
            ylab("") +
            labs(colour = "") +
            theme(legend.position = "bottom") +
            theme(plot.title = element_text(hjust = 0.5, face = "italic")) +
            scale_color_viridis_c(option = "plasma", na.value="red") +
            coord_fixed() +
            guides(colour = guide_colourbar(title.position = "top", title.hjust = 0.5)) +
            theme(legend.key.width = unit(0.5, "inches")) +
            # theme(plot.margin = margin(10,0,10,0)) +
            theme(plot.title = element_text(size = 20)) +
            theme(axis.title = element_text(size = 15)) +
            theme(legend.title = element_text(size = 15)) +
            labs(colour = labels[i]) +
            NULL
    })

    # grid.arrange(grobs=plot.vals, ncol=n)
    do.call("grid.arrange", c(plot.vals, ncol=n))
}

scatterPlot<-function(df){
    plotvar <-ggplot(df, aes(x=x, y=y)) +
        geom_point(size=2, shape=23) +
        labs(x="Cell Count " ,y= "Maximum Eigen value")
       


    plot(plotvar)


}
