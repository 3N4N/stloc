weightMatrix_nD = function(x, span = 0.5)
{

    ncells = nrow(x)
    coords = as.matrix(x)
    d = as.matrix(dist(coords))

    # extract a weights vector per cell

    W_raw = sapply(seq_len(ncells), function(cell) {
                       dvec = d[cell,]
                       vals = rep(0, ncells)
                       vals[order(dvec)[1:ceiling(span*ncells)]] = seq(1, 0, length.out = ceiling(span*ncells))
                       # vals[order(dvec)[1:ceiling(span*ncells)]] = 1
                       return(vals)
            }, simplify = FALSE)

    W = do.call(rbind, W_raw)

    return(W)
}

weightMatrix_gaussian  = function(x, l=20)
{

    ncells = nrow(x)
    coords = as.matrix(x)
    d = as.matrix(dist(coords))

    # extract a weights vector per cell

    W_raw = sapply(seq_len(ncells), function(cell) {
        dvec = d[cell,]

        vals = rep(0, ncells)
        # l=20
        for(i in 1:length(vals)) {
            vals[i]= exp(-d[cell,i]^2/l^2)
        }
        return(vals)
    }, simplify = FALSE)

    W = do.call(rbind, W_raw)

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

maxSingVal = function(x, w=1) {
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
    zscore_for_single_gene = matrix(nrow = nrow(x), ncol = 1);

    sm = 0
    for(i in 1:(nrow(x))) {
        zscore_for_single_gene[i, 1] = (x[i, index] - meanOfZenes[i]) / sdOfZenes[i];
        sm = sm + zscore_for_single_gene[i, 1]
    }
    aggrZscore = sm / sqrt(nrow(x));

    return (aggrZscore)
}

plotdf = function(df, vals1, vals2, vals3, label1, label2, label3) {

    plot_vals1 <- ggplot(df, aes(x = x, y = -y)) +
        geom_point(aes(colour = vals1), size = 5) +
        theme_minimal() +
        theme(panel.grid = element_blank()) +
        theme(axis.text = element_blank()) +
        xlab("") +
        ylab("") +
        labs(colour = "") +
        theme(legend.position = "bottom") +
        theme(plot.title = element_text(hjust = 0.5, face = "italic")) +
        scale_color_viridis_c(option = "plasma") +
        coord_fixed() +
        guides(colour = guide_colourbar(title.position = "top",
                                        title.hjust = 0.5)) +
        theme(legend.key.width = unit(0.5, "inches")) +
        theme(plot.title = element_text(size = 20)) +
        theme(axis.title = element_text(size = 15)) +
        theme(legend.title = element_text(size = 15)) +
        labs(colour = label1) +
        NULL

    plot_vals2 <- ggplot(df, aes(x = x, y = -y)) +
        # geom_point(aes(colour = vals2), size = 5) +
        geom_point(aes(colour = -log10(vals2)), size = 5) +
        theme_minimal() +
        theme(panel.grid = element_blank()) +
        theme(axis.text = element_blank()) +
        xlab("") +
        ylab("") +
        labs(colour = "") +
        theme(legend.position = "bottom") +
        theme(plot.title = element_text(hjust = 0.5, face = "italic")) +
        scale_color_viridis_c(option = "plasma") +
        coord_fixed() +
        guides(colour = guide_colourbar(title.position = "top",
                                        title.hjust = 0.5)) +
        theme(legend.key.width = unit(0.5, "inches")) +
        theme(plot.title = element_text(size = 20)) +
        theme(axis.title = element_text(size = 15)) +
        theme(legend.title = element_text(size = 15)) +
        labs(colour = label2) +
        NULL

    plot_vals3 <- ggplot(df, aes(x = x, y = -y)) +
        geom_point(aes(colour = -log10(vals3)), size = 5) +
        theme_minimal() +
        theme(panel.grid = element_blank()) +
        theme(axis.text = element_blank()) +
        xlab("") +
        ylab("") +
        labs(colour = "") +
        theme(legend.position = "bottom") +
        theme(plot.title = element_text(hjust = 0.5, face = "italic")) +
        scale_color_viridis_c(option = "plasma") +
        coord_fixed() +
        guides(colour = guide_colourbar(title.position = "top",
                                        title.hjust = 0.5)) +
        theme(legend.key.width = unit(0.5, "inches")) +
        theme(plot.title = element_text(size = 20)) +
        theme(axis.title = element_text(size = 15)) +
        theme(legend.title = element_text(size = 15)) +
        labs(colour = label3) +
        NULL


    vals1_leg = as_ggplot(get_legend(plot_vals1))
    vals2_leg = as_ggplot(get_legend(plot_vals2))
    vals3_leg = as_ggplot(get_legend(plot_vals3))

    gridExtra::grid.arrange(plot_vals1 + theme(legend.position = "none")
                            + theme(plot.margin = margin(10,0,-10,0)),
                            plot_vals2 + theme(legend.position = "none")
                            + theme(plot.margin = margin(10,0,-10,0)),
                            plot_vals3 + theme(legend.position = "none")
                            + theme(plot.margin = margin(10,0,-10,0)),
                            vals1_leg, vals2_leg, vals3_leg,
                            layout_matrix = matrix(c(1,1,2,2,3,3,
                                                     1,1,2,2,3,3,
                                                     4,4,5,5,6,6),
                                                   ncol = 6,
                                                   byrow = TRUE))
}

ploteig = function(df.eig, vals, loc, label) {

    if (loc != 0) vals[loc] = NA

    plot.vals = ggplot(df.eig, aes(x = x, y = -y)) +
        geom_point(aes(colour = vals), size = 4) +
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
        theme(plot.title = element_text(size = 20)) +
        theme(axis.title = element_text(size = 15)) +
        theme(legend.title = element_text(size = 15)) +
        labs(colour = label) +
        NULL

    vals.leg = as_ggplot(get_legend(plot.vals))

    gridExtra::grid.arrange(plot.vals + theme(legend.position = "none")
                            + theme(plot.margin = margin(10,0,-10,0)),
                            vals.leg,
                            layout_matrix = rbind(c(1,1,1),
                                                  c(1,1,1),
                                                  c(2,2,2)))
}
