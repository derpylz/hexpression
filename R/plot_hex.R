plot_hex_minimal <- function(
    df,
    bins = 200,
    hex_cut = c(0.02, 0.98),
    aggr_fun = "mean",
    colorscale = rev(viridis::inferno(256))) {
  values <- ggplot2::ggplot_build(
    ggplot2::ggplot() +
      ggplot2::stat_summary_hex(
        data = df,
        ggplot2::aes(x = df[[1]], y = df[[2]], z = df[[3]]),
        fun = "mean",
        bins = bins,
        color = NA
      )
  )$data[[1]]$value
  if (!is.null(hex_cut)) {
    limits <- stats::quantile(values, c(min(hex_cut), max(hex_cut)),
      na.rm = TRUE
    )
  }
  p <- ggplot2::ggplot(df, ggplot2::aes(x = df[[1]], y = df[[2]], z = df[[3]])) +
    ggplot2::stat_summary_hex(fun = "mean", bins = bins, color = NA) +
    ggplot2::scale_fill_gradientn(
      colors = colorscale,
      limits = limits,
      oob = scales::squish
    )
  return(p)
}



plot_hex <- function(
    df,
    bins = 200,
    hex_cut = c(0.02, 0.98),
    aggr_fun = "mean",
    colorscale = rev(viridis::inferno(256)),
    color_label = "Expression",
    axis_labels = c("UMAP1", "UMAP2")) {
  p <- plot_hex_minimal(df, bins, hex_cut, aggr_fun, colorscale) +
    ggplot2::theme_classic() +
    ggplot2::xlab(axis_labels[1]) +
    ggplot2::ylab(axis_labels[2]) +
    ggplot2::guides(fill=ggplot2::guide_legend(title = color_label)) +
    ggplot2::coord_fixed()
  return(p)
}

#' @title HexPlot
#' @description Plot a hex bin plot of a feature in a Seurat object
#' @param srt Seurat object
#' @param feature Feature to plot
#' @param dims Dimensions to plot
#' @param aggr.fun Aggregation function
#' @param color.scale Color scale
#' @param bins Number of hex bins
#' @param min.cutoff Minimum expression cutoff
#' @param max.cutoff Maximum expression cutoff
#' @param reduction Which dimensionality reduction to use
#' @param split.by A factor in object metadata to split the feature plot by,
#' pass 'ident' to split by cell identity'. Used by Seurat::SplitObject
#' @param slot Which slot to use for the feature expression
#' @param coord.fixed Whether to use a fixed aspect ratio
#' @return A ggplot2 object if split.by is NULL, otherwise a patchwork object
#' @export
#' @examples
#' library(Seurat)
#' srt <- pbmc_small
#' HexPlot(srt, "CD3D", reduction = "tsne")
HexPlot <- function(
    srt, feature, dims = c(1, 2), aggr.fun = "mean",
    color.scale = rev(viridis::inferno(256)), bins = 200,
    min.cutoff = NA, max.cutoff = NA, reduction = "umap",
    split.by = NULL, slot = "data", coord.fixed = FALSE) {
  # Check if Seurat is installed, if not, warn the user and exit
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop(
      "This is a convenience function to use with Seurat objects.",
      "Please install Seurat first or use the plot_hex function."
    )
  }
  if (length(dims) != 2) {
    stop("dims must be a vector of length 2")
  }
  if (!is.null(split.by)) {
    srt_list <- Seurat::SplitObject(srt, split.by)
    plot_list <- lapply(srt_list, function(srt) {
      p <- HexPlot(
        srt, feature, dims, aggr.fun, color.scale, bins,
        min.cutoff, max.cutoff, reduction, NULL, slot, coord.fixed
      )
      p <- p +
        ggplot2::theme(legend.position = "none")
      return(p)
    })
    return(
      patchwork::wrap_plots(plot_list) +
        patchwork::plot_annotation(title = feature)
    )
  }
  coordinates <- Seurat::Embeddings(srt, reduction)[, dims]
  exprs <- Seurat::FetchData(srt, feature, slot = slot)[[feature]]
  df <- data.frame(coordinates, exprs)
  p <- plot_hex_minimal(
    df,
    bins,
    c(min.cutoff, max.cutoff),
     aggr.fun,
     color.scale
    ) +
    ggplot2::theme_classic() +
    ggplot2::labs(fill = "") +
    ggplot2::xlab(colnames(coordinates)[1]) +
    ggplot2::ylab(colnames(coordinates)[2]) +
    ggplot2::ggtitle(feature)
  return(p)
}
