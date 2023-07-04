#' @title Minimal hex plot from a data frame
#' @description Plot a hex bin plot from a data frame. Good starting point for
#' customizing the plot.
#' @param df Data frame with three columns, the first two are the coordinates
#' and the third is the expression
#' @param bins Number of hex bins
#' @param hex_cut Quantiles to use for the color scale
#' @param aggr_fun Aggregation function (e.g. mean, median)
#' @param colorscale Color scale
#' @return A ggplot2 object
#' @export
plot_hex_minimal <- function(
    df,
    bins = 200,
    hex_cut = c(0.02, 0.98),
    aggr_fun = "mean",
    colorscale = rev(viridis::inferno(256))) {
  dfnames <- colnames(df)
  values <- ggplot2::ggplot_build(
    ggplot2::ggplot() +
      ggplot2::stat_summary_hex(
        data = df,
        ggplot2::aes(x = .data[[dfnames[1]]], y = .data[[dfnames[2]]], z = .data[[dfnames[3]]]),
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
  p <- ggplot2::ggplot(
    df, ggplot2::aes(x = .data[[dfnames[1]]], y = .data[[dfnames[2]]], z = .data[[dfnames[3]]])) +
    ggplot2::stat_summary_hex(fun = "mean", bins = bins, color = NA) +
    ggplot2::scale_fill_gradientn(
      colors = colorscale,
      limits = limits,
      oob = scales::squish
    )
  return(p)
}

#' @title Hex plot from a data frame
#' @description Plot a hex bin plot from a data frame with some default
#' parameters for a nice looking plot, when not using a Seurat object.
#' @param df Data frame with three columns, the first two are the coordinates
#' and the third is the expression
#' @param bins Number of hex bins
#' @param hex_cut Quantiles to use for the color scale
#' @param aggr_fun Aggregation function (e.g. mean, median)
#' @param colorscale Color scale
#' @param color_label Color scale label
#' @param axis_labels Axis labels
#' @return A ggplot2 object
#' @export
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
    ggplot2::labs(fill = color_label) +
    ggplot2::coord_fixed()
  return(p)
}

#' @title HexPlot
#' @description Plot a hex bin plot of a feature in a Seurat object
#' @param srt Seurat object
#' @param features Feature(s) to plot. Either a single feature or a vector of
#' features
#' @param dims Dimensions to plot
#' @param aggr.fun Aggregation function (e.g. mean, median)
#' @param color.scale Color scale
#' @param bins Number of hex bins
#' @param min.cutoff Minimum expression cutoff
#' @param max.cutoff Maximum expression cutoff
#' @param reduction Which dimensionality reduction to use
#' @param split.by A factor in object metadata to split the feature plot by,
#' pass 'ident' to split by cell identity'. Used by Seurat::SplitObject
#' @param slot Which slot to use for the feature expression
#' @param coord.fixed Whether to use a fixed aspect ratio
#' @param combine Whether to combine the plots into a single plot using
#' patchwork or return a list of plots. Only applies when split.by is not NULL
#' or features is a vector of length > 1
#' @return A ggplot2 object if split.by is NULL, otherwise a patchwork object
#' @export
#' @importFrom rlang .data
#' @examples
#' library(Seurat)
#' srt <- pbmc_small
#' HexPlot(srt, "CD3D", reduction = "tsne")
HexPlot <- function(
    srt, features, dims = c(1, 2), aggr.fun = "mean",
    color.scale = rev(viridis::inferno(256)), bins = 200,
    min.cutoff = NA, max.cutoff = NA, reduction = "umap",
    split.by = NULL, slot = "data", coord.fixed = FALSE, combine = TRUE) {
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
  coordinates <- Seurat::Embeddings(srt, reduction)[, dims]
  if (!is.null(split.by)) {
    xlimits <- range(coordinates[, 1])
    ylimits <- range(coordinates[, 2])
    srt_list <- Seurat::SplitObject(srt, split.by)
    plot_list <- lapply(srt_list, function(srt) {
      p <- HexPlot(
        srt, features, dims, aggr.fun, color.scale, bins,
        min.cutoff, max.cutoff, reduction, NULL, slot, coord.fixed, FALSE
      )
      if (length(features) > 1) {
        p <- lapply(p, function(p_i) {
          p_i <- p_i +
            ggplot2::labs(subtitle = as.character(unique(srt@meta.data[[split.by]]))) +
            ggplot2::theme(legend.position = "none") +
            ggplot2::xlim(floor(xlimits[1]), ceiling(xlimits[2])) +
            ggplot2::ylim(floor(ylimits[1]), ceiling(ylimits[2]))
          return(p_i)
        })
        return(p)
      }
      p <- p +
        ggplot2::ggtitle(as.character(unique(srt@meta.data[[split.by]]))) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::xlim(floor(xlimits[1]), ceiling(xlimits[2])) +
        ggplot2::ylim(floor(ylimits[1]), ceiling(ylimits[2]))
      return(p)
    })
    if (length(features) > 1) {
      plot_list <- unlist(plot_list, recursive = FALSE)
    }
    if (combine) {
      if (length(features) > 1) {
        return(
          patchwork::wrap_plots(plot_list, ncol = length(srt_list), byrow = F) +
            patchwork::plot_annotation(title = split.by)
        )
      } else {
        return(
          patchwork::wrap_plots(plot_list) +
            patchwork::plot_annotation(title = features)
        )
      }
    } else {
      return(plot_list)
    }
  }
  if (length(features) > 1) {
    plot_list <- lapply(features, function(feature) {
      p <- HexPlot(
        srt, feature, dims, aggr.fun, color.scale, bins,
        min.cutoff, max.cutoff, reduction, NULL, slot, coord.fixed, FALSE
      )
      p <- p +
        ggplot2::ggtitle(feature) +
        ggplot2::theme(legend.position = "none")
      return(p)
    })
    if (combine) {
      return(
        patchwork::wrap_plots(plot_list)
      )
    } else {
      return(plot_list)
    }
  }
  exprs <- Seurat::FetchData(srt, features, slot = slot)[[features]]
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
    ggplot2::ggtitle(features)
  return(p)
}
