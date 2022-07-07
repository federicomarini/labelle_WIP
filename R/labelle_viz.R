annotation_alluvial <- function(sce,
                                colData1,
                                colData2,
                                color_by) {

  # function checks on
  sce
  colData1
  colData2
  color_by

  cd_anno <- colData(sce)

  anno_left <- colData1
  anno_right <- colData2
  anno_colorfill <- color_by

  tbl_anno <- table(
    cd_anno[[anno_left]],
    cd_anno[[anno_right]]
  )
  names(dimnames(tbl_anno)) <-
    c(anno_left, anno_right)

  df_tbl_anno <- as.data.frame(tbl_anno)
  colnames(df_tbl_anno)[3] <- "Number of cells"

  p <- ggplot(df_tbl_anno,
              aes(y = `Number of cells`,
                  axis1 = .data[[anno_left]], axis2 = .data[[anno_right]])) +
    theme_bw() +
    geom_alluvium(aes(fill = .data[[anno_colorfill]]), width = 1/12) +
    geom_stratum(width = 1/12, fill = "black", color = "grey") +
    geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_x_discrete(limits = c(colData1, colData2), expand = c(.2, .2)) +
    scale_fill_brewer(type = "qual", palette = "Set1")

  p
}



annotation_chord <- function(sce,
                             colData1,
                             colData2,
                             color_by) {

  # function checks on
  sce
  colData1
  colData2
  color_by

  cd_anno <- colData(sce)

  anno_left <- colData1
  anno_right <- colData2
  anno_colorfill <- color_by

  tbl_anno <- table(
    cd_anno[[anno_left]],
    cd_anno[[anno_right]]
  )

  names(dimnames(tbl_anno)) <-
    c(anno_left, anno_right)

  df_tbl_anno <- as.data.frame(tbl_anno)

  m <- as.matrix(tbl_anno)

  message("here1")
  message(dim(m))
  row_cols <- hcl.colors(nrow(m), "Temps")
  message("here2")

  if(nrow(m) > 1) {
    cols_grid_rows <- hcl.colors(nrow(m), "Temps")
  } else {
    cols_grid_rows <- "lightgrey"
  }

  if(ncol(m) > 1) {
    cols_grid_columns <- hcl.colors(ncol(m), "Temps")
  } else {
    cols_grid_columns <- "lightgrey"
  }


  grid_cols <- setNames(
    object = c(
      cols_grid_rows,
      cols_grid_columns
    ),
    nm = c(rownames(m), colnames(m))
  )

  message("here")

  circlize::chordDiagram(m,
                         col = NULL,
                         row.col = row_cols,
                         grid.col = grid_cols,
                         transparency = 0.5,
                         link.lwd = 1,    # Line width
                         link.lty = 1,    # Line type
                         link.border = 1) # Border color)


}

annotation_chord_interactive <- function(sce,
                             colData1,
                             colData2,
                             color_by) {

  # function checks on
  sce
  colData1
  colData2
  color_by

  cd_anno <- colData(sce)

  anno_left <- colData1
  anno_right <- colData2
  anno_colorfill <- color_by

  tbl_anno <- table(
    cd_anno[[anno_left]],
    cd_anno[[anno_right]]
  )

  names(dimnames(tbl_anno)) <-
    c(anno_left, anno_right)

  df_tbl_anno <- as.data.frame(tbl_anno)

  m <- as.matrix(tbl_anno)

  message("here1")
  message(dim(m))
  row_cols <- hcl.colors(nrow(m), "Temps")
  message("here2")

  if(nrow(m) > 1) {
    cols_grid_rows <- hcl.colors(nrow(m), "Temps")
  } else {
    cols_grid_rows <- "lightgrey"
  }

  if(ncol(m) > 1) {
    cols_grid_columns <- hcl.colors(ncol(m), "Temps")
  } else {
    cols_grid_columns <- "lightgrey"
  }


  grid_cols <- setNames(
    object = c(
      cols_grid_rows,
      cols_grid_columns
    ),
    nm = c(rownames(m), colnames(m))
  )

  message("here")

  chorddiag::chorddiag(m, type = "bipartite")
}



annotation_donut <- function(sce,
                             colData1) {

  cd_anno <- colData(sce)

  df_anno <- as.data.frame(
    table(cd_anno[[colData1]])
  )

  category <- colData1
  colnames(df_anno) <- c(category, "Number of cells")

  # Compute percentages
  df_anno$fraction <- df_anno$`Number of cells` / sum(df_anno$`Number of cells`)

  # Compute the cumulative percentages (top of each rectangle)
  df_anno$ymax <- cumsum(df_anno$fraction)

  # Compute the bottom of each rectangle
  df_anno$ymin <- c(0, head(df_anno$ymax, n=-1))

  # Compute label position
  df_anno$labelPosition <- (df_anno$ymax + df_anno$ymin) / 2

  # Compute a good label
  df_anno$label <- paste0(df_anno[[category]], "\n value: ", df_anno$`Number of cells`)

  # Make the plot
  p <- ggplot(df_anno, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=.data[[category]])) +
    geom_rect() +
    geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
    scale_fill_brewer(type = "qual", palette = "Set1") +
    coord_polar(theta="y") +
    xlim(c(2, 4)) +
    theme_void() +
    theme(legend.position = "none")

  p
}
