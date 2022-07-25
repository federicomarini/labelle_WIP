.actionbutton_biocstyle <- "color: #ffffff; background-color: #0092AC"
.helpbutton_biocstyle <- "color: #0092AC; background-color: #FFFFFF; border-color: #FFFFFF"

editor_to_vector_sanitized <- function(txt) {
  rn <- strsplit(txt, split="\n")[[1]]
  rn <- sub("#.*", "", rn)
  rn <- sub("^ +", "", rn)
  sub(" +$", "", rn)
}


celltype_2_html <- function(celltype,
                            source = NULL) {
  celltypist_button <- .link2celltypist(celltype)

  mycontent <- paste0(
    tags$b(celltype), tags$br(),
    "Link to the celltypist: ", celltypist_button, tags$br(),
    "Link to the GeneCards database: ", "some gene ", tags$br(),
    tags$code("Something about the cell type")
  )
  return(HTML(mycontent))
}

.link2celltypist <- function(val) {
  sprintf(
    '<a href = "https://www.celltypist.org/encyclopedia/?celltype%s" target = "_blank" class = "btn btn-primary" style = "%s">%s</a>',
    val,
    .actionbutton_biocstyle,
    val
  )
}
