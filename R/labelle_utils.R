.actionbutton_biocstyle <- "color: #ffffff; background-color: #0092AC"
.helpbutton_biocstyle <- "color: #0092AC; background-color: #FFFFFF; border-color: #FFFFFF"

editor_to_vector_sanitized <- function(txt) {
  rn <- strsplit(txt, split="\n")[[1]]
  rn <- sub("#.*", "", rn)
  rn <- sub("^ +", "", rn)
  sub(" +$", "", rn)
}

