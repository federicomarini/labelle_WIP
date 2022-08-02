library(shiny)
# library(scRNAseq)
library(scater)
library(scran)
library(dplyr)
library(shinyAce)
library(shinydashboard)
library(ggalluvial)

# sce <- HCATonsilData::HCATonsilData(cellType = "epithelial")


cell_dictionary <- c(
  "cell A",
  "cell B",
  "cell D"
)


#' Title
#'
#' @param sce
#' @param cell_dictionary
#'
#' @return
#' @export
#'
#' @examples
labelle <- function(sce, cell_dictionary = NULL) {
  rv <- reactiveValues()
  rv$anno_sce <- sce


  rv$cell_dictionary <- cell_dictionary

  rv$cell_selected <- NULL

  # initializing the log of the operations
  rv$log_entries <- data.frame(
    timestamp = vector(mode = "character"),
    label_assigned = vector(mode = "character"),
    which_cells = vector(mode = "character"),
    rationale = vector(mode = "character"),
    user = vector(mode = "character")
  )

  # here: some sanitization on the cell_dictionary?


  #
  isolate({
    colData(rv$anno_sce)[["new_anno"]] <- "BLANK"
    colData(rv$anno_sce)[["timestamp_annotation"]] <- Sys.time()
    colData(rv$anno_sce)[["anno_rationale"]] <- ""
    # coldata_names <- c("new_anno", colnames(colData(sce)))
    coldata_names <- colnames(colData(rv$anno_sce))
  })

  labelle_ui <- shinydashboard::dashboardPage(
    header = shinydashboard::dashboardHeader(disable = TRUE),
    sidebar = shinydashboard::dashboardSidebar(disable = TRUE),
    body = shinydashboard::dashboardBody(
      # activate rintrojs based tours
      rintrojs::introjsUI(),
      ## Define output size and style of error messages
      tags$head(
        tags$style(
          HTML(
            ".shiny-output-error-validation {
            font-size: 15px;
            color: forestgreen;
            text-align: center;
            }
            "
          )
        )
      ),

      fluidRow(
        column(
          width = 11
        ),
        column(
          width = 1,
          actionButton(
            "tour_labelle",
            label = "", icon = icon("question-circle"),
            style = .helpbutton_biocstyle
          )
        )
      ),

      fluidRow(
        shinydashboard::box(
          title = "Provided data",
          width = 12,
          collapsible = TRUE,
          collapsed = FALSE,
          fluidRow(
            column(
              width = 6,
              plotOutput("reddim_plot"),

            ),
            column(
              width = 6,
              DT::dataTableOutput("dt_coldata"),
              hr()
              # ,
              # DT::dataTableOutput("dt_anno")
            )
          )
        )
      ),
      hr(),
      fluidRow(
        shinydashboard::box(
          title = "Annotating the data",
          width = 12,
          collapsible = TRUE,
          collapsed = FALSE,
          fluidRow(
            column(
              width = 6,
              textInput(inputId = "which_cells",
                        label = "which cells to update"),
              shinyAce::aceEditor(
                "editor_cells",
                theme = "solarized_light",
                height = "200px",
                readOnly = FALSE,
                wordWrap = TRUE,
                placeholder = paste(
                  c(
                    "Enter some cell identifiers",
                    "For example:",
                    head(colnames(sce), 3)
                  ),
                  collapse = "\n"
                )
              ),
              actionButton("load_cells",
                           label = "Select cells subset",
                           icon = icon("upload"),
                           style = .actionbutton_biocstyle),

              actionButton("btn_randomizer",
                           label = "Select some random cells",
                           icon = icon("dice"),
                           style = .actionbutton_biocstyle),

              valueBoxOutput("vb_cells_subset",
                             width = 12),

              selectizeInput(
                "id_label",
                label = "Label to apply",
                choices = cell_dictionary,
                multiple = TRUE,
                options = list(create = TRUE)
              ),

              actionButton(
                inputId = "btn_import_dictionary",
                label = "import the dictionary",
                icon = icon("upload"),
                style = .actionbutton_biocstyle
              ),

              textAreaInput(
                inputId = "anno_rationale",
                label = "Enter the rationale for annotating the selected cells",
                placeholder =
                  "You can specify for example: 'Expression of marker genes CD3D, TRBC2'"
              ),



              actionButton(inputId = "btn_label",
                           label = "labelle!",
                           icon = icon("edit"),
                           style = .actionbutton_biocstyle),

              actionButton(inputId = "btn_save_dictionary",
                           label = "Save the current dictionary!",
                           icon = icon("database"),
                           style = .actionbutton_biocstyle)
            ),
            column(
              width = 6,
              DT::dataTableOutput("dt_anno")
              # DT::dataTableOutput("dt_anno")
            )

          )
        )
      ),

      fluidRow(
        shinydashboard::box(
          title = "Annotation summary",
          width = 12,
          collapsible = TRUE,
          collapsed = FALSE,
          fluidRow(
            column(
              width = 4,
              selectizeInput(
                "cd_anno",
                label = "col data to color by",
                choices = coldata_names,
                selected = "new_anno"
              )
            ),
            column(
              width = 4,
              plotOutput("anno_plot")
            ),
            column(
              width = 4,
              verbatimTextOutput("overview_anno"),
              actionButton(inputId = "btn_download_anno_sce",
                           label = "Download the annotated sce",
                           icon = icon("download"),
                           style = .actionbutton_biocstyle),
              verbatimTextOutput("overview_log"),
              actionButton(inputId = "btn_export_log",
                           label = "Download the log of labelle's usage",
                           icon = icon("pen"),
                           style = .actionbutton_biocstyle),
              actionButton(inputId = "btn_store_log",
                           label = "Store the log as metadata",
                           icon = icon("pen"),
                           style = .actionbutton_biocstyle),
              verbatimTextOutput("log_as_metadata")

            )
          ),
          fluidRow(
            column(
              width = 6,
              uiOutput("ui_anno_analytics")

            ),
            column(
              width = 6,
              h2("som'thin' to think about!"),
              uiOutput("ui_ideas_collection"),

              hr(),
              uiOutput("ui_cl_matching")
            )
          )
        )

      ),

      fluidRow(
        shinydashboard::box(
          title = "Exploring the cell types",
          width = 12,
          collapsible = TRUE,
          collapsed = TRUE,
          fluidRow(
            column(
              width = 6,
              uiOutput("ui_celltypedia")
            ),
            column(
              width = 6,
              "some links or some other info",
              HTML(
                paste0(
                  tags$b("Even something"), tags$br(),
                  "Like this: ", tags$br(),
                  tags$code(fortunes::fortune())
                )
              )
            )
          )
        )
      )
    )
  )


  labelle_server <- function(input, output, session) {
    output$reddim_plot <- renderPlot({
      plotUMAP(sce, colour_by = "annotation_20220215")
    })

    output$dt_coldata <- DT::renderDataTable({
      dt <- DT::datatable(
        as.data.frame(colData(sce)),
        options = list(scrollX = TRUE)
      )
    })

    output$vb_cells_subset <- renderValueBox({
      valueBox(
        value = "Selected cells to annotate",
        subtitle = length(rv$cell_selected),
        icon = icon("glyphicon-option-vertical", lib = "glyphicon")
      )
    })

    output$dt_anno <- DT::renderDataTable({
      my_df <- as.data.frame(colData(rv$anno_sce)) %>%
        relocate(c("new_anno", "timestamp_annotation", "anno_rationale"))

      dt <- DT::datatable(
        my_df,
        options = list(scrollX = TRUE)
      )
    })

    output$anno_plot <- renderPlot({
      plotUMAP(rv$anno_sce,
               colour_by = input$cd_anno)
    })

    output$overview_anno <- renderPrint({
      table(rv$anno_sce[["new_anno"]])
    })

    output$overview_log <- renderPrint({
      list(
        head(rv$log_entries),
        dim(rv$log_entries)
      )
    })

    output$log_as_metadata <- renderPrint({
      md <- metadata(rv$anno_sce)
      md_labelle <- md[["labelle_log"]]

      md_labelle
    })

    output$ui_anno_analytics <- renderUI({
      tagList(
        selectInput("start_annotation", "select the annotation (left)",
                    choices = coldata_names,
                    selected = "annotation_20220215"),
        selectInput("end_annotation", "select the annotation (right)",
                    choices = coldata_names,
                    selected = "new_anno"),
        selectInput("color_annotation", "select the annotation (color fill)",
                    choices = coldata_names,
                    selected = "annotation_20220215"),
        plotOutput("anno_alluvial"),
        plotOutput("anno_chord"),
        plotOutput("anno_donut") # make interactive?
      )
    })

    output$anno_alluvial <- renderPlot({
      annotation_alluvial(
        sce = rv$anno_sce,
        colData1 = input$start_annotation,
        colData2 = input$end_annotation,
        color_by = input$color_annotation)
    })

    output$anno_chord <- renderPlot({
      annotation_chord(
        sce = rv$anno_sce,
        colData1 = input$start_annotation,
        colData2 = input$end_annotation,
        color_by = input$color_annotation)
    })

    output$anno_donut <- renderPlot({
      annotation_donut(rv$anno_sce, input$end_annotation)
    })

    output$ui_ideas_collection <- renderUI({
      tagList(
        h4("Interface in some way to celltypist?"),
        h4("Is there a way to communicate properly with cellxgene?"),
        h4("Can we think of a clever way to resolve conflicts in annotation?"),
        h4("Some interfacing to the API/content of the Cell Ontology project?"),
        h4("Some form of reporting a small-mid set of summary representations for the annotated features?"),
        h4("Some form of celltypepedia, especially if retrieved from other sources where the infos are in?")
      )
    })

    output$ui_cl_matching <- renderUI({
      tagList(
        DT::dataTableOutput("tbl_cl_matching")
      )
    })

    output$tbl_cl_matching <- DT::renderDataTable({
      originally_provided <- sort(unique(sce[[input$start_annotation]]))
      cl_matched <- data.frame(
        provided_names = originally_provided
      )

      library("ontoProc")
      library("stringdist")
      cl <- getCellOnto()

      # store and update the reactive!
      all_cl_terms <- cl$name
      all_cl_ids <- cl$id

      sorted_cl_terms <- sort(all_cl_terms)

      # IT IS IMPERFECT, I KNOW. but still returns some hits with very low effort...
      cl_matched$matched_names <-
        all_cl_terms[amatch(originally_provided, all_cl_terms, maxDist = 7)]

      # idea: could return the top 3,4,5 hits for each level
      DT::datatable(cl_matched,
                    options = list(scrollX = TRUE))

    })


    output$ui_celltypedia <- renderUI({
      tagList(
        selectInput(
          inputId = "celltypedia_select",
          label = "Select a cell type to retrieve more info from",
          choices = cell_dictionary
        ),
        htmlOutput("selected_celltype")
      )
    })

    output$selected_celltype <- renderUI({
      celltype_2_html(input$celltypedia_select)
    })

    observeEvent(input$tour_labelle, {
      tour <- read.delim(system.file("extdata", "tour_welcome.txt", package = "labelle"),
                         sep = ";", stringsAsFactors = FALSE,
                         row.names = NULL, quote = ""
      )
      rintrojs::introjs(session, options = list(steps = tour))
    })

    observeEvent(input$btn_label, {
      # read in which cells to rename
      # cells_to_edit <- input$which_cells
      cells_to_edit <- rv$cell_selected
      # To do: take it from a

      cells_valid <- intersect(cells_to_edit, colnames(rv$anno_sce))
      cells_id <- match(cells_valid, colnames(rv$anno_sce))


      # read in which label to apply
      label_to_apply <- input$id_label
      if (length(label_to_apply) > 1) {
        label_to_apply <- paste0(label_to_apply, collapse = "|")
      }
      # to do: expand what to do - can be it is a combi of more

      showNotification(
        ui = sprintf("Applying label %s to the %d selected cells...",
                     label_to_apply,
                     length(cells_id)),
        type = "message"
      )
      # apply the label to that subset
      right_now <- Sys.time()
      colData(rv$anno_sce)[["new_anno"]][cells_id] <- label_to_apply
      # add a timestamp while we are at it
      colData(rv$anno_sce)[["timestamp_annotation"]][cells_id] <- right_now
      # storing also the rationale behind the annotation given
      colData(rv$anno_sce)[["anno_rationale"]][cells_id] <- input$anno_rationale

      # also, logging this operation in the logging df
      rv$log_entries <- rbind(
        rv$log_entries,
        data.frame(
          timestamp = as.character(right_now),
          label_assigned = label_to_apply,
          which_cells = paste0(cells_id, collapse = ","),
          rationale = input$anno_rationale,
          user = "Fede"  # can be changed and read by the person running it (or also given as param)
        )
      )

    })

    observeEvent(input$load_cells, {
      showNotification("Loading cells from editor...", type = "default")

      cells_to_select <- editor_to_vector_sanitized(input$editor_cells)
      message("cells:")
      message(cells_to_select[1])
      cells_to_select_checked <- intersect(cells_to_select,
                                     colnames(rv$anno_sce))
      cells_to_select_not_in_sce <- setdiff(cells_to_select,
                                       colnames(rv$anno_sce))

      invalid_idx <- !cells_to_select %in% colnames(rv$anno_sce) & nzchar(cells_to_select)
      cells_to_select[invalid_idx] <- paste0("# ", cells_to_select[invalid_idx])
      newcontent_editor_cells <- paste0(cells_to_select, collapse = "\n")

      message("invalid")
      message(invalid_idx)

      if (sum(invalid_idx) > 0)
        showNotification("Some cells identifiers provided are not found in the sce object, please review the content of the editor",
                         type = "warning")


      message("checked")
      message(length(cells_to_select_checked))

      if (length(cells_to_select_checked) > 0)
        showNotification(
          ui = sprintf("Found %d valid cell identifiers in the editor", length(cells_to_select_checked)),
          type = "message"
        )

      rv$cell_selected <- cells_to_select_checked

      shinyAce::updateAceEditor(session, editorId = "editor_cells", value = newcontent_editor_cells)

    })

    observeEvent(input$btn_randomizer, {
      showNotification("Selecting some random cells from the provided sce...", type = "default")

      how_many <- c(5, 10, 20)

      cells_to_select <- sample(colnames(rv$anno_sce),
                                size = sample(how_many, 1))

      newcontent_editor_cells <- paste0(cells_to_select, collapse = "\n")
      shinyAce::updateAceEditor(session, editorId = "editor_cells", value = newcontent_editor_cells)
    })

    observeEvent(input$btn_save_dictionary, {
      showNotification("TODO - saving the dictionary to rds?...", type = "default")
    })

    observeEvent(input$btn_export_log, {
      showNotification("TODO - exporting the logfile to rds...", type = "default")
      message("things about the log")
      message(nrow(rv$log_entries), "entries")
      message("latest label set", tail(rv$log_entries$label_assigned, 1))
    })

    observeEvent(input$btn_store_log, {
      showNotification("TODO - storing the logs as metadata...", type = "default")
      metadata(rv$anno_sce)[["labelle_log"]] <- rv$log_entries
      message("stored ", nrow(rv$log_entries), " log entries as metadata")
    })

    observeEvent(input$btn_import_dictionary, {
      showNotification("TODO - importing the dictionary on the fly?...", type = "default")
      showNotification("We do get something via ontoProc!", type = "message")

      # load ontoProc and co.
      library("ontoProc")
      cl <- getCellOnto()

      # store and update the reactive!
      all_cl_terms <- cl$name
      all_cl_ids <- cl$id

      sorted_cl_terms <- sort(all_cl_terms)

      # keeping only the cl ones?
      sorted_cl_terms <- sorted_cl_terms[grepl("^CL:", names(sorted_cl_terms))]

      message("found ", length(sorted_cl_terms), " entries")
      rv$cell_dictionary <- unname(sorted_cl_terms)

      message("found ", length(rv$cell_dictionary), " entries")
      updateSelectizeInput(session,
                           inputId = "id_label",
                           choices = rv$cell_dictionary,
                           server = TRUE)
      showNotification("updated!!!", type = "message")
    })

  }

  shinyApp(labelle_ui, labelle_server)

}

