# Copyright 2016 by Katie Baker
#
# https://ics.hutton.ac.uk/blastmap/
# https://github.com/kb-bioinf/BLASTmap
#
# Distributed under the MIT licence
# Please see LICENSE.txt for terms of distribution

library(shiny)
library(d3heatmap)
library(reshape2)
library(htmlwidgets)
library(RColorBrewer)
library(gplots)
library(grDevices)

options(expressions=10000)
options(shiny.maxRequestSize=10*1024^2) 
`%then%` <- shiny:::`%OR%`


shinyServer(function(input, output, session) {

  selected_hits <- "All"
  selected_queries <- "All"
  selected_hits_pf <- "All"
  selected_queries_pf <- "All"
  file_selected_hits_pf <- NULL
  file_selected_queries_pf <- NULL
  file_selected_hits <- NULL
  file_selected_queries <- NULL
  input_blastn <- data.frame()
  filtered_blastn <- data.frame()
  plot_blastn <- data.frame()
  matrix_blastn <- as.matrix(data.frame())
  
  values <- reactiveValues()
  values$blastnFile <- NULL
  
  ##################################################################################################################
  ##################################################################################################################
  #
  # render UI
  #
  # tabs:
  # A. import
  # B. interactive heat map
  # C. table
  # D. export and static heat map
  # E. help
  
  
  #########################################################
  #########################################################
  # A. import
  
  #########################################################
  # 1. renderImport - ui.R file
  output$renderImport <- renderUI({
    renderImport_tc <- tryCatch({
      
      values$go <- paste0("go",runif(1))
      
      fluidRow(
        column(10,
               h3("Import data"),
               p("Default BLAST output is assumed with the order of the columns as follows: Query name, hit name, percentage identity, alignment length, mismatches, gaps, query start, query end, hit start, hit end, e-value, bitscore."),
               fileInput('blastnFile','Choose BLAST output'),
               h4("Test data"),
               p("Nucleotide BLAST of cloned and predicted Solanaceae R genes against self."),
               checkboxInput('testData','Use test data?'),
               br(),
               prefilterOpts()
        )
      )
      
      
    }, error = function(err) {
      msg <- paste0(err)
      print(msg)
    }, warning = function(warning) {
      msg <- paste0("WARNING: ",warning)
      print(msg)
    })
  })
  
  prefilterOpts <- reactive({
    uiOutput('renderPrefilter')
  })
  
  #########################################################
  # 2. renderPrefilter - ui.R file
  output$renderPrefilter <- renderUI({
    renderPrefilter_tc <- tryCatch({
      
      if (is.null(values$blastnFile))
        return(NULL) 
      
      if (nrow(input_blastn) == 0)
        return(NULL) 
      
      if (is.null(values$mat_size))
        return(NULL)
      
      if (values$mat_size <= 50000) {
        values$prefilter <- FALSE
        mainPanel(width=14,
                  h4("Matrix size prediction"),
                  p(paste0("Predicted matrix size: ",as.character(values$mat_size))),
                  p("Matrix size OK")
        )
      } else {
        values$prefilter <- TRUE
        mainPanel(width=14,
                  h4("Matrix size prediction"),
                  p(paste0("Predicted matrix size: ",as.character(values$mat_size))),
                  p(paste0("Matrix size not OK - import smaller data set or filter current data set to \u2264 50,000 matrix cells")), 
                  br(),
                  checkboxInput('filterOptions','Filter current data set?',value=FALSE),
                  br(),
                  uiOutput('renderFilterOptions')
        )
      }
    }, error = function(err) {
      msg <- paste0(err)
      print(msg)
    }, warning = function(warning) {
      msg <- paste0("WARNING: ",warning)
      print(msg)
    })
  })
  
  
  #########################################################
  # 3. renderFilterOptions - from output$renderPrefilter
  output$renderFilterOptions <- renderUI({
    renderFilterOptions_tc <- tryCatch({
      if (is.null(values$blastnFile))
        return(NULL) 
      
      values$filterOptions <- input$filterOptions
      
      if (is.null(values$filterOptions))
        return(NULL) 
      
      if (!values$filterOptions)
        return(NULL)
      
      if (values$filterOptions) {
        uiOutput('renderOptions')
      }
      
    }, error = function(err) {
      msg <- paste0(err)
      print(msg)
    }, warning = function(warning) {
      msg <- paste0("WARNING: ",warning)
      print(msg)
    })
    
  })
  
  
  output$sidebarOptions <- renderUI({
    if (is.null(values$blastnFile))
      return(NULL) 
    
    if (is.null(values$filterOptions) || !values$filterOptions)
      return(NULL) 
    
    if (values$mat_size < 50000) {
      return(NULL) 
      
    } else {
      ############################
      # update slider inputs with minimum and maximum values
      min_hitRange <- min(input_blastn$hit_start,input_blastn$hit_end) 
      max_hitRange <- max(input_blastn$hit_start,input_blastn$hit_end) 
      min_queryRange <- min(input_blastn$query_start,input_blastn$query_end) 
      max_queryRange <- max(input_blastn$query_start,input_blastn$query_end) 
      
      mainPanel(width = 14,
                style = "overflow-y:scroll; max-height: 700px",
                h3("Filter options"),
                selectInput('query_select_pf', 'Select query names', choices = c("All",sort(unique(input_blastn$query))), selected = "", multiple = TRUE,width="90%"),
                uiOutput('query_file_pf_resettable'),
                selectInput('hit_select_pf', 'Select hit names', choices = c("All",sort(unique(input_blastn$hit))), selected = "", multiple = TRUE,width="90%"),
                uiOutput('hit_file_pf_resettable'),
                actionButton('query_file_pf_rm','Clear queries',width="50%"),
                actionButton('hit_file_pf_rm','Clear hits',width="40%"),
                br(),br(),
                sliderInput('queryRange_pf', 'Filter on query location', min=min_queryRange, max=max_queryRange, value=c(min_queryRange,max_queryRange),step=1,width="90%",round=TRUE),
                sliderInput('hitRange_pf', 'Filter on hit location', min=min_hitRange, max=max_hitRange, value=c(min_hitRange,max_hitRange),step=1,width="90%",round=TRUE),
                sliderInput('id_filter_pf', 'Filter on percentage identity', min=min(input_blastn$id), max=max(input_blastn$id), value=c(min(input_blastn$id),max(input_blastn$id)),step=0.1,width="90%"),
                sliderInput('alnRange_pf', 'Filter on alignment lengths', min=min(input_blastn$aln.len), max=max(input_blastn$aln.len), value=c(min(input_blastn$aln.len),max(input_blastn$aln.len)),step=1,width="90%",round=TRUE),
                sliderInput('bitscore_filter_pf', 'Filter on bitscore', min=min(input_blastn$bitscore), max=max(input_blastn$bitscore), value=c(min(input_blastn$bitscore),max(input_blastn$bitscore)),step=1,width="90%",round=TRUE),
                sliderInput('mismatch_filter_pf', 'Filter on mismatches', min=min(input_blastn$mismatch), max=max(input_blastn$mismatch), value=c(min(input_blastn$mismatch),max(input_blastn$mismatch)),step=1,width="90%",round=TRUE),
                sliderInput('gaps_filter_pf', 'Filter on gaps', min=min(input_blastn$gaps), max=max(input_blastn$gaps), value=c(min(input_blastn$gaps),max(input_blastn$gaps)),step=1,width="90%",round=TRUE),
                sliderInput('evalue_filter_pf', 'Maximum e-value', min=min(input_blastn$evalue), max=max(input_blastn$evalue), value=max(input_blastn$evalue),step=max(input_blastn$evalue)/100,width="90%"),
                br()
      )
    }
  })
  
  #########################################################
  # 4. renderOptions - from output$renderFilterOptions
  output$renderOptions <- renderUI({
    renderOptions_tc <- tryCatch({
      if (is.null(values$blastnFile))
        return(NULL) 
      
      if (is.null(values$filterOptions))
        return(NULL) 
      
      mainPanel(width=12,
                h4("New matrix summary"),
                br(),
                htmlOutput('matrixConfirmation'),
                tags$head(tags$style("#matrixConfirmation{color: red;
                                     font-size: 18px;
                                     font-style: bold;
    }"
               )
                ),
               br(),
               uiOutput('matrixSummary'),
               br()
                )
      
    }, error = function(err) {
      msg <- paste0(err)
      print(msg)
    }, warning = function(warning) {
      msg <- paste0("WARNING: ",warning)
      print(msg)
    })
    
  })
  
  
  #########################################################
  # 5. resettable inputs - from output$renderOptions
  output$query_file_pf_resettable <- renderUI({
    input$query_file_pf_rm
    
    fileInput('query_file_pf','Upload file with query names')
  })
  
  output$hit_file_pf_resettable <- renderUI({
    input$hit_file_pf_rm
    
    fileInput('hit_file_pf','Upload file with hit names')
  })
  
  
  #########################################################
  # 6. matrix summaries - from output$renderOptions
  output$matrixSummary <- renderUI({
    if (is.null(values$blastnFile))
      return(NULL) 
    
    if (is.null(values$filterOptions))
      return(NULL) 
    
    dPrefilter()
    tableOutput('matrixSummaryTable')
  })
  
  output$matrixSummaryTable <- renderTable({
    matrixSummary_tc <- tryCatch({
      if (is.null(values$blastnFile))
        return(NULL) 
      
      if (is.null(values$filterOptions))
        return(NULL) 
      
      if (nrow(filtered_blastn) == 0) 
        return(NULL)
      
      queries <- toString(unique(filtered_blastn$query))
      num_queries <- length(unique(filtered_blastn$query))
      hits <- toString(unique(filtered_blastn$hit))
      num_hits <- length(unique(filtered_blastn$hit))
      
      if (num_queries > 1000)
        output_queries <- ">1000 queries"
      else {
        output_queries <- queries
      }
      
      if (num_hits > 1000)
        output_hits <- ">1000 hits"
      else {
        output_hits <- hits
      }
      
      table_summary <- data.frame(c(length(filtered_blastn$query),values$new_mat_size,num_queries,num_hits),c("Total number of BLAST alignments","Matrix size calculated from: query sequences x hit sequences",output_queries,output_hits),row.names=c("BLAST alignments","Matrix size","Query sequences","Hit sequences"))
      names(table_summary) <- c("Number","Information/Sequences")
      
      return(table_summary)
      
    }, error = function(err) {
      msg <- paste0(err)
      print(msg)
    }, warning = function(warning) {
      msg <- paste0("WARNING: ",warning)
      print(msg)
    })
  })
  
  
  output$matrixConfirmation <- renderText({
    matrixConfirmation_tc <- tryCatch({
      if (is.null(values$blastnFile))
        return(NULL) 
      
      if (is.null(values$filterOptions))
        return(NULL) 
      
      if (nrow(filtered_blastn) == 0) 
        return(NULL)
      
      if (values$new_mat_size < 4) {
        "Matrix size too small"
      } else if (values$new_mat_size >= 50000) {
        "Matrix size too big"
      } else {
        "Matrix size OK"
      }
    }, error = function(err) {
      msg <- paste0(err)
      print(msg)
    }, warning = function(warning) {
      msg <- paste0("WARNING: ",warning)
      print(msg)
    })
  })
  
  
  
  #########################################################
  #########################################################
  # B. interactive heat map
  
  #########################################################
  # 1. heatMapOptions - ui.R file
  output$heatMapOptions <- renderUI({
    heatMapOptions_tc <- tryCatch({    
      
      if (is.null(values$blastnFile))
        return(NULL) 
      
      if (nrow(input_blastn) == 0)
        return(NULL)
      
      dPrefilter()
      
      if (values$new_mat_size > 50000) {
        return(NULL)
        
      } else {
        
        selected_hits <<- "All"
        selected_queries <<- "All"
        
        ############################
        # update slider inputs with minimum and maximum values
        min_hitRange <- min(filtered_blastn$hit_start,filtered_blastn$hit_end) 
        max_hitRange <- max(filtered_blastn$hit_start,filtered_blastn$hit_end) 
        min_queryRange <- min(filtered_blastn$query_start,filtered_blastn$query_end) 
        max_queryRange <- max(filtered_blastn$query_start,filtered_blastn$query_end) 
        
        mainPanel(width=14,style = "overflow-y:scroll; max-height: 400px",
                  h3("Heat map options"),
                  checkboxInput('bestBlast','Best BLAST only?', value = FALSE),
                  selectInput('display', 'BLAST values to visualise:', c("Percentage identity"="id", "Bit score"="bitscore", "Alignment length"="aln.len", "Query start"="query_start", "Query end"="query_end", "Hit start"="hit_start", "Hit end"="hit_end"), 'id',width="90%"),
                  sliderInput('queryRange', 'Filter on query location', min=min_queryRange, max=max_queryRange, value=c(min_queryRange,max_queryRange),step=1,width="90%",round=TRUE),
                  sliderInput('hitRange', 'Filter on hit location', min=min_hitRange, max=max_hitRange, value=c(min_hitRange,max_hitRange),step=1,width="90%",round=TRUE),
                  sliderInput('id_filter', 'Filter on percentage identity', min=min(filtered_blastn$id), max=max(filtered_blastn$id), value=c(min(filtered_blastn$id),max(filtered_blastn$id)),step=0.1,width="90%"),
                  sliderInput('alnRange', 'Filter on alignment lengths', min=min(filtered_blastn$aln.len), max=max(filtered_blastn$aln.len), value=c(min(filtered_blastn$aln.len),max(filtered_blastn$aln.len)),step=1,width="90%",round=TRUE),
                  sliderInput('bitscore_filter', 'Filter on bitscore', min=min(filtered_blastn$bitscore), max=max(filtered_blastn$bitscore), value=c(min(filtered_blastn$bitscore),max(filtered_blastn$bitscore)),step=1,width="90%",round=TRUE),
                  sliderInput('mismatch_filter', 'Filter on mismatches', min=min(filtered_blastn$mismatch), max=max(filtered_blastn$mismatch), value=c(min(filtered_blastn$mismatch),max(filtered_blastn$mismatch)),step=1,width="90%",round=TRUE),
                  sliderInput('gaps_filter', 'Filter on gaps', min=min(filtered_blastn$gaps), max=max(filtered_blastn$gaps), value=c(min(filtered_blastn$gaps),max(filtered_blastn$gaps)),step=1,width="90%",round=TRUE),
                  sliderInput('evalue_filter', 'Maximum e-value', min=min(filtered_blastn$evalue), max=max(filtered_blastn$evalue), value=max(filtered_blastn$evalue),step=max(filtered_blastn$evalue)/100,width="90%"),
                  selectInput('query_select', 'Select query names', choices = c("All",sort(unique(filtered_blastn$query))), selected = "", multiple = TRUE,width="90%"),
                  uiOutput('query_file_resettable'),
                  selectInput('hit_select', 'Select hit names', choices = c("All",sort(unique(filtered_blastn$hit))), selected = "", multiple = TRUE,width="90%"),
                  uiOutput('hit_file_resettable'),
                  actionButton('query_file_rm','Clear queries',width="50%"),actionButton('hit_file_rm','Clear hits',width="40%")
        )
      }
    }, error = function(err) {
      msg <- paste0("output$heatMapOptions: ",err)
      print(msg)
    }, warning = function(warning) {
      msg <- paste0("output$heatMapOptions: ",warning)
      print(msg)
    })
  })
  
  
  output$displayOptions <- renderUI({
    
    if (is.null(values$blastnFile))
      return(NULL) 
    
    if (nrow(input_blastn) == 0)
      return(NULL)
    
    if (values$new_mat_size > 50000) {
      return(NULL)
      
    } else {
      mainPanel(width=14,style = "overflow-y:scroll; max-height: 300px",
                h3("Display options"),
                selectInput('dendrogram', 'Dendrogram', c(None='none', Row='row', Column='column', Both='both'), 'both',width="90%"),
                sliderInput('plotHeight','Plot height (px)',min=400,max=3000,value=800,step=100,width="90%",round=TRUE),
                sliderInput('width','Margin width (px)',min=60,max=400,value=180,step=20,width="90%",round=TRUE),
                sliderInput('height','Margin height (px)',min=60,max=400,value=120,step=20,width="90%",round=TRUE),
                sliderInput('fontSize','Font size (pt)',min=2,max=40,value=10,step=2,width="90%",round=TRUE)
      )
      
    }
    
  })
  
  
  #########################################################
  # 2. resettable inputs - from output$heatmapOptions
  output$query_file_resettable <- renderUI({
    input$query_file_rm
    
    fileInput('query_file','Upload file with query names')
  })
  
  output$hit_file_resettable <- renderUI({
    input$hit_file_rm
    
    fileInput('hit_file','Upload file with hit names')
  })
  
  
  #########################################################
  # 3. renderInteractive - ui.R file
  output$renderInteractive <- renderUI({
    if (is.null(values$blastnFile))
      return(NULL) 
    
    dPrefilter()
    
    if (values$new_mat_size > 50000) {
      return(h3('Prefilter matrix'))
    } else {      
      mainPanel(width=14,
                uiOutput('heatmapPanel')
      )
    }
    
  })
  
  
  output$heatmapPanel <- renderUI({
    plotHeatmap()
  })
  
  
  #########################################################
  # 5. plotHeatmap - from output$heatmap
  plotHeatmap <- eventReactive(values$go, {
    d3Heatmap_tc <- tryCatch({
      
      if (is.null(values$blastnFile))
        return(NULL) 
      
      height <- renderHeight()
      
      if (is.null(height))
        return(NULL)
      
      d3heatmapOutput('heatmap',width="100%",height=height)
      
    }, error = function(err) {
      msg <- paste0(err)
      print(msg)
    }, warning = function(warning) {
      msg <- paste0("WARNING: ",warning)
      print(msg)
    })
  })
  
  
  #########################################################
  # 4. heatmap - from output$renderInteractive
  output$heatmap <- renderD3heatmap({
    
    if (is.null(values$blastnFile))
      return(NULL) 
    
    dDataframe()
    
    if (length(matrix_blastn) == 0)
      return(NULL)
    
    d3heatmap(matrix_blastn,xaxis_height=marginHeight(),yaxis_width=marginWidth(),xaxis_font_size=fontSize(),yaxis_font_size=fontSize(),dendrogram=dendrogram(),height=renderHeight(),colors=rev(brewer.pal(n = 9, "Spectral")),anim_duration=0)
  })
  
  marginHeight <- eventReactive(values$go, {
    input$height
  })
  
  marginWidth <- eventReactive(values$go, {
    input$width
  })
  
  fontSize <- eventReactive(values$go, {
    input$fontSize
  })
  
  dendrogram <- eventReactive(values$go, {
    input$dendrogram
  })
  
  
  
  #########################################################
  #########################################################
  # C. table
  
  
  #########################################################
  # 1. renderTable - ui.R file
  output$renderTable <- renderUI({
    
    renderTable_tc <- tryCatch({
      
      if (is.null(values$blastnFile))
        return(NULL) 
      
      mainPanel(width=14,
                uiOutput('preTable')
      )
      
    }, error = function(err) {
      msg <- paste0(err)
      print(msg)
    }, warning = function(warning) {
      msg <- paste0("WARNING: ",warning)
      print(msg)
    })
  })
  
  
  #########################################################
  # 2. preTable - from output$renderTable
  output$preTable <- renderUI({
    if (is.null(values$blastnFile))
      return(NULL) 
    
    if (is.null(values$blastnFile))
      return(NULL) 
    
    dataTableOutput('table')
    
  })
  
  
  #########################################################
  # 4. table - from output$renderTable
  output$table <- renderDataTable({
    outputTable()
  }, options = list(
    pageLength = 10
  ))
  
  
  #########################################################
  # 5. outputTable - from out$table
  outputTable <- eventReactive(values$tableGo, {
    if (is.null(values$blastnFile))
      return(NULL) 
    
    if (nrow(plot_blastn) > 0) {
      df <- plot_blastn
    } else {
      if (nrow(filtered_blastn) > 0) {
        df <- filtered_blastn
      } else {
        df <- input_blastn
      }
    }
    
    df <- df[c(1:12)]
    df <- df[order(df$query,df$hit),]
    return(df)
  })
  
  
  
  #########################################################
  #########################################################
  # D. export and static heat map
  
  
  #########################################################
  # 1. renderExport - ui.R file
  output$renderExport <- renderUI({
    if (is.null(values$blastnFile))
      return(NULL) 
    
    if (values$new_mat_size >= 50000) {
      mainPanel(width=14,
                h3("Download options"),
                downloadButton('downloadData','Download data in TSV or CSV format'),
                br(),
                br(),
                p("Prefilter matrix before accessing heat map download options.")
      )
      
    } else {
      
      mainPanel(width=14,
                h3("Download options"),
                downloadButton('downloadData','Download data in TSV or CSV format'),
                br(),
                br(),
                uiOutput('renderPreview'),
                br()
      )
      
    }
    
  })
  
  
  #########################################################
  # 2. renderPreview - from output$renderExport
  output$renderPreview <- renderUI({
    if (is.null(values$blastnFile))
      return(NULL) 
    
    if (is.null(values$bestBlast)) {
      textOutput('renderFail')
    } else {
      uiOutput('renderExportOptions')
    }
  })
  
  
  #########################################################
  # 3. renderFail - from output$renderPreview
  output$renderFail <- renderText({
    if (is.null(values$blastnFile))
      return(NULL) 
    
    "Navigate to interactive heat map before previewing static heat map and accessing heat map download options."
  })
  
  
  #########################################################
  # 4. renderExportOptions - from output$renderPreview
  output$renderExportOptions <- renderUI({
    if (is.null(values$blastnFile))
      return(NULL) 
    
    mainPanel(width=14,
              downloadButton('downloadHeatmap','Download interactive heat map in HTML format'),
              br(),
              br(),
              selectInput("format", "File type for static heat map", choices=c("jpeg","png","tiff"), selected = "png", multiple = FALSE, selectize = TRUE),
              downloadButton('downloadHeatmap2','Download static heat map'),
              br(),
              br(),
              br(),
              h3("Static heat map preview"),
              br(),
              imageOutput('heatmap2')
    )
  })
  
  
  #########################################################
  # 5. heatmap2 - from output$renderExportOptions
  output$heatmap2 <- renderImage({
    if (is.null(values$blastnFile))
      return(NULL) 
    
    renderHeatmap()
  },deleteFile=TRUE)
  
  
  
  #########################################################
  # 6. renderHeatmap - from output$heatmap2 
  renderHeatmap <- reactive({
    if (is.null(values$blastnFile))
      return(NULL)
    
    file <- tempfile(fileext=".png")
    png(file, height=renderHeightOut(), width=renderWidth(), res=input$res, units="in")
    par(mai=c(input$height/input$res,input$width/input$res,input$height/input$res,input$width/input$res))
    par(pin=c(renderWidth()-input$width/(input$res/2),renderHeightOut()-input$height/(input$res/2)))
    colPalette <- colorRampPalette(rev(brewer.pal(n = 9, "Spectral")))
    heatmap_plot <- tryCatch({
      par(cex.lab=input$titleSize_static)
      heatmap.2(matrix_blastn,dendrogram=input$dendrogram_static,col=colPalette(256), sepcolor = "white", sepwidth = c(0.0001,0.0001),colsep=1:ncol(matrix_blastn),rowsep=1:nrow(matrix_blastn),trace="none",margins=c(input$height_static,input$width_static),cexRow=input$fontSize_static/10,cexCol=input$fontSize_static/10,density.info = "none",key.xlab=paste("BLAST attribute:",input$display), srtCol=45, key.title="",lwid=c(input$keysize_width,10),lhei=c(input$keysize_height,10))
      
    }, error = function(err) {
      msg <- paste0("MY_ERROR: ",err)
      return(msg)
    })
    
    dev.off()
    
    list(src=file,width=renderWidthPx(),height=renderHeightOutPx())
  })
  
  
  
  #########################################################
  #########################################################
  # E. help
  
  
  #########################################################
  # 1. renderHelp - from ui.R
  output$renderHelp <- renderUI({
    mainPanel(width=14,
              h3("General help", id="generalhelp"),
              p("BLASTmap is a Shiny web app designed to visualise BLAST data as interactive heat maps. The only input required is a BLAST data set (additional details given below). The heat map can be navigated and filtered to display query-hit pairs of interest. There are a number of export options: data table, interactive heat map and static heat map. For optimum BLASTmap display and performance, Chrome, Vivaldi and Opera are the recommended internet browsers. Please send any queries, issues or requests to", a("Katie Baker", href="mailto:katie.baker@hutton.ac.uk"),"."),
              hr(),
              
              h3("Import help", id="importhelp"),
              p("Data is imported in the 'Import' tab. The input required is tab- or column-delimited BLAST output in the form: Query name, hit name, percentage identity, alignment length, mismatches, gaps, query start, query end, hit start, hit end, e-value, bitscore. Import data using the 'Choose File' button or try out BLASTmap using test data (select the test data checkbox)."),
              #br(),
              p("Once uploaded, the number of matrix cells in the heat map is calculated by multiplying the number of query sequences by the number of hit sequences. If the data set meets the matrix size limitation (50,000 cells), a heat map may be generated with the data. However, if the data exceeds this, you have the option to filter the current data set. The data may be filtered on query/hit names by selecting names from a dropdown box and/or importing files with a list of query/hit names (one name per line). In addition to filtering by query/hit name, BLAST output attributes may be used to filter the dataset by using sliders to select thresholds, e.g. only visualising BLAST pairs which have fewer than 3 mismatches. As filtering options are applied, the number of query and hit sequences, and therefore the matrix size, is updated and once the matrix size decreases below 50,000 a heat map can be generated."),
              hr(),
              
              
              h3("Heat map help", id="heatmaphelp"),
              p("The heat map is generated in the 'Interactive heat map' tab and is rendered after clicking 'Plot heat map', each time a filtering parameter is applied this button must be pressed. The queries are on the y-axis and hits are on the x-axis. The heat map is interactive, so can be zoomed into by drawing a box around the cells using the mouse. Selecting a hit or query highlights that row or column. Mousing over each cell gives the query and hit name, as well as the value of the BLAST attribute being visualised."),
              h4("Heat map options"),
              p("You can select from a dropdown box which BLAST attribute to visualise, the default is percentage identity and the other options are bit score, alignment length, query start, query end, hit start and hit end. There is an option to visualise only best BLAST pairs, the default is to visualise all BLAST hits for each query. Each BLAST attribute can be filtered on, using the sliders to select ranges. Only BLAST pairs within these ranges will be visualised, for example, if you are only interested in the longest alignment lengths the slider can be used to exclude the shortest query-hit alignments. Query and hit sequences can also be selected by name, either using dropdown box or uploading a files with a list of query/hit names (one name per line). To remove the name filtering, just click the 'Clear queries' or 'Clear hits' button. When filtering on query/hit names, at least two queries/hits must be selected otherwise the heat map cannot be generated."),
              h4("Display options"),
              p("By default the data set is clustered and dendrograms accompany the heat map for both rows and columns. The dendrogram display can be changed by using a dropdown box to select only rows, only columns, no dendrogram or both (default). The plot height can be altering using a slider, the plot width cannot be changed manually and automatically extends to the width of the window. The query/hit label font size can be increased or decreased using a slider. When changing the font size, the query/hit names may extend beyond the plot boundaries. The margin width/height sliders can be used to adjust the dimensions of the margins to fully display the query/hit names."),
              hr(),
              
              
              h3("Table help", id="tablehelp"),
              p("The 'BLAST output' tab gives the BLAST data of all the query-hit pairs currently being visualised in the interactive heat map. If an interactive heat map has not been generated yet, for instance if the data set is larger than 50,000 cells and requires pre-filtering, all of the BLAST data is displayed in this tab. The data can be sorted by each column and filtered using the text boxes at the bottom of the table or by using the search box in the top right corner."),
              hr(),
              
              
              h3("Export help", id="exporthelp"),
              p("Data can be exported in the 'Export' tab. Firstly, the file name can be set using the text input box at the top of the left sidebar, this file name is applied to all of the export options. The data table (which can be explored in the 'BLAST output' tab) can be exported in tab-delimited (TSV) or comma-delimited (CSV) format. The format is selected using the radio buttons on the left sidebar under 'Download options for data table'."),
              p("The heat map can be exported in two forms. As an interactive heat map in HTML format which appears exactly as in the 'Interactive heat map' tab or it can also be exported as a static heat map (JPEG, PNG or TIFF format) and appears exactly as given in the preview. The static heat map also has additional output options which changes the heat map appearance. The resolution and the plot dimensions can be altered using sliders. The dendrogram can be set to both axes, only rows, only columns or neither axes. The font sizes of the margins and the key can be changed, as can the margin widths of the plot. Finally, the dendrogram and key dimensions can be changed using sliders."),
              hr(),
              
              
              h3("Links", id="links"),
              p(a("BLASTmap on GitHub", href="https://www.github.com/katie-baker/BLASTmap")),
              p(a("The James Hutton Institute ICS homepage", href="https://ics.hutton.ac.uk")),
              p(a("The James Hutton Institute ICS twitter", href="https://twitter.com/HuttonICS")),
              p(a("The James Hutton Institute ICS GitHub", href="https://github.com/HuttonICS")),
              p(a("Katie Baker's twitter", href="https://twitter.com/kb_bioinf")),
              hr()
    )
  })
  
  
  
  ##################################################################################################################
  ##################################################################################################################
  #
  # reactive functions
  #
  # create dataframes
  # dParse()
  # dPrefilter()
  # dDataframe()
  
  
  #########################################################
  # import data
  dParse <- reactive({
    
    if (input$testData) {
      
      blastn <- read.table("test.blastn", header=FALSE, sep="\t", stringsAsFactors=FALSE, comment.char="")
      
    } else {
        
      if (is.null(values$blastnFile))
        return(NULL)
      
      if (input$sep == "Tab") {
        blastn <- read.table(input$blastnFile$datapath, header=FALSE, sep="\t", stringsAsFactors=FALSE, comment.char="")
      } else if (input$sep == "Comma") {
        blastn <- read.table(input$blastnFile$datapath, header=FALSE, sep=",", stringsAsFactors=FALSE, comment.char="") 
      }
    
    }
    
    ############################
    # validate input
    validate(
      need(ncol(blastn) == 12, "BLAST file must be standard tab- or comma-delimited format with 12 columns corresponding to: Query name, hit name, percentage identity, alignment length, mismatches, gaps, query start, query end, hit start, hit end, e-value, bit score.")
    )
    
    ############################
    # rename columns
    names(blastn) <- c("query","hit","id","aln.len","mismatch","gaps","query_start","query_end","hit_start","hit_end","evalue","bitscore")
    
    validate(
      need(is.numeric(blastn$id), "Column 3 should be numeric and correspond to percentage identity."),
      need(is.numeric(blastn$aln.len), "Column 4 should be numeric and correspond to alignment length."),
      need(is.numeric(blastn$mismatch), "Column 5 should be numeric and correspond to mismatches."),
      need(is.numeric(blastn$gaps), "Column 6 should be numeric and correspond to gaps."),
      need(is.numeric(blastn$query_start), "Column 7 should be numeric and correspond to query start."),
      need(is.numeric(blastn$query_end), "Column 8 should be numeric and correspond to query end."),
      need(is.numeric(blastn$hit_start), "Column 9 should be numeric and correspond to hit start."),
      need(is.numeric(blastn$hit_end), "Column 10 should be numeric and correspond to hit end."),
      need(is.numeric(blastn$evalue), "Column 11 should be numeric and correspond to e-value."),
      need(is.numeric(blastn$bitscore), "Column 12 should be numeric and correspond to bit score.")
    )
    
    blastn$pair <- paste(blastn$query,blastn$hit,sep="-")
    blastn <- blastn[order(blastn$pair,-blastn$bitscore),]
    blastn <- blastn[!duplicated(blastn$pair),]
    
    ############################
    # check matrix size
    num_q <- length(unique(blastn$query))
    num_h <- length(unique(blastn$hit))
    values$mat_size <- num_q*num_h
    values$new_mat_size <- num_q*num_h
    
    input_blastn <<- blastn

    return(blastn)
    
  })
  
  
  #########################################################
  # prefilter the data based on import options
  dPrefilter <- reactive({
    if (is.null(values$blastnFile))
      return(NULL)
    
    if (nrow(input_blastn) == 0)
      return(NULL)
    
    if (is.null(values$filterOptions) || values$mat_size < 50000) {
      filtered_blastn <<- input_blastn
    } else if (values$filterOptions) {
      ############################
      # get input variables
      hit_min <- min(input$hitRange_pf)
      hit_max <- max(input$hitRange_pf)
      query_min <- min(input$queryRange_pf)
      query_max <- max(input$queryRange_pf)
      aln_min <- min(input$alnRange_pf)
      aln_max <- max(input$alnRange_pf)
      id_min <- min(input$id_filter_pf)
      id_max <- max(input$id_filter_pf)
      bitscore_min <- min(input$bitscore_filter_pf)
      bitscore_max <- max(input$bitscore_filter_pf)
      mismatch_min <- min(input$mismatch_filter_pf)
      mismatch_max <- max(input$mismatch_filter_pf)
      gaps_min <- min(input$gaps_filter_pf)
      gaps_max <- max(input$gaps_filter_pf)
      evalue_max <- max(input$evalue_filter_pf)
      selected_hits_pf <<- input$hit_select_pf
      selected_queries_pf <<- input$query_select_pf
      
      ############################
      # set up hit selector
      # check for file listing gene names first
      hit.file <- values$hit_file_pf
      if (!is.null(hit.file)) {
        hit_table <- read.table(hit.file, header=FALSE, stringsAsFactors=FALSE, comment.char="")
        file_selected_hits_pf <<- as.list(hit_table$V1)
      } else {
        
        if ( length(selected_hits_pf) == 0 ) {
          selected_hits_pf <<- unique(input_blastn$hit)
          
        } else if ("All" %in% selected_hits_pf) {
          selected_hits_pf <<- unique(input_blastn$hit)
          updateSelectInput(session, "hit_select_pf", choices = c("All",sort(unique(input_blastn$hit))), selected = "" )
          
        } else if ( length(selected_hits_pf) == 1 ) {
          
          if (selected_hits_pf == "All") {
            selected_hits_pf <<- unique(input_blastn$hit)
            updateSelectInput(session, "hit_select_pf", choices = c("All",sort(unique(input_blastn$hit))), selected = "") 
            
          } else {
            validate(
              need(input$hit_select_pf <= 2, "Select at least two hits.") 
            )
          }
          
        } else {
          selected_hits_pf <<- input$hit_select_pf
          updateSelectInput(session, "hit_select_pf", choices = c("All",sort(unique(input_blastn$hit))), selected = input$hit_select_pf) 
        }
        
      }
      
      if (!is.null(file_selected_hits_pf)) {
        if (!is.null(selected_hits_pf)) {
          selected_hits_pf <<- c(selected_hits_pf,file_selected_hits_pf)
        } else {
          selected_hits_pf <<- file_selected_hits_pf
        }
      } 
      
      
      ############################
      # set up query selector
      # check for file listing gene names first
      query.file <- values$query_file_pf
      if (!is.null(query.file)) {
        query_table <- read.table(query.file, header=FALSE, stringsAsFactors=FALSE, comment.char="")
        file_selected_queries_pf <<- as.list(query_table$V1)
      } else {
        
        if ( length(selected_queries_pf) == 0 ) {
          selected_queries_pf <<- unique(input_blastn$query)
          
        } else if ("All" %in% selected_queries_pf) {
          selected_queries_pf <<- unique(input_blastn$query)
          updateSelectInput(session, "query_select_pf", choices = c("All",sort(unique(input_blastn$query))), selected = "" )
          
        } else if ( length(selected_queries_pf) == 1 ) {
          
          if ( selected_queries_pf == "All" ) {
            selected_queries_pf <<- unique(input_blastn$query)
            updateSelectInput(session, "query_select_pf", choices = c("All",sort(unique(input_blastn$query))), selected = "") 
            
          } else {
            validate(
              need(input$query_select_pf <= 2, "Select at least two queries.") 
            )
          }
          
        } else {
          selected_queries_pf <<- input$query_select_pf
          updateSelectInput(session, "query_select_pf", choices = c("All",sort(unique(input_blastn$query))), selected = input$query_select_pf) 
          
        }
      }
      
      if (!is.null(query.file)) {
        if (!is.null(selected_queries_pf)) {
          selected_queries_pf <<- c(selected_queries_pf,file_selected_queries_pf)
        } else {
          selected_queries_pf <<- file_selected_queries_pf
        }
      }
      
      
      ############################
      # filter data frame based on input filtering
      filtered_blastn <<- input_blastn[input_blastn$query %in% selected_queries_pf
                                       & input_blastn$hit %in% selected_hits_pf
                                       & input_blastn$bitscore>=bitscore_min
                                       & input_blastn$bitscore<=bitscore_max
                                       & input_blastn$mismatch>=mismatch_min
                                       & input_blastn$mismatch<=mismatch_max
                                       & input_blastn$gaps>=gaps_min
                                       & input_blastn$gaps<=gaps_max
                                       & input_blastn$evalue<=evalue_max
                                       & input_blastn$id>=id_min
                                       & input_blastn$id<=id_max
                                       & input_blastn$hit_start>=hit_min
                                       & input_blastn$hit_end<=hit_max
                                       & input_blastn$query_start>=query_min
                                       & input_blastn$query_end<=query_max
                                       & input_blastn$aln.len>=aln_min
                                       & input_blastn$aln.len<=aln_max
                                       ,]
      
      
    } else {
      return(NULL)
    }
    
    num_q <- length(unique(filtered_blastn$query))
    num_h <- length(unique(filtered_blastn$hit))
    values$new_mat_size <- num_q*num_h
    
    return(filtered_blastn)
    
  })
  
  
  #########################################################
  # filter data for interactive heat map
  dDataframe <- eventReactive(values$go, {
    if (is.null(values$blastnFile))
      return(NULL)
    
    if (nrow(filtered_blastn) == 0)
      return(NULL)
    
    if (is.null(input$bestBlast))
      return(h3("Prefilter matrix"))
    
    ############################
    # get input variables
    hit_min <- min(input$hitRange)
    hit_max <- max(input$hitRange)
    query_min <- min(input$queryRange)
    query_max <- max(input$queryRange)
    aln_min <- min(input$alnRange)
    aln_max <- max(input$alnRange)
    id_min <- min(input$id_filter)
    id_max <- max(input$id_filter)
    bitscore_min <- min(input$bitscore_filter)
    bitscore_max <- max(input$bitscore_filter)
    mismatch_min <- min(input$mismatch_filter)
    mismatch_max <- max(input$mismatch_filter)
    gaps_min <- min(input$gaps_filter)
    gaps_max <- max(input$gaps_filter)
    evalue_max <- max(input$evalue_filter)
    selected_hits <<- input$hit_select
    selected_queries <<- input$query_select
    values$bestBlast <- input$bestBlast
    
    if (values$bestBlast) {
      heatmap_df <- filtered_blastn[order(filtered_blastn$query,-filtered_blastn$bitscore),]
      heatmap_df <- heatmap_df[!duplicated(heatmap_df$query),]
    } else {
      heatmap_df <- filtered_blastn
    }
    
    
    ############################
    # set up hit selector
    # check for file listing gene names first
    hit.file <- values$hit_file
    if (!is.null(hit.file)) {
      hit_table <- read.table(hit.file, header=FALSE, stringsAsFactors=FALSE, comment.char="")
      file_selected_hits <<- as.list(hit_table$V1)
    } else {
      
      if ( length(selected_hits) == 0 ) {
        selected_hits <<- unique(heatmap_df$hit)
        
      } else if ("All" %in% selected_hits) {
        selected_hits <<- unique(heatmap_df$hit)
        updateSelectInput(session, "hit_select", choices = c("All",sort(unique(heatmap_df$hit))), selected = "" )
        
      } else if ( length(selected_hits) == 1 ) {
        
        if ( selected_hits == "All" ) {
          selected_hits <<- unique(heatmap_df$hit)
          updateSelectInput(session, "hit_select", choices = c("All",sort(unique(heatmap_df$hit))), selected = "") 
          
        } else {
          validate(
            need(input$hit_select <= 2, "Select at least two hits.") 
          )
        }
        
      } else {
        selected_hits <<- input$hit_select
        updateSelectInput(session, "hit_select", choices = c("All",sort(unique(heatmap_df$hit))), selected = input$hit_select) 
      }
    }
    
    if (!is.null(hit.file)) {
      if (!is.null(selected_hits)) {
        selected_hits <<- c(selected_hits,file_selected_hits)
      } else {
        selected_hits <<- file_selected_hits
      }
    }
    
    ############################
    # set up query selector
    # check for file listing gene names first
    query.file <- values$query_file
    if (!is.null(query.file)) {
      query_table <- read.table(query.file, header=FALSE, stringsAsFactors=FALSE, comment.char="")
      file_selected_queries <<- as.list(query_table$V1)
    } else {
      
      if ( length(selected_queries) == 0 ) {
        selected_queries <<- unique(heatmap_df$query)
        
      } else if ("All" %in% selected_queries) {
        selected_queries <<- unique(heatmap_df$query)
        updateSelectInput(session, "query_select", choices = c("All",sort(unique(heatmap_df$query))), selected = "" )
        
      } else if ( length(selected_queries) == 1 ) {
        
        if ( selected_queries == "All" ) {
          selected_queries <<- unique(heatmap_df$query)
          updateSelectInput(session, "query_select", choices = c("All",sort(unique(heatmap_df$query))), selected = "") 
          
        } else {
          validate(
            need(input$query_select <= 2, "Select at least two queries.") 
          )
        }
        
      } else {
        selected_queries <<- input$query_select
        updateSelectInput(session, "query_select", choices = c("All",sort(unique(heatmap_df$query))), selected = input$query_select) 
      }
    }
    
    if (!is.null(query.file)) {
      if (!is.null(selected_queries)) {
        selected_queries <<- c(selected_queries,file_selected_queries)
      } else {
        selected_queries <<- file_selected_queries
      }
    }
    
    ############################
    # filter data frame based on input filtering
    plot_blastn <<- heatmap_df[heatmap_df$query %in% selected_queries 
                               & heatmap_df$hit %in% selected_hits
                               & heatmap_df$bitscore>=bitscore_min
                               & heatmap_df$bitscore<=bitscore_max
                               & heatmap_df$mismatch>=mismatch_min
                               & heatmap_df$mismatch<=mismatch_max
                               & heatmap_df$gaps>=gaps_min
                               & heatmap_df$gaps<=gaps_max
                               & heatmap_df$evalue<=evalue_max
                               & heatmap_df$id>=id_min
                               & heatmap_df$id<=id_max
                               & heatmap_df$hit_start>=hit_min
                               & heatmap_df$hit_end<=hit_max
                               & heatmap_df$query_start>=query_min
                               & heatmap_df$query_end<=query_max
                               & heatmap_df$aln.len>=aln_min 
                               & heatmap_df$aln.len<=aln_max
                               ,]
    
    ############################
    # check there are enough data points for heatmap plotting
    query_num <- length(unique(plot_blastn$query))
    hit_num <- length(unique(plot_blastn$hit))
    values$new_mat_size <- query_num*hit_num
    
    validate(
      need(query_num > 1 & hit_num > 1, "Too few hits and queries, cannot produce heat map. Reduce filtering stringency or make sure at least two query-hit pairs are selected.") %then%
        need(query_num > 1, "Need more data points to produce heat map, too few queries. Reduce filtering stringency.\n") %then%
        need(hit_num > 1, "Need more data points to produce heat map, too few hits. Reduce filtering stringency.\n")
    )
    
    heatmap_matrix <- acast(plot_blastn, query ~ hit, value.var=input$display)
    heatmap_matrix[is.na(heatmap_matrix)] <- 0
    
    # save in matrix
    matrix_blastn <<- heatmap_matrix
    
    ############################
    # return statement
    return(plot_blastn)
    
  })
  
  
  
  ##################################################################################################################
  ##################################################################################################################
  #
  # reactive functions
  #
  # image rendering parameters 
  
  
  #########################################################
  # get user inputted plot height - output
  renderHeightOut <- reactive({
    height <- (input$plotHeight_out/input$res) + (input$height/(input$res/2))
    return(height)
  })
  
  
  #########################################################
  # get user inputted plot height in pixels - output
  renderHeightOutPx <- reactive({
    height <- input$plotHeight_out+input$height
    return(height)
  })
  
  
  #########################################################
  # get user inputted plot width
  renderWidth <- reactive({
    width <- (input$plotWidth/input$res)+(input$width/(input$res/2))
    return(width)
  })
  
  
  #########################################################
  # get user inputted plot width in pixels
  renderWidthPx <- reactive({
    width <- input$plotWidth+input$width
    return(width)
  })
  
  
  #########################################################
  # get user inputted plot height
  renderHeight <- eventReactive(values$go, {
    input$plotHeight
  })
  
  
  
  ##################################################################################################################
  ##################################################################################################################
  #
  # observe events
  #
  
  
  #########################################################
  # button presses
  observeEvent(input$go, {
    values$go <- paste0("go",runif(1))
    values$tableGo <- paste0("go",runif(1))
  })
  
  observeEvent(values$go, {
    updateCheckboxInput(session, 'bestBlast', value = input$bestBlast)
    updateSelectInput(session, 'display', selected = input$display)
    updateSliderInput(session, 'queryRange', value=c(min(input$queryRange),max(input$queryRange)))
    updateSliderInput(session, 'hitRange', value=c(min(input$hitRange),max(input$hitRange)))
    updateSliderInput(session, 'id_filter', value=c(min(input$id_filter),max(input$id_filter)))
    updateSliderInput(session, 'alnRange', value=c(min(input$alnRange),max(input$alnRange)))
    updateSliderInput(session, 'bitscore_filter', value=c(min(input$bitscore_filter),max(input$bitscore_filter)))
    updateSliderInput(session, 'mismatch_filter', value=c(min(input$mismatch_filter),max(input$mismatch_filter)))
    updateSliderInput(session, 'gaps_filter', value=c(min(input$gaps_filter),max(input$gaps_filter)))
    updateSliderInput(session, 'evalue_filter', value=max(input$evalue_filter))
  })
  
  observeEvent(input$hit_file_pf_rm, {
    if (!is.null(input$hit_file_pf)) {
      file.remove(input$hit_file_pf$datapath)
    }
    values$hit_file_pf <- NULL
    updateSelectInput(session,'hit_select_pf', choices = c("All",sort(unique(input_blastn$hit))), selected = "")
  })
  
  observeEvent(input$query_file_pf_rm, {
    if (!is.null(input$query_file_pf)) {
      file.remove(input$query_file_pf$datapath)
    }
    values$query_file_pf <- NULL
    updateSelectInput(session,'query_select_pf', choices = c("All",sort(unique(input_blastn$query))), selected = "")
  })
  
  observeEvent(input$hit_file_rm, {
    if (!is.null(input$hit_file)) {
      file.remove(input$hit_file$datapath)
    }
    values$hit_file <- NULL
    updateSelectInput(session,'hit_select', choices = c("All",sort(unique(filtered_blastn$hit))), selected = "")
  })
  
  observeEvent(input$query_file_rm, {
    if (!is.null(input$query_file)) {
      file.remove(input$query_file$datapath)
    }
    values$query_file <- NULL
    updateSelectInput(session,'query_select', choices = c("All",sort(unique(filtered_blastn$query))), selected = "")
  })
  
  observeEvent(input_blastn, {
    values$tableGo <- paste0("go",runif(1))
  })
  
  observeEvent(filtered_blastn, {
    values$tableGo <- paste0("go",runif(1))
  })
  
  observeEvent(plot_blastn, {
    values$tableGo <- paste0("go",runif(1))
  })
  
  
  
  #########################################################
  # new file uploads
  observeEvent(input$blastnFile$name, {
    dParse()
    filtered_blastn <<- data.frame()
    plot_blastn <<- data.frame()
    matrix_blastn <<- as.matrix(data.frame())
    selected_hits_pf <<- "All"
    selected_queries_pf <<- "All"
    file_selected_queries_pf <<- NULL
    file_selected_hits_pf <<- NULL
    values$hit_file_pf <- NULL
    values$query_file_pf <- NULL
    values$hit_file <- NULL
    values$query_file <- NULL
    values$bestBlast <- NULL
    values$prefilter <- FALSE
    values$go <- NULL
    values$tableGo <- paste0("go",runif(1))
    values$filterOptions <- NULL
    values$blastnFile <- paste0("dataImported",runif(1))
    updateCheckboxInput(session, 'testData', value = FALSE)
  })
  
  observeEvent(input$hit_file_pf$name, {
    values$hit_file_pf <- input$hit_file_pf$datapath
    values$go <- NULL
  })
  
  observeEvent(input$query_file_pf$name, {
    values$go <- NULL
    values$query_file_pf <- input$query_file_pf$datapath
  })
  
  observeEvent(input$hit_file$name, {
    values$hit_file <- input$hit_file$datapath
  })
  
  observeEvent(input$query_file$name, {
    values$query_file <- input$query_file$datapath
  })
  
  observeEvent(input$testData, {
    dParse()
    filtered_blastn <<- data.frame()
    plot_blastn <<- data.frame()
    matrix_blastn <<- as.matrix(data.frame())
    selected_hits_pf <<- "All"
    selected_queries_pf <<- "All"
    file_selected_queries_pf <<- NULL
    file_selected_hits_pf <<- NULL
    values$hit_file_pf <- NULL
    values$query_file_pf <- NULL
    values$hit_file <- NULL
    values$query_file <- NULL
    values$bestBlast <- NULL
    values$prefilter <- FALSE
    values$go <- NULL
    values$tableGo <- paste0("go",runif(1))
    values$filterOptions <- NULL
    values$blastnFile <- paste0("dataImported",runif(1))
  })
  
  ##################################################################################################################
  ##################################################################################################################
  #
  # output functions
  #
  
  
  #########################################################
  # download statements
  output$downloadData <- downloadHandler(
    filename <- function(){
      if (input$tableFormat == "Tab") {
        filename <- paste0(input$fileName, ".tsv", sep="")
      } else {
        filename <- paste0(input$fileName, ".csv", sep="")
      }
      return(filename)
    },
    content <- function(file) {
      if (input$tableFormat == "Tab") {
        write.table(outputTable(), file, sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
      } else {
        write.table(outputTable(), file, sep=",",row.names=FALSE,col.names=FALSE,quote=FALSE)
      }
    }
  )
  
  
  output$downloadHeatmap <- downloadHandler(
    filename <- function() {
      paste0(input$fileName, ".html", sep="")
    },
    content <- function(file) {
      saveWidget(plotHeatmap(), file)
    }
  )
  
  
  output$downloadHeatmap2 <- downloadHandler(
    filename <- function() {paste0(input$fileName, ".", input$format, sep="")},
    content <- function(file) {
      #dev.off()
      if(input$format=="png") {
        png(file, height=renderHeightOut(), width=renderWidth(), res=input$res, units="in")
        par(mai=c(input$height/input$res,input$width/input$res,input$height/input$res,input$width/input$res))
        par(pin=c(renderWidth()-input$width/(input$res/2),renderHeightOut()-input$height/(input$res/2)))
        plotHeatmap2_png()
        dev.off()
      } else if(input$format=="tiff") {
        tiff(file, height=renderHeightOut(), width=renderWidth(), res=input$res, units="in", compression="lzw")
        par(mai=c(input$height/input$res,input$width/input$res,input$height/input$res,input$width/input$res))
        par(pin=c(renderWidth()-input$width/(input$res/2),renderHeightOut()-input$height/(input$res/2)))
        plotHeatmap2_tiff()
        dev.off()
      } else if(input$format=="jpeg") {
        jpeg(file, height=renderHeightOut(), width=renderWidth(), res=input$res, units="in", quality=100)
        par(mai=c(input$height/input$res,input$width/input$res,input$height/input$res,input$width/input$res))
        par(pin=c(renderWidth()-input$width/(input$res/2),renderHeightOut()-input$height/(input$res/2)))
        plotHeatmap2_jpeg()
        dev.off()
      }
    }
  )
  
  
  plotHeatmap2_png <- reactive({
    if (is.null(values$blastnFile))
      return(NULL)
    
    colPalette <- colorRampPalette(rev(brewer.pal(n = 9, "Spectral")))
    heatmap_plot <- tryCatch({
      par(cex.lab=input$titleSize_static)
      heatmap.2(matrix_blastn,dendrogram=input$dendrogram_static,col=colPalette(256), sepcolor = "white", sepwidth = c(0.0001,0.0001),colsep=1:ncol(matrix_blastn),rowsep=1:nrow(matrix_blastn),trace="none",margins=c(input$height_static,input$width_static),cexRow=input$fontSize_static/10,cexCol=input$fontSize_static/10,density.info = "none",key.xlab=paste("BLAST attribute:",input$display), srtCol=45, key.title="",lwid=c(input$keysize_width,10),lhei=c(input$keysize_height,10))
    }, error = function(err) {
      msg <- paste0("MY_ERROR: ",err)
      return(msg)
    })
  })
  
  
  plotHeatmap2_tiff <- reactive({
    if (is.null(values$blastnFile))
      return(NULL)
    
    colPalette <- colorRampPalette(rev(brewer.pal(n = 9, "Spectral")))
    heatmap_plot <- tryCatch({
      par(cex.lab=input$titleSize_static)
      heatmap.2(matrix_blastn,dendrogram=input$dendrogram_static,col=colPalette(256), sepcolor = "white", sepwidth = c(0.0001,0.0001),colsep=1:ncol(matrix_blastn),rowsep=1:nrow(matrix_blastn),trace="none",margins=c(input$height_static,input$width_static),cexRow=input$fontSize_static/10,cexCol=input$fontSize_static/10,density.info = "none",key.xlab=paste("BLAST attribute:",input$display), srtCol=45, key.title="",lwid=c(input$keysize_width,10),lhei=c(input$keysize_height,10))
    }, error = function(err) {
      msg <- paste0("MY_ERROR: ",err)
      return(msg)
    })
  })
  
  
  plotHeatmap2_jpeg <- reactive({
    if (is.null(values$blastnFile))
      return(NULL)
    
    colPalette <- colorRampPalette(rev(brewer.pal(n = 9, "Spectral")))
    heatmap_plot <- tryCatch({
      par(cex.lab=input$titleSize_static)
      heatmap.2(matrix_blastn,dendrogram=input$dendrogram_static,col=colPalette(256), sepcolor = "white", sepwidth = c(0.0001,0.0001),colsep=1:ncol(matrix_blastn),rowsep=1:nrow(matrix_blastn),trace="none",margins=c(input$height_static,input$width_static),cexRow=input$fontSize_static/10,cexCol=input$fontSize_static/10,density.info = "none",key.xlab=paste("BLAST attribute:",input$display), srtCol=45, key.title="",lwid=c(input$keysize_width,10),lhei=c(input$keysize_height,10))
    }, error = function(err) {
      msg <- paste0("MY_ERROR: ",err)
      return(msg)
    })
  })
  
  
})