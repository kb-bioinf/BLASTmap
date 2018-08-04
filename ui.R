# Copyright 2016 by Katie Baker
#
# https://ics.hutton.ac.uk/blastmap/
# https://github.com/kb-bioinf/BLASTmap
#
# Distributed under the MIT licence
# Please see LICENSE.txt for terms of distribution

library(shiny)
library(shinythemes)
library(d3heatmap)
library(reshape2)

shinyUI(fluidPage(id = "mainpage", theme = shinytheme("spacelab"),
                  titlePanel("BLASTmap"),
                  sidebarLayout(
                    sidebarPanel(width=2,
                                 tags$head(tags$style(type="text/css", "
                                                      #loadmessage {
                                                      position: fixed;
                                                      top: 0px;
                                                      left: 0px;
                                                      width: 15%;
                                                      padding: 0px 0px 0px 0px;
                                                      text-align: center;
                                                      font-weight: bold;
                                                      font-size: 100%;
                                                      color: #000000;
                                                      background-color: #8FBC8F;
                                                      z-index: 105;
                                                      }
                                                      ")),
                                 
                                 conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                  tags$div("Loading...",id="loadmessage")),
                                 
                                 conditionalPanel(condition = "input.tabs == 'Import'",
                                                  h3("Import options"),
                                                  radioButtons('sep', 'Column delimiter', choices=c("Tab","Comma"), selected = "Tab", inline = TRUE),hr(),
                                                  uiOutput('sidebarOptions')
                                 ),
                                 
                                 conditionalPanel(condition = "input.tabs == 'Interactive heat map'",
                                                  actionButton('go','Plot heat map',width='100%'),
                                                  hr(),
                                                  uiOutput('heatMapOptions'),hr(),
                                                  uiOutput('displayOptions')
                                 ), 
                                 
                                 conditionalPanel(condition = "input.tabs == 'BLAST output'",
                                                  h3("BLAST output"),
                                                  p("Change filtering parameters in the interactive heat map tab")
                                                  
                                 ), 
                                 
                                 conditionalPanel(condition = "input.tabs == 'Export'",
                                                  h3("Export options"),
                                                  textInput("fileName","File name (no extension)", value="BLASTmap-output"),
                                                  hr(),
                                                  h4("Download options for data table"),
                                                  radioButtons('tableFormat', 'Column delimiter', choices=c("Tab","Comma"), selected = "Tab", inline = TRUE),
                                                  hr(),
                                                  h4("Download options for static heatmap"),
                                                  sliderInput('res','Resolution (ppi)',min=72,max=300,value=120,step=12),
                                                  sliderInput('plotHeight_out','Plot height (px)',min=400,max=6000,value=800,step=100),
                                                  sliderInput('plotWidth','Plot width (px)',min=400,max=6000,value=1200,step=100),
                                                  selectInput('dendrogram_static', 'Dendrogram', c(None='none', Row='row', Column='column', Both='both'), 'both'),
                                                  sliderInput('fontSize_static','Axis font size (arbitary units)',min=0.1,max=10,value=5,step=0.1),
                                                  sliderInput('titleSize_static','Key font size (arbitary units)',min=0.1,max=2,value=1,step=0.1),
                                                  sliderInput('width_static','Y-axis label margin width (arbitary units)',min=0,max=20,value=10,step=1),
                                                  sliderInput('height_static','X-axis margin height (arbitary units)',min=0,max=20,value=10,step=1),
                                                  sliderInput('keysize_height','Key/dendrogram height (arbitary units)',min=0,max=5,value=2.5,step=0.1),
                                                  sliderInput('keysize_width','Key/dendrogram width (arbitary units)',min=0,max=5,value=2.5,step=0.1)
                                 ), 
                                 
                                 conditionalPanel(condition = "input.tabs == 'Help'",
                                                  h3("Help"),
                                                  a("General help", href = "#generalhelp"),
                                                  br(),
                                                  a("Import help", href = "#importhelp"),
                                                  br(),
                                                  a("Heat map help", href = "#heatmaphelp"),
                                                  br(),
                                                  a("Table help", href = "#tablehelp"),
                                                  br(),
                                                  a("Export help", href = "#exporthelp"),
                                                  br(),
                                                  a("Links", href = "#links")
                                 )
                                 ),
                    
                    
                    mainPanel(width=10,
                              tabsetPanel(id = "tabs",
                                          tabPanel("Import",
                                                   uiOutput('renderImport')
                                          ),
                                          tabPanel("Interactive heat map",
                                                   uiOutput('renderInteractive')
                                          ),
                                          tabPanel("BLAST output",
                                                   uiOutput('renderTable')
                                          ),
                                          tabPanel("Export",
                                                   uiOutput('renderExport')
                                          ),
                                          tabPanel("Help",
                                                   uiOutput('renderHelp')
                                          )
                              )
                    )
                    )
                  
                    ))