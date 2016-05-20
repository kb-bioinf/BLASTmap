# BLASTmap v1.0

Source code of [BLASTmap] Shiny web app.

BLASTmap is written in [R] (version 3.2.0 tested). To run BLASTmap locally, [RStudio] is recommended. Required R packages:
* *shiny*
* *shinythemes*
* *d3heatmap*
* *reshape2*
* *htmlwidgets*
* *RColorBrewer*
* *gplots*
* *grDevices*

The input required is tab- or column-delimited BLAST output in the form: Query name, hit name, percentage identity, alignment length, mismatches, gaps, query start, query end, hit start, hit end, e-value, bitscore.

BLASTmap is distributed under the MIT licence.

[//]: # 

   [BLASTmap]: <https://ics.hutton.ac.uk/blastmap/>
   [RStudio]: <https://www.rstudio.com>
   [R]: <https://www.r-project.org>

