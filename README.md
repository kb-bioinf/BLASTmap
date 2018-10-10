# BLASTmap v1.0

Source code of BLASTmap Shiny web app.

BLASTmap is written in R (version 3.2.0 tested). To run BLASTmap locally, RStudio is recommended. Required R packages:
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



### Citation

If you use BLASTmap, please cite:  

Baker K., Stephen G., Strachan S., Armstrong M., Hein I. (2018) BLASTmap: A Shiny-Based Application to Visualize BLAST Results as Interactive Heat Maps and a Tool to Design Gene-Specific Baits for Bespoke Target Enrichment Sequencing. In: Ma W., Wolpert T. (eds) Plant Pathogenic Fungi and Oomycetes. Methods in Molecular Biology, vol 1848. Humana Press, New York, NY

