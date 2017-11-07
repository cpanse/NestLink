#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(protViz)
library(ggplot2)

shinyServer(function(input, output, session) {
  loadNB <- reactive({
    progress <- shiny::Progress$new(session = session, min = 0, max = 1)
    progress$set(message = "load NBs ...")
    on.exit(progress$close())
    
    NB <- read.table(system.file("extdata/NB.tryptic", package = "NestLink"), col.names = "peptide")
    NB$cond <- "NB"
    NB$peptide <- (as.character(NB$peptide))
    NB$pim <- parentIonMass(NB$peptide)
    NB <- NB[nchar(NB$peptide) > 2, ]
    NB$ssrc <- sapply(NB$peptide, ssrc)
    NB$peptideLength <- nchar(as.character(NB$peptide))
    NB
    })
  
  getNB <- reactive({
    NB <- loadNB()
    filter <- (input$pimrange[1] < NB$pim 
      & NB$pim < input$pimrange[2] 
      & input$ssrcrange[1] < NB$ssrc 
      & NB$ssrc < input$ssrcrange[2])
    
    NB[filter,]
  })
  
  getUniqueNB <- reactive({
    uNB <- getNB()
    uNB$cond <- 'uNB'
    uNB <- unique(uNB)
    uNB
  })
  
  getUniqueFC <- reactive({
    uFC <- getFC()
    uFC$cond <- 'uFC'
    uFC <- unique(uFC)
    uFC
  })
  
  loadFC <- reactive({
    progress <- shiny::Progress$new(session = session, min = 0, max = 1)
    progress$set(message = "load FlyCodes ...")
    on.exit(progress$close())
    
    FC <- read.table(system.file("extdata/FC.tryptic", package = "NestLink"), col.names = "peptide")
    FC$peptide <- (as.character(FC$peptide))
    idx <- grep (input$FCPattern, FC$peptide)
   
    FC$cond <- "FC"
    FC$pim <- parentIonMass(FC$peptide)
    FC <- FC[nchar(FC$peptide) >2, ]
    FC$ssrc <- sapply(FC$peptide, ssrc)
    FC$peptideLength <- nchar(as.character(FC$peptide))
    FC[idx,]
})
    
  getFC <- reactive({
    FC<-loadFC()
    filter <- input$pimrange[1] < FC$pim & FC$pim < input$pimrange[2] & input$ssrcrange[1] < FC$ssrc & FC$ssrc < input$ssrcrange[2]
    FC[filter,]
  })
  
  output$hist2dNB <- renderPlot({
    progress <- shiny::Progress$new(session = session, min = 0, max = 1)
    progress$set(message = "render 2D histogram for NB ...")
    on.exit(progress$close())
    
  ggplot(getNB(), aes(ssrc, pim)) +
    theme_light() +
    stat_bin2d(bins = input$bins) +  
    labs(title = "NanoBody 2D histogram", 
         subtitle = "using AA sequence in description line of p1875 db10 FASTA") +
    labs(x = "hydrophobicity based on Sequence Specific Retention Calculator", 
         y= "patent ion mass")
  
  })
  
  getDat <- reactive({

    rr <- getFC()
    if (input$plotFC == FALSE){
      rr <- rr[FALSE, ]
          }
    if (input$plotuFC){
      rr <- rbind(rr, getUniqueFC())
    }
    
      if (input$plotNB){
      rr <- rbind(rr, getNB())
    }
    if (input$plotuNB){
      rr <- rbind(rr, getUniqueNB())
    }
    print(names(rr))
    rr
  })
  
 output$hist2dFC <- renderPlot({
   progress <- shiny::Progress$new(session = session, min = 0, max = 1)
   progress$set(message = "render 2D histogram for  FC ...")
   on.exit(progress$close())
   
   ggplot(getFC(), aes(x = ssrc, y= pim)) +
    theme_light() +
    stat_bin2d(bins = input$bins) + 
    labs(title = "FlyCodes 2D histogram", 
         subtitle = paste("using tryptic digested AA sequence of applied regex pattern filter", input$FCPattern, ".")) +
    labs(x = "hydrophobicity based on Sequence Specific Retention Calculator", 
         y= "patent ion mass")
 
  })
 
 output$histPim <- renderPlot({
   progress <- shiny::Progress$new(session = session, min = 0, max = 1)
   progress$set(message = "render parent ion histogram ...")
   on.exit(progress$close())
   
  p <- ggplot(getDat(), aes(x=pim, fill=cond)) +
     geom_histogram(bins=input$bins, alpha=input$alpha, position="identity") + 
     labs(x = "parent ion mass") +
     theme_light() 
     
   if (input$ggplot_facet_wrap){
     p <- p + facet_wrap(~ cond)
   }
  p
 })
 
 output$histSsrc <- renderPlot({
   progress <- shiny::Progress$new(session = session, min = 0, max = 1)
   progress$set(message = "render ssrc histogram ...")
   on.exit(progress$close())
   
   p <- ggplot(getDat(), aes(x=ssrc, fill=cond)) +
     geom_histogram(bins=input$bins, alpha=input$alpha, position="identity") +
     labs(x = "hydrophobicity based on Sequence Specific Retention Calculator") +
     theme_light() 
   
   if (input$ggplot_facet_wrap){
     p <- p + facet_wrap(~ cond)
   }
   
   p
 })
 
 output$FlyCodeTable<- DT::renderDataTable({
   getFC()
 })
 
 output$overview <- renderPrint({
   plot(table(unlist(strsplit(substr(getFC()$peptide, 3, 9), ""))))
 })
 
 # Downloadable csv of selected dataset ----
 output$downloadFC <- downloadHandler(
   filename = function() {
     paste("NestLink_p1875_FlyCodes", ".csv", sep = "")
   },
   content = function(file) {
     write.csv(getFC(), file, row.names = FALSE)
   }
 )
 
 output$sessionInfo <- renderPrint({
   capture.output(sessionInfo())
 })
})
