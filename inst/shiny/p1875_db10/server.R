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
# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
  
  
 
  
  loadNB <- reactive({
    progress <- shiny::Progress$new(session = session, min = 0, max = 1)
    progress$set(message = "load NBs ...")
    on.exit(progress$close())
    
    NB <- read.table(system.file("extdata/NB.tryptic", package = "NestLink"),
                     col.names = c("peptide", "ESP_Prediction"), header = TRUE)
    NB$cond <- "NB"
    NB$peptide <- (as.character(NB$peptide))
    NB$pim <- parentIonMass(NB$peptide)
    NB <- NB[nchar(NB$peptide) >2, ]
    NB$ssrc <- sapply(NB$peptide, ssrc)
    NB$peptideLength <- nchar(as.character(NB$peptide))
    NB
    })
  
  getNB <- reactive({
    NB <- loadNB()
    filter <- input$pimrange[1] < NB$pim & NB$pim < input$pimrange[2] & input$ssrcrange[1] < NB$ssrc & NB$ssrc < input$ssrcrange[2]
    
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
    
    FC <- read.table(system.file("extdata/FC.tryptic", package = "NestLink"),
                     col.names = c("peptide", "ESP_Prediction"), header = TRUE)
    
    FC$peptide <- (as.character(FC$peptide))
    idx <- grep (input$FCPattern, FC$peptide)
   
    FC$cond <- "FC"
    #FC$peptideLength <- nchar(as.character(FC$peptide))
    #FC$peptideLength <- nchar(FC$peptide)
    # FC$peptide <- (as.character(FC$peptide))
    FC$pim <- parentIonMass(FC$peptide)
    FC <- FC[nchar(FC$peptide) >2, ]
    FC$ssrc <- sapply(FC$peptide, ssrc)
    FC$peptideLength <- nchar(as.character(FC$peptide))
    FC[idx,]
    #FC
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
    
    p <- ggplot(getNB(), aes(ssrc, pim)) 
    p <- p + stat_bin2d(bins = input$bins) +  labs(title = "NB", subtitle = "plot1")
    p
    
  })
  
  getDat <- reactive({
    print(input$plotFC)
    print(input$plotuFC)
    print(input$plotNB)
    print(input$plotuNB)
    # rr <- data.frame(colnames = c("peptide", "cond",   "pim",     "ssrc") )
   #rr <- rbind(getUniqueNB(), getFC())
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
   
    p <- ggplot(getFC(), aes(ssrc, pim)) 
    p <- p + stat_bin2d(bins = input$bins) +  labs(title = "FlyCodes", subtitle = "plot2")+ labs(x = "eine achse", y= "eine andere achse")
   p
  })
 
 output$histPim <- renderPlot({
   progress <- shiny::Progress$new(session = session, min = 0, max = 1)
   progress$set(message = "render parent ion histogram ...")
   on.exit(progress$close())
   
   ggplot(getDat(), aes(x=pim, fill=cond)) +
     geom_histogram(bins=input$bins, alpha=input$alpha, position="identity")
 })
 
 output$histSsrc <- renderPlot({
   progress <- shiny::Progress$new(session = session, min = 0, max = 1)
   progress$set(message = "render ssrc histogram ...")
   on.exit(progress$close())
   
   ggplot(getDat(), aes(x=ssrc, fill=cond)) +
     geom_histogram(bins=input$bins, alpha=input$alpha, position="identity")
 })
 
 output$histESP_Prediction <- renderPlot({
   progress <- shiny::Progress$new(session = session, min = 0, max = 1)
   progress$set(message = "render ESP_Prediction histogram ...")
   on.exit(progress$close())
   
   ggplot(getDat(), aes(x=ESP_Prediction, fill=cond)) +
     geom_histogram(bins=input$bins, alpha=input$alpha, position="identity")
 })
 
 output$FlyCodeTable<- DT::renderDataTable({
   getFC()
 })
 
 output$overview <- renderPrint({
   plot(table(unlist(strsplit(substr(FC$peptide, 3, 9), ""))))
   #capture.output(table(nchar(as.character(getFC()$peptide))))
 })
 output$sessionInfo <- renderPrint({
   capture.output(sessionInfo())
 })
})
