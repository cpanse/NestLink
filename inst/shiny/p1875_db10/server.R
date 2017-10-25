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
    
    NB <- read.table(system.file("extdata/NB.tryptic", package = "NestLink"), col.names = "peptide")
    NB$cond <- "NB"
    NB$peptide <- (as.character(NB$peptide))
    NB$pim <- parentIonMass(NB$peptide)
    NB <- NB[nchar(NB$peptide) >2, ]
    NB$ssrc <- sapply(NB$peptide, ssrc)
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
    
    FC <- read.table(system.file("extdata/FC.tryptic", package = "NestLink"), col.names = "peptide")
    FC$cond <- "FC"
    FC$peptide <- (as.character(FC$peptide))
    FC$pim <- parentIonMass(FC$peptide)
    FC <- FC[nchar(FC$peptide) >2, ]
    FC$ssrc <- sapply(FC$peptide, ssrc)
    FC
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
    
    rr <- rbind(getUniqueNB(), getFC())
    rr <- rbind(rr, getNB())
    rbind(rr, getUniqueFC())
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
     geom_histogram(bins=input$bins, alpha=.5, position="identity")
 })
 output$histSsrc <- renderPlot({
   progress <- shiny::Progress$new(session = session, min = 0, max = 1)
   progress$set(message = "render ssrc histogram ...")
   on.exit(progress$close())
   
   ggplot(getDat(), aes(x=ssrc, fill=cond)) +
     geom_histogram(bins=input$bins, alpha=.5, position="identity")
 })
})
