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
shinyServer(function(input, output) {
  
  getNB <- reactive({
    NB <- read.table("~/data/NB.tryptic", col.names = "peptide")
    NB$cond <- "NB"
    NB$peptide <- (as.character(NB$peptide))
    NB$pim <- parentIonMass(NB$peptide)
    NB <- NB[nchar(NB$peptide) >2, ]
    NB$ssrc <- sapply(NB$peptide, ssrc)
    NB
  })
  
  getFC <- reactive({
    FC <- read.table("~/data/FC.tryptic", col.names = "peptide")
    FC$cond <- "FC"
    FC$peptide <- (as.character(FC$peptide))
    FC$pim <- parentIonMass(FC$peptide)
    FC <- FC[nchar(FC$peptide) >2, ]
    FC$ssrc <- sapply(FC$peptide, ssrc)
    FC
  })
  
  output$hist2dNB <- renderPlot({
    
    p <- ggplot(getNB(), aes(ssrc, pim)) 
    p <- p + stat_bin2d(bins = input$bins) +  labs(title = "NB", subtitle = "plot1")
    p
    
  })
  
  getDat <- reactive({
    rbind(getNB(), getFC())
  })
 output$hist2dFC <- renderPlot({
    p <- ggplot(getFC(), aes(ssrc, pim)) 
    p <- p + stat_bin2d(bins = input$bins) +  labs(title = "FC", subtitle = "plot2")+ labs(x = "eine achse", y= "eine andere achse")
   p
  })
 output$histPim <- renderPlot({
   
   
   ggplot(getDat(), aes(x=pim, fill=cond)) +
     geom_histogram(bins=input$bins, alpha=.5, position="identity")
   
  

 })
 output$histSsrc <- renderPlot({
   
   
   ggplot(getDat(), aes(x=ssrc, fill=cond)) +
     geom_histogram(bins=input$bins, alpha=.5, position="identity")
   
   
   
 })
})
