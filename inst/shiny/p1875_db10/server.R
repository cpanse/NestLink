#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(NestLink)

shinyServer(function(input, output, session) {

  getNBr <- reactive({
    NB <- getNB()
    filter <- input$pimrange[1] < NB$pim & NB$pim < input$pimrange[2] & input$ssrcrange[1] < NB$ssrc & NB$ssrc < input$ssrcrange[2]
    
    NB[filter,]
  })
  
  getUniqueNB <- reactive({
    uNB <- getNBr()
    uNB$cond <- 'uNB'
    uNB <- unique(uNB)
    uNB
  })
  
  getFCr <- reactive({
    FC <- getFC()
    filter <- input$pimrange[1] < FC$pim & FC$pim < input$pimrange[2] & input$ssrcrange[1] < FC$ssrc & FC$ssrc < input$ssrcrange[2]
    FC[filter,]
  })
  
  getUniqueFC <- reactive({
    uFC <- getFCr()
    uFC$cond <- 'uFC'
    uFC <- unique(uFC)
    uFC
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
   
   p <- ggplot(getDat(), aes(x=pim, fill=cond)) +
     geom_histogram(bins=input$bins, alpha=input$alpha, position="identity") +
     theme_light()
   
   
   if (input$ggplot_facet_wrap){
     p <- p + facet_wrap(~cond)
   }
   p
 })
 
 output$histSsrc <- renderPlot({
   progress <- shiny::Progress$new(session = session, min = 0, max = 1)
   progress$set(message = "render ssrc histogram ...")
   on.exit(progress$close())
   
   p <- ggplot(getDat(), aes(x=ssrc, fill=cond)) +
     geom_histogram(bins=input$bins, alpha=input$alpha, position="identity") +
     theme_light()
   
   
   if (input$ggplot_facet_wrap){
     p <- p + facet_wrap(~cond)
   }
   p
 })
 
 output$histESP_Prediction <- renderPlot({
   progress <- shiny::Progress$new(session = session, min = 0, max = 1)
   progress$set(message = "render ESP_Prediction histogram ...")
   on.exit(progress$close())
   
   p<- ggplot(getDat(), aes(x=ESP_Prediction, fill=cond)) +
     geom_histogram(bins=input$bins, alpha=input$alpha, position="identity") +
     theme_light()
   
   
   if (input$ggplot_facet_wrap){
     p <- p + facet_wrap(~cond)
   }
   p
 })
 
 output$FlyCodeTable<- DT::renderDataTable({
   getFC()
 })
 
 output$overview <- renderPrint({
   plot(table(unlist(strsplit(substr(FC$peptide, 3, 9), ""))))
 })
 
 output$sessionInfo <- renderPrint({
   capture.output(sessionInfo())
 })
 
 output$downloadData <- downloadHandler(
   filename = function() {
     paste("NestLink-FGCZ-p1875-DB10-", Sys.Date(), ".csv", sep="")
   },
   content = function(file) {
     write.csv(getDat(), file, row.names = FALSE)
   }
 )
})
