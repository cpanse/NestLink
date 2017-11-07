#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("p1875 DB10 tryptic peptides NB verus FC"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      textInput("FCPattern", "FC pattern", 
                value = "^GS[ASTNQDEVLYWGP]{7}(WR|WLTVR|WQEGGR|WLR|WQSR)$", 
                width = NULL, placeholder = NULL),
       sliderInput("bins",
                   "Number of bins:",
                   min = 1,
                   max = 250,
                   value = 50),
       
       sliderInput("pimrange", "parent ion mass - range:",
                   min = 0, max = 6000, value = c(0, 6000)),
       
       sliderInput("ssrcrange", "ssrc - range:",
                   min = 0, max = 70, value = c(0, 70)),
      checkboxInput("ggplot_facet_wrap", "facet_wrap", FALSE), 
       sliderInput("alpha", "alpha blending:",
                   min = 0, max = 1, value = 0.3),
       HTML('<hr>Histogram only filters:'),
       checkboxInput("plotFC", "Plot FC", TRUE),
       checkboxInput("plotuFC", "Plot Unique FC", TRUE),
       checkboxInput("plotNB", "Plot NB", TRUE),
       checkboxInput("plotuNB", "Plot Unique NB", TRUE),
      HTML("<hr>"),
      # Button
      downloadButton("downloadFC", "Download FCs")
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(
        tabPanel("histogram",
      list(
        plotOutput("histPim"),
        plotOutput("histSsrc"),
        uiOutput("overview")
  #     plotOutput("hist2dNB"),
  #     plotOutput("hist2dFC")
      
       )),
  tabPanel("2Dhistogram",
           list(plotOutput("hist2dNB"), plotOutput("hist2dFC"))
  ),
  
  tabPanel("FlyCode table",
           DT::dataTableOutput("FlyCodeTable")
  ),
  tabPanel("sessionInfo", verbatimTextOutput("sessionInfo"))
  )
    )
  ),HTML("<hr>source: <a href='https://github.com/cpanse/NestLink'>NestLink</a>"),
  HTML("cite: <a href='https://CRAN.R-project.org/package=protViz'>https://CRAN.R-project.org/package=protViz</a>")
  
))
