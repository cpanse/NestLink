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
       sliderInput("bins",
                   "Number of bins:",
                   min = 1,
                   max = 250,
                   value = 50)
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      list(
       plotOutput("hist2dNB"),
       plotOutput("hist2dFC"),
       plotOutput("histPim"),
       plotOutput("histSsrc"))
    )
  )
))
