#
library(shiny)
require(conmolfields)

n <- 200
# Use a fluid Bootstrap layout
fluidPage(

  # Give the page a title
  titlePanel("KRR Model"),

  # Generate a row with a sidebar
  sidebarLayout(

    # Define the sidebar with one input

    ##
    sidebarPanel(
      
      #
      selectInput("activity", "Activity File:",
                  list(`Activity` = c("activity-train.txt", "activity-train-1.txt", "activity-train-2.txt")
                       #`West Coast` = c("WA", "OR", "CA"),
                       #`Midwest` = c("MN", "WI", "IA")))
      #textOutput("result")
     )), numericInput("act_colnum", "Column number:", 1, min = 1, max = 3),
     actionButton("goButton", "Go!"),
     p("Click the button to update the value displayed in the main panel.")
     ),
   
    ##
   mainPanel(
     tags$style(type="text/css",
                ".shiny-output-error { visibility: hidden; }",
                ".shiny-output-error:before { visibility: hidden; }"),
     #tableOutput("values"),
     #textOutput("text1"),
     plotOutput("myPlot"),
     textOutput("text1")
     
     
     
     )
  )
)

