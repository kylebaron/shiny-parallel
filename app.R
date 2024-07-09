#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(mrgsolve)
library(dplyr)
library(shiny)
library(mrgsim.parallel)
library(ggplot2)
library(future)
library(future.callr)
options(future.fork.enable = TRUE)

#mod <- modlib("1005", end = 28*24, delta = 1)

mod <- mcode("shiny",'
$PLUGIN autodec, nm-vars
$PARAM THETA1 = 1, THETA2 = 20, THETA3 = 0.9
$CMT @number 2
$PK
CL = THETA1*exp(ETA(1));
V = THETA2*exp(ETA(2));
$OMEGA 1
$DES
dxdt_A1 = -THETA3 * A1;
dxdt_A2 = THETA3 * A1 - (CL/V)* A2;
$ERROR
capture IPRED = A2/V;
', end = 28*24, delta = 1)

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("Cmin distribution"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      numericInput("ids",
                   "Number of IDs; ramp this up to 3000",
                   min = 10,
                   max = 25000,
                   value = 30),
    
    numericInput("chunks", "Chunks/cores; once IDs at 3000, reduce to 3, then 2, then 1 ", 
                 min = 1, max = 8, value = 4)
    ),
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("distPlot")
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  output$distPlot <- renderPlot({
    data <- expand.ev(amt = c(100, 300, 1000, 3000), ii = 24, addl = 27, ID = 1:input$ids)
    data <- mutate(data, dose = amt)
    #plan(callr, workers = 5L)
    options(mc.cores = input$chunks)
     x <- system.time({
    out <- mc_mrgsim_d(mod, data, nchunk = input$chunks, 
                       carry_out = "dose", ) %>% as_tibble()
    })
    cmin <- filter(out, time==240)
    ggplot(data = cmin, aes(x = IPRED)) + geom_histogram() + scale_x_log10() + 
      facet_wrap(~dose, ncol = 4) + ggtitle(glue::glue("Time {signif(x[3], 3)} (seconds)"))
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
