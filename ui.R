library(rCharts) 
require(rjson)

shinyUI(pageWithSidebar(
  headerPanel("Sophia task #2"),
  sidebarPanel(
    h1('Draw menu'),
    selectInput("inputY", "observable Y:", choices = c('count'='count', 'average quality'='meanQual','standard deviation'='sdQual','accuracy'='acc')),

    sliderInput('xlow', 'Minimum X',value = 0,   min = 0, max = 100, step = 1,),
    sliderInput('xhigh','Maximum X',value = 100, min = 0, max = 100, step = 1,),

    sliderInput('ylow', 'Minimum Y',value = 0,   min = 0, max = 100, step = 1,),
    sliderInput('yhigh','Maximum Y',value = 100, min = 0, max = 100, step = 1,)

#    numericInput('id1', 'Numeric input, labeled id1', 0, min = 0, max = 10, step = 1),
#    checkboxGroupInput("id2", "Checkbox",
#          c("Value 1" = "1",
#            "Value 2" = "2",
#            "Value 3" = "3")),
#    dateInput("date", "Date:")
  ),
  mainPanel(
    p('Interactive visualization engine. Analyze patterns of length=10. Details are on the ', a("github page", href="https://github.com/koskot77/sophia", target="_blank"),'.'),
    plotOutput('myHist'),
    verbatimTextOutput("events"),
    showOutput("myChart", "dimple")
  )
))
