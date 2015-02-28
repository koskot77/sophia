library(UsingR)
library(ggplot2)
library(rCharts)

an <- read.csv(file="output_10.csv", header=T, sep=',')

shinyServer(
  function(input, output) {
    nEvents <- 0
    output$myHist <- renderPlot({
      x <- 'ID' #input$inputX
      y <- input$inputY
      xlow  <- input$xlow
      xhigh <- input$xhigh
      ylow  <- input$ylow
      yhigh <- input$yhigh

#      if( x=='ID' )
          xlab <- "pattern ID"
#      if( x=='label' )
#          xlab <- "pattern"
      if( y=='count' )
          ylab <- "# of found patterns"
      if( y=='meanQual' )
          ylab <- "average quality"
      if( y=='sdQual' )
          ylab <- "standard deviation of qualities"
      if( y=='acc' )
          ylab <- "accuracy (probability of errorless read)"

        rangeX = max(an[[x]]) - min(an[[x]])
        rangeY = max(an[[y]]) - min(an[[y]])
        a <- subset(an, an[[x]] > min(an[[x]]) + xlow/100.*rangeX & an[[x]] < min(an[[x]]) + xhigh/100.*rangeX &
                        an[[y]] > min(an[[y]]) + ylow/100.*rangeY & an[[y]] < min(an[[y]]) + yhigh/100.*rangeY)
        ggplot(data=a, aes_string(x=x,y=y)) + geom_point(alpha=0.1,color="blue") + 
           theme(title = element_text(size = 15), axis.title.x = element_text(size = 15)) +
           labs(x=xlab,y=ylab,title="")
    })
    output$events <- renderPrint({
      x <- 'ID' #input$inputX
      y <- input$inputY
      xlow  <- input$xlow
      xhigh <- input$xhigh
      ylow  <- input$ylow
      yhigh <- input$yhigh

      nEvents <- 0

        rangeX = max(an[[x]]) - min(an[[x]])
        rangeY = max(an[[y]]) - min(an[[y]])
        a <- subset(an, an[[x]] > min(an[[x]]) + xlow/100.*rangeX & an[[x]] < min(an[[x]]) + xhigh/100.*rangeX &
                        an[[y]] > min(an[[y]]) + ylow/100.*rangeY & an[[y]] < min(an[[y]]) + yhigh/100.*rangeY)
        nEvents <- dim(a)[1]

      if( nEvents > 100 )
         paste(nEvents," events captured, sampling 100 random events below")
      else
         paste(nEvents," events captured, all are shown below")
    })

    output$myChart = renderChart({

      x <- 'ID' #input$inputX
      y <- input$inputY
      xlow  <- input$xlow
      xhigh <- input$xhigh
      ylow  <- input$ylow
      yhigh <- input$yhigh

      nEvents <- 0

        rangeX = max(an[[x]]) - min(an[[x]])
        rangeY = max(an[[y]]) - min(an[[y]])
        a <- subset(an, an[[x]] > min(an[[x]]) + xlow/100.*rangeX & an[[x]] < min(an[[x]]) + xhigh/100.*rangeX &
                        an[[y]] > min(an[[y]]) + ylow/100.*rangeY & an[[y]] < min(an[[y]]) + yhigh/100.*rangeY)
        nEvents <- dim(a)[1]

      if( nEvents>100 ){
          sample <- data.frame( sel=rbinom(dim(a)[[1]], 1, 100./nEvents) )
          a <- a[as.numeric( rownames( subset(sample,sel==1) ) ),]
      }

        d1 <- dPlot(x=x, y=y, groups = c("label"), data = a, type = "bubble") #, height=800, width=1000)

        d1$addParams(dom = 'myChart')

        return(d1)

    })

  }
)
