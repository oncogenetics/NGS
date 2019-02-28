

# https://github.com/rstudio/DT/issues/267#issuecomment-458078711

library(shiny)
library(DT)

myModal <- function() {
  div(id = "test",
      modalDialog(downloadButton("download1","Download iris as csv"),
                  br(),
                  br(),
                  downloadButton("download2","Download iris as csv2"),
                  easyClose = TRUE, title = "Download Table")
  )
}

ui <- basicPage(
  DTOutput("dtable")
)

server <- function(input, output, session){
  output$dtable <- renderDT(
    datatable(iris,
              extensions = 'Buttons',
              options = list(
                dom = 'Bfrtip',
                buttons = list(
                  "copy",
                  list(
                    extend = "collection",
                    text = 'download entire dataset',
                    action = DT::JS("function ( e, dt, node, config ) {
                                    Shiny.setInputValue('test', true, {priority: 'event'});
}")
                  )
                    )
                  )
                )
              )
  
  observeEvent(input$test, {
    print("hello")
    showModal(myModal())
  })
  
  
  output$download1 <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(iris, file)
    }
  )
  
  output$download2 <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv2(iris, file)
    }
  )
}

shinyApp(ui, server)