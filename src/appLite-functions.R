finish_screen <-function()
{
  renderUI({
    fluidRow(column(5,h4("Results")),
             column(5, downloadButton("downloadData", label = "Download")))
  })
}