library(shiny)
library(bsicons)
library(bslib)

# Define UI ----
ui <- page_sidebar(
  title = "tiiiitle",
  sidebar = sidebar("Sidebar"),
  card(
    card_header("limbxlimb"),
    "yeh smeeerking",
    card_image("../dog.png"),
    card_footer("dem doge")
  ),
  value_box(
    title = "price for de doog",
    value = 100000000,
    showcase = bsicons::bs_icon("bar-chart"),
    theme = "teal"
  )
  
)

# Define server logic ----
server <- function(input, output) {

}

shinyApp(ui = ui, server = server)
