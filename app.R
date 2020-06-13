library(shiny)
library(shinydashboard)

load("./data-wrangling/diff_tx_corrected.rda")

ui <- dashboardPage(
    dashboardHeader(title = "LITCovid19 Text Analysis"),
    dashboardSidebar(
        sidebarMenu(
            menuItem("Entities", tabName = "Entities", icon = icon("bar-chart")),
            menuItem("Network", tabName = "Network", icon = icon("arrows-alt"))
        )
    ),
    dashboardBody(
        tabItems(
            tabItem(tabName = "Entities",
                    fluidRow(
                        box(title = "Selecione o gene", status = "warning", 
                            collapsible = TRUE, solidHeader = TRUE,
                            textInput("genebox", label = NULL, placeholder = NULL)),
                        box(title = "Selecione a região cerebral", status = "warning", 
                            collapsible = TRUE, solidHeader = TRUE,
                            selectInput(inputId = "regionbox", label = NULL, choices = NULL))
                    ),
                    fluidRow(
                        box()
                    )
            ),

            tabItem(tabName = "Network",
                    fluidRow(
                        box(),
                        box()
                    ),
                    fluidRow(
                        box()
                    )
            )
        )
    )
)
server <- function(input, output, session) {

    set.seed(122)
    
    choices <- reactive({
        choices <- df_res_padj %>%
            separate(group, c("region", "sex"), "_") %>% 
            filter(geneID == input$genebox) %>%
            # filter(geneID == "ENSG00000153310") %>%
            distinct(region) %>% 
            pull(region)
    })
    
    observe({
        updateSelectInput(session = session, inputId = "regionbox", choices = choices())
    })

}

# Run the application
shinyApp(ui = ui, server = server)
