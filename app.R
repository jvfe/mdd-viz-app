library(shiny)
library(shinydashboard)
library(tidyverse)
library(plotly)
library(DT)

df_genes_with_symbols <- readRDS("./data-wrangling/df_genes_with_symbols.rds")

ui <- dashboardPage(skin = "purple",
    dashboardHeader(title = "MDD"),
    dashboardSidebar(
        sidebarMenu(
            menuItem("Gráficos", tabName = "graphs", icon = icon("bar-chart")),
            menuItem("Dados", tabName = "data", icon = icon("arrows-alt"))
        )
    ),
    dashboardBody(
        tabItems(
            tabItem(tabName = "graphs",
                    fluidRow(
                        box(title = "Selecione o gene", status = "warning", 
                            collapsible = TRUE, solidHeader = TRUE,
                            textInput("genebox", label = NULL, placeholder = NULL)),
                        box(title = "Selecione a região cerebral", status = "warning", 
                            collapsible = TRUE, solidHeader = TRUE,
                            selectInput(inputId = "regionbox", label = NULL, choices = NULL))
                    ),
                    fluidRow(
                        box(title = "DGE", status = "primary", solidHeader = TRUE,
                            plotlyOutput("dge_plot")),
                        box(title = "DTE", status = "primary", solidHeader = TRUE,
                            HTML("<b>Apenas com padj<0.05</b></br>"), br(),
                            plotlyOutput("dte_plot"))
                    )
            ),

            tabItem(tabName = "data",
                    fluidRow(
                        box(title = "Tabela de Dados", status = "success", solidHeader = TRUE,
                            dataTableOutput("tabledata"), width = "90%")
                    )
            )
        )
    )
)
server <- function(input, output, session) {

    set.seed(112358)
    
    choices <- reactive({
        choices <- df_genes_with_symbols %>%
            filter(hgnc_symbol == input$genebox) %>%
            distinct(region) %>% 
            pull(region)
    })
    
    observe({
        updateSelectInput(session = session, inputId = "regionbox", choices = choices())
    })
    
    output$dge_plot <- renderPlotly({
        req(input$genebox)
        req(input$regionbox)
        ggplotly(df_genes_with_symbols %>% 
            filter(hgnc_symbol == input$genebox & region == input$regionbox) %>% 
            ggplot(aes(x = str_to_title(sex), y = gene, 
                       text = paste("Region:", region, "\np-adj:", gene))) +
            geom_point(alpha = 0.7, color = "#296e6b") +
            guides(color = FALSE) +
            expand_limits(y = 0) +
            labs(
                x = NULL,
                y = "p-value (adjusted)"
            ), tooltip = "text")
        

    })
    
    output$dte_plot <- renderPlotly({
        req(input$genebox)
        req(input$regionbox)
        ggplotly(df_genes_with_symbols %>% 
                     filter(hgnc_symbol == input$genebox & region == input$regionbox & transcript < 0.05) %>% 
                     ggplot(aes(x = str_to_title(sex), y = transcript,
                                text = paste("Region:", region, "\nTranscript ID:", txID, "\np-adj:", transcript))) +
                     geom_jitter(alpha = 0.7, color = "#296e6b", 
                                 height = 0, width = 0.4) +
                     guides(color = FALSE) +
                     labs(
                         x = NULL,
                         y = "p-value (adjusted)"
                     ), tooltip = "text")
        
        
    })
    
    output$tabledata <- renderDataTable({
        df_genes_with_symbols %>% 
            setNames(c("ENSG", "ENST", "Gene Exp. padj", "Transcript Exp. padj", "Region", "Gender", "Gene Name")) %>% 
            datatable()
    })

}

# Run the application
shinyApp(ui = ui, server = server)
