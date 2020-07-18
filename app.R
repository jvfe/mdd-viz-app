library(shiny)
library(shinydashboard)
library(dplyr)
library(plotly)
library(DT)
library(igraph)
library(visNetwork)

make_gene_plot <- function(data, ...) {
    ggplotly(data %>% 
        ggplot(aes(...)) +
            geom_point(alpha = 0.7, color = "#296e6b") +
            guides(color = FALSE) +
            expand_limits(y = 0) +
            labs(
                x = NULL,
                y = "p-value (adjusted)"
            ), tooltip = "text")
}


all_nets <- readRDS("./data-wrangling/dge_nets.rds")
df_genes_with_symbols <- readRDS("./data-wrangling/df_genes_with_symbols.rds")
theme_set(theme_bw())

ui <- dashboardPage(skin = "purple",
    dashboardHeader(title = "MDD"),
    dashboardSidebar(
        sidebarMenu(
            menuItem("Gráficos", tabName = "graphs", icon = icon("bar-chart")),
            menuItem("Redes PPI", tabName = "nets", icon = icon("cloudsmith")),
            menuItem("Dados", tabName = "data", icon = icon("table"))
        )
    ),
    dashboardBody(
        tabItems(
            tabItem(tabName = "graphs",
                    fluidRow(
                        box(title = "Escolha o gene", status = "warning", 
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
            
            tabItem(tabName = "nets",
                    fluidRow(
                        box(title = "Selecione a região cerebral", status = "warning", 
                            collapsible = TRUE, solidHeader = TRUE,
                            selectInput(inputId = "regionet", label = NULL, 
                                        choices = unique(all_nets$region), selected = "Nac"))
                    ),
                    fluidRow(
                        box(status = "primary", visNetworkOutput("network"), width = "90%")
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
        data_gen <- df_genes_with_symbols %>% 
            filter(hgnc_symbol == input$genebox & region == input$regionbox)
        
        make_gene_plot(data_gen, x = stringr::str_to_title(sex), y = gene, 
                       text = paste("Region:", region, "\np-adj:", gene))
        

    })
    
    output$dte_plot <- renderPlotly({
        req(input$genebox)
        req(input$regionbox)
        data_trans <- df_genes_with_symbols %>% 
            filter(hgnc_symbol == input$genebox & region == input$regionbox & transcript < 0.05)
        
        make_gene_plot(data_trans, x = stringr::str_to_title(sex), y = transcript,
                       text = paste("Region:", region, "\nTranscript ID:", txID, "\np-adj:", transcript))
        
        
    })
    
    output$network <- renderVisNetwork({
        net <- all_nets %>% 
            select(-c(stringId_A, stringId_B)) %>% 
            filter(region == input$regionet) %>% 
            graph_from_data_frame(directed = FALSE)
        
        visIgraph(net) %>%
            visIgraphLayout(layout = "layout_nicely") %>%
            visNodes(size = 10) %>%
            visOptions(highlightNearest = list(enabled = T, hover = T), 
                       nodesIdSelection = T)
    })
    
    output$tabledata <- renderDataTable({
        df_genes_with_symbols %>% 
            setNames(c("ENSG", "ENST", "Gene Exp. padj", "Transcript Exp. padj", "Region", "Gender", "Gene Name")) %>% 
            datatable()
    })

}

# Run the application
shinyApp(ui = ui, server = server)
