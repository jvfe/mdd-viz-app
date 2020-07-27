library(shiny)
library(shinydashboard)
library(dplyr)
library(plotly)
library(DT)
library(igraph)
library(visNetwork)

# For heroku deployment
# port <- Sys.getenv('PORT')

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

handle_plotly_exception <- function(data, ...) {
  if (nrow(data) == 0) {
    ggplotly(ggplot() +
      theme_void() +
      geom_text(aes(0, 0, label = "Não há dados para as requisições enviadas")) +
      xlab(NULL))
  } else {
    make_gene_plot(data, ...)
  }
}

make_isoform_plot <- function(isa_data, ...) {
  if (nrow(isa_data) == 0) {
    ggplotly(ggplot() +
               theme_void() +
               geom_text(aes(0, 0, label = "Não há dados para as requisições enviadas")) +
               xlab(NULL))
  } else {
    ggplotly(ggplot(isa_filter, aes(x = phenotype, y = val, 
                                    group = enst, colour = enst, ...)) + 
               geom_line() +
               geom_point() +
               scale_colour_discrete(guide = 'none') +
               theme(axis.text.x = element_text(face = "italic", size = 10),
                     legend.position='none') +
               labs(x="", y="Isoform Switch Value"))
  }
}

all_nets <- readRDS("./data-wrangling/dge_nets.rds")
df_genes_with_symbols <- readRDS("./data-wrangling/df_genes_with_symbols.rds")
isoforms_results <- read.csv("./data-wrangling/total_results_isa_long.csv")
theme_set(theme_bw())

ui <- dashboardPage(
  skin = "purple",
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
      tabItem(
        tabName = "graphs",
        fluidRow(
          box(
            title = "Escolha o gene", status = "warning",
            collapsible = TRUE, solidHeader = TRUE,
            HTML("<b>Ex.:</b> MAPK10</br>"),
            textInput("genebox", label = NULL, placeholder = NULL)
          ),

          box(
            title = "Selecione a região cerebral", status = "warning",
            collapsible = TRUE, solidHeader = TRUE,
            selectInput(inputId = "regionbox", label = NULL, choices = NULL)
          )
        ),
        fluidRow(
          box(
            title = "DGE", status = "primary", solidHeader = TRUE,
            plotlyOutput("dge_plot")
          ),
          box(
            title = "DTE", status = "primary", solidHeader = TRUE,
            HTML("<b>Apenas com padj<0.05</b></br>"), br(),
            plotlyOutput("dte_plot")
          )
        ),
        fluidRow(
          box(
            title = "ISA - Mulher", status = "primary", solidHeader = TRUE,
            plotlyOutput("isa_fem")
          ),
          box(
            title = "ISA - Homem", status = "primary", solidHeader = TRUE,
            plotlyOutput("isa_mal")
          )
        ),
      ),

      tabItem(
        tabName = "nets",
        fluidRow(
          box(
            title = "Selecione a região cerebral", status = "warning",
            collapsible = TRUE, solidHeader = TRUE,
            selectInput(
              inputId = "regionet", label = NULL,
              choices = unique(all_nets$region), selected = "Nac"
            )
          )
        ),
        fluidRow(
          box(status = "primary", visNetworkOutput("network"), width = "90%")
        )
      ),

      tabItem(
        tabName = "data",
        fluidRow(
          box(
            title = "Tabela de Dados", status = "success", solidHeader = TRUE,
            dataTableOutput("tabledata"), width = "90%"
          )
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

  curr_data <- reactive({
    req(input$genebox)
    req(input$regionbox)
    df_genes_with_symbols %>%
      filter(hgnc_symbol == input$genebox & region == input$regionbox)
  })

  trans_filter <- reactive({
    curr_data() %>%
      filter(transcript < 0.05)
  })

  output$dge_plot <- renderPlotly({
    handle_plotly_exception(curr_data(),
      x = stringr::str_to_title(sex), y = gene,
      text = paste("Region:", region, "\np-adj:", gene)
    )
  })

  output$dte_plot <- renderPlotly({
    handle_plotly_exception(trans_filter(),
      x = stringr::str_to_title(sex), y = transcript,
      text = paste("Region:", region, "\nTranscript ID:", txID, "\np-adj:", transcript)
    )
  })
  
  output$isa_fem <- renderPlotly({
    req(input$regionbox)
    isa_filter <- reactive({isoforms_results %>% 
      dplyr::filter(region == input$regionbox & sex == 'female')})
    make_isoform_plot(isa_filter())
  })

  output$isa_mal <- renderPlotly({
    req(input$regionbox)
    isa_filter <- reactive({isoforms_results %>% 
        dplyr::filter(region == input$regionbox & sex == 'male')})
    make_isoform_plot(isa_filter())
  })
  
  output$network <- renderVisNetwork({
    net <- all_nets %>%
      select(-c(stringId_A, stringId_B)) %>%
      filter(region == input$regionet) %>%
      graph_from_data_frame(directed = FALSE)

    visIgraph(net) %>%
      visIgraphLayout(layout = "layout_nicely") %>%
      visNodes(size = 10) %>%
      visOptions(
        highlightNearest = list(enabled = T, hover = T),
        nodesIdSelection = T
      )
  })

  output$tabledata <- renderDataTable({
    df_genes_with_symbols %>%
      setNames(c("ENSG", "ENST", "Gene Exp. padj", "Transcript Exp. padj", "Region", "Gender", "Gene Name")) %>%
      datatable()
  })
}

# Run the application
app <- shinyApp(ui = ui, server = server)

# For heroku deployment

# runApp(
#     appDir = app,
#     host = '0.0.0.0',
#     port = as.numeric(port)
# )
