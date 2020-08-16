library(shiny)
library(shinydashboard)
library(dplyr)
library(plotly)
library(ggrepel)
library(DT)

# For heroku deployment
# port <- Sys.getenv('PORT')

handle_plotly_exception <- function(data, func) {
  if (nrow(data) == 0) {
    ggplotly(ggplot() +
      theme_void() +
      geom_text(aes(0, 0, label = "Não há dados para as requisições enviadas")) +
      xlab(NULL))
  } else {
    func
  }
}

make_gene_plot <- function(data, ...) {
  handle_plotly_exception(
    data,
    ggplotly(data %>%
      ggplot(aes(...)) +
      geom_point(alpha = 0.7, color = "#296e6b") +
      guides(color = FALSE) +
      expand_limits(y = 0) +
      labs(
        x = NULL,
        y = "p-value (adjusted)"
      ), tooltip = "text")
  )
}

make_isoform_plot <- function(isa_data) {
  handle_plotly_exception(
    isa_data,
    ggplotly(ggplot(isa_data, aes(
      x = phenotype, y = val,
      group = enst, colour = enst,
      text = paste(
        "Region:", region,
        "\nTranscript ID:", enst,
        "\nSymbol:", symbol,
        "\nIF:", val
      )
    )) +
      geom_line() +
      geom_point() +
      scale_colour_discrete(guide = "none") +
      theme(
        axis.text.x = element_text(face = "italic", size = 10),
        legend.position = "none"
      ) +
      labs(x = "", y = "Isoform Switch Value"), tooltip = "text")
  )
}

load("./data-wrangling/plot_tables.rda")
isoforms_results <- read.csv("./data-wrangling/total_results_isa_long.csv")
theme_set(theme_bw())

ui <- dashboardPage(
  skin = "purple",
  dashboardHeader(title = "MDD"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Expressão Diferencial", tabName = "graphs", icon = icon("bar-chart")),
      menuItem("ISAnalyzeR", tabName = "isanalyzer", icon = icon("bar-chart")),
      menuItem(" Isoformas", tabName = "isoforms", icon = icon("cloudsmith")),
      menuItem("Dados", tabName = "data", icon = icon("table"))
    )
  ),
  dashboardBody(
    tabItems(
      
      ## UI - DGE e DTE
      tabItem(
        tabName = "graphs",
        fluidRow(
          box(
            title = "Escolha o gene", status = "warning",
            collapsible = TRUE, solidHeader = TRUE,
            HTML("<b>Ex.:</b> BDNF</br>"),
            textInput("genebox", label = NULL, placeholder = NULL)
          ),

          box(
            title = "Selecione a região cerebral", status = "warning",
            collapsible = TRUE, solidHeader = TRUE,
            selectInput(inputId = "regionbox", label = NULL, choices = NULL)
          ),
        ),
        fluidRow(
          box(width=12,
              status = "warning",
              textOutput(outputId ="gene_descr")
          ),
        ),
        fluidRow(
          box(width=12,
            status = "primary",
            plotOutput("dte_plot_iara")
          ),
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
          ),
        ),
      ),
      
      ## UI - IsanalyzeR pdfs
      tabItem(
        tabName = "isanalyzer",
        fluidRow(
          box(width=12,
              status = "warning",
            selectInput("groups", "Select a group", choices = unique(isa_df$group), 
                        selected = "aINS_female", multiple = F),
            selectInput("id_dtu", "Select a gene (name + ensG)", choices = unique(isa_df$id)),
            tags$h5("Description:"),
            textOutput(outputId = "description_dtu")
            )
          ),
        fluidRow(
          box(width=12,
              status = "primary", htmlOutput("file")
          )
        )
        ),
      
      ## UI - Isoformas - Slope plots
      tabItem(
        tabName = "isoforms",
        fluidRow(
          box(
            title = "Selecione a região cerebral", status = "warning",
            collapsible = TRUE, solidHeader = TRUE,
            selectInput(
              inputId = "regionisa", label = NULL,
              choices = unique(isoforms_results$region), selected = "Nac"
            )
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
      
      # UI - dados
      tabItem(
        tabName = "data",
        fluidRow(
          box(width=12,
            title = "Tabela de Dados DGE/DTE", status = "success", solidHeader = TRUE,
            dataTableOutput("tabledata")
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  set.seed(112358)

  # Server - DGE/DTE 
  choices <- reactive({
    choices <- plot_table %>%
      filter(hgnc_symbol == input$genebox) %>%
      distinct(region) %>%
      pull(region)
  })

  observe({
    updateSelectInput(session = session, inputId = "regionbox", choices = choices())
  })
  
  iara_dte_plot <- reactive({
    plot_table %>%
      filter(hgnc_symbol %in% input$genebox)
  })
  
  curr_data <- reactive({
    req(input$genebox)
    req(input$regionbox)
    plot_table %>%
      filter(hgnc_symbol == input$genebox & region == input$regionbox)
  })

  trans_filter <- reactive({
    curr_data() %>%
      filter(padj_tx < 0.05)
  })

  output$dge_plot <- renderPlotly({
    make_gene_plot(curr_data(),
      x = stringr::str_to_title(gender), y = padj_gene,
      text = paste("Region:", region, "\np-adj:", padj_gene)
    )
  })

  output$dte_plot <- renderPlotly({
    make_gene_plot(trans_filter(),
      x = stringr::str_to_title(gender), y = padj_tx,
      text = paste("Region:", region, "\nTranscript ID:", tx, "\np-adj:", padj_tx)
    )
  })
  
  output$gene_descr <- renderText({
    req(input$genebox)
    desc_list[input$genebox]
  })
  
  output$dte_plot_iara <- renderPlot({
    if (nrow(iara_dte_plot()) != 0){
    iara_dte_plot() %>% 
    ggplot() +
      facet_wrap(region ~ gender, drop = T) +
      geom_hline(aes(yintercept = logFC_gene), lty = 2, color = "#0000b1ff") +
      geom_hline(yintercept = 0, color = "#00000a8a") +
      geom_text(aes(x = 0.2, y = logFC_gene+.1, label = "gene logFC"), color = "#0000b1ff") + 
      scale_x_continuous(limits = c(0, 1), breaks = seq(0,1, 0.1)) +
      geom_point(aes(x = 0.5, y = logFC_tx, col = signif, size = prop)) + 
      scale_y_continuous(name = "logFC") +
      scale_size_continuous(name = "Transcript proportion") + 
      labs(title = paste0("DTE analysis for ", input$genebox, " gene"), x = "", y = "logFC") + 
      geom_label_repel(aes(x = 0.5, y = logFC_tx, label = tx, color = signif),
                       direction = "both",
                       segment.colour = "grey20",
                       nudge_x = 0.3,
                       segment.size = 0.3, box.padding = 0.5, show.legend = F) +
      scale_color_manual(name = "Transcript significance (padj <= 0.05)", values = c("S" = "red", "NS" = alpha("black", 0.7)),
                         labels = c("Not significant", "Significant")) +
      theme_bw() + 
      theme(axis.ticks.x = element_blank(), 
            axis.text.x = element_blank())}
    else{
      ggplot() +
        theme_void() +
        geom_text(aes(0, 0, label = "Não há dados para as requisições enviadas")) +
        xlab(NULL)
    }
  })
  
  # DTU server ----
  
  # Create dynamic selectInput for genes in each group
  observeEvent(input$groups, {
    x <- input$groups
    y <- isa_df$id[isa_df$group == x]
    updateSelectInput(session, inputId = "id_dtu", label = "Select a gene (name + ensG)",
                      choices = y,
                      selected = tail(y, 1))
  })
  
  # List all files in www/
  files <- reactive({
    list.files("www", full.names = T)
  })
  
  # Get the gene
  gene_ens <- reactive({
    gene <- as.character(unlist(strsplit(input$id_dtu, split = "_")))[1]
    ens <- as.character(unlist(strsplit(input$id_dtu, split = "_")))[2]
    return(c(gene, ens))
  })

  # Get the file corresponding to chosen gene
  graph <- reactive({
    x <- files()[( grepl(gene_ens()[1], files()) & grepl(gene_ens()[2], files()) ) & grepl(input$groups, files())]
    x <- gsub("www/", "", x)
    return(x)
  })
  
  # Outputs
  output$description_dtu <- renderText({
    desc_list[gene_ens()[1]]
  })
  
  observeEvent(input$id_dtu, {
    output$file <- renderUI({
      tags$iframe(src = graph(), style = "height:800px; width:100%;scrolling=yes")
    })
  })

  # observeEvent(input$id_dtu, {
  #   print(isa_df$entrezgene_id[match(gene_ens()[1], isa_df$gene_name)])
  # })
  
  ## Server - ISA slope Plots
  
  output$isa_fem <- renderPlotly({
    req(input$regionisa)
    isa_filter <- reactive({
      isoforms_results %>%
        dplyr::filter(region == input$regionisa & sex == "female")
    })
    make_isoform_plot(isa_filter())
  })

  output$isa_mal <- renderPlotly({
    req(input$regionisa)
    isa_filter <- reactive({
      isoforms_results %>%
        dplyr::filter(region == input$regionisa & sex == "male")
    })
    make_isoform_plot(isa_filter())
  })
  
  # SERVER - Data tables
  output$tabledata <- renderDataTable(
    plot_table,
    options = list(searchHighlight = TRUE,
                                scrollX = T), 
      filter = 'top'
  )
}

# Run the application
app <- shinyApp(ui = ui, server = server)

# For heroku deployment

# runApp(
#     appDir = app,
#     host = '0.0.0.0',
#     port = as.numeric(port)
# )
