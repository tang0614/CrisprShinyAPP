#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#


library(dplyr)
library(tidyr)
library(tibble)
library(gtools)
library(shiny)
library(shinyWidgets)
library(shinythemes)
library(shinycssloaders)
library(scales)
library(treemap)
library(gtools)
library(highcharter)
library(stringr)
library(lexicon)
library(visNetwork)
library(rsconnect)
library(networkD3)
library(magrittr)

if (FALSE) {
  library(RSQLite)
  library(dbplyr)
}
library(shiny)

#CSS

df <- read.csv("corr_rank_within70_new.txt", header = TRUE,stringsAsFactors = FALSE)
#default_gene = c('ATG14','ATG101','ULK1','WIPI2')


# Define UI for dataset viewer app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Network Analysis in Crispr2020Q1 Dataset"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      #Input Panel Title
      titlePanel("Input Panel"),
      hr(),
      helpText("Select the sign of pearson correlation."),
      
      fluidRow(
        column(11, 
               radioButtons("corr", "Correlation:",
                            c("Positive" = "positive",
                              "Negative" = "negative",
                              "All" = "all"),
                            selected = "positive"),
               
        )
      ),
      
      
      
      fluidRow(
        column(11, 
               helpText("Enter capitalized gene symbol separatd by space/new lines."),
               textAreaInput("caption", "Gene Symbol", "ATG14 ATG101 ULK1")
        ),
        
      ),
      
      fluidRow(
        column(11,
               helpText("Select the top n correlated genes for each selected gene above."),
               sliderInput("rank", "Rank:",
                           min = 1, max = 140,
                           value = 10)
               
               
        )
      ),
      
      
      fluidRow(
        column(11,
               helpText("Select the numeric value indicating either the strength of the node repulsion (negative value) or attraction (positive value)."),
               sliderInput("charge", "Charge:",
                           min = -500, max = 5,
                           value = -300)
               
               
        )
      ),
      
      
      
      fluidRow(
        column(11,
               helpText("Refresh and submit the value you entered above."),
               submitButton("Submit", icon("refresh"))
               
        )
      ),
      
      
      
      
      
      br(),
      hr(),
      br(),
      
      fluidRow(
        column(11, 
               
               tags$p(HTML("<b>Dataset Description</b> :This <a href=\"https://depmap.org/portal/download/\">DepMap release</a> contains data from CRISPR knockout screens from project Achilles which contains the results of genome-scale CRISPR knockout screens for 18,333 genes in 739 cell lines. These networks enable people to discover novel clusters of genes that might share functional
                            associations and connections with wach other based on their essentiality profiles from screen data across multiple cell lines.")),
               
        ))
      
      
      
      
      
    ),
    
    mainPanel(
      
      fluidRow(
        
        column(11,offset = 1,
               h3('Network Visualization'),
               helpText("The width of links indicates the strengh of gene correlations. Thick line represents greater positive/negative correlations."),
               helpText("Use mouse to zoom in/out to change the size of the networks. ")
        )
      ),
      
      fluidRow(
        
        
        column(3, offset = 8,
               uiOutput("render_ui_genes"),
               downloadButton("downloadNet", "Download network below")
               #conditionalPanel("output.show",  uiOutput("render_ui_genes"))
        )
      ),
      
      fluidRow(
        
        column(11,align="center",
               tabPanel("Force Network", forceNetworkOutput("force"))
        ),
        
        
      ),
      
      
      
      
      
      
      
      br(),
      br(),
      
      br(),
      hr(),
      
      
      
      fluidRow(
        column(4,
               
               h3('Summary Table I'),
               helpText("Genes that are showing in the plot."),
               downloadButton("downloadData", "Download Table")
               
        ),
        
        
        
        column(5,offset = 1,
               
               h3('Summary Table II'),
               helpText("Top correlated genes from the entire dataset."),
               downloadButton("downloadDataAll", "Download Table")
               
        )
        
        
      ),
      br(),
      br(),
      
      
      
      
      fluidRow(
        column(5,
               h3('Additional Info'),
               helpText("Genes that are not showing in the plot."),
               downloadButton("downloadMissingData", "Download Missing Genes")
               
               
               
        )
        
      ),
      fluidRow(
        column(5, 
               verbatimTextOutput("value", placeholder = TRUE)
               
               
               
        )
        
      )
      
      
    )
    
    
    
  )
)

# Define server logic to summarize and view selected dataset ----
server <- function(input, output,session) {
  
  
  input_genes<- reactive({
    g<-str_replace_all(input$caption, "[\r\n]" , " ")
    l<-unique(unlist(strsplit(g, split=" ")))
    l
    
  })
  
  
  
  missed_genes <- reactive({
    g<-input_genes()
    missed_genes<-g[!(g %in% df$Gene)]
    a<-str_replace_all(missed_genes, "[\r\n]" , " ")
    b<-unlist(strsplit(a, split=" "))
    my_vector =paste(b, collapse=";")
    my_vector
    
    
  })
  
  genes_in_db <-reactive({
    
    g<-input_genes()
    filted_genes<-g[(g %in% df$Gene)]
    a<-str_replace_all(filted_genes, "[\r\n]" , " ")
    b<-unlist(strsplit(a, split=" "))
    b
    
  })
  
  output$render_ui_genes<- renderUI({
    selectInput("gene", label = "Selected genes", 
                choices=mixedsort(as.vector(unique(df$Gene))), selected=genes_in_db(),multiple=TRUE)
  })
  
  
  
  output$value <- renderText({ missed_genes() })
  
  #select button
  
  
  dataInput1 <- reactive({
    m<-df
    if(input$corr == 'positive'){
      m<-filter(m, Correlation >= 0)
      return(m)
      
    }
    else if(input$corr == 'negative'){
      m<-filter(m, Correlation < 0)
      return(m)
      
    }
    else{return(m)}
    
  })
  
  
  dataInput2 <- reactive({
    m<-dataInput1()
    list_1<-m$Gene
    sub_list1<- list_1 %in% input$gene
    
    if (is.null(input$gene)){
      req(input$gene)
    } else {
      m <- subset(m, sub_list1)
    } 
    m
  })
  
  dataSort <- reactive({
    m<-dataInput2()
    r <- input$rank
    
    if (is.null(input$rank)){
      req(input$rank)
    } else {
      m<-m %>%group_by(Gene)%>%
        arrange(desc(Correlation, .group_by=TRUE))%>%
        slice(1:r)
      
    } 
    
    
    
    m<-subset(m, select = -c(X))
    m
    
  })
  

  
  dataSort_entire <- reactive({
    m<-dataInput1()
    r <- input$rank
    
    if (is.null(input$rank)){
      req(input$rank)
    } else {
      m<-m %>%group_by(Gene)%>%
        arrange(desc(Correlation, .group_by=TRUE))%>%
        slice(1:r)
      
    } 
    
    
    
    m<-subset(m, select = -c(X))
    m
    
  })
  
 
  
  
  
  MisNodes<- reactive({
    src <-input$gene
    
    m<-dataSort()
    m<-data.frame(unique_gene=union(m$Gene, m$Gene2))
    m$ID <- 0:(nrow(m)-1)
    m$size <- 15
    m$group <- ifelse(m$unique_gene %in% src, "lions", "tigers")
    m<-subset(m, select = c(ID,unique_gene,size,group))
    m
    
    
    
    
  })
  
  MisLinks<- reactive({
    m<-dataSort()
    nodes<-MisNodes()
    
    m<-m%>% 
      left_join(nodes, by = c("Gene" = "unique_gene")) %>% rename(source = ID)
    m<- m%>% 
      left_join(nodes, by = c("Gene2" = "unique_gene")) %>% rename(target = ID)
    
    m<-subset(m, select = c(source,target,Correlation))
    m
  })
  
  output$force <- renderForceNetwork({
    l<-MisLinks()
    n<-MisNodes()
    
    ColourScale <- 'd3.scaleOrdinal()
            .domain(["lions", "tigers"])
           .range(["#D93434","#3FBFBF"]);'
    
    
    MisLinks<-as.data.frame(l)
    MisNodes<-as.data.frame(n)
    forceNetwork(Links = MisLinks, Nodes = MisNodes, Source = "source",
                 Target = "target", 
                 Value = "Correlation", 
                 NodeID = "unique_gene",
                 Nodesize = "size",
                 Group = "group",
                 fontSize = 10, # Font size
                 linkDistance = networkD3::JS("function(d) { return 20*(1/Math.abs(d.value));}"),
                 linkWidth = networkD3::JS("function(d) { return 15*Math.abs(d.value); }"),
                 opacity = 0.8,
                 opacityNoHover = 0.8,
                 zoom = TRUE,
                 bounded = TRUE,
                 charge = input$charge,
                 colourScale = JS(ColourScale))
  })
  
  vis2<-reactive({
    ColourScale <- 'd3.scaleOrdinal()
            .domain(["lions", "tigers"])
           .range(["#D93434","#3FBFBF"]);'
    
    l<-MisLinks()
    n<-MisNodes()
    
    
    MisLinks<-as.data.frame(l)
    MisNodes<-as.data.frame(n)
    
    forceNetwork(Links = MisLinks, Nodes = MisNodes, Source = "source",
                 Target = "target", 
                 Value = "Correlation", 
                 NodeID = "unique_gene",
                 Nodesize = "size",
                 Group = "group",
                 fontSize = 10, # Font size
                 linkDistance = networkD3::JS("function(d) { return 20*(1/Math.abs(d.value));}"),
                 linkWidth = networkD3::JS("function(d) { return 15*Math.abs(d.value); }"),
                 opacity = 0.8,
                 opacityNoHover = 1,
                 zoom = TRUE,
                 charge = input$charge,
                 colourScale = JS(ColourScale)) 
  })
  
  

  
  
  # Downloadable csv of selected dataset ----
  output$downloadData <- downloadHandler(
    filename = function() {
      paste('table', ".csv", sep = "")
    },
    content = function(file) {
      write.csv(dataSort(), file, row.names = FALSE)
    }
  )
  
  output$downloadDataAll <- downloadHandler(
    filename = function() {
      paste('table', ".csv", sep = "")
    },
    content = function(file) {
      write.csv(dataSort_entire(), file, row.names = FALSE)
    }
  )
  
  output$downloadMissingData <- downloadHandler(
    filename = function() {
      paste('missingData', ".csv", sep = "")
    },
    content = function(file) {
      write.csv(missed_genes(), file, row.names = FALSE)
      
    }
  )
  
  output$downloadNet <- downloadHandler(
    filename = function() {
      paste('network', ".html", sep = "")
    },
    content = function(file) {
      saveNetwork(vis2(),file )
      
    }
  )
  
  
  
}

shinyApp(ui, server)