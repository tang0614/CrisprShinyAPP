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

if (FALSE) {
    library(RSQLite)
    library(dbplyr)
}
library(shiny)

df <- read.csv("corr_rank_within70.txt", header = TRUE,stringsAsFactors = FALSE)
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
                column(8, 
                       radioButtons("corr", "Correlation:",
                                    c("Positive" = "positive",
                                      "Negative" = "negative",
                                      "All" = "all"),
                                    selected = "positive"),
                       
                )
            ),
            
            
            
            br(),
            fluidRow(
              column(11, 
                     helpText("Enter capitalized gene symbol separatd by space/new lines."),
                     textAreaInput("caption", "Gene Symbol", "ATG14 ATG101 ULK1")
              ),
            
            ),
            

    
            
            br(),
            fluidRow(
                column(11,
                       helpText("Select the top n correlated genes for each selected gene above."),
                       sliderInput("rank", "Rank:",
                                   min = 1, max = 140,
                                   value = 10)
                     
                       
                )
            ),
            
           
            
            br(),
            hr(),
            br(),
            
            fluidRow(
                column(11, 
               
                tags$p(HTML("<b>Dataset Description</b> :This DepMap release contains data from CRISPR knockout screens from project Achilles, as well as genomic characterization data from the CCLE project.
                This Achilles dataset contains the results of genome-scale CRISPR knockout screens for 18,333 genes in 739 cell lines.<a href=\"https://depmap.org/portal/download/\">References</a>")),
                
            ))
            
     
            
        ),
        
        mainPanel(
            
            fluidRow(
               
                column(11,offset = 1,
                       h3('Network Visualization for Gene Correlations'),
                      
                       helpText("Click and drag nodes below to see in detail."),
                       helpText("The width of links indicates the strengh of gene correlations. Thick line represents greater positive/negative correlations.")
                )
            ),
  
            fluidRow(
               
                column(9,align="center",
                       tabPanel("Force Network", forceNetworkOutput("force"))
                ),
                
                
               
               column(2, 
                         br(),
                         br(),
                         uiOutput("render_ui_genes")
                         #conditionalPanel("output.show",  uiOutput("render_ui_genes"))
                )
            ),
            
            fluidRow(
                column(9,align="center",
                       downloadButton("downloadV2", "Download above network")
                       
                )
            ),
            
            
            br(),
            br(),
            
            br(),
            br(),
            
            br(),
            br(),
            
           
      

            fluidRow(
                column(4,offset = 1,
                       h3('Summary Table'),
                       helpText("Click button to download the table.")
                ),
                
            
                
                column(4,offset = 3,
                       h3('Missing Genes'),
                       helpText("Genes that are not showing in the plot.")
                )
              
                
            ),
            br(),
            
           
           
            fluidRow(
                column(5, offset = 1,
                       tableOutput('table1'),
                       downloadButton("downloadData", "Download Table")
                ),
             
                column(3,offset = 2,
                       br(),
                       verbatimTextOutput("value", placeholder = TRUE),
                       downloadButton("downloadMissingData", "Download Missing Genes")
                ),
                
                
              
            ),
            
          
    
            
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
    
   
    
    MisNodes<- reactive({
        src <-input$gene
      
        m<-dataSort()
        m<-data.frame(unique_gene=union(m$Gene, m$Gene2))
        m$ID <- 0:(nrow(m)-1)
        m$size <- 40
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
                     fontSize = 15, # Font size
                     linkDistance = networkD3::JS("function(d) { return Math.abs(d.value)/0.004;}"),
                     linkWidth = networkD3::JS("function(d) { return 15*Math.abs(d.value); }"),
                     opacity = 0.8,
                     opacityNoHover = 0.8,
                     colourScale = JS(ColourScale))
    })
    
    vis2<-reactive({
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
                     fontSize = 15, # Font size
                     linkDistance = networkD3::JS("function(d) { return Math.abs(d.value)/0.004;}"),
                     linkWidth = networkD3::JS("function(d) { return 15*Math.abs(d.value); }"),
                     opacity = 0.8,
                     opacityNoHover = 1,
                     colourScale = JS(ColourScale))
    })
    
    
    
    

    
    vis1<-reactive({
        m<- dataSort()
        src <- m$Gene
        target <- m$Gene2
        networkData <- data.frame(src, target)
        simpleNetwork(networkData, opacity = input$opacity)
    })
    
    output$table1 <- renderTable(dataSort())
    
    
    # Downloadable csv of selected dataset ----
    output$downloadData <- downloadHandler(
        filename = function() {
            paste(dataSort(), ".csv", sep = "")
        },
        content = function(file) {
            write.csv(dataSort(), file, row.names = FALSE)
        }
    )
    
    output$downloadMissingData <- downloadHandler(
      filename = function() {
        paste(dataSort(), ".csv", sep = "")
      },
      content = function(file) {
        write.csv(missed_genes(), file, row.names = FALSE)
      }
    )
    

    
    output$downloadV1 <- downloadHandler(
        filename = function() {
            'Net1.html'
        },
        content = function(file) {
            saveNetwork(vis1(), file)
        }
    )
    
    output$downloadV2 <- downloadHandler(
        filename = function() {
            'Net2.html'
        },
        content = function(file) {
            saveNetwork(vis2(), file)
        }
    )
    
    


}

# Create Shiny app ----
shinyApp(ui, server)