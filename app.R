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
df <- read.csv("autophagy.txt", header = TRUE,stringsAsFactors = FALSE)

default_gene = c('ATG14','ATG101','ULK1','WIPI2')


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
                column(8, 
                       
                       
                       selectInput("gene", label = "Gene Symbol", 
                                   choices=mixedsort(as.vector(unique(df$Gene))), selected=default_gene,multiple=TRUE)
                       
                ),
                br(),
                br(),
                column(4,prettyCheckbox("gene_selectall", "Select All"))
            ),
            
            
            br(),
            fluidRow(
                column(8,
                       helpText("Select the top n correlated genes for each selected gene above."),
                       sliderInput("rank", "Rank:",
                                   min = 1, max = 50,
                                   value = 10)
                     
                       
                )
            ),
            
            br(),
           
            fluidRow(
                column(8,
                       sliderInput("opacity", "Opacity for network", 0.8, min = 0.1,
                                   max = 1, step = .1)
                 
                       
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
               
                column(8,offset = 1,
                       h3('Network Visualization for Gene Correlations'),
                      
                       helpText("Click and drag nodes below to see in detail."),
                       helpText("Both networks represent same data in different designs"),
                       helpText("For the left network, the width of links indicates the strengh of gene correlations. Thick line represents greater positive/negative correlations."),
                )
            ),
  
            fluidRow(
                column(5,align="center",
                       tabPanel("Simple Network", simpleNetworkOutput("simple"))
                )
                ,
                column(7,align="center",
                       tabPanel("Force Network", forceNetworkOutput("force"))
                )
            ),
            fluidRow(
                column(5,offset = 1,align="center",
                       downloadButton("downloadV1", "Download above network")
                       
                )
                ,
                column(5,align="center",
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
                column(4,
                       h3('Summary Table'),
                       helpText("Click button to download the table.")
                ),
                column(3,
                       br(),
                       br(),
                       downloadButton("downloadData", "Download Table")
                )
                
               
            ),
            br(),
            
            fluidRow(
                column(6,
                       tableOutput('table1')
                )
              
            )
            
            
            
            
            
            
        )
        
       
    )
)

# Define server logic to summarize and view selected dataset ----
server <- function(input, output,session) {
    
    #select button

    
    observe({
        if (input$gene_selectall){
            updateSelectInput(session, "gene", selected = as.vector(unique(df$Gene)))
        }else{
            updateSelectInput(session, "gene", selected =  default_gene)
        }
        
    })
    

    
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
        m<-dataSort()
        m<-data.frame(unique_gene=union(m$Gene, m$Gene2))
        m$ID <- 0:(nrow(m)-1)
        m$size <- 40
        m<-subset(m, select = c(ID,unique_gene,size))
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
       
        
        MisLinks<-as.data.frame(l)
        MisNodes<-as.data.frame(n)
        forceNetwork(Links = MisLinks, Nodes = MisNodes, Source = "source",
                     Target = "target", 
                     Value = "Correlation", 
                     NodeID = "unique_gene",
                     Nodesize = "size",
                     Group = 1,
                     fontSize = 30, # Font size
                     linkDistance = networkD3::JS("function(d) { return Math.abs(d.value)/0.004;}"),
                     linkWidth = networkD3::JS("function(d) { return 15*Math.abs(d.value); }"),
                     opacity = input$opacity)
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
                     Group = 1,
                     fontSize = 30, # Font size
                     linkDistance = networkD3::JS("function(d) { return Math.abs(d.value)/0.004;}"),
                     linkWidth = networkD3::JS("function(d) { return 15*Math.abs(d.value); }"),
                     opacity = input$opacity)
    })
    
    
    
    
    output$simple <- renderSimpleNetwork({
        m<- dataSort()
        src <- m$Gene
        target <- m$Gene2
        networkData <- data.frame(src, target)
        simpleNetwork(networkData, opacity = input$opacity)
          
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