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
df <- read.csv("autophagy.txt", header = TRUE,stringsAsFactors = FALSE)

default_gene = c('ATG14','ATG101','ULK1')


# Define UI for dataset viewer app ----
ui <- fluidPage(
    
    # App title ----
    titlePanel("Autophagy Gene Correlation Network Analysis"),

    # Sidebar layout with input and output definitions ----
    sidebarLayout(
        
        # Sidebar panel for inputs ----
        sidebarPanel(
            
            #Input Panel Title
            titlePanel("Title"),
            hr(),
            helpText("helpText"),
            
           
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
                       
                       sliderInput("rank", "Rank:",
                                   min = 0, max = 50,
                                   value = 10)
                      # selectInput("rank", label = "Rank", 
                                   #choices=c(1:1000), selected=3,multiple=FALSE)
                       
                )
            ),
            
            br(),
           
            fluidRow(
                column(8,
                       sliderInput("opacity", "Opacity for network", 0.6, min = 0.1,
                                   max = 1, step = .1)
                 
                       
                )
            )
            
     
            
        ),
        
        mainPanel(
            fluidRow(
                column(11,offset = 1,
                       helpText("For Heatmap: NaN values (missing values) are represented as empty boxes."),
                       helpText("“To generate new violin plots plots without changing the heatmaps, click ‘Keep All Heatmaps”")
                )
            ),
            
            fluidRow(
                column(12,align="center",
                       tabPanel("Simple Network", simpleNetworkOutput("simple"))
                )
            ),
            
            fluidRow(
                column(11,offset = 1,
                       tableOutput('table')
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
        
    })
    
    output$table <- renderTable(dataSort())
    
    output$simple <- renderSimpleNetwork({
        m<- dataSort()
        src <- m$Gene
        target <- m$Gene2
        networkData <- data.frame(src, target)
        simpleNetwork(networkData, opacity = input$opacity)
           # saveNetwork(file = 'Net1.html')
    })
    
 
   
}

# Create Shiny app ----
shinyApp(ui, server)