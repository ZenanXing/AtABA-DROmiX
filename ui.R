
# Load the required packages ----------------------------------------------

library(shiny)
library(shinythemes)
library(shinycssloaders)
library(shinyjs)
library(bslib)
library(openxlsx)
library(tidyverse)
library(heatmaply)
library(plotly)
library(htmlwidgets)
library(DT)
library(ggpp)
library(ggpubr)

# Define UI for application that draws a histogram ------------------------

navbarPage(
  
  title = HTML("<b>AtABA-DROmiX</b> <small>- <b>A</b>rabidopsis <b>t</b>haliana <b>ABA</b> <b>D</b>ose-<b>R</b>esponse Transcript<b>omi</b>cs E<b>x</b>plorer</small>"),
  
  # Theme from bslib
  theme = bs_theme(bootswatch = "cerulean"), 
  useShinyjs(),
  
  # SidebarLayout
  tabPanel(title = "",
    sidebarLayout(
      # Siderbar
      sidebarPanel(
        h5("Upload the Gene List:"),
        
        # Select the way to input list
        div(style = "margin-top: -20px;"), 
        selectInput(inputId = "input_select",
                    label = "",
                    choices = c("Paste the list" = "pst",
                                "Upload the list" = "upld"),
                    selected = "pst"),
        
        
        # Paste the list
        conditionalPanel(condition = "input.input_select == 'pst'", 
                         
                         # TextArea
                         tags$b(h6("Paste list below:")), 
                         div(style = "margin-top: -20px;"), 
                         textAreaInput(inputId = "text1", 
                                       value = "AT3G11410\nAT5G52310\nAT3G61430\nAT2G38310", 
                                       height = "160px", label = ""), 
                         
                         # Note
                         p(style = "text-align: left; margin-top: 10px;", 
                           tags$b("Note: "), "Please enter AGI locus codes (Arabidopsis Genome Initiative, e.g., AT5G52310) on separate lines. The examples provided above in the textbox include four well-known ABA-responsive genes."),
                         
                         # Clear the list
                         actionButton(inputId = "clearText_Butn", 
                                      label = "Clear list")
                         
        ),
        
        # Upload the list
        conditionalPanel(condition = "input.input_select == 'upld'", 
                         tags$b(h6("Choose the File:")), 
                         div(style = "margin-top: -20px;"), 
                         fileInput(inputId = "file1", label = ""), 
                         div(style = "margin-top: 20px;"), 
                         div(style = "text-align: right;", 
                             downloadButton(outputId = "dl_smp", label = "Download Sample List")), 
                         p(style = "text-align: left; margin-top: 10px;", 
                           tags$b("Note: "), "Please provide the AGI locus code (AGI, Arabidopsis Genome Initiative, e.g. AT5G52310). You may upload files in TSV, CSV, or Excel formats. The sample list includes four well-known ABA-responsive genes.")
        ), 
        
        # Confirm the upload
        div(style = "text-align: right; margin-top: 10px;", 
            actionButton(inputId = "upldData_Butn", 
                         label = "Go"))
      ),
      
      # Main
      mainPanel(
        tabsetPanel(
          tabPanel(title = "Heatmap",
                   uiOutput(outputId = "hm")
          ),
          tabPanel(title = "Sensitivity",
                   uiOutput(outputId = "ed_table")
          ),
          tabPanel(title = "Dose-Response Curves",
                   # Demo
                   div(style = "margin-top: 20px"),
                   h5("Demo"),
                   #div(style = "margin-top: -10px"),
                   #hr(),
                   div(style = "margin-top: 20px"),
                   selectizeInput(inputId = "agi_slt",
                                  label = "Select the genes (≤ 2):",
                                  choices = NULL, multiple = TRUE,
                                  options = list(maxItems = 2)),
                   div(style = "margin-top: 20px;"),  
                   "Show the ED-related values and the corresponding responses：", 
                   div(),
                   div(style = "vertical-align: top;", 
                       checkboxInput(inputId = "plot_ed50_ck", label = HTML(paste0("ED", tags$sub("50"))), value = FALSE)), 
                   div(style = "vertical-align: top; margin-top: -15px;", 
                       checkboxInput(inputId = "plot_bmd_ck", label = HTML(paste0("BMD", tags$sub("1SD"))), value = FALSE)),
                   uiOutput(outputId = "lds_m"), 
                   div(), 
                   div(style = "display: inline-block; vertical-align: top; margin-top: -15px;", 
                       checkboxInput(inputId = "plot_resline_ck", label = "Max & Min Responses", value = FALSE)), 
                   div(style = "display: inline-block; vertical-align: top; margin-top: -15px;", 
                       checkboxInput(inputId = "plot_ci_ck", label = "Confidence Intervals", value = FALSE)),
                   plotlyOutput("dr_curve", width = "100%", height = "400px") %>% shinycssloaders::withSpinner(), 
                   
                   # Download
                   div(style = "margin-top: 20px"),
                   h5("Download"),
                   div(style = "margin-top: -10px"),
                   hr(),
                   "Show Info. of the genes：", 
                   div(style = "display: inline-block; vertical-align: top;", 
                       checkboxInput(inputId = "gene_info", label = "Gene Info.", value = FALSE)),
                   div(), 
                   div(style = "display: inline-block; vertical-align: top;", 
                       textInput(inputId = "pdf_name", label = "Enter a file name: ", value = Sys.time())
                   ),
                   div(style = "width: 150px;",
                       downloadButton(outputId = "dl_pdf", label = "Download")),
                   div(style = "margin-top: 20px")
          )
        )
      )
    )
  )
)
