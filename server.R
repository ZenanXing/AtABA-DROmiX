# Define server logic
function(input, output, session) {
  
# Gene List ---------------------------------------------------------------
  data <- eventReactive(input$upldData_Butn,
                        {if (input$input_select == "upld") {
                            req(input$file1)
                            inFile <- input$file1
                            if (grepl("csv", inFile$datapath)){
                              data <- read.table(inFile$datapath, sep = ",", header = FALSE)
                            } else {
                              if (grepl("tsv", inFile$datapath)){
                                data <- read.table(inFile$datapath, sep = "\t", header = FALSE)
                              } else {
                                data <- read.xlsx(inFile$datapath, colNames = FALSE)
                              }
                            }
                          } else {
                            if (input$input_select == "slt") {
                              df_all <- readRDS("data/ED_related.rds") %>% dplyr::select(AGI, Cluster, Trend, Sensitivity)
                              if (input$class_tp == "clst") {
                                data <- df_all %>% filter(Cluster == as.integer(input$class_1))
                              } else {
                                data_temp <- df_all %>% filter(Trend == input$class_2)
                                if (input$class_2 == "Bell") {
                                  data <- data_temp %>% filter(Sensitivity == input$class_3_1)
                                } else {
                                  if (input$class_2 == "Up") {
                                    data <- data_temp %>% filter(Sensitivity == input$class_3_2)
                                  } else {
                                    data <- data_temp %>% filter(Sensitivity == input$class_3_3)
                                  }
                                }
                              }
                              data <- data.frame(AGI = data$AGI)
                            } else {
                              req(input$text1)
                              # Input the string from textArea
                              tmp <- matrix(strsplit(input$text1, "\n")[[1]])
                              data <- data.frame(AGI = tmp) 
                            }
                          }
                          return(data)   
                        })
  
# Reactive Variables - data_values ----------------------------------------
  data_values <- reactiveValues()
  observe({
    data_values$n <- nrow(data())
    data_values$agi_list <- data()[, 1]
  })

# Gene_Description --------------------------------------------------------
  df_desc <- eventReactive(input$upldData_Butn, {
    desc <- readRDS("data/ATGeneDescription.rds") %>% filter(AGI %in% data_values$agi_list)
    return(desc)
  })

# Pairwise Comparison -----------------------------------------------------
  df_pc <- eventReactive(input$upldData_Butn, {
    pc <- readRDS("data/PC_ABAvsMock_Col.rds") %>% filter(AGI %in% data_values$agi_list)
    return(pc)
  })
  

# ED-related Table --------------------------------------------------------
  
  df_ed <- eventReactive(input$upldData_Butn, {
    ed <- readRDS("data/ED_related.rds") %>% filter(AGI %in% data_values$agi_list)
    return(ed)
  })
  
# Point Data --------------------------------------------------------------
  df_point <- eventReactive(input$upldData_Butn, {
    df_point <- readRDS("data/DF_Point.rds") %>% filter(AGI %in% data_values$agi_list)
    df_point$ABA_nM[df_point$ABA_nM == 0] <- 0.03
    df_point <- df_point %>% group_by(AGI) %>% nest()
    return(df_point)
  })
  df_point_mean <- eventReactive(input$upldData_Butn, {
    df_point_mean <- df_point() %>% 
      unnest() %>% 
      group_by(AGI, ABA_nM) %>% 
      dplyr::summarise(TPM_Mean = mean(TPM), 
                       TPM_SD = sd(TPM))
    return(df_point_mean)
  })
  
# Curve_Data --------------------------------------------------------------
  df_curve <- eventReactive(input$upldData_Butn, {
    df_curve_4_8 <- readRDS("data/DF_Curve_2.rds")
    df_curve_all <- readRDS("data/DF_Curve_1.rds") %>% rbind(df_curve_4_8)
    df_curve_1 <- df_curve_all %>% filter(AGI %in% data_values$agi_list) %>% unnest()
    df_curve_2 <- df_point_mean() %>% filter(AGI %in% setdiff(data_values$agi_list, df_curve_1$AGI)) %>% 
      mutate(Lower = NA, 
             Upper = NA,
             Prediction = TPM_Mean) %>% 
      dplyr::select(AGI, Prediction, Lower, Upper, ABA_nM)
    if (!is.null(df_curve_1)) {
      df_curve <- rbind(df_curve_1, df_curve_2)
    } else {
      df_curve <- df_curve_2
    }
    df_curve <- df_curve %>% group_by(AGI) %>% nest()
    return(df_curve)
  })
  
# Download the sample list ------------------------------------------------
  output$dl_smp <- downloadHandler(
    filename = function(){"Sample_List.csv"},
    content = function(file) {
      smp <- read.csv("data/Sample_List.csv", header = FALSE,
                      stringsAsFactors = FALSE)
      write.table(smp, file = file, sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
    }
  )
  
  
# Clear Data Button -------------------------------------------------------
  observeEvent(input$clearText_Butn, {
    updateTextAreaInput(session, inputId = "text1", label = "", value = "")
  })
  
# Heatmap -----------------------------------------------------------------
  
  # Dataframe_all
  df_hm_all <- reactive({
    req(nrow(df_pc()) != 0)
    ## Change the label
    lab <- c("0.003", "0.01", "0.03",  "0.09", "0.27", "0.8", 
             "2", "7", "22", "102", "160", "200")
    df_lab <- data.frame(
      ABA_Conc = c(3.33, 9.99, 29.97, 89.90, 269.70, 809.09, 
                   2427.26, 7281.78, 21845.33, 102400.00, 160000.00, 200000.00),
      ABA = factor(lab, levels = lab)
    )
    df_hm_all <- df_pc() %>% 
      mutate(ABA_Conc = round(ABA_nM, 2)) %>% 
      left_join(df_lab, by = "ABA_Conc") %>% 
      dplyr::select(-ABA_Conc)
    return(df_hm_all)
  })
  
  # Exported dataframe
  df_exp_hm <- reactive({
    req(nrow(df_pc()) != 0)
    df_exp_hm <- df_hm_all() %>% 
      dplyr::select(1:4) %>% 
      arrange(factor(Variable, level = c("logFC", "FDR", "logCPM", "F", "PValue"))) %>% 
      pivot_wider(names_from = Variable, values_from = Value) %>% 
      left_join(df_desc(), by = "AGI") %>% 
      dplyr::select(1, 8, 2:7, 9) %>% 
      arrange(factor(AGI, levels = data_values$agi_list), ABA_nM)
    return(df_exp_hm)
  })
  
  # Heatmap
  Hmap <- reactive({
    req(nrow(df_pc()) != 0)
    ## FC dataframe
    df_hm_temp <- df_hm_all() %>% 
      filter(Variable == "logFC") %>% 
      dplyr::select("AGI", "ABA", "Value") %>% 
      spread(ABA, Value) %>% 
      arrange(factor(AGI, levels = data_values$agi_list))
    df_hm_fc <- df_hm_temp[2:13]
    rownames(df_hm_fc) <- df_hm_temp$AGI
    
    ## FDR annotation
    df_hm_temp <- df_hm_all() %>% 
      filter(Variable == "FDR") %>% 
      left_join(df_desc(), by = "AGI") %>% 
      mutate(Lab = paste0("FDR: ", formatC(as.numeric(Value), format = "E", digit = 2), "\n",
                          "Gene: ", tair_symbol)) %>% 
      dplyr::select("AGI", "ABA", "Lab") %>% 
      spread(ABA, Lab) %>% 
      arrange(factor(AGI, levels = data_values$agi_list))
    df_hm_fdr <- df_hm_temp[2:13]
    rownames(df_hm_fdr) <- df_hm_temp$AGI
    
    ## Heatmap
    p <- heatmaply(
      df_hm_fc,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      column_text_angle = 0,
      #grid_color = "white", grid_size = 0.1*5/data_values$n,
      dendrogram = "none",
      seriate = "none",
      main = "ABA vs Mock",
      key.title = "logFC",
      xlab = "ABA (μM)",
      label_names = c("AGI", "ABA (μM)", "logFC (ABAvsMock)"),
      custom_hovertext = df_hm_fdr,
      scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
        low = "#2171B5", 
        high = "#CB181D", 
        midpoint = 0
      )
    ) %>%
      plotly::layout(
        margin = list(l = 50, r = 50, t = 50, b = 100)
      )
    
    return(p)
    
  })
  
  output$heatmap <- renderPlotly({
    return(Hmap())
  })
  
  output$hm <- renderUI({
    req(nrow(df_pc()) != 0)
    tagList(
      # Table
      div(style = "margin-top: 30px"),
      plotlyOutput("heatmap", width = "100%", height = "600px") %>% shinycssloaders::withSpinner(), 
      p(tags$b("Note:"), "Only genes considered as expressed in our experiment are displayed."),
      # Download
      div(style = "margin-top: 20px"),
      h5("Download"),
      div(style = "margin-top: -10px"),
      hr(),
      div(style = "margin-top: -10px"),
      div(style = "vertical-align: top;", 
          textInput(inputId = "file_name_hm", label = "Enter a file name: ", value = Sys.time())
      ),
      div(style = "vertical-align: top;", 
          selectInput(inputId = "file_type_hm", label = "Select file type for the dataframe:", 
                      choices = list("EXCEL", "CSV", "TSV", "TXT"), selected = "EXCEL")
      ),
      div(),
      div(style = "display: inline-block; vertical-align: top; width: 200px;",
          downloadButton(outputId = "dl_df_hm", label = "Download Dataframe")),
      div(style = "display: inline-block; vertical-align: top; width: 200px;",
          downloadButton(outputId = "dl_hm", label = "Download Heatmap")),
      div(style = "margin-top: 20px")
    )
  })
  
  # Download the dataframe
  output$dl_df_hm <- downloadHandler(
    filename = function() {
      if (input$file_type_hm == "EXCEL") { ext <- ".xlsx" } else { ext <- paste0(".", tolower(input$file_type_hm))}
      paste0(input$file_name_hm, ext)
    },
    
    content = function(file) {
      if (input$file_type_hm == "EXCEL") {
        openxlsx::write.xlsx(df_exp_hm(), file)
      } else {
        sep <- switch(input$file_type_hm, "TXT" = " ", "CSV" = ",", "TSV" = "\t" )
        write.table(x = df_exp_hm(), file = file, sep = sep, quote = FALSE, row.names = FALSE)
      }
    }
  )
  
  # Download the heatmap
  output$dl_hm <- downloadHandler(
    filename = function() {
      paste0(input$file_name_hm, ".html")
    },
    content = function(file) {
      htmlwidgets::saveWidget(Hmap(), file, selfcontained = TRUE)
    }
  )
  
# Sensitivity -------------------------------------------------------------

## Reactive Tables ---------------------------------------------------------
  
  ## exported data frame
  df_ed_exp <- reactive({
    req(df_ed())
    df_temp <- df_desc() %>% 
      dplyr::select(AGI, tair_symbol) %>% 
      left_join(df_ed(), by = "AGI") %>% 
      arrange(factor(AGI, levels = data_values$agi_list))
    colnames(df_temp) <- c("AGI", "TAIR_Symbol", "Cluster", "Trend", "Sensitivity", "Membership", "Maximum_Response", "Minimum_Response",
                           "Response_at_BMD", "BMD", "BMD_LowerBound", "BMD_UpperBound", 
                           "Response_at_Low_ED50", "Low_ED50", "Low_ED50_LowerBound", "Low_ED50_UpperBound", 
                           "Response_at_High_ED50", "High_ED50", "High_ED50_LowerBound", "High_ED50_UpperBound", 
                           "Response_at_LDS", "LDS", "LDS_LowerBound", "LDS_UpperBound", 
                           "Response_at_M", "M", 
                           "Model", "Neill's test", "No effet test")
    return(df_temp)
  })
  
  ## tables for display
  # low ED50s
  ED50_table <- reactive({
    
    req(df_ed_exp())
    df_temp <- df_ed_exp() %>% 
      dplyr::select(1:5, 14:16, 10:12) %>% 
      filter(!is.na(Low_ED50)|!is.na(BMD))
    colnames(df_temp) <- c("AGI", "TAIR_Symbol", "Cluster", "Trend", "Sensitivity", 
                           "ED\u2085\u2080", "ED\u2085\u2080\nLowerBound", "ED\u2085\u2080\nUpperBound", 
                           "BMD", "BMD\nLowerBound", "BMD\nUpperBound")
    df_temp <- df_temp %>% mutate(across(6:ncol(df_temp), ~ map_chr(.x, display_format)))
    return(df_temp)
    
  })
  
## Output UI ---------------------------------------------------------------
  
  output$tb_ed <- DT::renderDataTable({
    DT::datatable(
      ED50_table(),
      escape = FALSE,
      options = list(pageLength = 10,
                     lengthMenu = c(5, 10, 15, 20), scrollX = T)
    )
  })
  
  output$ed_table <- renderUI({
    req(df_ed())
    tagList(
      # Table
      div(style = "margin-top: 30px"),
      h5(HTML(paste0("ED", tags$sub("50"), " & BMD Estimation Table")), align = 'center'),
      div(style = "margin-top: -10px"),
      #hr(),
      #div(style = "margin-top: -10px"),
      DT::dataTableOutput("tb_ed") %>% shinycssloaders::withSpinner(),
      p(tags$b("Note:"), HTML(paste0("Only genes with valid ED", tags$sub("50"), " or BMD values are displayed."))),
      # Download
      div(style = "margin-top: 20px"),
      h5("Download"),
      div(style = "margin-top: -10px"),
      hr(),
      div(style = "margin-top: -10px"),
      div(style = "display: inline-block; vertical-align: top;", 
          textInput(inputId = "file_name", label = "Enter a file name: ", value = Sys.time())
      ),
      div(style = "display: inline-block; vertical-align: top; width: 150px;", 
          selectInput(inputId = "file_type", label = "Select file type:", 
                      choices = list("EXCEL", "CSV", "TSV", "TXT"), selected = "EXCEL")
      ),
      div(style = "width: 150px;",
          downloadButton(outputId = "dl_df", label = "Download")),
      div(style = "margin-top: 10px"), 
      tags$b("References:"),
      p(em("Ritz C, Baty F, Streibig JC, Gerhard D (2015) Dose-Response Analysis Using R. PLoS One. 10(12)")), 
      div(style = "margin-top: -15px"), 
      p(em("Serra A. Et al. (2020) BMDx: a graphical Shiny application to perform Benchmark Dose analysis for transcriptomics data. Bioinformatics 36: 2932–2933"))
    )
  })
  
  output$dl_df <- downloadHandler(
    filename = function() {
      if (input$file_type == "EXCEL") { ext <- ".xlsx" } else { ext <- paste0(".", tolower(input$file_type))}
      paste0(input$file_name, ext)
    },
    
    content = function(file) {
      if (input$file_type == "EXCEL") {
        openxlsx::write.xlsx(df_ed_exp(), file)
      } else {
        sep <- switch(input$file_type, "TXT" = " ", "CSV" = ",", "TSV" = "\t" )
        write.table(x = df_ed_exp(), file = file, sep = sep, quote = FALSE, row.names = FALSE)
      }
    }
  )
  
# Dose-Response Curves ----------------------------------------------------
  
  # Update the selectize 
  observeEvent(data_values$agi_list, {
    updateSelectizeInput(session, "agi_slt", choices = data_values$agi_list, 
                         selected = head(data_values$agi_list, 1), server = TRUE)
  })
  
  # Selected dataframe
  df_ed_slt <- reactive({
    req(input$agi_slt)
    agi_slt <- isolate({unname(input$agi_slt)})
    df_ed_slt <- df_ed_exp() %>% filter(AGI %in% agi_slt)
  })
  
  # Biphasic LDS & M
  output$lds_m <- renderUI({
    req(df_ed_exp())
    if (any(df_ed_exp()$Cluster == 1) && all(!is.na(df_ed_exp()$Cluster))) {
      tagList(
        div(), 
        div(style = "display: inline-block; vertical-align: top; margin-top: -15px;", 
            checkboxInput(inputId = "plot_lds_m_ck", label = "LDS & M", value = FALSE)
        )
      )
    }
  })
  
## Demo Plot ---------------------------------------------------------------
  
  demo_p <- reactive({
    req(input$agi_slt, data())
    agi_slt <- unname(input$agi_slt)
    df_curve_slt <- df_curve() %>% filter(AGI %in% agi_slt) %>% unnest()
    df_point_slt <- df_point() %>% filter(AGI %in% agi_slt) %>% unnest()
    df_point_mean_slt <- df_point_mean() %>% filter(AGI %in% agi_slt)
    df_ed_slt <- df_ed_slt()
    
    df_curve_slt$AGI <- factor(df_curve_slt$AGI, levels = agi_slt)
    df_point_slt$AGI <- factor(df_point_slt$AGI, levels = agi_slt)
    df_point_mean_slt$AGI <- factor(df_point_mean_slt$AGI, levels = agi_slt)
    df_ed_slt$AGI <- factor(df_ed_slt$AGI, levels = agi_slt)
    
    p_slt <- ggplot() + 
      # the major plot
      geom_point(data = df_point_slt, aes(x = ABA_nM, y = TPM, 
                                          text = paste("ABA [nM]:", round(ABA_nM, 2), "<br>Expression (TPM):", round(TPM, 2))), 
                 color = "#D55E00", alpha = 1/3) + 
      geom_line(data = df_curve_slt, aes(x = ABA_nM, y = Prediction), color = "#D55E00") +
      facet_wrap(. ~ ordered(AGI), ncol = 2, scales = "free_y") +
      scale_x_log10() + 
      labs(x = "ABA [nM]",  y = "Expression (TPM)") + 
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) + 
      theme(axis.title = element_text(size = 14),
            axis.text = element_text(size = 10),
            strip.text.x = element_text(size = 12),
            legend.title = element_text(size = 14), 
            legend.text = element_text(size = 10)) +
      theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "in"))
    
    # ed50 lines
    if (input$plot_ed50_ck == TRUE) {
      p_slt <- p_slt +
        geom_hline(data = df_ed_slt, aes(yintercept = Response_at_Low_ED50), linetype = "longdash", color = "#0072B2", alpha = 0.5) + 
        geom_vline(data = df_ed_slt, aes(xintercept = Low_ED50, 
                                         text = paste("ED50:", round(Low_ED50, 2),"<br>Response_at_ED50:", round(Response_at_Low_ED50, 2))),
                   linetype = "longdash", color = "#0072B2", alpha = 0.5)
      if (input$plot_ci_ck == TRUE) {
        p_slt <- p_slt + 
          geom_vline(data = df_ed_slt, aes(xintercept = Low_ED50_LowerBound), linetype = "dotted", color = "#0072B2", alpha = 0.5) + 
          geom_vline(data = df_ed_slt, aes(xintercept = Low_ED50_UpperBound), linetype = "dotted", color = "#0072B2", alpha = 0.5)
      }
    }
    
    # bmd lines
    if (input$plot_bmd_ck == TRUE) {
      p_slt <- p_slt +
        geom_hline(data = df_ed_slt, aes(yintercept = Response_at_BMD), linetype = "dashed", color = "#56B4E9", alpha = 0.5) + 
        geom_vline(data = df_ed_slt, aes(xintercept = BMD, 
                                         text = paste("BMD:", round(BMD, 2),"<br>Response_at_BMD:", round(Response_at_BMD, 2))), 
                   linetype = "dashed", color = "#56B4E9", alpha = 0.5)
      if (input$plot_ci_ck == TRUE) {
        p_slt <- p_slt +
          geom_vline(data = df_ed_slt, aes(xintercept = BMD_LowerBound), linetype = "dotted", color = "#56B4E9", alpha = 0.5) + 
          geom_vline(data = df_ed_slt, aes(xintercept = BMD_UpperBound), linetype = "dotted", color = "#56B4E9", alpha = 0.5)
      }
    }
    
    ## max & min response
    if (input$plot_resline_ck == TRUE) {
      p_slt <- p_slt +
        geom_hline(data = df_ed_slt, aes(yintercept = Maximum_Response, 
                                         text = paste("Max:", round(Maximum_Response, 2))), linetype = "dotted", color = "#009E73", alpha = 0.5) + 
        geom_hline(data = df_ed_slt, aes(yintercept = Minimum_Response, 
                                         text = paste("Min:", round(Minimum_Response, 2))), linetype = "dotted", color = "#009E73", alpha = 0.5)
    }
    
    # lds & m
    if (any(df_ed_slt$Cluster == 1) && all(!is.na(df_ed_slt$Cluster))) {
      if (input$plot_lds_m_ck == TRUE) {
        p_slt <- p_slt +
          # lds lines
          geom_hline(data = df_ed_slt, aes(yintercept = Response_at_LDS), linetype = "dashed", color = "#009E73", alpha = 0.5) + 
          geom_vline(data = df_ed_slt, aes(xintercept = LDS, 
                                           text = paste("LDS:", round(LDS, 2),"<br>Response_at_LDS:", round(Response_at_LDS, 2))),
                     linetype = "dashed", color = "#009E73", alpha = 0.5) + 
          # m lines
          geom_hline(data = df_ed_slt, aes(yintercept = Response_at_M), linetype = "dashed", color = "gold2", alpha = 0.5) + 
          geom_vline(data = df_ed_slt, aes(xintercept = M, 
                                           text = paste("M:", round(M, 2),"<br>Response_at_M:", round(Response_at_M, 2))), 
                     linetype = "dashed", color = "gold2", alpha = 0.5)
        
        if (input$plot_ci_ck == TRUE) {
          p_slt <- p_slt +
            geom_vline(data = df_ed_slt, aes(xintercept = LDS_LowerBound), linetype = "dotted", color = "#009E73", alpha = 0.5) + 
            geom_vline(data = df_ed_slt, aes(xintercept = LDS_UpperBound), linetype = "dotted", color = "#009E73", alpha = 0.5)
        }
      }
    }
    
    ## plotly
    p <- ggplotly(p_slt, tooltip = "text")
    
    return(p)
    
  })
  
  output$dr_curve <- renderPlotly({
    return(demo_p())
  })

## Download pdf ------------------------------------------------------------
  
  output$dl_pdf <- downloadHandler(
    
    filename = function() {
      paste0(input$pdf_name, ".pdf")
    },
    
    content = function(file) {
      
      df_ed_exp_val <- df_ed_exp() %>%
        arrange(factor(AGI, levels = data_values$agi_list))
      df_point_val <- df_point() %>%
        arrange(factor(AGI, levels = data_values$agi_list))
      df_curve_val <- df_curve() %>%
        arrange(factor(AGI, levels = data_values$agi_list))
      
      # Open the PDF device using 'file' as the destination
      pdf(file = file, width = 1000/72, height = 1250/72)
      
      N <- 8
      agi_list_val <- data_values$agi_list
      mx <- ceiling(length(agi_list_val)/N)
      
      for (j in seq_len(mx)) {
        a <- (j - 1) * N + 1
        b <- if (j < mx) j*N else length(agi_list_val)
        n <- b - a + 1
        
        # Subset data for pages
        df_point_slt <- df_point_val[a:b, ] %>% unnest() 
        df_curve_slt <- df_curve_val[a:b, ] %>% unnest() 
        df_ed_slt<- df_ed_exp_val[a:b, ]
        
        df_curve_slt$AGI <- factor(df_curve_slt$AGI, levels = unique(df_curve_slt$AGI))
        df_point_slt$AGI <- factor(df_point_slt$AGI, levels = unique(df_point_slt$AGI))
        df_ed_slt$AGI <- factor(df_ed_slt$AGI, levels = unique(df_ed_slt$AGI))
        
        # plot
        x <- ceiling(n/2)
        p_slt <- ggplot() + 
          # the major plot
          geom_point(data = df_point_slt, aes(x = ABA_nM, y = TPM, 
                                              text = paste("ABA [nM]:", round(ABA_nM, 2), "<br>Expression (TPM):", round(TPM, 2))), 
                     color = "#D55E00", alpha = 1/3) + 
          geom_line(data = df_curve_slt, aes(x = ABA_nM, y = Prediction), color = "#D55E00") +
          facet_wrap(. ~ ordered(AGI), ncol = 2, scales = "free_y") +
          scale_x_log10() + 
          labs(x = "ABA [nM]",  y = "Expression (TPM)") + 
          theme_bw() +
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()) + 
          theme(axis.title = element_text(size = 14),
                axis.text = element_text(size = 10),
                strip.text.x = element_text(size = 12),
                legend.title = element_text(size = 14), 
                legend.text = element_text(size = 10)) +
          theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "in"))
        # ed50 lines
        if (input$plot_ed50_ck == TRUE) {
          p_slt <- p_slt +
            geom_hline(data = df_ed_slt, aes(yintercept = Response_at_Low_ED50), linetype = "longdash", color = "#0072B2", alpha = 0.5) + 
            geom_vline(data = df_ed_slt, aes(xintercept = Low_ED50, 
                                             text = paste("ED50:", round(Low_ED50, 2),"<br>Response_at_ED50:", round(Response_at_Low_ED50, 2))),
                       linetype = "longdash", color = "#0072B2", alpha = 0.5)
          if (input$plot_ci_ck == TRUE) {
            p_slt <- p_slt + 
              geom_vline(data = df_ed_slt, aes(xintercept = Low_ED50_LowerBound), linetype = "dotted", color = "#0072B2", alpha = 0.5) + 
              geom_vline(data = df_ed_slt, aes(xintercept = Low_ED50_UpperBound), linetype = "dotted", color = "#0072B2", alpha = 0.5)
          }
        }
        
        # bmd lines
        if (input$plot_bmd_ck == TRUE) {
          p_slt <- p_slt +
            geom_hline(data = df_ed_slt, aes(yintercept = Response_at_BMD), linetype = "dashed", color = "#56B4E9", alpha = 0.5) + 
            geom_vline(data = df_ed_slt, aes(xintercept = BMD, 
                                             text = paste("BMD:", round(BMD, 2),"<br>Response_at_BMD:", round(Response_at_BMD, 2))), 
                       linetype = "dashed", color = "#56B4E9", alpha = 0.5)
          if (input$plot_ci_ck == TRUE) {
            p_slt <- p_slt +
              geom_vline(data = df_ed_slt, aes(xintercept = BMD_LowerBound), linetype = "dotted", color = "#56B4E9", alpha = 0.5) + 
              geom_vline(data = df_ed_slt, aes(xintercept = BMD_UpperBound), linetype = "dotted", color = "#56B4E9", alpha = 0.5)
          }
        }
        
        ## max & min response
        if (input$plot_resline_ck == TRUE) {
          p_slt <- p_slt +
            geom_hline(data = df_ed_slt, aes(yintercept = Maximum_Response, 
                                             text = paste("Max:", round(Maximum_Response, 2))), linetype = "dotted", color = "#009E73", alpha = 0.5) + 
            geom_hline(data = df_ed_slt, aes(yintercept = Minimum_Response, 
                                             text = paste("Min:", round(Minimum_Response, 2))), linetype = "dotted", color = "#009E73", alpha = 0.5)
        }
        
        # lds & m
        if (any(df_ed_slt$Cluster == 1) && all(!is.na(df_ed_slt$Cluster))) {
          if (input$plot_lds_m_ck == TRUE) {
            p_slt <- p_slt +
              # lds lines
              geom_hline(data = df_ed_slt, aes(yintercept = Response_at_LDS), linetype = "dashed", color = "#009E73", alpha = 0.5) + 
              geom_vline(data = df_ed_slt, aes(xintercept = LDS, 
                                               text = paste("LDS:", round(LDS, 2),"<br>Response_at_LDS:", round(Response_at_LDS, 2))),
                         linetype = "dashed", color = "#009E73", alpha = 0.5) + 
              # m lines
              geom_hline(data = df_ed_slt, aes(yintercept = Response_at_M), linetype = "dashed", color = "gold2", alpha = 0.5) + 
              geom_vline(data = df_ed_slt, aes(xintercept = M, 
                                               text = paste("M:", round(M, 2),"<br>Response_at_M:", round(Response_at_M, 2))), 
                         linetype = "dashed", color = "gold2", alpha = 0.5)
            
            if (input$plot_ci_ck == TRUE) {
              p_slt <- p_slt +
                geom_vline(data = df_ed_slt, aes(xintercept = LDS_LowerBound), linetype = "dotted", color = "#009E73", alpha = 0.5) + 
                geom_vline(data = df_ed_slt, aes(xintercept = LDS_UpperBound), linetype = "dotted", color = "#009E73", alpha = 0.5)
            }
          }
        }
        
        # gene info.
        if (input$gene_info == TRUE){
          p_slt <- p_slt +
            ggpp::geom_text_npc(data = df_ed_slt, aes(npcx = "left", npcy = "top", label = TAIR_Symbol), size = 6, color = "black")
        }
        
        # Print the figure(s)
        if (n == 1) {
          print(
            ggarrange(p_slt, NULL, NULL,
                      ncol = 2, nrow = 2,
                      widths = c(7, 4.8),
                      heights = c(x/4 + 1/28, (4 - x)/4))
          )
        } else {
          print(
            ggarrange(p_slt, NULL,
                      ncol = 1, nrow = 2,
                      heights = c(x/4 + 1/28, (4 - x)/4))
          )
        }
      }
      
      dev.off()
    }
  )
  
}
