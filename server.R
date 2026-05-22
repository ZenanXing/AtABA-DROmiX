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
    # genotype
    df_geno <- data.frame(
      ckgrp = c("col", "i", "ii", "iii"),
      Genotype = c("Col-0", "sfkoI", "sfkoII", "sfkoIII")
    )
    df_geno_slt <- data.frame(
      ckgrp = input$genotype
    ) %>% left_join(df_geno, by = "ckgrp")
    data_values$genotype_list <- df_geno_slt$Genotype
    
  })

# Gene_Description --------------------------------------------------------
  df_desc <- eventReactive(input$upldData_Butn, {
    desc <- readRDS("data/ATGeneDescription.rds") %>% filter(AGI %in% data_values$agi_list)
    return(desc)
  })

# Pairwise Comparison -----------------------------------------------------
  
  # PC - ABA vs Mock
  df_pc_aba <- eventReactive(input$upldData_Butn, {
    pc <- readRDS("data/PC_ABAvsMock.rds") %>% filter(AGI %in% data_values$agi_list)
    return(pc)
  })
  # PC - Mutant vs Col-0
  df_pc_mut <- eventReactive(input$upldData_Butn, {
    pc <- readRDS("data/PC_MutantvsWildtype.rds") %>% filter(AGI %in% data_values$agi_list)
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
    df_point <- df_point %>% group_by(Genotype, AGI) %>% nest()
    return(df_point)
  })
  df_point_mean <- eventReactive(input$upldData_Butn, {
    df_point_mean <- df_point() %>% 
      unnest() %>% 
      group_by(Genotype, AGI, ABA_nM) %>% 
      dplyr::summarise(TPM_Mean = mean(TPM), 
                       TPM_SD = sd(TPM))
    return(df_point_mean)
  })
  
# Curve_Data --------------------------------------------------------------
  
  # wildtype - col-0
  df_curve_wt <- eventReactive(input$upldData_Butn, {
    df_curve_4_8 <- readRDS("data/DF_Curve_wt_2.rds")
    df_curve_all <- readRDS("data/DF_Curve_wt_1.rds") %>% rbind(df_curve_4_8)
    df_curve_1 <- df_curve_all %>% filter(AGI %in% data_values$agi_list) %>% unnest()
    df_curve_2 <- df_point_mean() %>% filter(AGI %in% setdiff(data_values$agi_list, df_curve_1$AGI)) %>% 
      mutate(Lower = NA, 
             Upper = NA,
             Prediction = TPM_Mean) %>% 
      dplyr::select(AGI, Genotype, Prediction, Lower, Upper, ABA_nM)
    if (!is.null(df_curve_1)) {
      df_curve <- rbind(df_curve_1, df_curve_2)
    } else {
      df_curve <- df_curve_2
    }
    df_curve <- df_curve %>% group_by(Genotype, AGI) %>% nest()
    return(df_curve)
  })
  
  # sfkoI
  df_curve_sfko1 <- eventReactive(input$upldData_Butn, {
    df_curve_4_8 <- readRDS("data/DF_Curve_sfkoI_2.rds")
    df_curve_all <- readRDS("data/DF_Curve_sfkoI_1.rds") %>% rbind(df_curve_4_8)
    df_curve_1 <- df_curve_all %>% filter(AGI %in% data_values$agi_list) %>% unnest()
    df_curve_2 <- df_point_mean() %>% filter(AGI %in% setdiff(data_values$agi_list, df_curve_1$AGI)) %>% 
      mutate(Lower = NA, 
             Upper = NA,
             Prediction = TPM_Mean) %>% 
      dplyr::select(AGI, Genotype, Prediction, Lower, Upper, ABA_nM)
    if (!is.null(df_curve_1)) {
      df_curve <- rbind(df_curve_1, df_curve_2)
    } else {
      df_curve <- df_curve_2
    }
    df_curve <- df_curve %>% group_by(Genotype, AGI) %>% nest()
    return(df_curve)
  })
  
  # sfkoII
  df_curve_sfko2 <- eventReactive(input$upldData_Butn, {
    df_curve_4_8 <- readRDS("data/DF_Curve_sfkoII_2.rds")
    df_curve_all <- readRDS("data/DF_Curve_sfkoII_1.rds") %>% rbind(df_curve_4_8)
    df_curve_1 <- df_curve_all %>% filter(AGI %in% data_values$agi_list) %>% unnest()
    df_curve_2 <- df_point_mean() %>% filter(AGI %in% setdiff(data_values$agi_list, df_curve_1$AGI)) %>% 
      mutate(Lower = NA, 
             Upper = NA,
             Prediction = TPM_Mean) %>% 
      dplyr::select(AGI, Genotype, Prediction, Lower, Upper, ABA_nM)
    if (!is.null(df_curve_1)) {
      df_curve <- rbind(df_curve_1, df_curve_2)
    } else {
      df_curve <- df_curve_2
    }
    df_curve <- df_curve %>% group_by(Genotype, AGI) %>% nest()
    return(df_curve)
  })
  
  # sfkoIII
  df_curve_sfko3 <- eventReactive(input$upldData_Butn, {
    df_curve_4_8 <- readRDS("data/DF_Curve_sfkoIII_2.rds")
    df_curve_all <- readRDS("data/DF_Curve_sfkoIII_1.rds") %>% rbind(df_curve_4_8)
    df_curve_1 <- df_curve_all %>% filter(AGI %in% data_values$agi_list) %>% unnest()
    df_curve_2 <- df_point_mean() %>% filter(AGI %in% setdiff(data_values$agi_list, df_curve_1$AGI)) %>% 
      mutate(Lower = NA, 
             Upper = NA,
             Prediction = TPM_Mean) %>% 
      dplyr::select(AGI, Genotype, Prediction, Lower, Upper, ABA_nM)
    if (!is.null(df_curve_1)) {
      df_curve <- rbind(df_curve_1, df_curve_2)
    } else {
      df_curve <- df_curve_2
    }
    df_curve <- df_curve %>% group_by(Genotype, AGI) %>% nest()
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
  
# Heatmap - ABA vs Mock -----------------------------------------------------------------
  
  # Dataframe_all
  df_hm_aba_all <- reactive({
    req(nrow(df_pc_aba()) != 0)
    ## Change the label
    lab <- c("0.003", "0.01", "0.03",  "0.09", "0.27", "0.8", 
             "2", "7", "22", "102", "160", "200")
    df_lab <- data.frame(
      ABA_Conc = c(3.33, 9.99, 29.97, 89.90, 269.70, 809.09, 
                   2427.26, 7281.78, 21845.33, 102400.00, 160000.00, 200000.00),
      ABA = factor(lab, levels = lab)
    )
    df_hm_aba_all <- df_pc_aba() %>% 
      filter(Genotype %in% data_values$genotype_list) %>% 
      mutate(ABA_Conc = round(ABA_nM, 2)) %>% 
      left_join(df_lab, by = "ABA_Conc") %>% 
      dplyr::select(-ABA_Conc)
    return(df_hm_aba_all)
  })
  
  # Exported dataframe
  df_exp_aba_hm <- reactive({
    req(nrow(df_pc_aba()) != 0)
    
    # exported data frame
    df_exp_hm <- df_hm_aba_all() %>% 
      dplyr::select(AGI, Genotype, ABA_nM, Variable, Value) %>% 
      arrange(factor(Variable, level = c("logFC", "FDR", "logCPM", "F", "PValue"))) %>% 
      pivot_wider(names_from = Variable, values_from = Value) %>% 
      left_join(df_desc(), by = "AGI") %>% 
      mutate(Comparison = "ABA vs Mock") %>% 
      dplyr::select(Comparison, AGI, tair_symbol, Genotype, ABA_nM, logFC, FDR, logCPM, F, PValue, entrezgene_description) %>% 
      arrange(factor(AGI, levels = data_values$agi_list), factor(Genotype, levels = data_values$genotype_list), ABA_nM)
    
    return(df_exp_hm)
  })
  
  # Heatmap
  Hmap_aba <- reactive({
    req(nrow(df_pc_aba()) != 0)
    ## FC dataframe
    df_hm_temp <- df_hm_aba_all() %>% 
      filter(Variable == "logFC") %>% 
      unite(Genotype, ABA, col = "comb", sep = ", ABA(μM): ", remove = FALSE) %>% 
      arrange(factor(AGI, levels = data_values$agi_list), factor(Genotype, levels = data_values$genotype_list), ABA) %>% 
      dplyr::select(AGI, comb, Value) %>% 
      pivot_wider(names_from = comb, values_from = Value)
    df_hm_fc <- df_hm_temp[2:ncol(df_hm_temp)]
    rownames(df_hm_fc) <- df_hm_temp$AGI
    
    ## FDR annotation
    df_hm_temp <- df_hm_aba_all() %>% 
      filter(Variable == "FDR") %>% 
      left_join(df_desc(), by = "AGI") %>% 
      mutate(Lab = paste0("FDR: ", formatC(as.numeric(Value), format = "E", digit = 2), "\n",
                          "Gene: ", tair_symbol)) %>% 
      unite(Genotype, ABA, col = "comb", sep = ", ABA(μM): ", remove = FALSE) %>% 
      arrange(factor(AGI, levels = data_values$agi_list), factor(Genotype, levels = data_values$genotype_list), ABA) %>% 
      dplyr::select(AGI, comb, Lab) %>% 
      pivot_wider(names_from = comb, values_from = Lab)
    df_hm_fdr <- df_hm_temp[2:ncol(df_hm_temp)]
    rownames(df_hm_fdr) <- df_hm_temp$AGI
    
    # Genotype annotation
    df_hm_temp <- df_hm_aba_all() %>% 
      filter(Variable == "FDR") %>% 
      unite(Genotype, ABA, col = "comb", sep = ", ABA(μM): ", remove = FALSE) %>% 
      dplyr::select(comb, Genotype, ABA) %>% unique() %>% 
      arrange(factor(Genotype, levels = data_values$genotype_list), ABA)
    gntp <- factor(df_hm_temp$Genotype, levels = data_values$genotype_list)
    color_map <- c(
      "Col-0" = "#ffac6a", "sfkoI" = "#f2cf80",
      "sfkoII" = "#abd9f4", "sfkoIII" = "#80ceb9"
    )
    color_map_slt <- color_map[data_values$genotype_list]
    gntp_char <- as.character(df_hm_temp$Genotype)
    names(gntp_char) <- "Genotype"
    col_side_palette <- color_map_slt[levels(gntp)]
    col_side_colors_df <- data.frame(Genotype = gntp)
    
    ## Heatmap
    p <- heatmaply(
      df_hm_fc,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      col_side_colors = col_side_colors_df,
      col_side_palette = col_side_palette,
      grid_color = "white", grid_size = 0.1,
      dendrogram = "none",
      seriate = "none",
      key.title = "log2FC",
      xlab = "ABA (Low ⇢ High)",
      main = "ABA vs Mock",
      label_names = c("AGI", "Genotype", "log2FC (ABA vs Mock)"),
      custom_hovertext = df_hm_fdr,
      scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
        low = "#2171B5", 
        high = "#CB181D", 
        midpoint = 0
      )
    ) %>%
      plotly::layout(
        xaxis = list(showticklabels = FALSE, ticks = ""), 
        margin = list(l = 50, r = 50, t = 100, b = 100)
      )
    
    return(p)
    
  })
  
  output$heatmap_aba <- renderPlotly({
    return(Hmap_aba())
  })
  
  output$hm_aba <- renderUI({
    req(nrow(df_pc_aba()) != 0)
    tagList(
      # Heatmap
      div(style = "margin-top: 10px"),
      plotlyOutput("heatmap_aba", width = "100%", height = "600px") %>% shinycssloaders::withSpinner(), 
      p(tags$b("Note:"), "Only genes considered as expressed in our experiment are displayed."),
      # Download
      div(style = "margin-top: 20px"),
      h5("Download"),
      div(style = "margin-top: -10px"),
      hr(),
      div(style = "margin-top: -10px"),
      div(style = "vertical-align: top;", 
          textInput(inputId = "file_name_hm_aba", label = "Enter a file name: ", value = paste0("ABA_vs_Mock_", Sys.Date()))
      ),
      div(style = "vertical-align: top;", 
          selectInput(inputId = "file_type_hm_aba", label = "Select file type for the dataframe:", 
                      choices = list("EXCEL", "CSV", "TSV", "TXT"), selected = "EXCEL")
      ),
      div(),
      div(style = "display: inline-block; vertical-align: top; width: 200px;",
          downloadButton(outputId = "dl_df_hm_aba", label = "Download Dataframe")),
      div(style = "display: inline-block; vertical-align: top; width: 200px;",
          downloadButton(outputId = "dl_hm_aba", label = "Download Heatmap")),
      div(style = "margin-top: 20px")
    )
  })
  
  # Download the dataframe
  output$dl_df_hm_aba <- downloadHandler(
    filename = function() {
      if (input$file_type_hm_aba == "EXCEL") { ext <- ".xlsx" } else { ext <- paste0(".", tolower(input$file_type_hm_aba))}
      paste0(input$file_name_hm_aba, ext)
    },
    
    content = function(file) {
      if (input$file_type_hm_aba == "EXCEL") {
        openxlsx::write.xlsx(df_exp_aba_hm(), file)
      } else {
        sep <- switch(input$file_type_hm_aba, "TXT" = " ", "CSV" = ",", "TSV" = "\t" )
        write.table(x = df_exp_aba_hm(), file = file, sep = sep, quote = FALSE, row.names = FALSE)
      }
    }
  )
  
  # Download the heatmap
  output$dl_hm_aba <- downloadHandler(
    filename = function() {
      paste0(input$file_name_hm_aba, ".html")
    },
    content = function(file) {
      htmlwidgets::saveWidget(Hmap_aba(), file, selfcontained = TRUE)
    }
  )

# Heatmap - Mutant vs Col-0 -----------------------------------------------

  # Show/Hide this tab
  observe({
    if (length(input$genotype) == 1 && input$genotype == "col") {
      hideTab(inputId = "tabs1", target = "Heatmap - Mutant vs Col-0")
    } else {
      showTab(inputId = "tabs1", target = "Heatmap - Mutant vs Col-0")
    }
  })

  # Dataframe_all
  df_hm_mut_all <- reactive({
    req(nrow(df_pc_mut()) != 0)
    ## Change the label
    lab <- c("0", "0.003", "0.01", "0.03",  "0.09", "0.27", "0.8", 
             "2", "7", "22", "102", "160", "200")
    df_lab <- data.frame(
      ABA_Conc = c(0, 3.33, 9.99, 29.97, 89.90, 269.70, 809.09, 
                   2427.26, 7281.78, 21845.33, 102400.00, 160000.00, 200000.00),
      ABA = factor(lab, levels = lab)
    )
    mutant_list <- data_values$genotype_list[!data_values$genotype_list %in% "Col-0"]
    df_hm_mut_all <- df_pc_mut() %>% 
      filter(Genotype %in% mutant_list) %>% 
      mutate(ABA_Conc = round(ABA_nM, 2)) %>% 
      left_join(df_lab, by = "ABA_Conc") %>% 
      dplyr::select(-ABA_Conc)
    return(df_hm_mut_all)
  })
  
  # Exported dataframe
  df_exp_mut_hm <- reactive({
    req(nrow(df_pc_mut()) != 0)
    
    # exported data frame
    mutant_list <- data_values$genotype_list[!data_values$genotype_list %in% "Col-0"]
    df_exp_mut_hm <- df_hm_mut_all() %>% 
      dplyr::select(AGI, Genotype, ABA_nM, Variable, Value) %>% 
      arrange(factor(Variable, level = c("logFC", "FDR", "logCPM", "F", "PValue"))) %>% 
      pivot_wider(names_from = Variable, values_from = Value) %>% 
      left_join(df_desc(), by = "AGI") %>% 
      mutate(Comparison = "Mutant vs Col-0") %>% 
      dplyr::select(Comparison, AGI, tair_symbol, Genotype, ABA_nM, logFC, FDR, logCPM, F, PValue, entrezgene_description) %>% 
      arrange(factor(AGI, levels = data_values$agi_list), factor(Genotype, levels = mutant_list), ABA_nM)
    
    return(df_exp_mut_hm)
  })
  
  # Heatmap
  Hmap_mut <- reactive({
    req(nrow(df_pc_mut()) != 0)
    ## Mutant list
    mutant_list <- data_values$genotype_list[!data_values$genotype_list %in% "Col-0"]
    
    ## FC dataframe
    df_hm_temp <- df_hm_mut_all() %>% 
      filter(Variable == "logFC") %>% 
      unite(Genotype, ABA, col = "comb", sep = ", ABA(μM): ", remove = FALSE) %>% 
      arrange(factor(AGI, levels = data_values$agi_list), factor(Genotype, levels = mutant_list), ABA) %>% 
      dplyr::select(AGI, comb, Value) %>% 
      pivot_wider(names_from = comb, values_from = Value)
    df_hm_fc <- df_hm_temp[2:ncol(df_hm_temp)]
    rownames(df_hm_fc) <- df_hm_temp$AGI
    
    ## FDR annotation
    df_hm_temp <- df_hm_mut_all() %>% 
      filter(Variable == "FDR") %>% 
      left_join(df_desc(), by = "AGI") %>% 
      mutate(Lab = paste0("FDR: ", formatC(as.numeric(Value), format = "E", digit = 2), "\n",
                          "Gene: ", tair_symbol)) %>% 
      unite(Genotype, ABA, col = "comb", sep = ", ABA(μM): ", remove = FALSE) %>% 
      arrange(factor(AGI, levels = data_values$agi_list), factor(Genotype, levels = mutant_list), ABA) %>% 
      dplyr::select(AGI, comb, Lab) %>% 
      pivot_wider(names_from = comb, values_from = Lab)
    df_hm_fdr <- df_hm_temp[2:ncol(df_hm_temp)]
    rownames(df_hm_fdr) <- df_hm_temp$AGI
    
    # Genotype annotation
    df_hm_temp <- df_hm_mut_all() %>% 
      filter(Variable == "FDR") %>% 
      unite(Genotype, ABA, col = "comb", sep = ", ABA(μM): ", remove = FALSE) %>% 
      dplyr::select(comb, Genotype, ABA) %>% unique() %>% 
      arrange(factor(Genotype, levels = mutant_list), ABA)
    gntp <- factor(df_hm_temp$Genotype, levels = mutant_list)
    color_map <- c(
      "sfkoI" = "#f2cf80",
      "sfkoII" = "#abd9f4",
      "sfkoIII" = "#80ceb9"
    )
    color_map_slt <- color_map[mutant_list]
    gntp_char <- as.character(df_hm_temp$Genotype)
    names(gntp_char) <- "Genotype"
    col_side_palette <- color_map_slt[levels(gntp)]
    col_side_colors_df <- data.frame(Genotype = gntp)
    
    ## Heatmap
    p <- heatmaply(
      df_hm_fc,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      col_side_colors = col_side_colors_df,
      col_side_palette = col_side_palette,
      grid_color = "white", grid_size = 0.1,
      dendrogram = "none",
      seriate = "none",
      key.title = "log2FC",
      xlab = "ABA (Low ⇢ High)",
      main = "Mutant vs Col-0",
      label_names = c("AGI", "Genotype", "log2FC (Mutant vs Col-0)"),
      custom_hovertext = df_hm_fdr,
      scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
        low = "#2171B5", 
        high = "#CB181D", 
        midpoint = 0
      )
    ) %>%
      plotly::layout(
        xaxis = list(showticklabels = FALSE, ticks = ""), 
        margin = list(l = 50, r = 50, t = 100, b = 100)
      )
    
    return(p)
    
  })
  
  output$heatmap_mut <- renderPlotly({
    return(Hmap_mut())
  })
  
  output$hm_mut <- renderUI({
    req(nrow(df_pc_mut()) != 0)
    tagList(
      # Heatmap
      div(style = "margin-top: 10px"),
      plotlyOutput("heatmap_mut", width = "100%", height = "600px") %>% shinycssloaders::withSpinner(), 
      p(tags$b("Note:"), "Only genes considered as expressed in our experiment are displayed."),
      # Download
      div(style = "margin-top: 20px"),
      h5("Download"),
      div(style = "margin-top: -10px"),
      hr(),
      div(style = "margin-top: -10px"),
      div(style = "vertical-align: top;", 
          textInput(inputId = "file_name_hm_mut", label = "Enter a file name: ", value = paste0("Mutant_vs_Wildtype_", Sys.Date()))
      ),
      div(style = "vertical-align: top;", 
          selectInput(inputId = "file_type_hm_mut", label = "Select file type for the dataframe:", 
                      choices = list("EXCEL", "CSV", "TSV", "TXT"), selected = "EXCEL")
      ),
      div(),
      div(style = "display: inline-block; vertical-align: top; width: 200px;",
          downloadButton(outputId = "dl_df_hm_mut", label = "Download Dataframe")),
      div(style = "display: inline-block; vertical-align: top; width: 200px;",
          downloadButton(outputId = "dl_hm_mut", label = "Download Heatmap")),
      div(style = "margin-top: 20px")
    )
  })
  
  # Download the dataframe
  output$dl_df_hm_mut <- downloadHandler(
    filename = function() {
      if (input$file_type_hm_mut == "EXCEL") { ext <- ".xlsx" } else { ext <- paste0(".", tolower(input$file_type_hm_mut))}
      paste0(input$file_name_hm_mut, ext)
    },
    
    content = function(file) {
      if (input$file_type_hm_mut == "EXCEL") {
        openxlsx::write.xlsx(df_exp_mut_hm(), file)
      } else {
        sep <- switch(input$file_type_hm_mut, "TXT" = " ", "CSV" = ",", "TSV" = "\t" )
        write.table(x = df_exp_mut_hm(), file = file, sep = sep, quote = FALSE, row.names = FALSE)
      }
    }
  )
  
  # Download the heatmap
  output$dl_hm_mut <- downloadHandler(
    filename = function() {
      paste0(input$file_name_hm_mut, ".html")
    },
    content = function(file) {
      htmlwidgets::saveWidget(Hmap_mut(), file, selfcontained = TRUE)
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
    colnames(df_temp) <- c("AGI", "TAIR_Symbol", "Genotype", "Cluster", "Trend", "Sensitivity", "Membership", "Maximum_Response", "Minimum_Response",
                           "Response_at_BMD", "BMD", "BMD_LowerBound", "BMD_UpperBound", 
                           "Response_at_Low_EC50", "Low_EC50", "Low_EC50_LowerBound", "Low_EC50_UpperBound", 
                           "Response_at_High_EC50", "High_EC50", "High_EC50_LowerBound", "High_EC50_UpperBound", 
                           "Response_at_LDS", "LDS", "LDS_LowerBound", "LDS_UpperBound", 
                           "Response_at_M", "M", 
                           "Model", "Neill's test", "No effet test")
    return(df_temp)
  })
  
  ## tables for display
  # low EC50s
  ED50_table <- reactive({
    
    req(df_ed_exp())
    df_temp <- df_ed_exp() %>% 
      dplyr::select(1:6, 15:17, 11:13) %>% 
      filter(!is.na(Low_EC50)|!is.na(BMD))
    colnames(df_temp) <- c("AGI", "TAIR_Symbol", "Genotype", "Cluster", "Trend", "Sensitivity", 
                           "EC\u2085\u2080", "EC\u2085\u2080\nLowerBound", "EC\u2085\u2080\nUpperBound", 
                           "BMD", "BMD\nLowerBound", "BMD\nUpperBound")
    df_temp <- df_temp %>% mutate(across(7:ncol(df_temp), ~ map_chr(.x, display_format)))
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
      h5(HTML(paste0("EC", tags$sub("50"), " & BMD Estimation Table")), align = 'center'),
      div(style = "margin-top: -10px"),
      #hr(),
      #div(style = "margin-top: -10px"),
      DT::dataTableOutput("tb_ed") %>% shinycssloaders::withSpinner(),
      p(tags$b("Note:"), HTML(paste0("Only genes with valid EC", tags$sub("50"), " or BMD values are displayed."))),
      # Download
      div(style = "margin-top: 20px"),
      h5("Download"),
      div(style = "margin-top: -10px"),
      hr(),
      div(style = "margin-top: -10px"),
      div(style = "display: inline-block; vertical-align: top;", 
          textInput(inputId = "file_name", label = "Enter a file name: ", value = paste0("Sensitivity_Table_", Sys.Date()))
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
  
  # UI for selecting rep genes
  output$agi_slt_ui <- renderUI({
    req(data_values$agi_list)
    if (isTRUE(input$mut_drc)) {
      max <- 1
    } else {
      max <- 2
    }
    selectizeInput(inputId = "agi_slt",
                   label = "Select the genes:",
                   choices = data_values$agi_list, 
                   selected = head(data_values$agi_list, 1),
                   multiple = TRUE,
                   options = list(maxItems = max))
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
  
  # Include the mutants
  output$mut_drc_ck <- renderUI({
    if (!(length(input$genotype) == 1 && input$genotype == "col")) {
      checkboxInput(inputId = "mut_drc",
                    label = "Include the Mutants",
                    value = FALSE)
    }
  })
  
## Demo Plot ---------------------------------------------------------------

### wt DRC ------------------------------------------------------------------

  demo_p <- reactive({
    req(input$agi_slt, data())
    agi_slt <- unname(input$agi_slt)
    
    ## dataframe
    df_curve_slt <- df_curve_wt() %>% filter(AGI %in% agi_slt) %>% unnest()
    df_point_slt <- df_point() %>% filter(Genotype == "Col-0") %>% filter(AGI %in% agi_slt) %>% unnest()
    df_point_mean_slt <- df_point_mean() %>% filter(AGI %in% agi_slt)
    df_ed_slt <- df_ed_slt()
    
    df_curve_slt$AGI <- factor(df_curve_slt$AGI, levels = agi_slt)
    df_point_slt$AGI <- factor(df_point_slt$AGI, levels = agi_slt)
    df_point_mean_slt$AGI <- factor(df_point_mean_slt$AGI, levels = agi_slt)
    df_ed_slt$AGI <- factor(df_ed_slt$AGI, levels = agi_slt)
    
    
    ## plot
    source(file.path("src/drc_wt.R"), local = TRUE)$value
    # remove the legend
    p_slt <- p_slt + theme(legend.position = "none")
    p <- ggplotly(p_slt, tooltip = "text")
    
    return(p)
    
  })
  
### mutant DRC --------------------------------------------------------------

  demo_p_mut <- reactive({
    req(input$agi_slt, data())
    agi_slt <- unname(input$agi_slt)
    genotype_list <- data_values$genotype_list
    df_curve_slt <- rbind(df_curve_wt(), df_curve_sfko1(), df_curve_sfko2(), df_curve_sfko3()) %>% 
      filter(AGI %in% agi_slt) %>% 
      filter(Genotype %in% genotype_list) %>% 
      unnest()
    df_point_slt <- df_point() %>% 
      filter(AGI %in% agi_slt) %>% 
      filter(Genotype %in% genotype_list) %>% 
      unnest()
    df_point_mean_slt <- df_point_mean() %>% 
      filter(Genotype %in% genotype_list) %>% 
      filter(AGI %in% agi_slt)
    df_ed_slt <- df_ed_slt()
    
    ## AGI
    df_curve_slt$AGI <- factor(df_curve_slt$AGI, levels = agi_slt)
    df_point_slt$AGI <- factor(df_point_slt$AGI, levels = agi_slt)
    df_point_mean_slt$AGI <- factor(df_point_mean_slt$AGI, levels = agi_slt)
    df_ed_slt$AGI <- factor(df_ed_slt$AGI, levels = agi_slt)
    ## Genotype
    df_curve_slt$Genotype <- factor(df_curve_slt$Genotype, levels = genotype_list)
    df_point_slt$Genotype <- factor(df_point_slt$Genotype, levels = genotype_list)
    df_point_mean_slt$Genotype <- factor(df_point_mean_slt$Genotype, levels = genotype_list)
    df_ed_slt$Genotype <- factor(df_ed_slt$Genotype, levels = genotype_list)
    
    ## plot on the left -- Col-0
    
    source(file.path("src/drc_wt.R"), local = TRUE)$value
    # remove the legend
    p_slt <- p_slt + theme(legend.position = "none")
    
    ## plot on the right - with all selected genotype
    # palette
    pal <- c("#D55E00")
    if ("sfkoI" %in% genotype_list) {
      pal <- c(pal, "#E69F00")
    }
    if ("sfkoII" %in% genotype_list) {
      pal <- c(pal, "#56B4E9")
    }
    if ("sfkoIII" %in% genotype_list) {
      pal <- c(pal, "#009E73")
    }
    
    p2 <- ggplot() + 
      # the major plot
      geom_point(data = df_point_slt, 
                 aes(x = ABA_nM, y = TPM, color = Genotype, 
                     text = paste("ABA [nM]:", round(ABA_nM, 2), "<br>Expression (TPM):", round(TPM, 2))), 
                 alpha = 1/3) + 
      geom_line(data = df_curve_slt, aes(x = ABA_nM, y = Prediction, color = Genotype)) +
      facet_wrap(ordered(AGI) ~ ., ncol = 1, scales = "free_y") +
      scale_x_log10() + 
      # labs(x = "ABA [nM]",  y = "Expression (TPM)") + 
      scale_color_manual(values = pal) +
      scale_fill_manual(values = pal) +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) + 
      theme(axis.title.x = element_text(size = 14),
            axis.text.x = element_text(size = 10),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            strip.text.x = element_text(size = 12),
            legend.position = "none") + 
      theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0), "in"))
    
    ## plotly
    p <- subplot(
      ggplotly(p_slt),
      ggplotly(p2),
      nrows = 1,
      shareX = FALSE,
      shareY = FALSE
    ) %>%
      layout(
        margin = list(b = 80, l = 80),
        
        annotations = list(
          list(
            text = "ABA [nM]",
            x = 0.5,
            y = -0.18,
            xref = "paper",
            yref = "paper",
            xanchor = "center",
            yanchor = "top",
            showarrow = FALSE,
            font = list(size = 20)
          ),
          list(
            text = "Expression (TPM)",
            x = -0.08,
            y = 0.5,
            xref = "paper",
            yref = "paper",
            textangle = -90,
            xanchor = "center",
            yanchor = "middle",
            showarrow = FALSE,
            font = list(size = 20)
          )
        )
      )
    
    return(p)
    
  })
  
  output$dr_curve <- renderPlotly({
    if (isTRUE(input$mut_drc)) {
      return(demo_p_mut())
    } else {
      return(demo_p())
    }
    
  })

## Download pdf ------------------------------------------------------------
  
  output$dl_pdf <- downloadHandler(
    
    filename = function() {
      paste0(input$pdf_name, ".pdf")
    },
    
    content = function(file) {
      
      genotype_list <- data_values$genotype_list
      agi_list_val <- data_values$agi_list
      
      # ed-related dataframe
      df_ed_exp_val <- df_ed_exp() %>%
        arrange(factor(AGI, levels = data_values$agi_list))
      
      
      if (isTRUE(input$mut_drc)) {
        
        # dataframe
        df_point_val <- df_point() %>% 
          filter(Genotype %in% genotype_list) %>% 
          unnest() %>% 
          group_by(AGI) %>% nest() %>% 
          arrange(factor(AGI, levels = data_values$agi_list))
        df_curve_val <- rbind(df_curve_wt(), df_curve_sfko1(), df_curve_sfko2(), df_curve_sfko3()) %>% 
          filter(Genotype %in% genotype_list) %>% 
          unnest() %>% 
          group_by(AGI) %>% nest() %>% 
          arrange(factor(AGI, levels = data_values$agi_list))
        
        # export pdf
        pdf(file = file, width = 1200/72, height = 1500/72)
        N <- 4
        mx <- ceiling(length(agi_list_val)/N)
        
        for (j in 1: mx) {
          a <- (j - 1) * N + 1
          if (j < mx) {b <- j*N} else {b <- length(agi_list_val)}
          n <- b - a + 1
          
          # Subset data for pages
          df_point_slt <- df_point_val[a:b, ] %>% unnest() 
          df_curve_slt <- df_curve_val[a:b, ] %>% unnest() 
          df_ed_slt<- df_ed_exp_val[a:b, ]
          
          df_curve_slt$AGI <- factor(df_curve_slt$AGI, levels = unique(df_curve_slt$AGI))
          df_point_slt$AGI <- factor(df_point_slt$AGI, levels = unique(df_point_slt$AGI))
          df_ed_slt$AGI <- factor(df_ed_slt$AGI, levels = unique(df_ed_slt$AGI))
          
          # plot
          x <- n
          
          ## Col-0 only
          
          source(file.path("src/drc_wt.R"), local = TRUE)$value
          
          # change the facet order and font size
          p_slt <- p_slt + 
            facet_wrap(ordered(AGI) ~ ., ncol = 1, scales = "free_y") +
            theme(axis.title = element_text(size = 20),
                  axis.text = element_text(size = 16),
                  strip.text.x = element_text(size = 18),
                  legend.key = element_rect(fill = "white"),
                  legend.title = element_text(size = 20, color = "white"),
                  legend.text = element_text(size = 18, color = "white"),
                  legend.position = "bottom") + guides(fill = guide_legend(nrow = 1)) + guides(color = guide_legend(override.aes = list(color = NA))) + 
            theme(plot.margin = unit(c(0.5, 0.5, 0.5*floor(x/4), 0.5), "in"))
          
          # gene info.
          if (input$gene_info == TRUE){
            p_slt <- p_slt +
              ggpp::geom_text_npc(data = df_ed_slt, aes(npcx = "left", npcy = "top", label = TAIR_Symbol), size = 6, color = "black")
          }
          
          ## Mutants
          # palette
          pal <- c("#D55E00")
          if ("sfkoI" %in% genotype_list) {
            pal <- c(pal, "#E69F00")
          }
          if ("sfkoII" %in% genotype_list) {
            pal <- c(pal, "#56B4E9")
          }
          if ("sfkoIII" %in% genotype_list) {
            pal <- c(pal, "#009E73")
          }
          
          p2 <- ggplot() + 
            # the major plot
            geom_point(data = df_point_slt, 
                       aes(x = ABA_nM, y = TPM, color = Genotype, 
                           text = paste("ABA [nM]:", round(ABA_nM, 2), "<br>Expression (TPM):", round(TPM, 2))), 
                       alpha = 1/3) + 
            geom_line(data = df_curve_slt, aes(x = ABA_nM, y = Prediction, color = Genotype)) +
            facet_wrap(ordered(AGI) ~ ., ncol = 1, scales = "free_y") +
            scale_x_log10() + 
            labs(x = "ABA [nM]",  y = "Expression (TPM)") + 
            scale_color_manual(values = pal) +
            scale_fill_manual(values = pal) +
            theme_bw() +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()) + 
            theme(axis.title = element_text(size = 20),
                  axis.text = element_text(size = 16),
                  strip.text.x = element_text(size = 18),
                  legend.title = element_text(size = 20),
                  legend.text = element_text(size = 18),
                  legend.position = "bottom") + guides(fill = guide_legend(nrow = 1)) + 
            theme(plot.margin = unit(c(0.5, 0.5, 0.5*floor(x/4), 0.5), "in"))
          
          if (n == 4) {
            print(ggarrange(p_slt, p2, ncol = 2, nrow = 1))
          } else {
            print(ggarrange(p_slt, p2, NULL, NULL,  ncol = 2, nrow = 2, heights = c(n*0.25+0.05, 1-(n*0.25+0.05))))
          }
          
        }
        
        dev.off()
        
      } else {
        
        # dataframe
        df_point_val <- df_point() %>% 
          filter(Genotype == "Col-0") %>% 
          arrange(factor(AGI, levels = data_values$agi_list))
        df_curve_val <- df_curve_wt() %>% 
          arrange(factor(AGI, levels = data_values$agi_list))
        
        # export pdf
        pdf(file = file, width = 1000/72, height = 1250/72)
        
        N <- 8
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
          
          source(file.path("src/drc_wt.R"), local = TRUE)$value
          
          # change the font size
          p_slt <- p_slt + 
            theme(axis.title = element_text(size = 20),
                  axis.text = element_text(size = 16),
                  strip.text.x = element_text(size = 18),
                  legend.title = element_text(size = 20), 
                  legend.text = element_text(size = 18),
                  legned.position = "bottom")
          
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
    }
  )
  
}
