
p_slt <- ggplot() + 
  # the major plot
  geom_point(data = df_point_slt %>% filter(Genotype == "Col-0"), 
             aes(x = ABA_nM, y = TPM, color = Genotype, fill = Genotype,
                 text = paste("ABA [nM]:", round(ABA_nM, 2), "<br>Expression (TPM):", round(TPM, 2))), alpha = 1/3) + 
  geom_line(data = df_curve_slt %>% filter(Genotype == "Col-0"), aes(x = ABA_nM, y = Prediction, color = Genotype)) +
  facet_wrap(. ~ ordered(AGI) , scales = "free_y") + 
  scale_x_log10() + 
  labs(x = "ABA [nM]",  y = "Expression (TPM)") + 
  scale_color_manual(values = "#D55E00") +
  scale_fill_manual(values = "#D55E00") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 10),
        strip.text.x = element_text(size = 12)) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "in"))

# ed50 lines
if (input$plot_ed50_ck == TRUE) {
  p_slt <- p_slt +
    geom_hline(data = df_ed_slt, aes(yintercept = Response_at_Low_EC50), linetype = "longdash", color = "#CC79A7", alpha = 0.5) + 
    geom_vline(data = df_ed_slt, aes(xintercept = Low_EC50, 
                                     text = paste("EC50:", round(Low_EC50, 2),"<br>Response_at_EC50:", round(Response_at_Low_EC50, 2))),
               linetype = "longdash", color = "#CC79A7", alpha = 0.5)
  if (input$plot_ci_ck == TRUE) {
    p_slt <- p_slt + 
      geom_vline(data = df_ed_slt, aes(xintercept = Low_EC50_LowerBound), linetype = "dotted", color = "#CC79A7", alpha = 0.5) + 
      geom_vline(data = df_ed_slt, aes(xintercept = Low_EC50_UpperBound), linetype = "dotted", color = "#CC79A7", alpha = 0.5)
  }
}

# bmd lines
if (input$plot_bmd_ck == TRUE) {
  p_slt <- p_slt +
    geom_hline(data = df_ed_slt, aes(yintercept = Response_at_BMD), linetype = "dashed", color = "#CC79A7", alpha = 0.5) + 
    geom_vline(data = df_ed_slt, aes(xintercept = BMD, 
                                     text = paste("BMD:", round(BMD, 2),"<br>Response_at_BMD:", round(Response_at_BMD, 2))), 
               linetype = "dashed", color = "#CC79A7", alpha = 0.5)
  if (input$plot_ci_ck == TRUE) {
    p_slt <- p_slt +
      geom_vline(data = df_ed_slt, aes(xintercept = BMD_LowerBound), linetype = "dotted", color = "#CC79A7", alpha = 0.5) + 
      geom_vline(data = df_ed_slt, aes(xintercept = BMD_UpperBound), linetype = "dotted", color = "#CC79A7", alpha = 0.5)
  }
}

# max & min response
if (input$plot_resline_ck == TRUE) {
  p_slt <- p_slt +
    geom_hline(data = df_ed_slt, aes(yintercept = Maximum_Response, 
                                     text = paste("Max:", round(Maximum_Response, 2))), linetype = "dotted", color = "#999999", alpha = 0.5) + 
    geom_hline(data = df_ed_slt, aes(yintercept = Minimum_Response, 
                                     text = paste("Min:", round(Minimum_Response, 2))), linetype = "dotted", color = "#999999", alpha = 0.5)
}

# lds & m
if (any(df_ed_slt$Cluster == 1) && all(!is.na(df_ed_slt$Cluster))) {
  if (input$plot_lds_m_ck == TRUE) {
    p_slt <- p_slt +
      # lds lines
      geom_hline(data = df_ed_slt, aes(yintercept = Response_at_LDS), linetype = "dashed", color = "#999999", alpha = 0.5) + 
      geom_vline(data = df_ed_slt, aes(xintercept = LDS, 
                                       text = paste("LDS:", round(LDS, 2),"<br>Response_at_LDS:", round(Response_at_LDS, 2))),
                 linetype = "dashed", color = "#999999", alpha = 0.5) + 
      # m lines
      geom_hline(data = df_ed_slt, aes(yintercept = Response_at_M), linetype = "dashed", color = "#999999", alpha = 0.5) + 
      geom_vline(data = df_ed_slt, aes(xintercept = M, 
                                       text = paste("M:", round(M, 2),"<br>Response_at_M:", round(Response_at_M, 2))), 
                 linetype = "dashed", color = "#999999", alpha = 0.5)
    
    if (input$plot_ci_ck == TRUE) {
      p_slt <- p_slt +
        geom_vline(data = df_ed_slt, aes(xintercept = LDS_LowerBound), linetype = "dotted", color = "#999999", alpha = 0.5) + 
        geom_vline(data = df_ed_slt, aes(xintercept = LDS_UpperBound), linetype = "dotted", color = "#999999", alpha = 0.5)
    }
  }
}
