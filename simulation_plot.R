#######N<<M########

## nonlinear common and rare
data41 <- data.frame(  
  Method = factor(rep(c("iSVR", "VW_TOW_GE"), 2)), 
  level  = factor(c(rep(0.01, 2), rep(0.05, 2))),
  value  = c(
    0.411, 0.004,   # alpha = 0.01: iSVR power, VW_TOW_GE power
    0.698, 0.032    # alpha = 0.05: iSVR power, VW_TOW_GE power
  )
)


## nonlinear rare
data42 <- data.frame(
  Method = factor(rep(c("iSVR", "TOW_GE"), 2)),
  level  = factor(c(rep(0.01, 2), rep(0.05, 2))),
  value  = c(
    0.129, 0.050,   # alpha = 0.01: iSVR power, TOW_GE power
    0.310, 0.151    # alpha = 0.05: iSVR power, TOW_GE power
  )
)


library(ggplot2)
library(dplyr)


data41$Comparison <- "Common and rare varitants (iSVR vs VW_TOW_GE)"
data42$Comparison <- "Rare varitants (iSVR vs TOW_GE)"

# Combine data frames
combined_data <- bind_rows(data41, data42)
B <- 1000

combined_data <- combined_data %>%
  mutate(
    MCSE  = sqrt(value * (1 - value) / B),
    lower = pmax(0, value - 1.96 * MCSE),
    upper = pmin(1, value + 1.96 * MCSE),
    label_y = ifelse(value <= 0.04, upper + 0.02, value / 2)
    )

all_methods <- unique(c(data41$Method, data42$Method))


color_palette <- c(
  "iSVR" = "#3C9BC8",
  "VW_TOW_GE" = "#D9D69B",
  "TOW_GE" = "#8BA486"  # TOW_GE: a different color
)


color_palette <- color_palette[names(color_palette) %in% all_methods]

# Combined Diagram
p_combined <- ggplot(combined_data, aes(x = level, y = value, fill = Method)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.8),
           width = 0.8,
           color = "black") +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    position = position_dodge(width = 0.8),
    width = 0.25,
    linewidth = 0.6) +
  geom_text(
    aes(y = label_y, label = sprintf("%.3f", value)),
    size = 4,
    position = position_dodge(width = 0.8),
    vjust = 0.5) +
  labs(
    x = "Significance Level",
    y = "Power"
  ) +
  scale_y_continuous(limits = c(0, 0.8)) +
  scale_fill_manual(values = color_palette) +
  facet_wrap(~ Comparison, nrow = 1) +  
  theme_bw(base_size = 18) +
  theme(
    axis.text = element_text(colour = "black"),
    strip.background = element_rect(fill = "lightgray"),  # 
    strip.text = element_text(size = 14, face = "bold"),   # 
    legend.position = "bottom"
  )

# Graphic
print(p_combined)

# save PDF
ggsave("NM_nonlinear_plot.pdf", p_combined,  dpi = 350, width = 12, height = 6, units = "in")







## POWER_fig
#####nonlinear common and rare
data1 = data.frame(
  Method = factor(rep(c("GESAT","iSKAT","VW_TOW_GE","CKLRT","iSVR"), 2)),
  level  = factor(c(rep(0.01, 5), rep(0.05, 5))),
  value  = c(
    0.048, 0.049, 0.071, 0.599, 1.000,
    0.178, 0.179, 0.265, 0.789, 1.000
  )
)

data2 = data.frame(
  Method = factor(rep(c("GESAT","iSKAT","VW_TOW_GE","CKLRT","iSVR"), 2)),
  level  = factor(c(rep(0.01, 5), rep(0.05, 5))),
  value  = c(
    0.037, 0.039, 0.066, 0.572, 0.909,
    0.161, 0.164, 0.208, 0.789, 0.970
  )
)

data3 = data.frame(
  Method = factor(rep(c("GESAT","iSKAT","VW_TOW_GE","CKLRT","iSVR"), 2)),
  level  = factor(c(rep(0.01, 5), rep(0.05, 5))),
  value  = c(
    0.066, 0.066, 0.143, 0.687, 0.997,
    0.219, 0.221, 0.360, 0.832, 0.997
  )
)


#########nonlinear rare
case1 = data.frame(
  Method = factor(rep(c("GESAT","iSKAT","TOW_GE","CKLRT","iSVR"), 2)),
  level  = factor(c(rep(0.01, 5), rep(0.05, 5))),
  value  = c(
    0.030, 0.030, 0.085, 0.090, 0.237,
    0.105, 0.107, 0.221, 0.226, 0.401
  )
)

case2 = data.frame(
  Method = factor(rep(c("GESAT","iSKAT","TOW_GE","CKLRT","iSVR"), 2)),
  level  = factor(c(rep(0.01, 5), rep(0.05, 5))),
  value  = c(
    0.056, 0.057, 0.254, 0.296, 0.649,
    0.164, 0.167, 0.469, 0.534, 0.826
  )
)

case3 = data.frame(
  Method = factor(rep(c("GESAT","iSKAT","TOW_GE","CKLRT","iSVR"), 2)),
  level  = factor(c(rep(0.01, 5), rep(0.05, 5))),
  value  = c(
    0.026, 0.026, 0.062, 0.069, 0.122,
    0.098, 0.100, 0.185, 0.193, 0.265
  )
)

############Graphic####
library(ggplot2)
library(dplyr)
library(cowplot)

# p_combined (Common and Rare Variants)-top

data1$Dataset <- "Case 1"
data2$Dataset <- "Case 2"
data3$Dataset <- "Case 3"

combined_data <- bind_rows(data1, data2, data3)
## Monte Carlo uncertainty for power
## MCSE = sqrt(p_hat * (1 - p_hat) / 1000)
## 95% CI = p_hat +/- 1.96 * MCSE

B <- 1000

combined_data <- combined_data %>%
  mutate(
    MCSE = sqrt(value * (1 - value) / B),
    lower = pmax(0, value - 1.96 * MCSE),
    upper = pmin(1, value + 1.96 * MCSE),
    label_y = ifelse(value <= 0.08, upper + 0.03, value / 2)
  )

method_order <- c("GESAT", "iSKAT", "VW_TOW_GE", "CKLRT", "iSVR")
combined_data$Method <- factor(combined_data$Method, levels = method_order)

color_palette <- c(
  "GESAT" = "#9384B4",
  "iSKAT" = "#B17F9F",
  "TOW_GE" = "#8BA486",
  "VW_TOW_GE" = "#D9D69B",
  "CKLRT" = "#C1E0DB",
  "iSVR" = "#3C9BC8"
)

# Common and Rare Variants-top
p1 <- ggplot(combined_data, aes(x = level, y = value, fill = Method)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.8),
           width = 0.8,
           color = "black",
           linewidth = 0.3) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    position = position_dodge(width = 0.8),
    width = 0.25,
    linewidth = 0.6
  ) +
  geom_text(
    aes(y = label_y, label = sprintf("%.3f", value)),
    size = 3.2,
    position = position_dodge(width = 0.8),
    vjust = 0.5
  ) +
  labs(
    x = "Significance Level",  
    y = "Power",
    title = "(a) Common and rare variants"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    expand = expansion(mult = c(0, 0.1))
  ) +
  scale_fill_manual(values = color_palette) +
  facet_wrap(~ Dataset, nrow = 1) +
  theme_bw(base_size = 14) +
  theme(
    axis.text = element_text(colour = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    strip.background = element_rect(fill = "#2C3E50"),
    strip.text = element_text(size = 12, face = "bold", color = "white"),
    legend.position = "none",  
    plot.title = element_text(hjust = 0, face = "bold", size = 16, 
                              margin = margin(b = 10)),  
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(t = 10, r = 10, b = 5, l = 10, unit = "pt")  #
  )

# p_combined1 (Rare Variants) - bottom

case1$Dataset <- "Case 1"
case2$Dataset <- "Case 2"
case3$Dataset <- "Case 3"

combined_data1 <- bind_rows(case1, case2, case3)
## Monte Carlo uncertainty for power
## MCSE = sqrt(p_hat * (1 - p_hat) / 1000)
## 95% CI = p_hat +/- 1.96 * MCSE

combined_data1 <- combined_data1 %>%
  mutate(
    MCSE = sqrt(value * (1 - value) / B),
    lower = pmax(0, value - 1.96 * MCSE),
    upper = pmin(1, value + 1.96 * MCSE),
    label_y = ifelse(value <= 0.08, upper + 0.03, value / 2)
  )
method_order1 <- c("GESAT", "iSKAT", "TOW_GE", "CKLRT", "iSVR")
combined_data1$Method <- factor(combined_data1$Method, levels = method_order1)

# (Rare Variants) - bottom
p2 <- ggplot(combined_data1, aes(x = level, y = value, fill = Method)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.8),
           width = 0.8,
           color = "black",
           linewidth = 0.3) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    position = position_dodge(width = 0.8),
    width = 0.25,
    linewidth = 0.6
  ) +
  geom_text(
    aes(y = label_y, label = sprintf("%.3f", value)),
    size = 3.2,
    position = position_dodge(width = 0.8),
    vjust = 0.5
  ) +
  labs(
    x = "Significance Level",
    y = "Power",
    title = "(b) Rare variants"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    expand = expansion(mult = c(0, 0.1))
  ) +
  scale_fill_manual(values = color_palette) +
  facet_wrap(~ Dataset, nrow = 1) +
  theme_bw(base_size = 14) +
  theme(
    axis.text = element_text(colour = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    strip.background = element_rect(fill = "#2C3E50"),
    strip.text = element_text(size = 12, face = "bold", color = "white"),
    legend.position = "none",  
    plot.title = element_text(hjust = 0, face = "bold", size = 16, 
                              margin = margin(b = 10)),  
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(t = 5, r = 10, b = 10, l = 10, unit = "pt")  
  )

##Legend
legend_data <- data.frame(
  Method = factor(c("GESAT", "iSKAT", "TOW_GE", "VW_TOW_GE", "CKLRT", "iSVR"), 
                  levels = c("GESAT", "iSKAT", "TOW_GE", "VW_TOW_GE", "CKLRT", "iSVR")),
  x = 1:6,
  y = 1:6
)

####Compact Legend
get_compact_legend <- function() {
  
  legend_plot_compact <- ggplot(legend_data, aes(x = x, y = y, fill = Method)) +
    geom_bar(stat = "identity", width = 0.5) +
    scale_fill_manual(values = color_palette) +
    guides(fill = guide_legend(
      nrow = 1,  
      byrow = TRUE,
      title.position = "top",
      label.position = "right",
      keywidth = unit(0.5, "cm"),  
      keyheight = unit(0.4, "cm"),  
      default.unit = "cm"
    )) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_blank(),
      legend.text = element_text(size = 10),  
      legend.spacing.x = unit(0.1, "cm"),  
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(-5, 0, 0, 0)  
    )
  
  return(get_legend(legend_plot_compact))
}

# Compact Legend Graphics
legend_compact <- get_compact_legend()

# Compact Legend Combination Graphics
final_plot_compact <- plot_grid(
  p1, p2, legend_compact,
  ncol = 1,
  rel_heights = c(2, 2, 0.1),  
  align = "v",
  axis = "lr"
)

#  Combination Graphics
print(final_plot_compact)


# save PDF
ggsave("combined_nonlinear_common_rare_variants.pdf", 
       final_plot_compact, dpi = 350,
       width = 16,
       height = 12,
       units = "in")









########type I error
####nonlinear#
##### common and rare

data11 = data.frame(
  Method = factor(rep(c("GESAT", "iSKAT", "VW_TOW_GE", "CKLRT", "iSVR"), 2)),
  level  = factor(c(rep(0.01, 5), rep(0.05, 5))),
  value  = c(
    0.006, 0.006, 0.008, 0.006, 0.015,
    0.037, 0.037, 0.052, 0.057, 0.058
  )
)


######### rare

case11 = data.frame(
  Method = factor(rep(c("GESAT", "iSKAT", "TOW_GE", "CKLRT", "iSVR"), 2)),
  level  = factor(c(rep(0.01, 5), rep(0.05, 5))),
  value  = c(
    0.007, 0.008, 0.015, 0.014, 0.002,
    0.045, 0.046, 0.058, 0.060, 0.039
  )
)
#########CASE 1##########
library(ggplot2)
library(dplyr)
library(cowplot)

case11$Dataset <- "Rare Variants"
data11$Dataset <- "Common and Rare Variants"


combined_data <- bind_rows(case11, data11)


method_order <- c("GESAT", "iSKAT", "TOW_GE", "VW_TOW_GE", "CKLRT", "iSVR")
combined_data$Method <- factor(combined_data$Method, levels = method_order)

## Monte Carlo uncertainty for type I error
## MCSE = sqrt(p_hat * (1 - p_hat) / 1000)
## 95% CI = p_hat +/- 1.96 * MCSE

B <- 1000

combined_data <- combined_data %>%
  mutate(
    MCSE  = sqrt(value * (1 - value) / B),
    lower = pmax(0, value - 1.96 * MCSE),
    upper = pmin(1, value + 1.96 * MCSE),
        label_y = ifelse(value <= 0.005, value + 0.004, value / 2)
  )

color_palette <- c(
  "GESAT" = "#9384B4",
  "iSKAT" = "#B17F9F",
  "TOW_GE" = "#8BA486",
  "VW_TOW_GE" = "#D9D69B",
  "CKLRT" = "#C1E0DB",
  "iSVR" = "#3C9BC8"
)


data_0_01 <- combined_data %>% filter(level == 0.01)
data_0_05 <- combined_data %>% filter(level == 0.05)

########### 0.01 ###########
p_0_01 <- ggplot(data_0_01, aes(x = Method, y = value, fill = Method)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.8),
           width = 0.8,
           color = "black",
           linewidth = 0.3) +
  
  #geom_hline(yintercept = 0.01, 
  #           linetype = "dashed", 
  #           color = "blue", 
  #           linewidth = 1.2,  
  #           alpha = 0.9) +
  
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    position = position_dodge(width = 0.8),
    width = 0.25,
    linewidth = 0.6
  ) +
  
  geom_text(
    aes(y = label_y, label = sprintf("%.3f", value)),
    size = 3.2,
    position = position_dodge(width = 0.8),
    vjust = 0.5
  ) +
  labs(
    x = "Method",
    y = "Type I Error Rate",
    title = "(a) Significance Level = 0.01"
  ) +
  scale_y_continuous(
    limits = c(0, 0.075),
    breaks = seq(0, 0.075, by = 0.01),
    expand = expansion(mult = c(0, 0.1))
  ) +
  scale_fill_manual(values = color_palette) +
  facet_wrap(~ Dataset, nrow = 1, scales = "free_x") +
  theme_bw(base_size = 14) +
  theme(
    axis.text = element_text(colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    strip.background = element_rect(fill = "#2C3E50"),
    strip.text = element_text(size = 12, face = "bold", color = "white"),
    legend.position = "none",  
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, 
                              margin = margin(b = 10)),  
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(t = 10, r = 10, b = 5, l = 10, unit = "pt")
  )

########### 0.05 ###########
p_0_05 <- ggplot(data_0_05, aes(x = Method, y = value, fill = Method)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.8),
           width = 0.8,
           color = "black",
           linewidth = 0.3) +
  #geom_hline(yintercept = 0.05, 
  #           linetype = "dashed", 
  #           color = "blue", 
  #           linewidth = 1.2,  
  #           alpha = 0.9) +
  
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    position = position_dodge(width = 0.8),
    width = 0.25,
    linewidth = 0.6
  ) +
  
  geom_text(
    aes(y = label_y, label = sprintf("%.3f", value)),
    size = 3.2,
    position = position_dodge(width = 0.8),
    vjust = 0.5
  ) +
  labs(
    x = "Method",
    y = "Type I Error Rate",
    title = "(b) Significance Level = 0.05"
  ) +
  scale_y_continuous(
    limits = c(0, 0.075),
    breaks = seq(0, 0.075, by = 0.01),
    expand = expansion(mult = c(0, 0.1))
  ) +
  scale_fill_manual(values = color_palette) +
  facet_wrap(~ Dataset, nrow = 1, scales = "free_x") +
  theme_bw(base_size = 14) +
  theme(
    axis.text = element_text(colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    strip.background = element_rect(fill = "#2C3E50"),
    strip.text = element_text(size = 12, face = "bold", color = "white"),
    legend.position = "none",  
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, 
                              margin = margin(b = 10)),  
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(t = 5, r = 10, b = 10, l = 10, unit = "pt")
  )

########### legend ###########
legend_data <- data.frame(
  Method = factor(c("GESAT", "iSKAT", "TOW_GE", "VW_TOW_GE", "CKLRT", "iSVR"), 
                  levels = c("GESAT", "iSKAT", "TOW_GE", "VW_TOW_GE", "CKLRT", "iSVR")),
  x = 1:6,
  y = 1:6
)


get_compact_legend <- function() {
  
  legend_plot_compact <- ggplot(legend_data, aes(x = x, y = y, fill = Method)) +
    geom_bar(stat = "identity", width = 0.5) +
    scale_fill_manual(values = color_palette) +
    guides(fill = guide_legend(
      nrow = 1,  
      byrow = TRUE,
      title.position = "top",
      label.position = "right",
      keywidth = unit(1, "cm"),  
      keyheight = unit(0.8, "cm"),  
      default.unit = "cm"
    )) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_blank(),
      legend.text = element_text(size = 10),  
      legend.spacing.x = unit(0.1, "cm"),  
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(-5, 0, 0, 0)  
    )
  
  return(get_legend(legend_plot_compact))
}


legend_compact <- get_compact_legend()


final_plot_horizontal <- plot_grid(
  p_0_01, p_0_05,
  ncol = 2,
  align = "h",
  axis = "tb"
)


final_plot_horizontal_with_legend <- plot_grid(
  final_plot_horizontal, legend_compact,
  ncol = 1,
  rel_heights = c(4, 1),  
  align = "v",
  axis = "lr"
)


print(final_plot_horizontal_with_legend)

# save PDF
ggsave("combined_case1_nonlinear_significance.pdf", 
       final_plot_horizontal_with_legend, dpi=350,
       width = 16, 
       height = 8,
       units = "in")










## POWER_fig
#####linear common and rare
data111 = data.frame(
  Method = factor(rep(c("GESAT","iSKAT","VW_TOW_GE","CKLRT","iSVR"), 2)),
  level  = factor(c(rep(0.01, 5), rep(0.05, 5))),
  value  = c(
    0.071, 0.073, 0.998, 0.340, 0.998,
    0.238, 0.238, 1.000, 0.579, 0.999
  )
)

data211 = data.frame(
  Method = factor(rep(c("GESAT","iSKAT","VW_TOW_GE","CKLRT","iSVR"), 2)),
  level  = factor(c(rep(0.01, 5), rep(0.05, 5))),
  value  = c(
    0.011, 0.011, 1.000, 0.125, 1.000,
    0.074, 0.074, 1.000, 0.318, 1.000
  )
)

data311 = data.frame(
  Method = factor(rep(c("GESAT","iSKAT","VW_TOW_GE","CKLRT","iSVR"), 2)),
  level  = factor(c(rep(0.01, 5), rep(0.05, 5))),
  value  = c(
    0.015, 0.015, 1.000, 0.121, 1.000,
    0.087, 0.088, 1.000, 0.310, 1.000
  )
)


#########linear rare
case111 = data.frame(
  Method = factor(rep(c("GESAT","iSKAT","TOW_GE","CKLRT","iSVR"), 2)),
  level  = factor(c(rep(0.01, 5), rep(0.05, 5))),
  value  = c(
    0.008, 0.008, 0.137, 0.078, 0.163,
    0.069, 0.070, 0.319, 0.230, 0.346
  )
)

case211 = data.frame(
  Method = factor(rep(c("GESAT","iSKAT","TOW_GE","CKLRT","iSVR"), 2)),
  level  = factor(c(rep(0.01, 5), rep(0.05, 5))),
  value  = c(
    0.006, 0.006, 0.061, 0.021, 0.210,
    0.041, 0.041, 0.180, 0.090, 0.423
  )
)

case311 = data.frame(
  Method = factor(rep(c("GESAT","iSKAT","TOW_GE","CKLRT","iSVR"), 2)),
  level  = factor(c(rep(0.01, 5), rep(0.05, 5))),
  value  = c(
    0.013, 0.013, 0.312, 0.086, 0.422,
    0.069, 0.070, 0.532, 0.247, 0.658
  )
)

############Graphic####
library(ggplot2)
library(dplyr)
library(cowplot)

# p_combined (Common and Rare Variants)-top

data111$Dataset <- "Case 1"
data211$Dataset <- "Case 2"
data311$Dataset <- "Case 3"

combined_data <- bind_rows(data111, data211, data311)
## Monte Carlo uncertainty for power
## MCSE = sqrt(p_hat * (1 - p_hat) / 1000)
## 95% CI = p_hat +/- 1.96 * MCSE

B <- 1000

combined_data <- combined_data %>%
  mutate(
    MCSE = sqrt(value * (1 - value) / B),
    lower = pmax(0, value - 1.96 * MCSE),
    upper = pmin(1, value + 1.96 * MCSE),
    label_y = ifelse(value <= 0.08, upper + 0.03, value / 2)
  )

method_order <- c("GESAT", "iSKAT", "VW_TOW_GE", "CKLRT", "iSVR")
combined_data$Method <- factor(combined_data$Method, levels = method_order)

color_palette <- c(
  "GESAT" = "#9384B4",
  "iSKAT" = "#B17F9F",
  "TOW_GE" = "#8BA486",
  "VW_TOW_GE" = "#D9D69B",
  "CKLRT" = "#C1E0DB",
  "iSVR" = "#3C9BC8"
)

# Common and Rare Variants-top
p1 <- ggplot(combined_data, aes(x = level, y = value, fill = Method)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.8),
           width = 0.8,
           color = "black",
           linewidth = 0.3) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    position = position_dodge(width = 0.8),
    width = 0.25,
    linewidth = 0.6
  ) +
  geom_text(
    aes(y = label_y, label = sprintf("%.3f", value)),
    size = 3.2,
    position = position_dodge(width = 0.8),
    vjust = 0.5
  ) +
  labs(
    x = "Significance Level",  
    y = "Power",
    title = "(a) Common and rare variants"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    expand = expansion(mult = c(0, 0.1))
  ) +
  scale_fill_manual(values = color_palette) +
  facet_wrap(~ Dataset, nrow = 1) +
  theme_bw(base_size = 14) +
  theme(
    axis.text = element_text(colour = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    strip.background = element_rect(fill = "#2C3E50"),
    strip.text = element_text(size = 12, face = "bold", color = "white"),
    legend.position = "none",  
    plot.title = element_text(hjust = 0, face = "bold", size = 16, 
                              margin = margin(b = 10)),  
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(t = 10, r = 10, b = 5, l = 10, unit = "pt")  #
  )

# p_combined1 (Rare Variants) - bottom

case111$Dataset <- "Case 1"
case211$Dataset <- "Case 2"
case311$Dataset <- "Case 3"

combined_data1 <- bind_rows(case111, case211, case311)
## Monte Carlo uncertainty for power
## MCSE = sqrt(p_hat * (1 - p_hat) / 1000)
## 95% CI = p_hat +/- 1.96 * MCSE

combined_data1 <- combined_data1 %>%
  mutate(
    MCSE = sqrt(value * (1 - value) / B),
    lower = pmax(0, value - 1.96 * MCSE),
    upper = pmin(1, value + 1.96 * MCSE),
    label_y = ifelse(value <= 0.08, upper + 0.03, value / 2)
  )
method_order1 <- c("GESAT", "iSKAT", "TOW_GE", "CKLRT", "iSVR")
combined_data1$Method <- factor(combined_data1$Method, levels = method_order1)

# (Rare Variants) - bottom
p2 <- ggplot(combined_data1, aes(x = level, y = value, fill = Method)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.8),
           width = 0.8,
           color = "black",
           linewidth = 0.3) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    position = position_dodge(width = 0.8),
    width = 0.25,
    linewidth = 0.6
  ) +
  geom_text(
    aes(y = label_y, label = sprintf("%.3f", value)),
    size = 3.2,
    position = position_dodge(width = 0.8),
    vjust = 0.5
  ) +
  labs(
    x = "Significance Level",
    y = "Power",
    title = "(b) Rare variants"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    expand = expansion(mult = c(0, 0.1))
  ) +
  scale_fill_manual(values = color_palette) +
  facet_wrap(~ Dataset, nrow = 1) +
  theme_bw(base_size = 14) +
  theme(
    axis.text = element_text(colour = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    strip.background = element_rect(fill = "#2C3E50"),
    strip.text = element_text(size = 12, face = "bold", color = "white"),
    legend.position = "none",  
    plot.title = element_text(hjust = 0, face = "bold", size = 16, 
                              margin = margin(b = 10)),  
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(t = 5, r = 10, b = 10, l = 10, unit = "pt")  
  )

##Legend
legend_data <- data.frame(
  Method = factor(c("GESAT", "iSKAT", "TOW_GE", "VW_TOW_GE", "CKLRT", "iSVR"), 
                  levels = c("GESAT", "iSKAT", "TOW_GE", "VW_TOW_GE", "CKLRT", "iSVR")),
  x = 1:6,
  y = 1:6
)

####Compact Legend
get_compact_legend <- function() {
  
  legend_plot_compact <- ggplot(legend_data, aes(x = x, y = y, fill = Method)) +
    geom_bar(stat = "identity", width = 0.5) +
    scale_fill_manual(values = color_palette) +
    guides(fill = guide_legend(
      nrow = 1,  
      byrow = TRUE,
      title.position = "top",
      label.position = "right",
      keywidth = unit(0.5, "cm"),  
      keyheight = unit(0.4, "cm"),  
      default.unit = "cm"
    )) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_blank(),
      legend.text = element_text(size = 10),  
      legend.spacing.x = unit(0.1, "cm"),  
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(-5, 0, 0, 0)  
    )
  
  return(get_legend(legend_plot_compact))
}

# Compact Legend Graphics
legend_compact <- get_compact_legend()

# Compact Legend Combination Graphics
final_plot_compact1 <- plot_grid(
  p1, p2, legend_compact,
  ncol = 1,
  rel_heights = c(2, 2, 0.1),  
  align = "v",
  axis = "lr"
)

#  Combination Graphics
print(final_plot_compact1)


# save PDF
ggsave("combined_linear_common_rare_variants.pdf", 
       final_plot_compact1, dpi = 350,
       width = 16,
       height = 12,
       units = "in")





##N<<M  linear common and rare
data43 <- data.frame(  
  Method = factor(rep(c("iSVR", "VW_TOW_GE"), 2)), 
  level  = factor(c(rep(0.01, 2), rep(0.05, 2))),
  value  = c(
    1.000, 1.000,   # alpha = 0.01: iSVR power, VW_TOW_GE power
    1.000, 1.000    # alpha = 0.05: iSVR power, VW_TOW_GE power
  )
)


## linear rare
data44 <- data.frame(
  Method = factor(rep(c("iSVR", "TOW_GE"), 2)),
  level  = factor(c(rep(0.01, 2), rep(0.05, 2))),
  value  = c(
    0.534, 0.342,   # alpha = 0.01: iSVR power, TOW_GE power
    0.807, 0.619    # alpha = 0.05: iSVR power, TOW_GE power
  )
)

library(ggplot2)
library(dplyr)


data43$Comparison <- "Common and rare varitants (iSVR vs VW_TOW_GE)"
data44$Comparison <- "Rare varitants (iSVR vs TOW_GE)"

# Combine data frames
combined_data1 <- bind_rows(data43, data44)

B <- 1000

combined_data1 <- combined_data1 %>%
  mutate(
    MCSE  = sqrt(value * (1 - value) / B),
    lower = pmax(0, value - 1.96 * MCSE),
    upper = pmin(1, value + 1.96 * MCSE),
    
  
    label_y = ifelse(value >= 0.12, value / 2, value * 0.65)
  )
all_methods <- unique(c(data43$Method, data44$Method))


color_palette <- c(
  "iSVR" = "#3C9BC8",
  "VW_TOW_GE" = "#D9D69B",
  "TOW_GE" = "#8BA486"  # TOW_GE: a different color
)


color_palette <- color_palette[names(color_palette) %in% all_methods]

# Combined Diagram
p_combined1 <- ggplot(combined_data1, aes(x = level, y = value, fill = Method)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.8),
           width = 0.8,
           color = "black") +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    position = position_dodge(width = 0.8),
    width = 0.25,
    linewidth = 0.6
  ) +
  geom_text(
    aes(y = label_y, label = sprintf("%.3f", value)),
    size = 4,
    position = position_dodge(width = 0.8),
    vjust = 0.5
  ) +
  labs(
    x = "Significance Level",
    y = "Power"
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_manual(values = color_palette) +
  facet_wrap(~ Comparison, nrow = 1) +  
  theme_bw(base_size = 18) +
  theme(
    axis.text = element_text(colour = "black"),
    strip.background = element_rect(fill = "lightgray"),  # 
    strip.text = element_text(size = 14, face = "bold"),   # 
    legend.position = "bottom"
  )

# Graphic
print(p_combined1)

# save PDF
ggsave("NM_linear_plot.pdf", p_combined1,  dpi = 350, width = 12, height = 6, units = "in")


