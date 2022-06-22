library(ggplot2)

df_plus <- read.table("plus_cov.tsv", header = F)
names(df_plus)[1] <- "chr"
names(df_plus)[2] <- "coordinate"
names(df_plus)[3] <- "coverage"
df_plus = subset(df_plus, select = -c(chr))

df_minus <- read.table("minus_cov.tsv", header = F)
names(df_minus)[1] <- "chr"
names(df_minus)[2] <- "coordinate"
names(df_minus)[3] <- "coverage"
df_minus = subset(df_minus, select = -c(chr))

av_cov <- read.table("average_cov.txt", header = F)
av_cov <- av_cov[1, 1]

plus_cov_norm <- df_plus$coverage / av_cov
minus_cov_norm <- -1 * df_minus$coverage / av_cov

df <- data.frame(df_plus$coordinate, plus_cov_norm, minus_cov_norm)
names(df)[1] <- "coordinate"


# Chi-sites

plot_chi <- function(start_coordinate, name) {
  
  end_coordinate <- start_coordinate + 7
  
  df_for_plot <- subset(df, coordinate >= start_coordinate - 100 & coordinate <= end_coordinate + 100)
  
  ggplot(df_for_plot, aes(x = coordinate)) +
    geom_col(aes(y = plus_cov_norm), fill = "#297372", width = 1) +
    geom_col(aes(y = minus_cov_norm), fill = "#cc1452", width = 1) +
    ylab("normalized coverage") +
    xlab("chromosome coordinate") +
    ggtitle(name) +
    scale_x_continuous(breaks = c(start_coordinate - 100, start_coordinate - 50, start_coordinate, end_coordinate, end_coordinate + 50, end_coordinate + 100),
                       labels = c("-100", "-50", "  Chi", "", "+50", "+100")) +
    geom_hline(yintercept = 0, linetype = "solid", colour = "black", size = 1) +
    geom_vline(xintercept = c(start_coordinate, end_coordinate), linetype = "dashed", colour = "grey50") +
    theme(text = element_text(size = 40),
          panel.background = element_rect(fill = "white"),
          axis.line = element_line(size = 2),
          axis.ticks = element_line(size = 1.5),
          axis.ticks.length = unit(0.3, "cm"),
          axis.text = element_text(size = 40, colour = 'black'),
          title = element_text(size = 40),
          plot.title = element_text(hjust = 0.5),
          panel.grid.major.y = element_line(color = "white",
                                            linetype = "longdash", 
                                            size = 0.3))
}

plot_chi(1524183, "first Chi-site left of TerC")
ggsave("chi_left_TerC.png", width = 15, height = 10, dpi = 400)

plot_chi(1626124, "first Chi-site left of TerB")
ggsave("chi_left_TerB.png", width = 15, height = 10, dpi = 400)

plot_chi(1371631, "first Chi-site right of TerA")
ggsave("chi_right_TerA.png", width = 15, height = 10, dpi = 400)



# Ter-sites

plot_ter <- function(start_coordinate, end_coordinate, ter_name) {
  
  df_for_plot <- subset(df, coordinate >= start_coordinate - 100 & coordinate <= end_coordinate + 100)
  
  ggplot(df_for_plot, aes(x = coordinate)) +
    geom_col(aes(y = plus_cov_norm), fill = "#297372", width = 1) +
    geom_col(aes(y = minus_cov_norm), fill = "#cc1452", width = 1) +
    ylab("normalized coverage") +
    xlab("chromosome coordinate") +
    ggtitle(ter_name) +
    scale_x_continuous(breaks = c(start_coordinate - 100, start_coordinate - 50, start_coordinate, end_coordinate, end_coordinate + 50, end_coordinate + 100),
                       labels = c("-100", "-50", "      Ter-site", "", "+50", "+100")) +
    geom_hline(yintercept = 0, linetype = "solid", colour = "black", size = 1) +
    geom_vline(xintercept = c(start_coordinate, end_coordinate), linetype = "dashed", colour = "grey50") +
    theme(text = element_text(size = 40),
          panel.background = element_rect(fill = "white"),
          axis.line = element_line(size = 2),
          axis.ticks = element_line(size = 1.5),
          axis.ticks.length = unit(0.3, "cm"),
          axis.text = element_text(size = 40, colour = 'black'),
          title = element_text(size = 40),
          plot.title = element_text(hjust = 0.5),
          panel.grid.major.y = element_line(color = "white",
                                            linetype = "longdash", 
                                            size = 0.3))
}


plot_ter(1328110, 1328136, "TerA")
ggsave("TerA.png", width = 15, height = 10, dpi = 400)

plot_ter(1629423, 1629445, "TerB")
ggsave("TerB.png", width = 15, height = 10, dpi = 400)

plot_ter(1557232, 1557254, "TerC")
ggsave("TerC.png", width = 15, height = 10, dpi = 400)

plot_ter(1267671, 1267697, "TerD")
ggsave("TerD.png", width = 15, height = 10, dpi = 400)
