library(ggplot2)

df <- read.table("coverage.tsv", header = T)
total_reads_number <- read.table("total_reads_aligned.txt", header = F)
total_reads_number <- total_reads_number[1, 1]
df$RPKM <- df$coverage/df$interval_length*1000/total_reads_number*1000000
df$interval_name <- df$interval_name * 1000


# Calculate the peaks
percent_left_term_peak <- sum(subset(df, interval_name >= 1328000 & interval_name <= 1350000)['coverage'])/total_reads_number*100
percent_right_term_peak <- sum(subset(df, interval_name >= 1524000 & interval_name <= 1557000)['coverage'])/total_reads_number*100


ggplot(df, aes(x = interval_name, y = RPKM)) +
  geom_col(fill = "black", width = 1000) +
  ylab("RPKM × 1000") +
  xlab("chromosome coordinate, Mb") +
  ggtitle("sDNA coverage") +
  scale_x_continuous(breaks = c(0, 500000, 1000000, 1500000, 2000000, 2500000, 3000000, 3500000, 4000000, 4500000),
                     labels = c(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5)) +
  scale_y_continuous(limits = c(0, 7000),
                     breaks = c(0, 2000, 4000, 6000),
                     labels = c("0", "2", "4", "6")) +
  geom_vline(xintercept = c(73000), linetype = "dashed", colour = "grey30", size = 1) +
  geom_vline(xintercept = c(1329000, 1558000), linetype = "dashed", colour = "grey30", size = 1) +
  geom_vline(xintercept = 3815000, linetype = "dashed", colour = "grey30", size = 1) +
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

ggsave("both.png", width = 15, height = 10, dpi = 400)

df_plus <- read.table("plus_coverage.tsv", header = T)
df_minus <- read.table("minus_coverage.tsv", header = T)
df_comb <- data.frame(df$interval_name, df$interval_length)
df_comb$RPKM_plus <- df_plus$coverage/df_plus$interval_length*1000/total_reads_number*1000000
df_comb$RPKM_minus <- df_minus$coverage/df_minus$interval_length*1000/total_reads_number*1000000*-1

ggplot(df_comb, aes(x = df.interval_name)) +
  geom_col(aes(y = RPKM_plus), fill = "turquoise4", width = 1000) +
  geom_col(aes(y = RPKM_minus), fill = "maroon4", width = 1000) +
  ylab("RPKM × 1000") +
  xlab("chromosome coordinate, Mb") +
  ggtitle("sDNA coverage") +
  scale_x_continuous(breaks = c(0, 500000, 1000000, 1500000, 2000000, 2500000, 3000000, 3500000, 4000000, 4500000),
                     labels = c(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5)) +
  scale_y_continuous(limits = c(-5000, 5500),
                     breaks = c(-4000, -2000, 0, 2000, 4000),
                     labels = c("-4", "-2", "0", "2", "4")) +
  geom_vline(xintercept = c(1329000, 1558000), linetype = "dashed", colour = "grey30", size = 1) +
  geom_vline(xintercept = 3815000, linetype = "dashed", colour = "grey30", size = 1) +
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

ggsave("plus_minus.png", width = 15, height = 10, dpi = 400)


chi_plus_df <- read.table("plus_chi.txt", header = F)
chi_plus_v <- chi_plus_df[,1]
chi_minus_df <- read.table("minus_chi.txt", header = F)
chi_minus_v <- chi_minus_df[,1]

ggplot(df_comb, aes(x = df.interval_name)) +
  geom_col(aes(y = RPKM_plus), fill = "turquoise4", width = 1000) +
  geom_col(aes(y = RPKM_minus), fill = "maroon4", width = 1000) +
  ylab("RPKM × 1000") +
  xlab("chromosome coordinate, Mb") +
  ggtitle("sDNA coverage") +
  scale_x_continuous(limits = c(1200000, 1700000),
                     breaks = c(1200000, 1300000, 1400000, 1500000, 1600000, 1700000),
                     labels = c(1.2, 1.3, 1.4, 1.5, 1.6, 1.7)) +
  scale_y_continuous(limits = c(-5000, 5500),
                     breaks = c(-4000, -2000, 0, 2000, 4000),
                     labels = c("-4", "-2", "0", "2", "4")) +
  geom_vline(xintercept = c(1267000, 1328000, 1557000, 1629000), linetype = "dashed", colour = "grey30", size = 1) +
  annotate("segment", x = chi_plus_v, xend = chi_plus_v, y = -5000, yend = -4200, col = "turquoise4", size = 1) +
  annotate("segment", x = chi_minus_v, xend = chi_minus_v, y = -5000, yend = -4200, col = "maroon4", size = 1) +
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

ggsave("Ter_region.png", width = 15, height = 10, dpi = 400)
