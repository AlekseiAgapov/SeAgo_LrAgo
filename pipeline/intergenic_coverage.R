library(ggplot2)

total_reads_number <- read.table("total_reads_aligned.txt", header = F)
total_reads_number <- total_reads_number[1, 1]

df_conv <- read.table("conv_coverage.tsv", header = T)
df_conv$RPKM <- df_conv$coverage/df_conv$interval_length*1000/total_reads_number*1000000
df_conv$type <- 'convergent'
conv_coverage <- mean(df_conv$RPKM)

df_div <- read.table("div_coverage.tsv", header = T)
df_div$RPKM <- df_div$coverage/df_div$interval_length*1000/total_reads_number*1000000
df_div$type <- 'divergent'
div_coverage <- mean(df_div$RPKM)

df_plus <- read.table("plus_coverage.tsv", header = T)
df_plus$RPKM <- df_plus$coverage/df_plus$interval_length*1000/total_reads_number*1000000
df_plus$type <- 'plus'
plus_coverage <- mean(df_plus$RPKM)

df_minus <- read.table("minus_coverage.tsv", header = T)
df_minus$RPKM <- df_minus$coverage/df_minus$interval_length*1000/total_reads_number*1000000
df_minus$type <- 'minus'
minus_coverage <- mean(df_minus$RPKM)

df <- rbind(df_div, df_conv, df_plus, df_minus)
df$type <- as.factor(df$type)
df <- subset(df, interval_name >= 1700000 | interval_name <= 1267000)

ggplot(df, aes(x = type, y = RPKM)) +
  geom_boxplot(aes(fill = factor(type))) +
  ylab("RPKM") +
  xlab("type") +
  ggtitle("coverage of intergenic regions") +
  scale_y_continuous(limits = c(0, 300)) +
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

ggsave("boxplot.png", width = 17, height = 10, dpi = 400)

df_sub <- rbind(subset(df_div, interval_name >= 1700000 | interval_name <= 1267000), subset(df_conv, interval_name >= 1700000 | interval_name <= 1267000))

ggplot(df_sub, aes(x = interval_name, y = RPKM)) +
  geom_point(aes(col = factor(type)), alpha = 0.5) +
  scale_color_manual(breaks = c('convergent', 'divergent'),
                     values=c("coral2", "chartreuse4")) +
  ylab("RPKM") +
  xlab("chromosome coordinate, Mb") +
  ggtitle("coverage of intergenic regions") +
  scale_x_continuous(breaks = c(0, 500000, 1000000, 1500000, 2000000, 2500000, 3000000, 3500000, 4000000, 4500000),
                     labels = c(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5)) +
  scale_y_continuous(limits = c(0, 1000)) +
  geom_hline(yintercept = conv_coverage, linetype = "dashed", colour = "coral2") +
  geom_hline(yintercept = div_coverage, linetype = "dashed", colour = "chartreuse4") +
  theme(text = element_text(size = 23),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        axis.text = element_text(size = 18, colour = 'black'),
        title = element_text(size = 20),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(color = "white",
                                          linetype = "longdash", 
                                          size = 0.3))

ggsave("scatter_plot.png", width = 15, height = 5, dpi = 600)
