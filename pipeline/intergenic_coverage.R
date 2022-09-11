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
  scale_y_continuous(limits = c(0, 400)) +
  theme(text = element_text(size = 40),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 2),
        axis.ticks = element_line(size = 1.5),
        axis.ticks.length = unit(0.3, "cm"),
        axis.text = element_text(size = 40, colour = 'black'))

ggsave("boxplot.png", width = 17, height = 10, dpi = 400)
