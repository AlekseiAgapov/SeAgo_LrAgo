library(ggplot2)

df <- read.table("./cipro/coverage_10000/coverage.tsv", header = T)
total_reads_number <- read.table("./cipro/alignment/total_reads_aligned.txt", header = F)
total_reads_number <- total_reads_number[1, 1]
df$RPKM <- df$coverage/df$interval_length*1000/total_reads_number*1000000
df$interval_name = df$interval_name * 10000

df_plus <- read.table("./cipro/coverage_10000/plus_coverage.tsv", header = T)
df_minus <- read.table("./cipro/coverage_10000/minus_coverage.tsv", header = T)

df_comb <- data.frame(df$interval_name, df$interval_length)
df_comb$RPKM_plus <- df_plus$coverage/df_plus$interval_length*1000/total_reads_number*1000000
df_comb$RPKM_minus <- df_minus$coverage/df_minus$interval_length*1000/total_reads_number*1000000*-1

df_comb$ratio <- df_comb$RPKM_plus/df_comb$RPKM_minus*-1

av_mean <- c()
for(i in 1:length(df_comb$ratio)){
  x <- (df_comb$ratio[(i-2)] + df_comb$ratio[(i-1)] + df_comb$ratio[i] + df_comb$ratio[(i+1)] + df_comb$ratio[(i+2)])/5
  av_mean[i] <- x
}

av_mean[1] <- (df_comb$ratio[(length(df_comb$ratio)-1)] + df_comb$ratio[length(df_comb$ratio)] + df_comb$ratio[1] + df_comb$ratio[2] + df_comb$ratio[3])/5
av_mean[2] <- (df_comb$ratio[length(df_comb$ratio)] + df_comb$ratio[1] + df_comb$ratio[2] + df_comb$ratio[3] + df_comb$ratio[4])/5
av_mean[(length(av_mean) - 1)] <- (df_comb$ratio[(length(df_comb$ratio)-3)] + df_comb$ratio[(length(df_comb$ratio)-2)] + df_comb$ratio[(length(df_comb$ratio)-1)] + df_comb$ratio[length(df_comb$ratio)] + df_comb$ratio[1])/5
av_mean[length(av_mean)] <- (df_comb$ratio[(length(df_comb$ratio)-2)] + df_comb$ratio[(length(df_comb$ratio)-1)] + df_comb$ratio[length(df_comb$ratio)] + df_comb$ratio[1] + df_comb$ratio[1])/5

df_comb$av_mean_ratio <- av_mean

df_ratio <- data.frame('interval_name' = df_comb$df.interval_name, 'cip' = df_comb$av_mean_ratio)

df_plus <- read.table("./no_cipro/coverage_10000/plus_coverage.tsv", header = T)
df_minus <- read.table("./no_cipro/coverage_10000/minus_coverage.tsv", header = T)

df_comb <- data.frame(df$interval_name, df$interval_length)
df_comb$RPKM_plus <- df_plus$coverage/df_plus$interval_length*1000/total_reads_number*1000000
df_comb$RPKM_minus <- df_minus$coverage/df_minus$interval_length*1000/total_reads_number*1000000*-1

df_comb$ratio <- df_comb$RPKM_plus/df_comb$RPKM_minus*-1

av_mean <- c()
for(i in 1:length(df_comb$ratio)){
  x <- (df_comb$ratio[(i-2)] + df_comb$ratio[(i-1)] + df_comb$ratio[i] + df_comb$ratio[(i+1)] + df_comb$ratio[(i+2)])/5
  av_mean[i] <- x
}

av_mean[1] <- (df_comb$ratio[(length(df_comb$ratio)-1)] + df_comb$ratio[length(df_comb$ratio)] + df_comb$ratio[1] + df_comb$ratio[2] + df_comb$ratio[3])/5
av_mean[2] <- (df_comb$ratio[length(df_comb$ratio)] + df_comb$ratio[1] + df_comb$ratio[2] + df_comb$ratio[3] + df_comb$ratio[4])/5
av_mean[(length(av_mean) - 1)] <- (df_comb$ratio[(length(df_comb$ratio)-3)] + df_comb$ratio[(length(df_comb$ratio)-2)] + df_comb$ratio[(length(df_comb$ratio)-1)] + df_comb$ratio[length(df_comb$ratio)] + df_comb$ratio[1])/5
av_mean[length(av_mean)] <- (df_comb$ratio[(length(df_comb$ratio)-2)] + df_comb$ratio[(length(df_comb$ratio)-1)] + df_comb$ratio[length(df_comb$ratio)] + df_comb$ratio[1] + df_comb$ratio[1])/5

df_comb$av_mean_ratio <- av_mean

df_ratio$no_cip <- df_comb$av_mean_ratio

df_ratio$ratio <- df_ratio$cip / df_ratio$no_cip

ggplot(df_ratio, aes(x = interval_name, y = ratio)) +
  geom_line(col = "black", size = 2.5) +
  ylab("+ cipro / - cipro") +
  xlab("chromosome coordinate, Mb") +
#  ggtitle("") +
  scale_x_continuous(breaks = c(0, 500000, 1000000, 1500000, 2000000, 2500000, 3000000, 3500000, 4000000, 4500000),
                     labels = c(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5)) +
  scale_y_continuous(limits = c(0.5, 2), breaks = c(0.5, 1, 2), trans = "log2") +
  geom_vline(xintercept = c(1320000, 1550000), linetype = "dashed", colour = "grey30", size = 1) +
  geom_vline(xintercept = 3810000, linetype = "dashed", colour = "grey30", size = 1) +
  geom_hline(yintercept = 1, colour = "black", size = 1) +
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

ggsave("strand_ratio.png", width = 12, height = 8, dpi = 400)
