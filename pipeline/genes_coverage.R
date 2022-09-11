library(ggplot2)

total_reads_number <- read.table("total_reads_aligned.txt", header = F)
total_reads_number <- total_reads_number[1, 1]

df_sense_co <- read.table("sense_co.tsv", header = T)
df_sense_co$RPKM <- df_sense_co$coverage/df_sense_co$interval_length*1000/total_reads_number*1000000
df_sense_co$rep <- "codirected"
df_sense_co$dir <- "sense"

df_sense_rev <- read.table("sense_rev.tsv", header = T)
df_sense_rev$RPKM <- df_sense_rev$coverage/df_sense_rev$interval_length*1000/total_reads_number*1000000
df_sense_rev$rep <- "reversely directed"
df_sense_rev$dir <- "sense"

df_antisense_co <- read.table("antisense_co.tsv", header = T)
df_antisense_co$RPKM <- df_antisense_co$coverage/df_antisense_co$interval_length*1000/total_reads_number*1000000
df_antisense_co$rep <- "codirected"
df_antisense_co$dir <- "antisense"

df_antisense_rev <- read.table("antisense_rev.tsv", header = T)
df_antisense_rev$RPKM <- df_antisense_rev$coverage/df_antisense_rev$interval_length*1000/total_reads_number*1000000
df_antisense_rev$rep <- "reversely directed"
df_antisense_rev$dir <- "antisense"

df <- rbind(df_sense_co, df_sense_rev, df_antisense_co, df_antisense_rev)
df$rep <- as.factor(df$rep)
df$dir <- as.factor(df$dir)

ggplot(df, aes(x = dir, y = RPKM)) +
  geom_boxplot(aes(fill = dir)) +
  ylab("RPKM") +
  scale_y_continuous(limits = c(0, 300)) +
#  scale_y_continuous(trans = "log10") +
  scale_x_discrete(name = "",
                   breaks = c("antisense", "sense"),
                   labels = c("", "")) +
  scale_fill_manual(breaks = c("antisense", "sense"),
                    values = c("salmon", "turquoise"),
                    name = "") +
  facet_wrap(~rep) +
  theme(text = element_text(size = 40),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        axis.ticks.length.y.left = unit(0.3, "cm"),
        axis.ticks.length.x = unit(0, "cm"),
        axis.text = element_text(size = 40, colour = 'black'))

ggsave("boxplot_rep_dir.png", width = 13, height = 7, dpi = 400)
