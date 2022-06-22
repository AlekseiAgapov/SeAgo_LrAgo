library(ggplot2)

table <- read.table("lengths.tsv", header = T)

ggplot(table, aes(x = length, y = percent)) +
  geom_col(fill = "turquoise4", alpha=0.7, width=0.65) +
  scale_x_continuous(breaks = c(14, 16, 18, 20, 22, 24)) +
  #  ggtitle("smDNA length distribution") +
  ylab("% of all aligned reads") +
  scale_y_continuous(limits = c(0, 30), 
                     breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  xlab("length, nt") +
  theme(text = element_text(size = 30),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 2),
        axis.ticks = element_line(size = 2),
        axis.ticks.length = unit(0.3, "cm"),
        axis.text = element_text(size = 30, colour = 'black'),
        axis.title = element_text(size = 30),
        title = element_text(size = 30),
        panel.grid.major.y = element_line(color = "grey66",
                                          linetype = "dashed", 
                                          size = 0.5))

ggsave("lengths_distribution.png", width = 10, height = 7, dpi = 300)
