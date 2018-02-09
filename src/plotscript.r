###############
library(ggplot2)

size <- c(240, 960, 2160, 4000)
# Mean
CPUseq <- cbind(size, "CPUSeq", c(0.47939565440028675, 7.482897269700091,37.21196495326694,123.72078290366699))
CPUOMP <- cbind(size, "CPU-OMP", c(0.159272048434165, 3.2126924417738336, 16.992475143164242, 56.72498899894029 ))
GPU <-  cbind(size, "GPU+OMP", c(0.1288769468703928,0.31956582149529517, 1.0928144277461493, 3.2652298854043087))
DF <- data.frame(rbind(CPUseq, CPUOMP, GPU))
names(DF) <- c("Size", "Version", "Mean")
DF$Size <- factor(DF$Size, levels = c("240", "960", "2160", "4000"))
DF$Mean <- as.numeric(as.character(DF$Mean))

Graph <- ggplot(DF, aes(x = Size, y = Mean, color = Version, group = Version)) + 
  geom_point(size = 2.5,mapping = aes(shape=Version)) + 
  geom_line(aes(linetype = Version), size=1.25) +
  theme_bw() + 
#  scale_colour_manual(values=c(cbbPalette)) +
  scale_fill_manual(values=c(cbbPalette)) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
              labels = scales::trans_format("log10", scales::math_format(10^.x)))  + 
  annotation_logticks(sides = "l")  +
  ylab("Mean (seconds)") + 
  xlab(expression(paste("Matrix Size"))) +
  theme(plot.title = element_text(family = "Times", face="bold", size=16)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title = element_text(family = "Times", face="bold", size=16)) +
  theme(axis.text  = element_text(family = "Times", face="bold", size=10, colour = "Black")) +
  theme(legend.title  = element_text(family = "Times", face="bold", size=0)) +
  theme(legend.text  = element_text(family = "Times", face="bold", size=12)) +
  theme(legend.key.size = unit(5, "cm")) +
  theme(legend.direction = "horizontal",
        legend.position = "bottom",
        legend.key=element_rect(size=0),
        legend.key.size = unit(1, "lines")) +
  guides(col = guide_legend(nrow = 1))
ggsave(paste("~/Giuliano-plots.pdf",sep=""), Graph,  height=10, width=15, units="cm")


# Median
CPUseq <- c(0.4790759129991784, 7.484063304498704, 37.15695717749986, 123.71510654349913 )
CPUOMP <- c(0.15920984800322913, 3.2040811575134285, 16.986593497509602, 56.69127115351148 )
GPU <-  c(0.12808945053257048, 0.3175889280100819, 1.0840124645328615,  3.266195518517634)
  
##################
