#!/usr/bin/python3
from statistics import *
import re
import os

Rcode = """
Graph <- ggplot(DF, aes(x = Size, y = Mean, color = Version, group = Version)) + 
  geom_point(size = 2.5,mapping = aes(shape=Version)) + 
  geom_line(aes(linetype = Version), size=1.25) +
  theme_bw() +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
              labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  annotation_logticks(sides = "l")  +
  geom_errorbar(aes(ymax = Mean + 1.96*DP/sqrt(30), ymin = Mean - 1.96*DP/sqrt(30)), width=0.4) +
  ylab("Mean (seconds)") + 
  xlab(expression(paste("Mesh size"))) +
  theme(plot.title = element_text(family = "Times", face="bold", size=12)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title = element_text(family = "Times", face="bold", size=12)) +
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
"""

float_regex = "([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)"
keys  = ("SHARED", "GHMATECE", "RIGID", "GHMATECD", "LINSOLVE", "INTEREC1")
mesh_numbers = (240, 960, 2160, 4000, 14400)
modes = ("cpu", "gpu") #, "gpu_sing")
threads = (1, 8, 24, 48)
executions = 30

def get_data_from_file(filename):
    result = {}
    f = open(filename, "r")
    regexes = (
        "Tempo gasto alocando os vetores compartilhados:\s*" + float_regex,
        "GHMATECE: Tempo na ...:\s*" + float_regex,
        "GHMATECE: Corpo rigido:\s*" + float_regex,
        "GHMATECD: Tempo na ...:\s*" + float_regex,
        "LINSOLVE: Tempo na ...:\s*" + float_regex,
        "INTEREC1: Tempo na ...:\s*" + float_regex
    )

    lines = f.read().split("\n")
    for line in lines:
        for i in range(len(keys)):
            key = keys[i]
            regex = regexes[i]

            match = re.search(regex, line)
            if (match):
                result[key] = float(match.group(1))
                break

    return result


def get_all_data(mode, mesh, threads, n):
    name = "results/results_" + mode + "_" + str(mesh) + "_" + str(threads) + "_"
    extention = ".txt"

    times = {}
    for key in keys:
        times[key] = []

    for i in range(1,n+1):
        filename = name + str(i) + extention
        time = get_data_from_file(filename)
        for key in time:
            times[key].append(time[key])
        
    return times

def generate_r_data(values, subroutine_name):
    r_format1 = "{} <- cbind(size, \"{}\", c("
    r_format2 = "))"

    DF1 = "data.frame(rbind("

    header_string = "library(ggplot2)\nsize <- c(240, 960, 2160, 4000, 14400)\n"

    for mode in modes:
        if (mode == "cpu"):
            threads=(1, 24, 48)
        if (mode == "gpu"):
            threads=[8]
        if (mode == "gpu_sing"):
            threads=[1]

        for thread in threads:
            string = ""
            for mesh in mesh_numbers:
                string += str(values[mode][thread][mesh][subroutine_name]["MEAN"]) + ","
       
            string = string[:-2]
            final_string = r_format1.format(mode + str(thread), mode + str(thread)) + string + r_format2
            header_string += final_string + "\n"

    string = ""
    for mode in modes:
        if (mode == "cpu"):
            threads=(1, 24, 48)
        if (mode == "gpu"):
            threads=[8]
        if (mode == "gpu_sing"):
            threads=[1]
        for thread in threads:
            for mesh in mesh_numbers:
                string += str(values[mode][thread][mesh][subroutine_name]["STDEV"]) + ","
       
    string = string[:-1]
    final_string = "DP <- c(" + string + ")"
    header_string += final_string + "\n"

    string = ""
    for mode in modes:
        if (mode == "cpu"):
            threads=(1, 24, 48)
        if (mode == "gpu"):
            threads=[8]
        if (mode == "gpu_sing"):
            threads=[1]
        for thread in threads:
            string += mode + str(thread) + ","
    string = string[:-1]
    final_string = "DF <- data.frame(rbind(" + string + "))"
    header_string += final_string + "\n"
    header_string += 'names(DF) <- c("Size", "Version", "Mean")' + "\n"
    header_string += 'DF$Size <- factor(DF$Size, levels = c("240", "960", "2160", "4000", "14400"))' + "\n"
    header_string += 'DF$Mean <- as.numeric(as.character(DF$Mean))' + "\n"
    header_string += 'DP <- as.numeric(as.character(DP))' + "\n"
    return header_string

def main():
    statistical_values = {}
    for mode in modes:
        statistical_values[mode] = {}
        for thread in threads:
            statistical_values[mode][thread] = {}
            for mesh in mesh_numbers:
                statistical_values[mode][thread][mesh] = {}
                t = get_all_data(mode, mesh, thread, executions)
                for key in t:
                    if (len(t[key]) == 0):
                        continue
                    m  = mean(t[key])
                    sd = stdev(t[key])
                   
                    if (mode != "cpu" and (key == "GHMATECD" or key == "GHMATECE")):
                        m = m + mean(t["SHARED"])

                    statistical_values[mode][thread][mesh][key] = {}
                    statistical_values[mode][thread][mesh][key]["MEAN"] = m
                    statistical_values[mode][thread][mesh][key]["STDEV"] = sd


    header = generate_r_data(statistical_values, "GHMATECE")
    temp_file = open("temp.r", "w")
    temp_file.write(header + Rcode)
    temp_file.close()

    os.system("Rscript temp.r")

if __name__ == "__main__":
    main()

