#!/usr/bin/python3
# Script used to plot graphics.
#

from statistics import *
import re
import os
import sys
from math import sqrt


float_regex = "([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)"
file_keys = ["SHARED", "GHMATECE", "GHMATECD", "LINSOLVE", "INTEREC1", "SIGMAEC"]
mesh_numbers = (510, 1820, 6840) #, 26480)
modes = ("cpu", "gpu") #, "gpu_sing")
all_threads = (1, 4, 8)
threads = {"cpu": all_threads, "gpu": [4, 8]}

executions = 30

keys = file_keys[:]
keys.append("ALL")

def get_data_from_file(filename):
    result = {}
    f = open(filename, "r")
    regexes = (
        "Tempo gasto alocando os vetores compartilhados:\s*" + float_regex,
        "GHMATECE: Tempo na ...:\s*" + float_regex,
        "GHMATECD: Tempo na ...:\s*" + float_regex,
        "LINSOLVE: Tempo na ...:\s*" + float_regex,
        "INTEREC1: Tempo na ...:\s*" + float_regex,
        "SIGMAEC: Tempo na ...:\s*"  + float_regex
    )

    lines = f.read().split("\n")
    for line in lines:
        for i in range(len(file_keys)):
            key = keys[i]
            regex = regexes[i]

            match = re.search(regex, line)
            if (match):
                result[key] = float(match.group(1))
                break
       
    acc = 0.0;
    for key in result:
       acc += result[key]
    
    result["ALL"] = acc

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


def generate_speedup_r_data(values, subroutine_name):
    
    Rcode = """
Graph <- ggplot(DF, aes(x = Size, y = Mean, color = Version, group = Version)) + 
  geom_point(size = 2.5,mapping = aes(shape=Version)) + 
  geom_line(aes(linetype = Version), size=1.25) +
  theme_bw() +
#  scale_y(breaks = scales::trans_breaks(function(x) x),
#         labels = scales::trans_format(scales::math_format(x))) + 
  annotation_logticks(sides = "l")  +
  ylab("Speedup") + 
  xlab(expression(paste("Number of Elements"))) +
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
ggsave(paste("~/speedup_double_prec.pdf",sep=""), Graph,  height=10, width=15, units="cm")
    """
    
    r_format1 = "{} <- cbind(size, \"{}\", c("
    r_format2 = "))"

    DF1 = "data.frame(rbind("

    header_string  = "library(ggplot2)\nsize <- c("
    
    for i in range(len(mesh_numbers)):
        header_string += str(mesh_numbers[i]) + ","
    header_string = header_string[:-1] + ")\n"

    for mode in modes:
        for thread in threads[mode]:
            
            if (thread == 1 and mode == "cpu"):
                continue

            string = ""
            for mesh in mesh_numbers:
                acc = values["cpu"][1][mesh][subroutine_name]["MEAN"]/values[mode][thread][mesh][subroutine_name]["MEAN"]
                string += str(acc) + ","

            string = string[:-2]
            final_string = r_format1.format(mode + str(thread), mode + str(thread)) + string + r_format2
            header_string += final_string + "\n"
    string = ""
    for mode in modes:
       for thread in threads[mode]:
            if (mode == "cpu" and thread == 1):
                continue
            string += mode + str(thread) + ","
    string = string[:-1]
    final_string = "DF <- data.frame(rbind(" + string + "))"
    header_string += final_string + "\n"
    header_string += 'names(DF) <- c("Size", "Version", "Mean")' + "\n"
    header_string += 'DF$Size <- factor(DF$Size, levels = c('
    
    for mesh in mesh_numbers:
        header_string += '"' + str(mesh) + '",'
    header_string = header_string[:-1]
    header_string += '))\n'
    header_string += 'DF$Mean <- as.numeric(as.character(DF$Mean))' + "\n"
#    header_string += 'DP <- as.numeric(as.character(DP))' + "\n"
    return header_string + Rcode

def generate_r_data(values, subroutine_name):
    
    Rcode = """
Graph <- ggplot(DF, aes(x = Size, y = Mean, color = Version, group = Version)) + 
  geom_point(size = 2.5,mapping = aes(shape=Version)) + 
  geom_line(aes(linetype = Version), size=1.25) +
  theme_bw() +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
              labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  annotation_logticks(sides = "l")  +
#  geom_errorbar(aes(ymax = Mean + DP, ymin = Mean - DP), width=0.4) +
  ylab("Mean (seconds)") + 
  xlab(expression(paste("Number of Elements"))) +
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
ggsave(paste("~/time_double_prec.pdf",sep=""), Graph,  height=10, width=15, units="cm")
    """
    
    r_format1 = "{} <- cbind(size, \"{}\", c("
    r_format2 = "))"

    DF1 = "data.frame(rbind("

    header_string  = "library(ggplot2)\nsize <- c("
    
    for i in range(len(mesh_numbers)):
        header_string += str(mesh_numbers[i]) + ","
    header_string = header_string[:-1] + ")\n"

    for mode in modes:
        for thread in threads[mode]:
            string = ""
            for mesh in mesh_numbers:
                acc = values[mode][thread][mesh][subroutine_name]["MEAN"]
                string += str(acc) + ","

            string = string[:-2]
            final_string = r_format1.format(mode + str(thread), mode + str(thread)) + string + r_format2
            header_string += final_string + "\n"

    string = ""
    for mode in modes:
        for thread in threads[mode]:
            for mesh in mesh_numbers:
                acc = values[mode][thread][mesh][subroutine_name]["STDEV"]
                string += str(acc) + ","
       
    string = string[:-1]
    final_string = "DP <- c(" + string + ")"
    header_string += final_string + "\n"

    string = ""
    for mode in modes:
       for thread in threads[mode]:
            string += mode + str(thread) + ","
    string = string[:-1]
    final_string = "DF <- data.frame(rbind(" + string + "))"
    header_string += final_string + "\n"
    header_string += 'names(DF) <- c("Size", "Version", "Mean")' + "\n"
    header_string += 'DF$Size <- factor(DF$Size, levels = c('
    
    for mesh in mesh_numbers:
        header_string += '"' + str(mesh) + '",'
    header_string = header_string[:-1]
    header_string += '))\n'
    header_string += 'DF$Mean <- as.numeric(as.character(DF$Mean))' + "\n"
    header_string += 'DP <- as.numeric(as.character(DP))' + "\n"
    return header_string + Rcode

def check_graphic_param(graphic):
    possible = ["time", "speedup"]
    for g in possible:
        if graphic == g:
            return True
    return False

def check_routine_param(routine):
    possible = file_keys
    for r in possible:
        if r == routine:
            return True
    return False

def main(argc, argv):
   
    if (argc != 3):
        print("""
Usage: analyser.py <ROUTINE> <GRAPHIC_TYPE>
Possible values to ROUTINE: SHARED GHMATECE GHMATECD LINSOLVE INTEREC1 SIGMAEC
Possible values to GRAPHIC_TYPE: time speedup 
"""     )
        sys.exit(1)

    routine = argv[1]
    graph = argv[2]
    
    statistical_values = {}
    for mode in modes:
        statistical_values[mode] = {}
        for thread in all_threads:
            statistical_values[mode][thread] = {}
            
            if (thread == 1):
                num_execs = 5
            else:
                num_execs = executions
            
            for mesh in mesh_numbers:
                statistical_values[mode][thread][mesh] = {}
                t = get_all_data(mode, mesh, thread, num_execs)
                for key in t:
                    if (len(t[key]) == 0):
                        continue
                    m  = mean(t[key])
                    sd = stdev(t[key])
                   
                    statistical_values[mode][thread][mesh][key] = {}
                    statistical_values[mode][thread][mesh][key]["MEAN"] = m
                    statistical_values[mode][thread][mesh][key]["STDEV"] = 1.96*sd/sqrt(num_execs)

    if (graph == "time"):
        header = generate_r_data(statistical_values, routine)
    elif (graph == "speedup"):
        header = generate_speedup_r_data(statistical_values, routine)

#    header = generate_r_data(statistical_values, "ALL")
    temp_file = open("temp.r", "w")
    temp_file.write(header)
    temp_file.close()

    os.system("Rscript temp.r")

if __name__ == "__main__":
    main(len(sys.argv), sys.argv)

