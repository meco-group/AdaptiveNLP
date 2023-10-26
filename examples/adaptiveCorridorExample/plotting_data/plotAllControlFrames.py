import matplotlib.pyplot as plt
import csv
import numpy as np
import pickle

import pathlib
folder = str(pathlib.Path(__file__).parent.resolve())
figure_folder = folder + "/../figures"


def latexify():
    fontsize = 15
    params = {'backend': 'ps',
              'axes.labelsize': fontsize,
              'axes.titlesize': fontsize,
              'legend.fontsize': fontsize,
              'xtick.labelsize': fontsize,
              'ytick.labelsize': fontsize,
              'text.usetex': True,
              'font.family': 'serif',
              'figure.figsize': [7,5],
              'text.latex.preamble': [r'\usepackage{bm}'],
              }

    plt.rcParams.update(params)

try:
    latexify()
except:
    print("Unable to latexify the figure.")

#####################
### read all data ###
#####################

# read "results_adaptive.csv"
timings_adaptive = []
total_timings_adaptive = []
nb_constraints_adaptive = []
nb_constraints_specific_adaptive = [[], [], [], []]
x0s_adaptive = []

with open(folder+"/results_adaptive.csv") as file:
    csv_reader = csv.reader(file)
    header = next(csv_reader)
    for row in csv_reader:
        timings_adaptive.append(float(row[0]))
        total_timings_adaptive.append(float(row[1]))
        nb_constraints_adaptive.append(float(row[2]))
        for j in range(4):
            nb_constraints_specific_adaptive[j].append(float(row[3+j]))
        x0s_adaptive.append([float(row[i]) for i in range(7,7+4)])

# read "solutions_adaptive.csv"
formatted_solutions_adaptive = \
    [[[] for j in range(6)] for k in range(len(timings_adaptive))]
with open(folder+"/solutions_adaptive.csv") as file:
    csv_reader = csv.reader(file)
    for i in range(len(timings_adaptive)):
        for j in range(6):
            nb_samples = int(next(csv_reader)[0])
            for k in range(nb_samples):
                formatted_solutions_adaptive[i][j].append(
                    float(next(csv_reader)[0]))

#################
### make plot ###
#################

for frame_nb in range(len(timings_adaptive)):
    print("(", frame_nb, "/", len(timings_adaptive) - 1, ")")
    fig, axs = plt.subplots(2,1,figsize=(7,6))

    axs[0].plot([formatted_solutions_adaptive[i][4][0] for i in range(frame_nb)],
                color='royalblue')
    axs[0].set_xlim(0, len(timings_adaptive))
    axs[0].set_ylim(min([formatted_solutions_adaptive[i][4][0] for i in 
                         range(len(formatted_solutions_adaptive))]) - 0.2,
                    max([formatted_solutions_adaptive[i][4][0] for i in 
                         range(len(formatted_solutions_adaptive))]) + 0.2)
    axs[0].set_ylabel("$a [$m$/$s$^2]$")

    axs[1].plot([formatted_solutions_adaptive[i][5][0] for i in range(frame_nb)],
                color='royalblue')
    axs[1].set_xlim(0, len(timings_adaptive))
    axs[1].set_ylim(min([formatted_solutions_adaptive[i][5][0] for i in 
                         range(len(formatted_solutions_adaptive))]) - 0.2,
                    max([formatted_solutions_adaptive[i][5][0] for i in 
                         range(len(formatted_solutions_adaptive))]) + 0.2)
    axs[1].set_ylabel("$\delta [$rad$^{-1}]$")

    axs[1].set_xlabel("iteration number")
    
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.2, left=0.15, right=0.85, top=0.9)

    plt.savefig(figure_folder + "/animation_controls/animation_frame_" + str(frame_nb) + ".png", dpi=200)
    
    plt.close()