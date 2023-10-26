import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
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

# read "circles.csv"
visible_constraints = []

with open(folder+"/circles.csv") as file:
    csv_reader = csv.reader(file)
    header = next(csv_reader)
    row = next(csv_reader)
    circles_x = [float(v) for v in row]
    row = next(csv_reader)
    circles_y = [float(v) for v in row]
    row = next(csv_reader)
    circles_r = [float(v) for v in row]
    row = next(csv_reader)
    circles_colors = [str(v) for v in row]
    
    for row in csv_reader:
        visible_constraints.append([bool(int(v)) for v in row])
# print("circles_x: ", circles_x)
# print("circles_y: ", circles_y)
# print("circles_r: ", circles_r)
# print("circles_colors: ", circles_colors)

# read "corridors.csv"
corridors = []
with open(folder+"/corridors.csv") as file:
    csv_reader = csv.reader(file)
    header = next(csv_reader)
    for row in csv_reader:
        corridors.append([float(v) for v in row])

corridor_1 = corridors[0]
corridor_2= corridors[1]

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
# print(timings_adaptive)
# print(total_timings_adaptive)
# print(nb_constraints_adaptive)
# print(nb_constraints_specific_adaptive)
# print(x0s_adaptive)

# read "results_ref.csv"
timings_ref = []
total_timings_ref = []
nb_constraints_ref = []
x0s_ref = []

with open(folder+"/results_ref.csv") as file:
    csv_reader = csv.reader(file)
    header = next(csv_reader)
    for row in csv_reader:
        timings_ref.append(float(row[0]))
        total_timings_ref.append(float(row[1]))
        nb_constraints_ref.append(float(row[2]))
        x0s_ref.append([float(row[i]) for i in range(3,3+4)])

# read "results_naive.csv"
timings_naive = []
total_timings_naive = []
nb_constraints_naive = []
x0s_naive = []

with open(folder+"/results_naive.csv") as file:
    csv_reader = csv.reader(file)
    header = next(csv_reader)
    for row in csv_reader:
        timings_naive.append(float(row[0]))
        total_timings_naive.append(float(row[1]))
        nb_constraints_naive.append(float(row[2]))
        x0s_naive.append([float(row[i]) for i in range(3,3+4)])

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

# read "solutions_ref.csv"
formatted_solutions_ref = \
    [[[] for j in range(6)] for k in range(len(timings_ref))]
with open(folder+"/solutions_ref.csv") as file:
    csv_reader = csv.reader(file)
    for i in range(len(timings_ref)):
        for j in range(6):
            nb_samples = int(next(csv_reader)[0])
            for k in range(nb_samples):
                formatted_solutions_ref[i][j].append(
                    float(next(csv_reader)[0]))

# read "solutions_naive.csv"
formatted_solutions_naive = \
    [[[] for j in range(6)] for k in range(len(timings_naive))]
with open(folder+"/solutions_naive.csv") as file:
    csv_reader = csv.reader(file)
    for i in range(len(timings_naive)):
        for j in range(6):
            nb_samples = int(next(csv_reader)[0])
            for k in range(nb_samples):
                formatted_solutions_naive[i][j].append(
                    float(next(csv_reader)[0]))
                
###############
### helpers ###
###############
def plotTravelledTrajectory(ax, x0s, c, frame_nb):
    ax.plot([x0s[i][0] for i in range(frame_nb+1)], 
            [x0s[i][1] for i in range(frame_nb+1)], color=c)

def plotCircle(ax, x, y, r, color):
    ax.add_patch(plt.Circle([x, y], r, color=color))

def plotCircleMargins(ax, xx, yy, rr, cc, v):
    for i in range(len(xx)):
        if v[i]:
            if cc[i] == 'r':
                ax.add_patch(plt.Circle((xx[i], yy[i]), rr[i]+4, color=[1.0, 0.8, 0.8]))
            elif cc[i] == 'b':
                ax.add_patch(plt.Circle((xx[i], yy[i]), rr[i], color=[0.8, 0.8, 1.0]))
            else:
                ax.add_patch(plt.Circle((xx[i], yy[i]), rr[i], color=str(cc[i])))

def plotCircles(ax, xx, yy, rr, cc, v):
    for i in range(len(xx)):
        if v[i]:
            if cc[i] == 'r':
                ax.add_patch(plt.Circle((xx[i], yy[i]), rr[i], color=str(cc[i])))
            elif cc[i] == 'b':
                ax.add_patch(plt.Circle((xx[i], yy[i]), 0.4, color=str(cc[i])))
            else:
                ax.add_patch(plt.Circle((xx[i], yy[i]), rr[i], color=str(cc[i])))

def plotCorridors(ax, corridors):
    for corridor in corridors:
        ax.plot([corridor[0], corridor[1], corridor[1], corridor[0], 
                  corridor[0]], 
                 [corridor[2], corridor[2], corridor[3], corridor[3], 
                  corridor[2]],
                  color="g")

#################
### make plot ###
#################
viewing_radius = 10
frame_nb = len(timings_adaptive) - 1
# fig_summary, axs_summary = plt.subplots(2,1)
gs = GridSpec(2, 2, height_ratios=[0.5, 0.5], width_ratios=[0.7, 0.3], hspace=0.2)
fig = plt.figure()
axs_summary_0 = fig.add_subplot(gs[0,0])

# axs_summary[0].plot(range(frame_nb), nb_constraints_naive[:frame_nb], 
#         color="orange", linewidth=2.5)
axs_summary_0.plot(range(frame_nb), nb_constraints_ref[:frame_nb], 
            color=[0.5,0.5,0.5], label="Maximal number")
axs_summary_0.plot(range(frame_nb), 
            nb_constraints_specific_adaptive[3][:frame_nb], "purple")
axs_summary_0.plot(range(frame_nb), 
            nb_constraints_specific_adaptive[2][:frame_nb], 
            "mediumslateblue")
axs_summary_0.plot(range(frame_nb), 
            nb_constraints_specific_adaptive[1][:frame_nb], "navy")
axs_summary_0.plot(range(frame_nb), 
            nb_constraints_specific_adaptive[0][:frame_nb], "blue")
axs_summary_0.set_ylim(min(nb_constraints_adaptive + nb_constraints_ref + 
                    nb_constraints_naive) - 10, 
                max(nb_constraints_adaptive + nb_constraints_ref +
                    nb_constraints_naive) + 10)

axs_summary_0.fill_between(range(frame_nb), 0, nb_constraints_specific_adaptive[0][:frame_nb], nb_constraints_ref[:frame_nb], color='blue', alpha=0.5, label='Basic constraints')
axs_summary_0.fill_between(range(frame_nb), nb_constraints_specific_adaptive[0][:frame_nb], nb_constraints_specific_adaptive[1][:frame_nb], color='navy', alpha=0.5, label='Corridor constraints')
axs_summary_0.fill_between(range(frame_nb), nb_constraints_specific_adaptive[1][:frame_nb], nb_constraints_specific_adaptive[2][:frame_nb], color='mediumslateblue', alpha=0.5, label='Obstacle constraints')
axs_summary_0.fill_between(range(frame_nb), nb_constraints_specific_adaptive[2][:frame_nb], nb_constraints_specific_adaptive[3][:frame_nb], color='purple', alpha=0.5, label='Safety constraints')

axs_summary_0.set_xlim(0, len(nb_constraints_adaptive))
axs_summary_0.set_ylabel('number of\nconstraints')

legend_ax = fig.add_subplot(gs[1,1])
legend_ax.axis('off')
legend = axs_summary_0.legend(bbox_to_anchor=(1.8, 0.5), loc="center right")

axs_summary_1 = fig.add_subplot(gs[1,0])

# plot computation times
axs_summary_1.plot(timings_ref[:frame_nb], color="gray", linewidth=0.4)
axs_summary_1.plot(total_timings_ref[:frame_nb], color="gray", label="CasADi Opti 1")
axs_summary_1.plot(timings_naive[:frame_nb], color="orange", linewidth=0.4)
axs_summary_1.plot(total_timings_naive[:frame_nb], color="orange", label="CasADi Opti 2")
axs_summary_1.plot(timings_adaptive[:frame_nb], color="royalblue", linewidth=0.4)
axs_summary_1.plot(total_timings_adaptive[:frame_nb], color="royalblue", label="AdaptiveNLP")

axs_summary_1.set_xlim(0, len(timings_adaptive))
axs_summary_1.set_ylim(0.0, max(total_timings_adaptive+total_timings_naive+total_timings_ref))

axs_summary_1.set_xlabel("iteration number")
axs_summary_1.set_ylabel("time [ms]")

legend_ax = fig.add_subplot(gs[1,1])
legend_ax.axis('off')
legend = axs_summary_1.legend(bbox_to_anchor=(1.7, 0.5), loc="center right")

plt.subplots_adjust(left=0.15, bottom=0.15)
# plt.tight_layout()
plt.savefig(figure_folder + "/animation_summary.pdf")

plt.close()


gs = GridSpec(2, 1, height_ratios=[0.6, 0.4], hspace=0.2)
fig = plt.figure(figsize=(7.4, 4.2))
axs_summary_0 = fig.add_subplot(gs[0,0])

iterations_to_show = [25, 65, 105]
axs_summary_0.vlines(iterations_to_show, 0, 400, "k", linestyle="--", 
                     alpha=0.2)

basic_color = "gray" # blue
corridor_color = "green" # navy
obstacle_color = "red" # mediumslateblue
safety_color = "blue" # purple

# axs_summary[0].plot(range(frame_nb), nb_constraints_naive[:frame_nb], 
#         color="orange", linewidth=2.5)
axs_summary_0.plot(range(frame_nb), nb_constraints_ref[:frame_nb], 
            color="k", linestyle="--",  linewidth=2.0,
            label="Maximal number of constraints")
axs_summary_0.plot(range(frame_nb), 
            nb_constraints_specific_adaptive[3][:frame_nb], safety_color)
axs_summary_0.plot(range(frame_nb), 
            nb_constraints_specific_adaptive[2][:frame_nb], 
            obstacle_color)
axs_summary_0.plot(range(frame_nb), 
            nb_constraints_specific_adaptive[1][:frame_nb], corridor_color)
axs_summary_0.plot(range(frame_nb), 
            nb_constraints_specific_adaptive[0][:frame_nb], basic_color)
axs_summary_0.set_ylim(min(nb_constraints_adaptive + nb_constraints_ref + 
                    nb_constraints_naive) - 10, 
                max(nb_constraints_adaptive + nb_constraints_ref +
                    nb_constraints_naive) + 10)

alpha = 0.3
axs_summary_0.fill_between(range(frame_nb), 0, nb_constraints_specific_adaptive[0][:frame_nb], nb_constraints_ref[:frame_nb], color=basic_color, alpha=alpha, label='Basic constraints')
axs_summary_0.fill_between(range(frame_nb), nb_constraints_specific_adaptive[0][:frame_nb], nb_constraints_specific_adaptive[1][:frame_nb], color=corridor_color, alpha=alpha, label='Corridor constraints')
axs_summary_0.fill_between(range(frame_nb), nb_constraints_specific_adaptive[1][:frame_nb], nb_constraints_specific_adaptive[2][:frame_nb], color=obstacle_color, alpha=alpha, label='Obstacle constraints')
axs_summary_0.fill_between(range(frame_nb), nb_constraints_specific_adaptive[2][:frame_nb], nb_constraints_specific_adaptive[3][:frame_nb], color=safety_color, alpha=alpha, label='Safety constraints')

axs_summary_0.set_xlim(0, len(nb_constraints_adaptive)-2)
axs_summary_0.set_xlabel('iteration')
axs_summary_0.set_ylabel('number of\nconstraints')

legend_ax = fig.add_subplot(gs[1,0])
legend_ax.axis('off')
legend = axs_summary_0.legend(bbox_to_anchor=(0.5, -0.4), loc="upper center", ncols=2)

axs_summary_0.set_yticks([136, 200, 300, 376])

plt.subplots_adjust(left=0.15, bottom=0.15)
# plt.tight_layout()
plt.savefig(figure_folder + "/animation_constraint_summary.pdf")

plt.close()
