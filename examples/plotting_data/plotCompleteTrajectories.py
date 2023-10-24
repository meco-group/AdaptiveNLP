import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import csv
import numpy as np
import pickle

import pathlib
folder = str(pathlib.Path(__file__).parent.resolve()) # "examples/plotting_data"
figure_folder = folder + "/../figures"

def latexify():
    fontsize = 9
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

latexify()

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
def plotTravelledTrajectory(ax, x0s, c, frame_nb, label, linewidth):
    ax.plot([x0s[i][0] for i in range(frame_nb+1)], 
            [x0s[i][1] for i in range(frame_nb+1)], color=c, label=label, 
            linewidth=linewidth)

def plotCircle(ax, x, y, r, color, label=None, alpha=1.0):
    if label is not None:
        circle = mpatches.Circle([x, y], r, color=color, label=label, 
                                 alpha=alpha)
    else:
        circle = mpatches.Circle([x, y], r, color=color, label=label,
                                 alpha=alpha)
    ax.add_patch(circle)

def plotCircleMargins(ax, xx, yy, rr, cc, v):
    label_person_added = False
    label_obstacle_added = False
    for i in range(len(xx)):
        if v[i]:
            if cc[i] == 'r':
                if label_obstacle_added:
                    circle = mpatches.Circle((xx[i], yy[i]), rr[i]+4, 
                                             color=[1.0, 0.8, 0.8])
                else:
                    circle = mpatches.Circle((xx[i], yy[i]), rr[i]+4, 
                                             color=[1.0, 0.8, 0.8], 
                                             label="Obstacle-aware zone")
                    label_obstacle_added = True
                ax.add_patch(circle)
            elif cc[i] == 'b':
                if label_person_added:
                    circle = mpatches.Circle((xx[i], yy[i]), rr[i], 
                                             color=[0.8, 0.8, 1.0])
                else:
                    circle = mpatches.Circle((xx[i], yy[i]), rr[i], 
                                             color=[0.8, 0.8, 1.0], 
                                             label="Safety zone")
                    label_person_added = True
                ax.add_patch(circle)
            else:
                ax.add_patch(plt.Circle((xx[i], yy[i]), rr[i], 
                                        color=str(cc[i])))

def plotCircles(ax, xx, yy, rr, cc, v):
    label_person_added = False
    label_obstacle_added = False
    for i in range(len(xx)):
        if v[i]:
            if cc[i] == 'r':
                if label_obstacle_added:
                    circle = mpatches.Circle((xx[i], yy[i]), rr[i], 
                                             color=str(cc[i]))
                else:
                    circle = mpatches.Circle((xx[i], yy[i]), rr[i], 
                                             color=str(cc[i]), 
                                             label="Obstacle")
                    label_obstacle_added = True
            elif cc[i] == 'b':
                if label_person_added:
                    circle = mpatches.Circle((xx[i], yy[i]), 0.4, 
                                             color=str(cc[i]))
                else:
                    circle = mpatches.Circle((xx[i], yy[i]), 0.4, 
                                             color=str(cc[i]), 
                                             label="Person")
                    label_person_added = True
            else:
                circle = mpatches.Circle((xx[i], yy[i]), rr[i], 
                                         color=str(cc[i]))
            ax.add_patch(circle)

def plotCorridors(ax, corridors):
    for corridor in corridors:
        if corridor == corridors[0]:
            ax.plot([corridor[0], corridor[1], corridor[1], corridor[0], 
                    corridor[0]], 
                    [corridor[2], corridor[2], corridor[3], corridor[3], 
                    corridor[2]],
                    color="g", label="Corridor")
        else:
            ax.plot([corridor[0], corridor[1], corridor[1], corridor[0], 
                    corridor[0]], 
                    [corridor[2], corridor[2], corridor[3], corridor[3], 
                    corridor[2]],
                    color="g")
            
def plotVehicle(ax, x0s, iteration, label=None):
    width = 1.5
    height = 0.8
    center = [x0s[iteration][0], x0s[iteration][1]]
    tilt = x0s[iteration][3]
    bottom_left = [center[0]-width/2*np.cos(tilt)-height/2*np.sin(tilt), 
                   center[1]-width/2*np.sin(tilt)-height/2*np.cos(tilt)]
    ax.add_patch(plt.Rectangle((bottom_left[0], bottom_left[1]), width,
                                  height, angle=tilt/(2*np.pi)*360,
                                  color=[0.0,0.0,0.0], zorder=2, 
                                  label=label))

#################
### make plot ###
#################
viewing_radius = 10
x0 = [-9.8, 1.0, 0.1, -0.1]

iterations_to_show = [25, 65, 105]
heights_to_print = [7, 7, -0.5]

gs = GridSpec(2, 1, height_ratios=[0.9, 0.1], hspace=0.2)
fig = plt.figure(figsize=(7, 3.6))
axs = fig.add_subplot(gs[0,0])

# plot viewing radii
plotCircle(axs, 0, 0, 10000, color=[0.9,0.9,0.9], label="Invisible region")
for i in iterations_to_show:
    plotCircle(axs, x0s_adaptive[i][0], x0s_adaptive[i][1], viewing_radius, 
               "white")

# plot corridors
plotCorridors(axs, corridors)

# plot obstacles and low-speed zones
plotCircleMargins(axs, circles_x, circles_y, circles_r, circles_colors, 
            [True for i in visible_constraints[0]])

for i in iterations_to_show:
    plotCircle(axs, x0s_adaptive[i][0], x0s_adaptive[i][1], viewing_radius,
               color="white", alpha=0.5)

plotCircles(axs, circles_x, circles_y, circles_r, circles_colors, 
            [True for i in visible_constraints[0]])

frame_nb = len(x0s_naive)-1
plotTravelledTrajectory(axs, x0s_ref, [0.5,0.5,0.5], frame_nb, 
                        # label="casadi::Opti 1",
                        label=None,
                        linewidth=1.5)
plotTravelledTrajectory(axs, x0s_naive, "orange", frame_nb, 
                        # label="casadi::Opti 2", 
                        label=None,
                        linewidth=1.5)
plotTravelledTrajectory(axs, x0s_adaptive, "royalblue", frame_nb, 
                        label="Travelled path", 
                        linewidth=1.5)

comp1 = sum([(x0s_naive[i][0] - x0s_ref[i][0])**2 
             for i in range(frame_nb+1)]) + \
        sum([(x0s_naive[i][1] - x0s_ref[i][1])**2 
             for i in range(frame_nb+1)])
comp2 = sum([(x0s_naive[i][0] - x0s_adaptive[i][0])**2 
             for i in range(frame_nb+1)]) + \
        sum([(x0s_naive[i][1] - x0s_adaptive[i][1])**2 
             for i in range(frame_nb+1)])
comp3 = sum([(x0s_adaptive[i][0] - x0s_ref[i][0])**2 
             for i in range(frame_nb+1)]) + \
        sum([(x0s_adaptive[i][1] - x0s_ref[i][1])**2 
             for i in range(frame_nb+1)])
print(comp1, comp2, comp3)


# plot vehicle
for i, h in zip(iterations_to_show, heights_to_print):
    if i == iterations_to_show[0]:
        plotVehicle(axs, x0s_adaptive, i, label="Vehicle")
    else:
        plotVehicle(axs, x0s_adaptive, i)
    axs.text(x0s_adaptive[i][0], h, str(i), ha="center", fontsize=10)

# plt.legend(bbox_to_anchor=(0.5, -0.57), loc='lower center', ncol=3)
axs.set_xlim(-18,68)
axs.set_ylim(-8,18)
axs.set_aspect('equal')


legend_ax = fig.add_subplot(gs[1, 0])
legend_ax.axis('off')
handles, labels = axs.get_legend_handles_labels()
# order = [0,1,2,3,4,5,6,7]
order = [7,0,1,6,2,4,3,5]
# order = [6, 0, 7, 2, 4, 1, 3, 5]
legend = axs.legend([handles[idx] for idx in order], 
                    [labels[idx] for idx in order], 
                    bbox_to_anchor=(0.5, -0.17), loc='upper center', ncol=4)

for t in legend.texts:
    if "\n" in t.get_text():
        t.set_verticalalignment("center")
    # t.set_multialignment('bottom')

plt.savefig(figure_folder + "/final_result.pdf")

plt.close()

fig, axs = plt.subplots(2, 1, figsize=(7, 3.8))
# fig, axs = plt.subplots(1, 2, figsize=(7, 2.0))
axs[0].plot([formatted_solutions_ref[i][4][0] for i in range(frame_nb)], color="gray")
axs[1].plot([formatted_solutions_ref[i][5][0] for i in range(frame_nb)], color="gray")
axs[0].plot([formatted_solutions_naive[i][4][0] for i in range(frame_nb)], color="orange")
axs[1].plot([formatted_solutions_naive[i][5][0] for i in range(frame_nb)], color="orange")
axs[0].plot([formatted_solutions_adaptive[i][4][0] for i in range(frame_nb)], color="royalblue")
axs[1].plot([formatted_solutions_adaptive[i][5][0] for i in range(frame_nb)], color="royalblue")

axs[0].set_ylabel("$a$")
axs[1].set_ylabel("$\delta$")
axs[1].set_xlabel("iteration number")

axs[0].set_xlim(0, frame_nb-1)
axs[1].set_xlim(0, frame_nb-1)

# plt.tight_layout()
plt.subplots_adjust(left=0.15, bottom=0.2, hspace=0.35)

# plt.savefig(figure_folder + "/final_result_controls.png", dpi=400)
plt.savefig(figure_folder + "/final_result_controls.pdf")