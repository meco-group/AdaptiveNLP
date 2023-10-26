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

def plotCircle(ax, x, y, r, color, alpha=1):
    ax.add_patch(plt.Circle([x, y], r, color=color, alpha=alpha))

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
for frame_nb in range(len(timings_adaptive)):
    print("(", frame_nb, "/", len(timings_adaptive) - 1, ")")
    fig, axs = plt.subplots(3,1,figsize=(7,6))

    # plot viewing radius
    plotCircle(axs[0], x0s_adaptive[frame_nb][0], x0s_adaptive[frame_nb][1], 
               1000*viewing_radius, [0.9,0.9,0.9])
    plotCircle(axs[0], x0s_adaptive[frame_nb][0], x0s_adaptive[frame_nb][1], 
               viewing_radius, [1.0,1.0,1.0])

    # plot travelled trajectory and computed solutions
    plotTravelledTrajectory(axs[0], x0s_naive, "orange", frame_nb)
    axs[0].plot(formatted_solutions_naive[frame_nb][0], 
                formatted_solutions_naive[frame_nb][1], color='orange')
    
    plotTravelledTrajectory(axs[0], x0s_ref, [0.5,0.5,0.5], frame_nb)
    axs[0].plot(formatted_solutions_ref[frame_nb][0], 
                formatted_solutions_ref[frame_nb][1], color=[0.5,0.5,0.5])
    
    plotTravelledTrajectory(axs[0], x0s_adaptive, "royalblue", frame_nb)
    axs[0].plot(formatted_solutions_adaptive[frame_nb][0], 
                formatted_solutions_adaptive[frame_nb][1], color='royalblue')
    axs[0].set_xlim(-15,65)
    axs[0].set_ylim(-5,12)
    axs[0].set_aspect('equal')

    # plot corridors
    plotCorridors(axs[0], corridors)

    # plot obstacles and low-speed zones
    plotCircleMargins(axs[0], circles_x, circles_y, circles_r, circles_colors, 
                visible_constraints[frame_nb])
    plotCircle(axs[0], x0s_adaptive[frame_nb][0], x0s_adaptive[frame_nb][1], 
               viewing_radius, [1.0,1.0,1.0], alpha=0.5)
    plotCircles(axs[0], circles_x, circles_y, circles_r, circles_colors, 
                visible_constraints[frame_nb])

    # plot vehicle
    width = 1.5
    height = 0.8
    center = [x0s_adaptive[frame_nb][0], x0s_adaptive[frame_nb][1]]
    tilt = x0s_adaptive[frame_nb][3]
    bottom_left = [center[0]-width/2*np.cos(tilt)-height/2*np.sin(tilt), 
                   center[1]-width/2*np.sin(tilt)-height/2*np.cos(tilt)]
    axs[0].add_patch(plt.Rectangle((bottom_left[0], bottom_left[1]), width,
                                  height, angle=tilt/(2*np.pi)*360,
                                  color=[0.0,0.0,0.0], zorder=2))
    

    # plot nunber of constraints
    basic_color = "gray"
    corridor_color = "green"
    obstacle_color = "red"
    safety_color = "blue" 

    # axs[1].plot(range(frame_nb), nb_constraints_naive[:frame_nb], 
    #         color="orange", linewidth=2.5)
    axs[1].plot(range(frame_nb), nb_constraints_ref[:frame_nb], 
                color="k", label="maximal number of constraints")
    axs[1].plot(range(frame_nb), 
                nb_constraints_specific_adaptive[3][:frame_nb], safety_color)
    axs[1].plot(range(frame_nb), 
                nb_constraints_specific_adaptive[2][:frame_nb], 
                obstacle_color)
    axs[1].plot(range(frame_nb), 
                nb_constraints_specific_adaptive[1][:frame_nb], corridor_color)
    axs[1].plot(range(frame_nb), 
                nb_constraints_specific_adaptive[0][:frame_nb], basic_color)
    axs[1].set_ylim(min(nb_constraints_adaptive + nb_constraints_ref + 
                        nb_constraints_naive) - 10, 
                    max(nb_constraints_adaptive + nb_constraints_ref +
                        nb_constraints_naive) + 10)

    alpha = 0.3
    axs[1].fill_between(range(frame_nb), 0, nb_constraints_specific_adaptive[0][:frame_nb], nb_constraints_ref[:frame_nb], color=basic_color, alpha=alpha, label='basic constraints')
    axs[1].fill_between(range(frame_nb), nb_constraints_specific_adaptive[0][:frame_nb], nb_constraints_specific_adaptive[1][:frame_nb], color=corridor_color, alpha=alpha, label='corridor constraints')
    axs[1].fill_between(range(frame_nb), nb_constraints_specific_adaptive[1][:frame_nb], nb_constraints_specific_adaptive[2][:frame_nb], color=obstacle_color, alpha=alpha, label='obstacle constraints')
    axs[1].fill_between(range(frame_nb), nb_constraints_specific_adaptive[2][:frame_nb], nb_constraints_specific_adaptive[3][:frame_nb], color=safety_color, alpha=alpha, label='safety constraints')
    axs[1].set_ylim(min(nb_constraints_adaptive + nb_constraints_ref + 
                         nb_constraints_naive) - 10, 
                    max(nb_constraints_adaptive + nb_constraints_ref +
                        nb_constraints_naive) + 10)
    axs[1].set_xlim(0, len(nb_constraints_adaptive)-2)
    axs[1].set_ylabel('number of\nconstraints')

    # plot computation times
    # axs[2].plot(timings_naive[:frame_nb], color="orange", linewidth=0.4)
    axs[2].plot(total_timings_naive[:frame_nb], color="orange")
    # axs[2].plot(timings_ref[:frame_nb], color="gray", linewidth=0.4)
    axs[2].plot(total_timings_ref[:frame_nb], color="gray")
    # axs[2].plot(timings_adaptive[:frame_nb], color="royalblue", linewidth=0.4)
    axs[2].plot(total_timings_adaptive[:frame_nb], color="royalblue")

    axs[2].set_xlim(0, len(timings_adaptive)-2)
    axs[2].set_ylim(0.0, max(total_timings_adaptive+total_timings_naive+total_timings_ref))

    axs[2].set_xlabel("iteration number")
    axs[2].set_ylabel("time [ms]")

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.2, left=0.15, right=0.85, top=0.9)

    # plt.show()
    plt.savefig(figure_folder + "/animation/animation_frame_" + str(frame_nb) + ".png", dpi=200)
    # plt.savefig("figures/animation/animation_frame_" + str(frame_nb) + ".pdf")

    if (frame_nb == len(timings_adaptive)-1):
        fig_summary, axs_summary = plt.subplots(2,1)
        # axs_summary[0] = pickle.loads(pickle.dumps(axs[1]))
        # axs_summary[1] = pickle.loads(pickle.dumps(axs[2]))
        
        # axs_summary[0].plot(axs[1].lines[0].get_xdata(), axs[1].lines[0].get_ydata())
        # axs_summary[1].plot(axs[2].lines[0].get_xdata(), axs[2].lines[0].get_ydata())
        # axs_summary[0].set_xlabel(axs[1].get_xlabel())
        # axs_summary[0].set_ylabel(axs[1].get_ylabel())
        # axs_summary[1].set_xlabel(axs[2].get_xlabel())
        # axs_summary[1].set_ylabel(axs[2].get_ylabel())

        # fig_summary = pickle.loads(pickle.dumps(fig))
        # fig_summary.delaxes(fig_summary.axes[0])
        
        # fig_summary, axs_summary = plt.subplots(2, 1)
        # axs_summary[0].plot(axs[1].lines[0].get_xdata(), axs[1].lines[0].get_ydata())
        # axs_summary[1].plot(axs[2].lines[0].get_xdata(), axs[2].lines[0].get_ydata())
        # axs_summary[0].set_xlabel(axs[1].get_xlabel())
        # axs_summary[0].set_ylabel(axs[1].get_ylabel())
        # axs_summary[1].set_xlabel(axs[2].get_xlabel())
        # axs_summary[1].set_ylabel(axs[2].get_ylabel())

        basic_color = "gray"
        corridor_color = "green"
        obstacle_color = "red"
        safety_color = "blue" 

        # axs_summary[0].plot(range(frame_nb), nb_constraints_naive[:frame_nb], 
        #         color="orange", linewidth=2.5)
        axs_summary[0].plot(range(frame_nb), nb_constraints_ref[:frame_nb], 
                    color="k", label="maximal number of constraints")
        axs_summary[0].plot(range(frame_nb), 
                    nb_constraints_specific_adaptive[3][:frame_nb], safety_color)
        axs_summary[0].plot(range(frame_nb), 
                    nb_constraints_specific_adaptive[2][:frame_nb], 
                    obstacle_color)
        axs_summary[0].plot(range(frame_nb), 
                    nb_constraints_specific_adaptive[1][:frame_nb], corridor_color)
        axs_summary[0].plot(range(frame_nb), 
                    nb_constraints_specific_adaptive[0][:frame_nb], basic_color)
        axs_summary[0].set_ylim(min(nb_constraints_adaptive + nb_constraints_ref + 
                            nb_constraints_naive) - 10, 
                        max(nb_constraints_adaptive + nb_constraints_ref +
                            nb_constraints_naive) + 10)

        alpha = 0.3
        axs_summary[0].fill_between(range(frame_nb), 0, nb_constraints_specific_adaptive[0][:frame_nb], nb_constraints_ref[:frame_nb], color=basic_color, alpha=alpha, label='basic constraints')
        axs_summary[0].fill_between(range(frame_nb), nb_constraints_specific_adaptive[0][:frame_nb], nb_constraints_specific_adaptive[1][:frame_nb], color=corridor_color, alpha=alpha, label='corridor constraints')
        axs_summary[0].fill_between(range(frame_nb), nb_constraints_specific_adaptive[1][:frame_nb], nb_constraints_specific_adaptive[2][:frame_nb], color=obstacle_color, alpha=alpha, label='obstacle constraints')
        axs_summary[0].fill_between(range(frame_nb), nb_constraints_specific_adaptive[2][:frame_nb], nb_constraints_specific_adaptive[3][:frame_nb], color=safety_color, alpha=alpha, label='safety constraints')

        # plot computation times
        # axs_summary[1].plot(timings_naive[:frame_nb], color="orange", linewidth=0.4)
        axs_summary[1].plot(total_timings_naive[:frame_nb], color="orange")
        # axs_summary[1].plot(timings_ref[:frame_nb], color="gray", linewidth=0.4)
        axs_summary[1].plot(total_timings_ref[:frame_nb], color="gray")
        # axs_summary[1].plot(timings_adaptive[:frame_nb], color="royalblue", linewidth=0.4)
        axs_summary[1].plot(total_timings_adaptive[:frame_nb], color="royalblue")

        axs_summary[1].set_xlim(0, len(timings_adaptive))
        axs_summary[1].set_ylim(0.0, max(total_timings_adaptive+total_timings_naive+total_timings_ref))

        axs_summary[1].set_xlabel("iteration number")
        axs_summary[1].set_ylabel("time [ms]")
        
        plt.savefig(figure_folder + "/animation_summary.png", dpi=400)

    plt.close()
