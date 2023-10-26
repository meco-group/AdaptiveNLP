import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
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

# read "corridors.csv"
corridors = []
with open(folder+"/corridors.csv") as file:
    csv_reader = csv.reader(file)
    header = next(csv_reader)
    for row in csv_reader:
        corridors.append([float(v) for v in row])

corridor_1 = corridors[0]
corridor_2 = corridors[1]
                
###############
### helpers ###
###############
def plotCircle(ax, x, y, r, color, label=None):
    if label is not None:
        circle = mpatches.Circle([x, y], r, color=color, label=label)
    else:
        circle = mpatches.Circle([x, y], r, color=color, label=label)
    ax.add_patch(circle)

def plotCircleMargins(ax, xx, yy, rr, cc, v):
    label_person_added = False
    label_obstacle_added = False
    for i in range(len(xx)):
        if v[i]:
            if cc[i] == 'r':
                if label_obstacle_added:
                    circle = mpatches.Circle((xx[i], yy[i]), rr[i]+4, color=[1.0, 0.8, 0.8])
                else:
                    circle = mpatches.Circle((xx[i], yy[i]), rr[i]+4, color=[1.0, 0.8, 0.8], label="Obstacle-aware zone")
                    label_obstacle_added = True
                ax.add_patch(circle)
            elif cc[i] == 'b':
                if label_person_added:
                    circle = mpatches.Circle((xx[i], yy[i]), rr[i], color=[0.8, 0.8, 1.0])
                else:
                    circle = mpatches.Circle((xx[i], yy[i]), rr[i], color=[0.8, 0.8, 1.0], label="Safety zone")
                    label_person_added = True
                ax.add_patch(circle)
            else:
                ax.add_patch(plt.Circle((xx[i], yy[i]), rr[i], color=str(cc[i])))

def plotCircles(ax, xx, yy, rr, cc, v):
    label_person_added = False
    label_obstacle_added = False
    for i in range(len(xx)):
        if v[i]:
            if cc[i] == 'r':
                if label_obstacle_added:
                    circle = mpatches.Circle((xx[i], yy[i]), rr[i], color=str(cc[i]))
                else:
                    circle = mpatches.Circle((xx[i], yy[i]), rr[i], color=str(cc[i]), label="Obstacle")
                    label_obstacle_added = True
            elif cc[i] == 'b':
                if label_person_added:
                    circle = mpatches.Circle((xx[i], yy[i]), 0.4, color=str(cc[i]))
                else:
                    circle = mpatches.Circle((xx[i], yy[i]), 0.4, color=str(cc[i]), label="Person")
                    label_person_added = True
            else:
                circle = mpatches.Circle((xx[i], yy[i]), rr[i], color=str(cc[i]))
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

#################
### make plot ###
#################
viewing_radius = 10
x0 = [-9.8, 1.0, 0.1, -0.1]

# fig = plt.figure(figsize=(6.4, 3.2))
# axs = plt.gca()
gs = GridSpec(2, 1, height_ratios=[0.7, 0.3], hspace=0.2)
fig = plt.figure(figsize=(6.4, 3.8))
axs = fig.add_subplot(gs[0,0])

# plot viewing radius
plotCircle(axs, x0[0], x0[1], 1000*viewing_radius, [0.9,0.9,0.9], label="Invisible region")
plotCircle(axs, x0[0], x0[1], viewing_radius, [1.0,1.0,1.0])

# plot corridors
plotCorridors(axs, corridors)

# plot obstacles and low-speed zones
plotCircleMargins(axs, circles_x, circles_y, circles_r, circles_colors, 
            [True for i in visible_constraints[0]])
plotCircles(axs, circles_x, circles_y, circles_r, circles_colors, 
            [True for i in visible_constraints[0]])

# plot vehicle
width = 1.5
height = 0.8
center = [x0[0], x0[1]]
tilt = x0[3]
bottom_left = [center[0]-width/2*np.cos(tilt)-height/2*np.sin(tilt), 
                center[1]-width/2*np.sin(tilt)-height/2*np.cos(tilt)]
vehicle = mpatches.Rectangle((bottom_left[0], bottom_left[1]), width,
                              height, angle=tilt/(2*np.pi)*360,
                              color=[0.0,0.0,0.0], zorder=2, label="Vehicle")
axs.add_patch(vehicle)
axs.set_xlim(-18,68)
axs.set_ylim(-8,18)
axs.set_aspect('equal')

# plt.legend(bbox_to_anchor=(0.5, -1.57), loc='lower center', ncol=2)
# Add the legend to the second row of the GridSpec
legend_ax = fig.add_subplot(gs[1, 0])
legend_ax.axis('off')  # Turn off the axes for the legend subplot
legend = axs.legend(bbox_to_anchor=(0.5, -0.3), loc='upper center', ncol=3)
# legend.set_in_layout(False)  # Prevent the legend from affecting layout

# plt.tight_layout()
plt.savefig(figure_folder + "/map.pdf")

# import tikzplotlib
# tikzplotlib.save(figure_folder + "/map.tex")

plt.close()
