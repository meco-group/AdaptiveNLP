import matplotlib.pyplot as plt
import csv
import numpy as np
from matplotlib.gridspec import GridSpec

folder = "plotting_data"

def latexify():
    fontsize = 15
    params = {'backend': 'ps',
              'font.size': fontsize,
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

Ns = []
t_solve_casadi1 = []
t_total_casadi1 = []
t_solve_casadi2 = []
t_total_casadi2 = []
t_solve_adaptive = []
t_total_adaptive = []

with open(folder+"/scaling_data.csv") as file:
    csv_reader = csv.reader(file)
    header = next(csv_reader)
    for row in csv_reader:
        Ns.append(int(row[0]))
        t_solve_casadi1.append(float(row[1]))
        t_total_casadi1.append(float(row[2]))
        t_solve_casadi2.append(float(row[3]))
        t_total_casadi2.append(float(row[4]))
        t_solve_adaptive.append(float(row[5]))
        t_total_adaptive.append(float(row[6]))
        

gs = GridSpec(2, 1, height_ratios=[0.8, 0.2], hspace=0.2)
fig = plt.figure()
ax = fig.add_subplot(gs[0,0])

ax.plot(Ns, t_solve_casadi1, color = "royalblue", marker="^", markersize=10,
         linewidth=3.0)
ax.plot(Ns, t_solve_casadi2, color = "royalblue", marker="o", markersize=10,
         linewidth=2.0)
ax.plot(Ns, t_solve_adaptive, color = "royalblue", marker="s", markersize=10,
         linewidth=3.0)

ax.plot(Ns, t_total_casadi1, color = "darkgray", marker="^", markersize=10,
         linewidth=2.0, linestyle="-")
ax.plot(Ns, t_total_casadi2, color = "darkgray", marker="o", markersize=10,
         linewidth=2.0, linestyle="-")
ax.plot(Ns, t_total_adaptive, color = "darkgray", marker="s", markersize=10,
         linewidth=2.0, linestyle="-")

ax.set_xlabel("$N$")
ax.set_ylabel("time [ms]")

time_categories = ["$t_{solve}$", "$t_{total}$"]
time_category_colors = ["royalblue", "darkgray"]
time_category_linestyles = ["-", "-"]

methods = ["CasADi Opti 1", "CasADi Opti 2", "AdaptiveNLP"]
method_markers = ["v", "o", "s"]

legend_labels = time_categories + methods
legend_handles = []
for label in legend_labels:
    if label in time_categories:
        handle = plt.Line2D([], [], 
                    color=time_category_colors[time_categories.index(label)],
                    linestyle=time_category_linestyles[time_categories.index(label)])
    else:
        handle = plt.Line2D([], [], color="k",
                    marker=method_markers[methods.index(label)],
                    markersize=10)
    legend_handles.append(handle)

legend_ax = fig.add_subplot(gs[1, 0])
legend_ax.axis('off')
legend = ax.legend(legend_handles, legend_labels, bbox_to_anchor=(0.5, -0.2), ncol=3, loc="upper center")

plt.savefig("figures/scaling_old.pdf")


time_categories = ["$t_{solve}$", "$t_{update}$", "$t_{total}$"]
time_category_colors = ["royalblue", "gray", "k"]
time_category_linestyles = ["-", "-", "-"]

methods = ["CasADi Opti 1", "CasADi Opti 2", "AdaptiveNLP"]
method_markers = ["v", "o", "s"]

gs = GridSpec(2, 3, height_ratios=[0.83, 0.17], width_ratios=[0.5, 0.5, 0.5], hspace=0.2)
fig = plt.figure(figsize=(6.8, 4))
ax1 = fig.add_subplot(gs[0,0])
alpha=1.0

# make plot
ax1.plot(Ns, t_solve_casadi1, color=time_category_colors[0], 
         marker=method_markers[0], markersize=10, linewidth=2.0, alpha=alpha)
ax1.plot(Ns, t_solve_casadi2, color = time_category_colors[0], 
         marker=method_markers[1], markersize=10, linewidth=2.0, alpha=alpha)
ax1.plot(Ns, t_solve_adaptive, color = time_category_colors[0], 
         marker=method_markers[2], markersize=8, linewidth=2.0, alpha=alpha)
ax1.set_xlabel("$N$")
ax1.set_ylabel("time [ms]")
ax1.set_title("$t_\mathrm{solve}$")

ax2 = fig.add_subplot(gs[0,1])
t_update_casadi1 = [t_total_casadi1[i]-t_solve_casadi1[i] for i in 
                    range(len(t_total_casadi1))]
t_update_casadi2 = [t_total_casadi2[i]-t_solve_casadi2[i] for i in 
                    range(len(t_total_casadi2))]
t_update_adaptive = [t_total_adaptive[i]-t_solve_adaptive[i] for i in 
                    range(len(t_total_adaptive))]
ax2.plot(Ns, t_update_casadi1, color = time_category_colors[1], 
         marker=method_markers[0], markersize=10, linewidth=2.0, linestyle="-", 
         alpha=alpha)
ax2.plot(Ns, t_update_casadi2, color = time_category_colors[1], 
         marker=method_markers[1], markersize=10, linewidth=2.0, linestyle="-", 
         alpha=alpha)
ax2.plot(Ns, t_update_adaptive, color = time_category_colors[1], 
         marker=method_markers[2], markersize=8, linewidth=2.0, linestyle="-", 
         alpha=alpha)
ax2.set_xlabel("$N$")
ax2.set_title("$t_\mathrm{update}$")

ax3 = fig.add_subplot(gs[0,2])
ax3.plot(Ns, t_total_casadi1, color = time_category_colors[2], 
         marker=method_markers[0], markersize=10, linewidth=2.0, linestyle="-")
ax3.plot(Ns, t_total_casadi2, color = time_category_colors[2], 
         marker=method_markers[1], markersize=10, linewidth=2.0, linestyle="-")
ax3.plot(Ns, t_total_adaptive, color = time_category_colors[2], 
         marker=method_markers[2], markersize=8, linewidth=2.0, linestyle="-")
ax3.set_xlabel("$N$")
# ax3.set_ylabel("time [ms]")
ax3.set_title("$t_\mathrm{total}$")

ylim1 = ax1.get_ylim()
ylim2 = ax2.get_ylim()
ylim3 = ax3.get_ylim()

ax1.set_ylim((min(ylim1[0], ylim2[0], ylim3[0]), 
              max(ylim1[1], ylim2[1], ylim3[1])))
ax2.set_ylim((min(ylim1[0], ylim2[0], ylim3[0]), 
              max(ylim1[1], ylim2[1], ylim3[1])))
ax3.set_ylim((min(ylim1[0], ylim2[0], ylim3[0]), 
              max(ylim1[1], ylim2[1], ylim3[1])))

# make legend
if len(time_categories) == len(methods):
    legend_labels = []
    for i in range(len(time_categories)):
        legend_labels.append(time_categories[i])
        legend_labels.append(methods[i])
    ncols = 3
else:
    legend_labels = time_categories + methods
    ncols = 2
legend_handles = []
for label in legend_labels:
    if label in time_categories:
        handle = plt.Line2D([], [], 
                    color=time_category_colors[time_categories.index(label)],
                    linestyle=time_category_linestyles[time_categories.index(label)])
    else:
        handle = plt.Line2D([], [], color="k",
                    marker=method_markers[methods.index(label)],
                    markersize=10, markerfacecolor='white')
    legend_handles.append(handle)

legend_ax = fig.add_subplot(gs[1, 0])
legend_ax.axis('off')

legend_handles_short = [legend_handles[i] for i in range(1, len(legend_handles), 2)]
legend_labels_short = [legend_labels[i] for i in range(1, len(legend_labels), 2)]

legend = ax1.legend(legend_handles_short, legend_labels_short, bbox_to_anchor=(1.75, -0.2), ncol=ncols, loc="upper center")

plt.subplots_adjust(wspace=0.3)

plt.savefig("figures/scaling.pdf")