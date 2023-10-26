import matplotlib.pyplot as plt
import csv
import numpy as np

import pathlib
folder = str(pathlib.Path(__file__).parent.resolve())
figure_folder = folder + "/../figures"

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

try:
    latexify()
except:
    print("Unable to latexify the figure.")

t_solve_adaptive = []
t_err_adaptive = []
t_update_adaptive = []
t_tot_adaptive = []
t_solve_casadi = []
t_err_casadi = []
t_update_casadi = []
t_tot_casadi = []

with open(folder+"/adaptive_gridding_timings.csv") as file:
    csv_reader = csv.reader(file)
    header = next(csv_reader)
    for row in csv_reader:
        t_solve_adaptive.append(float(row[0]))
        t_err_adaptive.append(float(row[1]))
        t_update_adaptive.append(float(row[2]))
        t_tot_adaptive.append(t_solve_adaptive[-1] + t_err_adaptive[-1] + 
                              t_update_adaptive[-1])
        t_solve_casadi.append(float(row[3]))
        t_err_casadi.append(float(row[4]))
        t_update_casadi.append(float(row[5]))
        t_tot_casadi.append(t_solve_casadi[-1] + t_err_casadi[-1] + 
                              t_update_casadi[-1])

max_t = max(max(t_tot_adaptive), max(t_tot_casadi))
medians = [np.median(t_solve_casadi), np.median(t_err_casadi), 
           np.median(t_update_casadi), np.median(t_tot_casadi),
           np.median(t_solve_adaptive), np.median(t_err_adaptive), 
           np.median(t_update_adaptive), np.median(t_tot_adaptive)]

# complete figure
plt.figure()
positions = [1, 2, 3, 4, 5, 6, 7, 8]
bplot = plt.boxplot([t_solve_casadi, t_err_casadi, t_update_casadi, 
                     t_tot_casadi, t_solve_adaptive, t_err_adaptive, 
                     t_update_adaptive, t_tot_adaptive],
                     patch_artist=True, widths=0.5, positions=positions)
for patch, c in zip(bplot['boxes'], ['gray', 'gray', 'gray', 'gray', 
                                     'royalblue', 'royalblue', 'royalblue', 
                                     'royalblue']):
    patch.set_facecolor(c)
    patch.set_color(c)

for median in bplot['medians']:
    median.set_color('black')

plt.xticks(positions, labels=["$t_{solve}$", "$t_{err}$", "$t_{update}$",
                              "$t_{total}$", "$t_{solve}$", "$t_{err}$", 
                              "$t_{update}$", "$t_{total}$"])
plt.xlim(positions[0]-0.4, positions[-1]+1.0)
plt.ylabel("$time$ $[ms]$")
plt.ylim(-1.0, max_t+5)
text_shift = 0.28
for j in range(len(medians)):
    plt.text(positions[j]+text_shift, medians[j], round(medians[j], 2))

plt.savefig(figure_folder + "/adaptive_gridding_timings.pdf")
# plt.savefig(figure_folder + "/adaptive_gridding_timings.png", dpi=400)



# figure casadi
plt.figure()
positions = [1, 2, 3, 4]
bplot = plt.boxplot([t_solve_casadi, t_err_casadi, t_update_casadi, 
                     t_tot_casadi],
                     patch_artist=True, widths=0.5, positions=positions)
for patch, c in zip(bplot['boxes'], ['gray', 'gray', 'gray', 'gray']):
    patch.set_facecolor(c)
    patch.set_color(c)

for median in bplot['medians']:
    median.set_color('black')

plt.xticks(positions, labels=["$t_{solve}$", "$t_{err}$", 
                              "$t_{update}$", "$t_{total}$"])
plt.xlim(positions[0]-0.4, positions[-1]+0.7)
plt.ylabel("$time$ $[ms]$")
plt.ylim(-1.0, max_t+5)
text_shift = 0.28
for j in range(4):
    plt.text(positions[j]+text_shift, medians[j], round(medians[j], 2))

plt.savefig(figure_folder + "/adaptive_gridding_timings_0.pdf")
# plt.savefig(figure_folder + "/adaptive_gridding_timings_0.png", dpi=400)

# figure adaptive
plt.figure()
positions = [1, 2, 3, 4]
bplot = plt.boxplot([t_solve_adaptive, t_err_adaptive, t_update_adaptive, 
                     t_tot_adaptive],
                     patch_artist=True, widths=0.5, positions=positions)
for patch, c in zip(bplot['boxes'], ['royalblue', 'royalblue', 'royalblue', 
                                     'royalblue']):
    patch.set_facecolor(c)
    patch.set_color(c)

for median in bplot['medians']:
    median.set_color('black')

plt.xticks(positions, labels=["$t_{solve}$", "$t_{err}$", 
                              "$t_{update}$", "$t_{total}$"])
plt.xlim(positions[0]-0.4, positions[-1]+0.7)
plt.ylabel("$time$ $[ms]$")
plt.ylim(-1.0, max_t+5)
text_shift = 0.28
for j in range(4):
    plt.text(positions[j]+text_shift, medians[4+j], round(medians[4+j], 2))

# plt.savefig(figure_folder + "/adaptive_gridding_timings_1.png", dpi=400)
plt.savefig(figure_folder + "/adaptive_gridding_timings_1.pdf")


from matplotlib.gridspec import GridSpec
width = 0.4
gs = GridSpec(2, 1, height_ratios=[0.9, 0.1], hspace=0.2)
fig = plt.figure(figsize=(7,3.8))
axs = fig.add_subplot(gs[0,0])
# axs = plt.gca()
bottoms = [0,0]
p = plt.bar(["CasADi Opti", "AdaptiveNLP"],
            [np.median(t_solve_casadi), np.median(t_solve_adaptive)],
            label="$t_{solve}$",
            bottom=bottoms,
            width=width,
            color = "royalblue")
plt.bar_label(p, label_type='center', fmt='%.2f', color="white")
bottoms[0] += np.median(t_solve_casadi)
bottoms[1] += np.median(t_solve_adaptive)

p = plt.bar(["CasADi Opti", "AdaptiveNLP"],
            [np.median(t_err_casadi), np.median(t_err_adaptive)],
            label="$t_{err}$",
            bottom=bottoms,
            width=width,
            color = "orange")
# plt.bar_label(p, label_type='center', fmt='%.2f', padding=9)
plt.text(p.get_children()[0].get_center()[0] + width/2 + 0.1, 
         bottoms[0] + 0.5*np.median(t_err_casadi) - 0.15,
         f'{np.median(t_err_casadi):.2f}',
         ha='center')
plt.text(p.get_children()[1].get_center()[0] + width/2 + 0.1, 
         bottoms[1] + 0.5*np.median(t_err_adaptive) - 0.15,
         f'{np.median(t_err_adaptive):.2f}',
         ha='center')
bottoms[0] += np.median(t_err_casadi)
bottoms[1] += np.median(t_err_adaptive)

p = plt.bar(["CasADi Opti", "AdaptiveNLP"],
            [np.median(t_update_casadi), np.median(t_update_adaptive)],
            label="$t_{update}$",
            bottom=bottoms,
            width=width,
            color = "gray")
plt.text(p.get_children()[0].get_center()[0], 
         bottoms[0] + 0.5*np.median(t_update_casadi),
         f'{np.median(t_update_casadi):.2f}',
         ha='center',
         color="white")
plt.text(p.get_children()[1].get_center()[0], 
         bottoms[1] + 0.5*np.median(t_update_adaptive) + 0.22,
         f'{np.median(t_update_adaptive):.2f}',
         ha='center')


dist = p.get_children()[1].get_center()[0] - \
    p.get_children()[0].get_center()[0]
plt.xlim([p.get_children()[0].get_center()[0] - dist/2, 
          p.get_children()[1].get_center()[0] + dist/2])


plt.ylabel("time [ms]")
plt.legend()

legend_ax = fig.add_subplot(gs[1, 0])
legend_ax.axis('off')  # Turn off the axes for the legend subplot
legend = axs.legend(bbox_to_anchor=(0.5, -0.12), loc='upper center', ncol=3)

plt.savefig(figure_folder + "/adpative_gridding_timings_bars.pdf")

plt.close()
