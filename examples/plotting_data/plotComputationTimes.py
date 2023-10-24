import matplotlib.pyplot as plt
import csv
import numpy as np

import pathlib
folder = str(pathlib.Path(__file__).parent.resolve()) # "examples/plotting_data"
figure_folder = folder + "/../figures"

def latexify():
    fontsize = 15
    params = {'backend': 'ps',
              'font.size': fontsize,
              'axes.labelsize': fontsize,
              'axes.titlesize': fontsize,
              'legend.fontsize': fontsize-5,
              'xtick.labelsize': fontsize,
              'ytick.labelsize': fontsize,
              'text.usetex': True,
              'font.family': 'serif',
              'figure.figsize': [7,5],
              'text.latex.preamble': [r'\usepackage{bm}'],
              }

    plt.rcParams.update(params)

latexify()

# read the data
nb_files = 3

t_solve = [[] for i in range(nb_files)]
t_update = [[] for i in range(nb_files)]
t_total = [[] for i in range(nb_files)]
medians = [[0 for i in range(3)] for i in range(nb_files)]
for i in range(nb_files):
    data = []
    with open(folder+"/timings_" + str(i) + ".csv") as file:
        csv_reader = csv.reader(file)
        header = next(csv_reader)
        for row in csv_reader:
            data.append([float(value) for value in row])
        
    t1, t2, t3 = zip(*data)
    t_solve[i] = list(t1[:-1])
    t_update[i] = list(t2[:-1])
    t_total[i] = list(t3[:-1])
    medians[i] = [t1[-1], t2[-1], t3[-1]]

# t_update = [[np.abs(t) for t in tt] for tt in t_update]
for i in range(nb_files):
    medians[i][1] = np.median(t_update[i])

# # compute max time so that all plots have the same axis limits
# max_t = np.max(t_total)

# # make the plots
# labels = [["$t_{solve}$", "$t_{update}$", "$t_{total}$"],
#           ["$t_{solve}$", "$t_{update}$", "$t_{total}$"],
#           ["$t_{solve}$", "$t_{update}$", "$t_{total}$"]]
# labels_combined = []
# for i in range(len(labels[0])):
#     labels_combined += [labels[0][i], labels[1][i], labels[2][i]]

# plt.rc('text', usetex=True)
# colors = [['gray', 'gray', 'gray'],
#           ['orange', 'orange', 'orange'],
#           ['royalblue', 'royalblue', 'royalblue']]
# colors_combined = []
# for i in range(len(colors[0])):
#     colors_combined += [colors[0][i], colors[1][i], colors[2][i]]

# # make separate boxplots
# for i in range(nb_files):
#     plt.figure()
#     bplot = plt.boxplot([t_solve[i], t_update[i], t_total[i]], widths = 0.6,
#                         patch_artist=True, boxprops=dict(color=colors[i][0]))
    
#     for patch, c in zip(bplot['boxes'], colors[i]):
#         patch.set_facecolor(c)

#     for median in bplot['medians']:
#         median.set_color('black')
    
#     plt.xticks([j+1 for j in range(3)], labels=labels[i])
#     plt.ylabel("$time$ $[ms]$")
#     plt.xlim(0, 4)
#     plt.ylim(-1.0, max_t+5)

#     text_shift = 0.4
#     for j in range(3):
#         plt.text(j+1+text_shift, medians[i][j], round(medians[i][j], 2))

#     # plt.savefig("figures/computation_times_" + str(i) + ".png", dpi=200)
#     plt.savefig("figures/computation_times_" + str(i) + ".pdf")
#     plt.close()


# # make combined figure
# plt.figure(figsize=[10.4, 4.8])
# my_positions = [1.05, 1.5, 1.95, 3.05, 3.5, 3.95, 5.05, 5.5, 5.95]
# bplot = plt.boxplot([t_solve[0], t_solve[1], t_solve[2], 
#                      t_update[0], t_update[1], t_update[2], 
#                      t_total[0], t_total[1], t_total[2]], widths=0.4,
#                      patch_artist=True, 
#                      positions=my_positions,
#                      boxprops=dict(color=colors_combined[0]))
# for patch, c in zip(bplot['boxes'], colors_combined):
#     patch.set_facecolor(c)
#     patch.set_color(c)

# for median in bplot['medians']:
#     median.set_color('black')

# # plt.xticks(my_positions, labels=labels_combined, fontsize=9)
# plt.xticks([my_positions[1], my_positions[4], my_positions[7]], 
#            labels=["$t_{solve}$", "$t_{update}$", "$t_{total}$"])
# plt.ylabel("$time$ $[ms]$")
# plt.xlim(my_positions[0]-0.3, my_positions[-1]+0.3)
# plt.ylim(-1.0, max_t+5)
# plt.legend([bplot['boxes'][0], bplot['boxes'][1], bplot['boxes'][2]],
#            ["CasADi::Opti 1", "CasADi::Opti 2", "AdaptiveNLP"])

# # text_shift = 0.15
# # for i in range(len(medians)):
# #     for j in range(3):
# #         plt.text(my_positions[3*j+i] + text_shift, medians[i][j], 
# #                  round(medians[i][j], 2),
# #                  fontsize=5)

# # plt.savefig("figures/computation_times_combined.png", dpi=400)
# plt.savefig("figures/computation_times_combined.pdf")


# plt.figure()
# my_positions = [1, 2, 3]
# bplot = plt.boxplot([t_total[0], t_total[1], t_total[2]], widths=0.5,
#                      patch_artist=True, 
#                      positions=my_positions)
# for patch, c in zip(bplot['boxes'], [colors[i][0] for i in range(3)]):
#     patch.set_facecolor(c)
#     patch.set_color(c)

# for median in bplot['medians']:
#     median.set_color('black')

# # plt.xticks(my_positions, labels=labels_combined, fontsize=9)
# plt.xticks(my_positions, 
#            labels=["CasADi::Opti 1", "CasADi::Opti 2", "AdaptiveNLP"])
# plt.ylabel("$time$ $[ms]$")
# plt.ylim(-1.0, max_t+5)
# text_shift = 0.28
# for j in range(3):
#     plt.text(my_positions[j]+text_shift, medians[j][2], round(medians[j][2], 2))

# # plt.savefig("figures/computation_times_combined-zoomed.png", dpi=400)
# plt.savefig("figures/computation_times_combined-zoomed.pdf")

# plt.close()


# Bar plot
plt.figure(figsize=(7,3.5))

medians_casadi_1 = medians[0]
medians_casadi_2 = medians[1]
medians_adaptiveNLP = medians[2]

width = 0.4

p = plt.bar(["CasADi Opti 1", "CasADi Opti 2", "AdaptiveNLP"],
        [medians_casadi_1[0], medians_casadi_2[0], medians_adaptiveNLP[0]],
        label="$t_{\mathrm{solve}}$", 
        bottom=[0,0,0],
        color="royalblue",
        width=width)
plt.bar_label(p, label_type='center', fmt='%.2f', color="white")

p = plt.bar(["CasADi Opti 1", "CasADi Opti 2", "AdaptiveNLP"],
        [medians_casadi_1[1], medians_casadi_2[1], medians_adaptiveNLP[1]],
        label="$t_{\mathrm{update}}$", 
        bottom=[medians_casadi_1[0], medians_casadi_2[0], medians_adaptiveNLP[0]],
        color="gray",
        width=width)
# plt.bar_label(p, label_type='center', fmt='%.2f', padding=7)
plt.text(p.get_children()[0].get_center()[0], 
         medians_casadi_1[0] + 0.5*np.median(medians_casadi_1[1] + 0.4),
         f'{np.median(medians_casadi_1[1]):.2f}',
         ha='center',
         color="black")
plt.text(p.get_children()[1].get_center()[0], 
         medians_casadi_2[0] + 0.5*np.median(medians_casadi_2[1]),
         f'{np.median(medians_casadi_2[1]):.2f}',
         ha='center',
         color="white")
plt.text(p.get_children()[2].get_center()[0], 
         medians_adaptiveNLP[0] + 0.5*np.median(medians_adaptiveNLP[1] + 0.6),
         f'{np.median(medians_adaptiveNLP[1]):.2f}',
         ha='center',
         color="black")

plt.ylabel("time [ms]")

plt.legend()

plt.savefig(figure_folder + "/computation_times_bars.pdf")