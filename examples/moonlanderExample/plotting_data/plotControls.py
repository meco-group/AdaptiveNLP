import matplotlib.pyplot as plt
import csv
import numpy as np

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

tt = []
uu = []

file_nb = 0
file_found = True

while (file_found):
    try:
        with open(folder+"/controls_" + str(file_nb) + ".csv") as file:
            csv_reader = csv.reader(file)
            header = next(csv_reader)
            nu = len(header)-1

            uu_to_add = [[] for i in range(nu)]
            tt_to_add = []
                
            for row in csv_reader:
                tt_to_add.append(float(row[0]))
                for i in range(nu):
                    uu_to_add[i].append(float(row[i+1]))

            uu.append(uu_to_add)
            tt.append(tt_to_add)

        file_nb += 1
    except Exception as err:
        # print(err)
        file_found = False

    
for i in range(file_nb):
    fig = plt.figure(figsize=(6.4, 2.8))
    for j in range(nu):
        plt.plot(tt[i], uu[i][j], marker='o', fillstyle='none', color='royalblue')
    plt.xlabel("time [s]")
    plt.ylabel("control input")
    plt.yticks([0, 4, 8])

    plt.legend([plt.Line2D([], [], color='royalblue', marker='o', 
                        linestyle='none', fillstyle='none'),
                plt.Line2D([], [], color='royalblue')],
            ["$u_k$", "interpolated values"])

    # plt.show()
    plt.subplots_adjust(bottom=0.3)
    plt.savefig(figure_folder + "/animation/controls_" + str(i) + ".png", dpi=400)

    if (i == file_nb-1):
        plt.savefig(figure_folder + "/controls.pdf", dpi=400)

    plt.close()