import matplotlib.pyplot as plt
import csv

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
xx = []
uu = []

file_nb = 0
file_found = True

nx = 4
nu = 2

while (file_found):
    try:
        with open(folder+"\solution_" + str(file_nb) + ".csv") as file:
            csv_reader = csv.reader(file)
            header = next(csv_reader)

            tt_to_add = []
            xx_to_add = [[] for i in range(nx)]
            uu_to_add = [[] for i in range(nu)]
                
            for row in csv_reader:
                tt_to_add.append(float(row[0]))
                for i in range(nx):
                    xx_to_add[i].append(float(row[1+i]))
                for i in range(nu):
                    uu_to_add[i].append(float(row[1+nx+i]))
            
            tt.append(tt_to_add)
            xx.append(xx_to_add)
            uu.append([uu_to_add[i][:-1] for i in range(nu)])
            
        file_nb += 1
    except Exception as err:
        # print(err)
        file_found = False

for i in range(file_nb):
    fig = plt.figure()

    if (i > 0):
        plt.gca().add_patch(plt.Circle([1.5, 0.4], 0.3, color="red"))
    plt.plot(xx[i][0], xx[i][1], marker='o')
    plt.gca().set_aspect("equal")
    plt.savefig(figure_folder + "/trajectory_" + str(i) + ".pdf")


    fig = plt.figure()
    for j in range(nu):
        plt.plot(tt[i][:-1], uu[i][j], marker='o')
    plt.xlabel("time [s]")
    plt.ylabel("control input")
    # plt.yticks([0, 4, 8])

    plt.legend(["acceleration", "steering angle"])

    plt.savefig(figure_folder + "/controls_" + str(i) + ".pdf")

    plt.close()