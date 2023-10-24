import matplotlib.pyplot as plt
import csv
import numpy as np

import pathlib
folder = str(pathlib.Path(__file__).parent.resolve()) # "examples/plotting_data"
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

# read the data
nb_obs = []
adaptive = []
casadi_1 = []
casadi_2 = []
t_comp_2 = []

data = []
with open(folder+"/virtualObstaclesScaling.csv") as file:
    csv_reader = csv.reader(file)
    header = next(csv_reader)
    for row in csv_reader:
        data.append([float(value) for value in row])
    
nb_obs, adaptive, casadi_1, casadi_2, t_comp_2= zip(*data)

plt.figure()
plt.xlabel("number of virtual obstacles")
plt.ylabel("computation time [ms]")
start_idx = 0
alpha = 1.0
for i in range(len(nb_obs)):
    if i == len(nb_obs)-1 or nb_obs[i+1] < nb_obs[i]:
        stop_idx = i+1
        # plt.plot(nb_obs[start_idx:stop_idx], adaptive[start_idx:stop_idx],
        #          linewidth=2, color="royalblue", alpha=alpha)
        plt.plot(nb_obs[start_idx:stop_idx], casadi_1[start_idx:stop_idx],
                 linewidth=2, color="gray", alpha=alpha)
        plt.plot(nb_obs[start_idx:stop_idx], casadi_2[start_idx:stop_idx],
                 linewidth=2, color="orange", alpha=alpha)
        plt.plot(nb_obs[start_idx:stop_idx], t_comp_2[start_idx:stop_idx],
                 linewidth=1, color="orange", alpha=alpha)
        start_idx=stop_idx
    


plt.savefig(figure_folder + "/virtual_obstacles_scaling.png", dpi=400)

plt.close()
