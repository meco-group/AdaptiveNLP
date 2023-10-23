import matplotlib.pyplot as plt
import csv
import numpy as np

folder = "plotting_data"

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

tt = []
uu = []

with open(folder+"/controls.csv") as file:
    csv_reader = csv.reader(file)
    header = next(csv_reader)
    nu = len(header)-1

    uu = [[] for i in range(nu)]
        
    for row in csv_reader:
        tt.append(float(row[0]))
        for i in range(nu):
            uu[i].append(float(row[i+1]))
    

# plt.subplots()
fig = plt.figure(figsize=(6.4, 2.8))
for i in range(nu):
    plt.plot(tt, uu[i], marker='o', fillstyle='none', color='royalblue')
plt.xlabel("time [s]")
plt.ylabel("control input")
plt.yticks([0, 4, 8])

plt.legend([plt.Line2D([], [], color='royalblue', marker='o', 
                       linestyle='none', fillstyle='none'),
            plt.Line2D([], [], color='royalblue')],
           ["$u_k$", "interpolated values"])

# plt.show()
plt.subplots_adjust(bottom=0.3)
plt.savefig("figures/controls.pdf")