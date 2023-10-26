import matplotlib.pyplot as plt
import csv
import numpy as np
from matplotlib.patches import Rectangle
from matplotlib.ticker import MultipleLocator

import pathlib
folder = str(pathlib.Path(__file__).parent.resolve())
figure_folder = folder + "/../figures"

def latexify():
    fontsize = 30
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

# read jacobian sparsities
jac_rows_adaptive = [[]]
jac_cols_adaptive = [[]]
jac_rows_casadi = [[]]
jac_cols_casadi = [[]]

iteration_ptr = 0

with open(folder+"/adaptive_gridding_sparsities_jac.csv") as file:
    csv_reader = csv.reader(file)
    header = next(csv_reader)
    for row in csv_reader:
        if (float(row[0]) != -1):
            jac_rows_adaptive[iteration_ptr].append(int(row[0]))
            jac_cols_adaptive[iteration_ptr].append(int(row[1]))
            jac_rows_casadi[iteration_ptr].append(int(row[2]))
            jac_cols_casadi[iteration_ptr].append(int(row[3]))
        else:
            iteration_ptr += 1
            jac_rows_adaptive.append([])
            jac_cols_adaptive.append([])
            jac_rows_casadi.append([])
            jac_cols_casadi.append([])

jac_rows_adaptive = jac_rows_adaptive[:-1]
jac_cols_adaptive = jac_cols_adaptive[:-1]
jac_rows_casadi = jac_rows_casadi[:-1]
jac_cols_casadi = jac_cols_casadi[:-1]

jac_max_rows_adaptive = max([max(rr) for rr in jac_rows_adaptive])
jac_max_cols_adaptive = max([max(rr) for rr in jac_cols_adaptive])
jac_max_rows_casadi = max([max(rr) for rr in jac_rows_casadi])
jac_max_cols_casadi = max([max(rr) for rr in jac_cols_casadi])


# read hessian sparsities
hess_rows_adaptive = [[]]
hess_cols_adaptive = [[]]
hess_rows_casadi = [[]]
hess_cols_casadi = [[]]

iteration_ptr = 0

with open(folder+"/adaptive_gridding_sparsities_hess.csv") as file:
    csv_reader = csv.reader(file)
    header = next(csv_reader)
    for row in csv_reader:
        if (float(row[0]) != -1):
            hess_rows_adaptive[iteration_ptr].append(int(row[0]))
            hess_cols_adaptive[iteration_ptr].append(int(row[1]))
            hess_rows_casadi[iteration_ptr].append(int(row[2]))
            hess_cols_casadi[iteration_ptr].append(int(row[3]))
        else:
            iteration_ptr += 1
            hess_rows_adaptive.append([])
            hess_cols_adaptive.append([])
            hess_rows_casadi.append([])
            hess_cols_casadi.append([])

hess_rows_adaptive = hess_rows_adaptive[:-1]
hess_cols_adaptive = hess_cols_adaptive[:-1]
hess_rows_casadi = hess_rows_casadi[:-1]
hess_cols_casadi = hess_cols_casadi[:-1]

hess_max_rows_adaptive = max([max(rr) for rr in hess_rows_adaptive])
hess_max_cols_adaptive = max([max(rr) for rr in hess_cols_adaptive])
hess_max_rows_casadi = max([max(rr) for rr in hess_rows_casadi])
hess_max_cols_casadi = max([max(rr) for rr in hess_cols_casadi])


def plotSparsity(rows, cols, max_row, max_col, hessian=False):
    f = plt.figure()

    max_col_local = max(cols)+1
    max_row_local = max(rows)+1
    if (hessian):
        max_col = max(max_row, max_col)
        max_row = max(max_row, max_col)
        max_col_local = max(max_row_local, max_col_local)
        max_row_local = max(max_row_local, max_col_local)

    m = np.zeros((max_row+1, max_col+1))
    for r, c in zip(rows, cols):
        m[r,c] = 1

    if (hessian):
        mt = np.transpose(m)
        plt.spy(mt, markersize=2.5, color='gray')

    plt.spy(m)#, markersize=10, marker='s', color='black')
    ax = plt.gca()
    ax.add_patch(Rectangle((max_col_local-0.5, -0.5), 
                           max_col-max_col_local+1, max_row+1,
                           edgecolor='lightgray',
                           facecolor='lightgray'))
    ax.add_patch(Rectangle((-0.5, max_row_local-0.5), 
                           max_col+1, max_row-max_row_local+1,
                           edgecolor='lightgray',
                           facecolor='lightgray'))
    min_tick_difference = 20
    ax.xaxis.set_major_locator(MultipleLocator(base=min_tick_difference))
    
    return f


for iteration in range(len(jac_rows_adaptive)):
# for iteration in range(1):
    plotSparsity(jac_rows_adaptive[iteration], jac_cols_adaptive[iteration], 
                 jac_max_rows_adaptive, jac_max_cols_adaptive, False)
    
    plt.subplots_adjust(top=0.8)
    plt.title(f"iteration {iteration}")
    plt.savefig(figure_folder + "/sparsities/jac_adaptive_"+str(iteration)+".png",
                dpi=400)
    # plt.savefig(figure_folder + "/sparsities/jac_adaptive_"+str(iteration)+".pdf")
    plt.close()
    
    plotSparsity(hess_rows_adaptive[iteration], hess_cols_adaptive[iteration], 
                 max(hess_max_rows_adaptive, hess_max_cols_adaptive),
                 max(hess_max_rows_adaptive, hess_max_cols_adaptive), True)
    
    plt.subplots_adjust(top=0.8)
    plt.title(f"iteration {iteration}")
    plt.savefig(figure_folder + "/sparsities/hess_adaptive_"+str(iteration)+".png",
                dpi=400)
    # plt.savefig(figure_folder + "/sparsities/hess_adaptive_"+str(iteration)+".pdf")
    
    plotSparsity(jac_rows_casadi[iteration], jac_cols_casadi[iteration], 
                 jac_max_rows_casadi, jac_max_cols_casadi, False)
    
    plt.subplots_adjust(top=0.8)
    plt.title(f"iteration {iteration}")
    plt.savefig(figure_folder + "/sparsities/jac_casadi_"+str(iteration)+".png",
                dpi=400)
    # plt.savefig(figure_folder + "/sparsities/jac_casadi_"+str(iteration)+".pdf")
    plt.close()
    
    plotSparsity(hess_rows_casadi[iteration], hess_cols_casadi[iteration], 
                 max(hess_max_rows_casadi, hess_max_cols_casadi),
                 max(hess_max_rows_casadi, hess_max_cols_casadi), True)
    
    plt.subplots_adjust(top=0.8)
    plt.title(f"iteration {iteration}")
    plt.savefig(figure_folder + "/sparsities/hess_casadi_"+str(iteration)+".png",
                dpi=400)
    # plt.savefig(figure_folder + "/sparsities/hess_casadi_"+str(iteration)+".pdf")
    plt.close()

