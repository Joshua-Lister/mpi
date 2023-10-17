import matplotlib.pyplot as plt
from matplotlib import animation
import pandas as pd
from glob import glob
import numpy as np 
from mpl_toolkits.mplot3d import Axes3D
def collate(file1, file2):
    with open(file1, "r") as f:
        fline = f.readline().rstrip().split("\t")
    m = int(fline[0])
    n = int(fline[1])
    with open(file2, "r") as f:
        for i in range(2):
            f.readline()
        fline = f.readline().rstrip().split("\t")
        ffline = f.readline().rstrip().split("\t")
        bc = f.readline().rstrip()
    xmax = int(fline[0])
    ymax = int(fline[1])
    c = float(fline[2])
    dt_out = float(ffline[1])
    tmax = int(ffline[0])
    list_of_files = glob('out/*.dat')
    list_of_files = sorted(list_of_files)
    df_from_each_file = [pd.read_csv(f, sep='\s+', header=None, skiprows=2) for f in list_of_files]
    df_dframe = []
    for i in range(0, len(list_of_files), m * n):
        df = pd.concat([df_from_each_file[j] for j in range(i, n + i)], axis = 1, ignore_index = True)
        for j in range(1, m):
            df1 = pd.concat([df_from_each_file[i + k + j * m] for k in range(0, n)], axis = 1, ignore_index = True)
            df = pd.concat([df, df1], axis = 0, ignore_index = True)
        df_dframe.append(df)

    return df_dframe, xmax, ymax, c, tmax, dt_out, bc

def plotwave(u, k):
    plt.clf()
    plt.axis('off')
    plt.title(f'{bc} Solution at t = {((k + 1) * dt_out):3f} seconds')
    plt.imshow(u, cmap = 'coolwarm')

    return plt

def animate(i):
    plotwave(df_dframe[i], i)


df_dframe, xmax, ymax, c, tmax, dt_out, bc = collate("out/output_1000_id0.dat","./input.txt")
dx = xmax / (df_dframe[0].shape[0] - 1)
dy = ymax / (df_dframe[0].shape[1] - 1)


frames = int(tmax / dt_out)
anim = animation.FuncAnimation(plt.figure(), animate, interval=1, frames=frames, repeat = False)
anim.save('gifs/wave_parallel.gif')
plt.tight_layout()
plt.show()
