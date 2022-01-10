import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd
import glob
import dask.dataframe as dd
from sklearn.metrics.pairwise import cosine_similarity
import math

if __name__ == "__main__":
    work_dir = '/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/examples/benchmark_rheological_yann_init/'
    
    df_lattice_average_energy = dd.read_csv(work_dir + 'outDir_5_dalco48/outDir/tmp/lattice_average_energy.csv', sep=',').compute().dropna()
    df_physical_average_velocity_x = dd.read_csv(work_dir + 'outDir_5_dalco48/outDir/tmp/lattice_average_velocity_x.csv', sep=',').compute().dropna()
    df_physical_average_velocity_gradient_x = dd.read_csv(work_dir + 'outDir_5_dalco48/outDir/tmp/lattice_average_velocity_gradient_x.csv', sep=',').compute().dropna()

    fig1 = plt.figure()
    fig1.suptitle('Lattice Average Energy Through the Iterations', fontsize=20)
    ay = fig1.add_subplot(111)
    ay.plot(df_lattice_average_energy['iT'], df_lattice_average_energy['average_energy'], 'k', label='getStoredAverageEnergy<T>(lattice)')
    ay.legend(loc='best', fontsize=20)
    axes = plt.gca()
    ay.set_xlabel('Lattice Iteration', fontsize=20)
    ay.set_ylabel('Lattice Average Energy', fontsize=20)
    #axes.set_xlim([min([min(df_polar_directly.z.values), min(df_polar_surface.z.values)]), max([max(df_polar_directly.z.values), max(df_polar_surface.z.values)])])
    #axes.set_ylim([min([min(df_polar_directly.fric.values), min(df_polar_surface.fric.values)]), max([max(df_polar_directly.fric.values), max(df_polar_surface.fric.values)])])

    df_physical_average_velocity_x_aux = (df_physical_average_velocity_x.iloc[len(df_physical_average_velocity_x)-9:len(df_physical_average_velocity_x)])[:].mean()[1:]
    fig1 = plt.figure()
    fig1.suptitle('Average Velocity X Through the Vertical Direction', fontsize=20)
    ay = fig1.add_subplot(111)
    ay.plot(df_physical_average_velocity_x_aux.to_numpy(), np.array(df_physical_average_velocity_x_aux.index).astype(np.int64), 'b')
    #ay.legend(loc='best', fontsize=10)
    axes = plt.gca()
    ay.set_xlabel('Physical Velocity (m/s)', fontsize=20)
    ay.set_ylabel('Lattice Position on the Vertical Direction', fontsize=20)

    #df_physical_average_velocity_gradient_x_aux = (df_physical_average_velocity_gradient_x.iloc[len(df_physical_average_velocity_gradient_x)-9:len(df_physical_average_velocity_gradient_x)])[:].mean()[1:]
    #df_physical_average_velocity_gradient_x_aux = (df_physical_average_velocity_gradient_x.iloc[0:len(df_physical_average_velocity_gradient_x)])[:].mean()[1:]
    df_physical_average_velocity_gradient_x_aux = df_physical_average_velocity_x_aux[::-1].diff()[1:][::-1]
    fig1 = plt.figure()
    fig1.suptitle('Average Velocity Gradient X Through the Vertical Direction', fontsize=20)
    ay = fig1.add_subplot(111)
    ay.plot(df_physical_average_velocity_gradient_x_aux.to_numpy(), np.array(df_physical_average_velocity_gradient_x_aux.index).astype(np.int64), 'r')
    #ay.legend(loc='best', fontsize=10)
    axes = plt.gca()
    ay.set_xlabel('Physical Velocity Gradient (1/s)', fontsize=20)
    ay.set_ylabel('Lattice Position on the Vertical Direction', fontsize=20)

    measured_strain_rate = df_physical_average_velocity_gradient_x_aux.to_numpy()[0]
    strain_rate = 100.0
    #strain_rate = 0.05
    relative_apparent_viscosity = measured_strain_rate / strain_rate
    print("relative_apparent_viscosity: " + str(relative_apparent_viscosity) )

    plt.show()
    import pdb
    pdb.set_trace()