import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd
import glob
import dask.dataframe as dd
from sklearn.metrics.pairwise import cosine_similarity
import math
from decimal import Decimal

if __name__ == "__main__":
    work_dir = "/home/pedro/singularity/singularity-ce-3.8.1/workspace/LBDEMcoupling-public/examples/benchmark_rheological_yann/outDir_dalco48_6/"
    log_file = open(work_dir + "benchmark_rheological_yann_outDir.log", "r")
    log_data = log_file.read()
    log_data_split_1 = log_data.split('\n')
    string_check = 'Step    Atoms         KinEng              1          ts[1]          ts[2]            CPU'
    index = 0
    bulk_data = []
    bulk_data_timestep = []
    bulk_data_kineng = []
    for log_line in log_data_split_1:
        if string_check in log_line:
            data_1 = log_data_split_1[index + 1].split(' ')
            aux_data = []
            for sub_string in data_1:
                try:
                    aux_data.append( float(Decimal(sub_string)) )
                except Exception as e:
                    pass
            
            bulk_data_timestep.append(aux_data[0])
            bulk_data_kineng.append(aux_data[2])
            bulk_data.append(aux_data)
          
        index += 1            

    fig1 = plt.figure()
    fig1.suptitle('Average Kinect Energy', fontsize=20)
    ay = fig1.add_subplot(111)
    ay.plot(bulk_data_timestep[0:-6], bulk_data_kineng[0:-6], 'k')
    #ay.legend(loc='best', fontsize=10)
    axes = plt.gca()
    ay.set_xlabel('Time Steps', fontsize=20)
    ay.set_ylabel('Kinect Energy', fontsize=20)
    ay.set_ylim([0.0, 3.3652345e-10])
    plt.show()

    import pdb
    pdb.set_trace()