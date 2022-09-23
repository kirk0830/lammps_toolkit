# q(z,t) version 2022.09.23
# at first i thought it was very easy to implement this time correlation function, however,
# large memory and CPU cost beat me. On some computers it is possible to run the very first
# version of q(z,t). Then program was modified to analyze data slice by slice. However this
# time the time-complexity becomes another fatal problem
# in this version, large amount of temporary files will be generated, but will be deleted at
# the end of run of program

import numpy as np

def writematrix2d21d(filename, matrix):
    '''
    Write Matrix in 2d to 1d to file\n
    Elements are seperated by whitespace, can be removed simply by split().
    '''
    nrow, ncol = np.shape(matrix)
    with open(filename, mode = 'w+', encoding = 'utf-8') as f:
        for irow in range(nrow):
            for icol in range(ncol):
                f.write(str(matrix[irow][icol]))
                if irow*icol != (nrow - 1)*(ncol - 1):
                    f.write(' ')
                    
def readmatrix1d22d(filename, shape, dtype = float):
    '''
    Read Matrix in 1d from file and recover its shape\n
    shape: tuple, say (8, 8) or etc., same as what used in numpy, or give -1 to have 1d list\n
    dtype: data type, to convert every element to certain data type using map function\n
    '''
    with open(filename, mode = 'r', encoding = 'utf-8') as f:
        B = list(map(dtype, f.readline().split()))
    return np.reshape(B, shape)

from lmp_readin import lmp2list
from coords_io import readcoords
from doi101038nphys2133 import ndigitsgen_2d
from os import remove
def lmptrj2ndigitsfiles(filename, elements, nx, ny, nz, dist):
    
    '''
    lmptrj2ndigitsfiles\n
    convert LAMMPS Trajectory file to multiple n-digits list files, marked as *_[t]_[z].dat\n
    filename: (str) name of LAMMPS Trajectory file\n
    elements: (str, list) elements in list, MUST have the same number of elements as LAMMPS Trajectory\n
    nx, ny, nz: (int) grid size of n-digits list, i.e., how many grid points are used to descretize present space\n
    dist: (float) distance interested in to calculate q(z,t) the distance dependent TCF
    '''
    lmp2list(filename = filename,
             elements = elements,
             iomode = 'files')
    idx_file = 0
    while True:
        try:
            icoords = readcoords(filename='{}_{}-pos-out.xyz'.format(filename, idx_file), convert=True)
            
            indigts2d = ndigitsgen_2d(icoords, nx, ny, nz, 'x', dist)
            writematrix2d21d('{}_{}_{}.dat'.format(filename, idx_file, dist), indigts2d)
            remove('{}_{}-pos-out.xyz'.format(filename, idx_file))
            idx_file += 1
        except FileNotFoundError:
            print('Try-Except report: FileNotFoundError raised, BREAK present loop')
            break

#lmptrj2ndigitsfiles('Pd_1000K.trj', ['Pd', 'Pd'], 100, 50, 50, 0.5)

def qzt(filename, z, max_idx_t, window_idx_t):
    
    A = 0
    B = 0
    for start_idx_t in range(max_idx_t - window_idx_t + 1):
        fname1 = '{}_{}_{}.dat'.format(filename, start_idx_t, z)
        n_vec_1 = readmatrix1d22d(fname1, -1)
        fname2 = '{}_{}_{}.dat'.format(filename, start_idx_t+window_idx_t, z)
        n_vec_2 = readmatrix1d22d(fname2, -1)
        A += np.dot(n_vec_1, n_vec_2)
        B += np.sum(n_vec_1)
    return A/B

def qz(filename, elements, nx, ny, nz, z, max_idx_t):
    
    lmptrj2ndigitsfiles(filename, elements, nx, ny, nz, z)
    qz_val = []
    for i_width_window in range(max_idx_t):
        A = 0
        B = 0
        for start_idx_t in range(max_idx_t - i_width_window + 1):
            fname1 = '{}_{}_{}.dat'.format(filename, start_idx_t, z)
            n_vec_1 = readmatrix1d22d(fname1, -1)
            fname2 = '{}_{}_{}.dat'.format(filename, start_idx_t+i_width_window, z)
            n_vec_2 = readmatrix1d22d(fname2, -1)
            A += np.dot(n_vec_1, n_vec_2)
            B += np.sum(n_vec_1)
        qzt_val = A/B
        qz_val.append(qzt_val)
    for idx_file in range(max_idx_t):
        remove('{}_{}_{}.dat'.format(filename, idx_file, z))
    return qz_val

qz_val = qz(filename = 'Pd_1000K.trj', elements = ['Pd', 'Pd'],
            nx = 100, ny = 50, nz = 50,
            z = 0.5, max_idx_t = 10)
from matplotlib import pyplot as plt
plt.plot(qz_val)
plt.show()