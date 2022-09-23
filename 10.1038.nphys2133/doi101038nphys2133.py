# this program is compatible with NPT ensemble.
# during simulation, the fractional coordinates are invariant with rescaling of cell dimension
# therefore program will directly make the [0,1]^3 space to grid
# but cell dimension is better to be provided to rationally make grid size
from lmp_readin import lmp2list
import numpy as np
import sys
VERBO = 3

def ndigitsgen(frame, nx, ny, nz):
    
    '''
    Generate a list in 3 dimension, from fractional xyz data, respect to certain level of griding\n
    For example dividing [0,1] along x-axis in 1001 slices, say, 0, 0.001, 0.002, 0.003, ..., 0.999, 1.000\n
    frame (2d list): in format of [['element1', x1, y1, z1], ['element2', x2, y2, z2], ...]\n
    nx, ny, nz (int): grid size
    '''
    dx = 1./nx # e.g., nx = 1000, dx = 0.001
    dy = 1./ny
    dz = 1./nz
    natom = len(frame)
    # xyz format: standard xyz format: element, x, y, z, without number of atoms and the second comment line
    n = np.zeros((nx, ny, nz), dtype=int)
    print('n-digit list generation information:\nGrid size: {}x{}x{}\nList information: {}'.format(
        nx, ny, nz, type(n)
    ))
    for iatom in frame:
        try:
            n[int(iatom[1]/dx)][int(iatom[2]/dy)][int(iatom[3]/dz)] += 1
        except IndexError:
            print('ERROR: index out of range:')
            print('i: {}\nj: {}\nk: {}'.format(int(iatom[1]/dx), int(iatom[2]/dy), int(iatom[3]/dz)))
            raise IndexError
        
    if np.max(n) > 1:
        
        print('ERROR: present precision of griding cannot generate required n-digits list, re-try with higher grid size!')
        raise ValueError
    
    # return n[nx+1, ny+1, nz+1]
    return n

def ndigitsgen_2d(frame, nx, ny, nz, axis, r):
    
    # due to limited memory, generate 2d slice instead of 3d
    
    n_2d = 0.
    dx = 1./nx # e.g., nx = 1000, dx = 0.001
    dy = 1./ny
    dz = 1./nz
    # subrountine configuration
    if VERBO > 3:
        print('q(z,t) | SUBROUTINE ndigitsgen_2d is called, n digit list generation info.:\n'
            +'         Frame has size: {}\n'.format(np.shape(frame))
            +'         Grid level on x-axis: {}\n'.format(nx)
            +'                    on y-axis: {}\n'.format(ny)
            +'                    on z-axis: {}\n'.format(nz)
            +'         Axis selected for slicing: {}\n'.format(axis)
            +'         Coodinate selected for slicing: {} = {}'.format(axis, r)
            )
    natom = len(frame)
    if (axis=='x') or (axis==0):
        __idx__ = int(r/dx)
        n_2d = np.zeros((ny, nz), dtype = int)
        for iatom in frame:
            if int(iatom[1]/dx) == __idx__:
                n_2d[int(iatom[2]/dy)][int(iatom[3]/dz)] += 1
    elif axis=='y' or (axis==1):
        __idx__ = int(r/dy)
        n_2d = np.zeros((nx, nz), dtype = int)
        for iatom in frame:
            if int(iatom[2]/dx) == __idx__:
                n_2d[int(iatom[1]/dx)][int(iatom[3]/dz)] += 1
    elif axis=='UNDEFINED':
        print('ERROR: axis is not properly defined!')
        raise ValueError
    else:
        __idx__ = int(r/dz)
        n_2d = np.zeros((nx, ny), dtype = int)
        for iatom in frame:
            if int(iatom[3]/dx) == __idx__:
                n_2d[int(iatom[1]/dx)][int(iatom[2]/dy)] += 1
    
    if np.max(n_2d) > 1:
        
        print('ERROR: present precision of griding cannot generate required n-digits list, re-try with higher grid size!')
        raise ValueError
    
    return n_2d
    
# from a frame of xyz return a frame of digits, n[nx+1, ny+1, nz+1]
def ndigits_t(trj, nx, ny, nz, mode = 'default', axis = 'UNDEFINED', r = 0):
    
    nframe = len(trj) # where trj is a Python dict
    ndigits_t_list = []
    for idx_frame in range(nframe):
        if mode == 'save_memory':
            ndigit_frame = ndigitsgen_2d(trj[idx_frame], nx, ny, nz, axis, r)
        else:
            ndigit_frame = ndigitsgen(trj[idx_frame], nx, ny, nz)
        ndigits_t_list.append(ndigit_frame)
        # if ndigit_frame:
        #     ndigits_t_list.append(ndigit_frame)
        # else:
        #     print('WARNING: increase nx, ny and nz to eliminate number larger than 1 in n-digits list!')
        #     raise TypeError
    return ndigits_t_list
# from a trajectory of xyz return a trajectory of digits, nlist[nt, nx, ny, nz]
def slicegen(list_txyz, axis, index, mode = 'default'):
    
    '''
    Slice generation\n
    list_txyz (4-d list): list_txyz[t, x, y, z]\n
    axis (str or int): 'x' or 0 will be interpreted as x-axis, 'y' or 1 corresponds to y-axis and otherwise z-axis by default\n
    index (int): the index along axis to slice 4d list
    '''
    if mode == 'save_memory':
        pass
    else:
        if (axis=='x') or (axis==0):
            # (t, y, z)
            #return np.transpose(list_txyz, (1,0,2,3))[index][:][:][:]
            return list_txyz[:,index,:,:]
        elif axis=='y' or (axis==1):
            # (t, x, z)
            #return np.transpose(list_txyz, (2,0,1,3))[index][:][:][:]
            return list_txyz[:,:,index,:]
        else:
            # (t, x, y)
            #return np.transpose(list_txyz, (3,0,1,2))[index][:][:][:]
            return list_txyz[:,:,:,index]
# slice one list of (t, x, y, z) to (t, x, y), (t, x, z) or (t, y, z)
# but it may be better only perform transposition for once instead of performing it every time

def qzt(list_t2d, idx_t_disp, mode = 'default'):
    
    # developer notes: here use indices to count frames instead of real time interval
    # list_t2d is a 3-dimension list
    len_t, len_dim1, len_dim2 = np.shape(list_t2d)
    qzt_val = 0.
    for idx_dim1 in range(len_dim1):
        for idx_dim2 in range(len_dim2):
            A = 0.
            B = 0.
            for idx_t in range(len_t - idx_t_disp + 1):
                # kernal loop
                A += list_t2d[idx_t][idx_dim1][idx_dim2]*list_t2d[idx_t+idx_t_disp][idx_dim1][idx_dim2]
                B += list_t2d[idx_t][idx_dim1][idx_dim2]
            qz_val += A/B
    return qzt_val

def qzt_2d(trj, nx, ny, nz, axis, r, idx_t_disp):

    if VERBO > 2:
        print('q(z,t)| SUBROUTINE qzt_2d is called, calculating qzt, iterating grids and time window...\n'
            +'         Present length of time window (param idx_t_disp): {}'.format(idx_t_disp))
    
    qzt_val = 0.
    len_t = len(trj)
    len_dim1 = 0
    len_dim2 = 0
    if (axis=='x') or (axis==0): len_dim1, len_dim2 = ny, nz
    elif (axis=='y') or (axis==1): len_dim1, len_dim2 = nx, nz
    else: len_dim1, len_dim2 = nx, ny
    ntask = len_dim1*len_dim2*(len_t - idx_t_disp)
    itask = 0
    for idx_dim1 in range(len_dim1):
        for idx_dim2 in range(len_dim2):
            # if DEBUG:
            #     print('Indices check:\nlen_dim1, idx_dim1, len_dim2, idx_dim2')
            #     print('{}, {}, {}, {}'.format(len_dim1, idx_dim1, len_dim2, idx_dim2))
            A = 0.
            B = 0.
            for idx_t in range(len_t - idx_t_disp):
                itask += 1
                protask = itask/ntask
                # if DEBUG:
                #     print('Indices check:\nlen_t, idx_t_disp, idx_t')
                #     print('{}, {}, {}'.format(len_t, idx_t_disp, idx_t))
                print("\r", end="")
                print(
                      '         Calculating...{}%'.format(round(protask*100, ndigits=2)), 
                      "▋" * (int(protask*100) // 2), 
                      end=""
                      )
                sys.stdout.flush()
                A += ndigitsgen_2d(
                    trj[idx_t], nx, ny, nz, axis, r
                    )[idx_dim1][idx_dim2]*ndigitsgen_2d(
                        trj[idx_t+idx_t_disp], nx, ny, nz, axis, r
                        )[idx_dim1][idx_dim2]
                B += ndigitsgen_2d(
                    trj[idx_t], nx, ny, nz, axis, r
                    )[idx_dim1][idx_dim2]
            try:
                if A == 0 and B == 0:
                    qzt_val += 1
                else:
                    qzt_val += A/B
            except RuntimeError:
                print('RuntimeWarning: invalid value encountered in double_scalars')
                print('A = {}'.format(A))
                print('B = {}'.format(B))
                exit()
    return qzt_val

# calculate one point of q(z, t)
def qz(list_t2d, idx_t_max, mode = 'default'):
    
    qz_val = []
    for idx_t in range(idx_t_max):
            
        qz_val.append(qzt(list_t2d, idx_t))
    return qz_val
# calculate a series of t of q(z,t)

def qz_2d(trj, nx, ny, nz, axis, r, idx_t_max):
    qz_val = []
    if VERBO > 1:
        print('q(z,t)| SUBROUTINE qz_2d is called, calculating qz(t), iterating length of time window')
    for idx_t in range(idx_t_max):
        if VERBO > 1:
            print('         Process...({}/{})'.format(idx_t+1, idx_t_max))
        qz_val.append(qzt_2d(trj, nx, ny, nz, axis, r, idx_t))
    return qz_val

def q(lmp_trjfile, elements, axis, rslist, idx_t_max, nx = 1000, ny = 1000, nz = 1000, mode = 'default'):

    if VERBO > 0:
        print(
            'q(z,t)| Many-body correlation function called, function from doi: 10.1038/nphys2133\n'
        +'        Visiting original paper via link: https://www.nature.com/articles/nphys2133 '
        +'\n'
        +'        Parameters information:\n'
        +'        Lammps trajectory file: {}\n'.format(lmp_trjfile)
        +'        Elements assigned for trj: {}\n'.format(elements)
        +'        Axis selected for calculating correlation: {}\n'.format(axis)
        +'        Scaled distances along {}-axis selected for calculating correlation: {}\n'.format(axis, rslist)
        +'        Max length of time window for this time correlation function: {}\n'.format(idx_t_max)
        +'        Grid level: {} x {} x {}\n'.format(nx, ny, nz)
        +'        Program running mode (normal or save_memory): {}'.format(mode)
        )
    q_val = []
    txyz = lmp2list(
        filename = lmp_trjfile,
        elements = elements,
        isort = True,
        vc = True
        )
    if mode == 'save_memory':
        idx_rs = 0
        for rs in rslist:
            print('q(z,t)| Calculate {} = {}, calling qz, process...({}/{})'.format(axis, rs, idx_rs+1, len(rslist)))
            q_val.append(
                         qz_2d(
                             trj = txyz, nx = nx, ny = ny, nz = nz,
                             axis = axis, r = rs, idx_t_max = idx_t_max)
                         )
            idx_rs += 1
    else:
        n_txyz = ndigits_t(trj = txyz, nx = nx, ny = ny, nz = nz)
        
        for rs in rslist:
            # rlist is also given by fractional z
            
            q_val.append(
                        qz(
                            list_t2d = slicegen(list_txyz = n_txyz, axis = axis, index = int(rs*1000)),
                            idx_t_max = idx_t_max
                        )
                        )
    return q_val
# calculate q, the spatial-time correlation function
