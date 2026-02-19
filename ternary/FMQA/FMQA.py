import fmqa
import numpy as np
import csv
import matplotlib.pyplot as plt
import sys
import copy
import time
import shutil
import os
import subprocess
from bit2POS import read_POSCAR, get_POSCAR

# Settings #########################################
reference_POSCAR_file='../condition/POSCAR'
file_initial_x='../initial_data/initial_bit_structures.csv'
file_initial_y='../initial_data/initial_TOTEN.csv'
alpha=0.5 # weight of penalty functions
N_iterations=30
mode = "FMQA"
#mode = "Random"
####################################################

# amplify settings #################################
import amplify
from amplify.client import FixstarsClient
from amplify import Solver
client = FixstarsClient() 
timeout_sec = 5 # sec 
client.parameters.timeout = 1000 * timeout_sec
client.token = "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
solver = Solver(client) 
####################################################

def calc_TOTEN(new_x,reference_POSCAR_file,iloop):
    get_POSCAR(reference_POSCAR_file,new_x,iloop)

    subprocess.run(['bash', 'VASP_exe.sh', str(iloop)])
    shutil.copy('./VASP_exe/OUTCAR', './data/OUTCAR_'+str(iloop))
    new_y =  np.asarray(np.loadtxt('./data/TOTEN_'+str(iloop)+'.csv',skiprows=0, delimiter=',') )

    return new_y

commentiline,scalingfactor,lattice,atom_name,atom_num,coordinateformat,coords=read_POSCAR(reference_POSCAR_file)
tot_atom_num = sum(atom_num) # Nuber of total atoms
atom_type = len(atom_num)

X = np.asarray(np.loadtxt(file_initial_x,skiprows=0, delimiter=',',dtype='int') )
y = np.asarray(np.loadtxt(file_initial_y,skiprows=0, delimiter=',') )

#judgement to generate random state
def contains(X, new_x):
    return np.any(np.all(X == new_x[None,:], axis=1))

samplingmethod=[]
for iloop in range(N_iterations):

    #training of FM model
    ynorm =  (y - np.min(y))/(np.max(y) - np.min(y))
    model = fmqa.FMBQM.from_data(X, ynorm)

    #calculation by amplify
    gen = amplify.SymbolGenerator(amplify.BinaryPoly)
    q = gen.array(model.num_variables)
    
    f =  sum([q[k] * model.linear[k] for k in model.linear])
    f += sum([q[k] * q[l] * model.quadratic[(k, l)] for (k, l) in model.quadratic])
    
    from amplify.constraint import equal_to
    from amplify.constraint import one_hot
    if len(atom_num)==2:
        g= alpha * equal_to(sum(q), atom_num[1]) #constraint on the number of Fe atoms
        constraint=g
    elif 2<len(atom_num):
        gs=[]
        for i in range(len(atom_num)):
            gs.append( alpha * equal_to(sum([q[k * atom_type + i] for k in range(tot_atom_num)]), atom_num[i]) )
        g=sum(gs) #constraint on the number of each atom
        h = alpha * sum(one_hot(q[k * atom_type:(k + 1) * atom_type]) for k in range(tot_atom_num)) #one-hot constraint
        constraint=g+h
    
    if mode == "FMQA":
        m = f + constraint
    
    if mode == "Random":
        m = constraint
    
    amplify_model = amplify.BinaryQuadraticModel(m)
    result = solver.solve(amplify_model)
    
    if len(result)!=0:
        samplingmethod.append([iloop,0])
        new_x = q.decode(result[0].values)
        new_x_int = [int(i) for i in new_x]
        new_x = np.array(new_x_int)
    else: # random
        samplingmethod.append([iloop,1])
        amplify_model = amplify.BinaryQuadraticModel(constraint)
        result = solver.solve(amplify_model)
        new_x = q.decode(result[0].values)
        new_x_int = [int(i) for i in new_x]
        new_x = np.array(new_x_int)

    while contains(X, new_x):
        samplingmethod.append([iloop,2])
        amplify_model = amplify.BinaryQuadraticModel(constraint)
        result = solver.solve(amplify_model)
        new_x = q.decode(result[0].values)
        new_x_int = [int(i) for i in new_x]
        new_x = np.array(new_x_int)
    new_y = calc_TOTEN(new_x,reference_POSCAR_file,iloop)
    X = np.r_[X, [new_x]]
    y = np.r_[y,[new_y]]

    with open('y_tmp.csv','w') as f:
        writer = csv.writer(f)
        writer.writerow(y)
    with open('X_tmp.csv','w') as f:
        writer = csv.writer(f)
        writer.writerows(X)
    with open('samplingmethod_tmp.csv','w') as f:
        writer = csv.writer(f)
        writer.writerows(samplingmethod)

with open('y_all.csv','w') as f:
    writer = csv.writer(f)
    writer.writerow(y)
with open('X_all.csv','w') as f:
    writer = csv.writer(f)
    writer.writerows(X)
with open('samplingmethod_all.csv','w') as f:
    writer = csv.writer(f)
    writer.writerows(samplingmethod)
