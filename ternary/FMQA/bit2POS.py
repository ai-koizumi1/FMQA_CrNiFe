import numpy as np
import random
import sys
import csv

def write_POSCAR(commentiline,scalingfactor,lattice,atom_name,atom_num,coordinateformat,coords,file_id):
    with open('./data/POSCAR_'+str(file_id),'w') as f:
        f.write(commentiline+"\n")
        f.write(f"{scalingfactor:.2f} \n")
        f.write(f"{lattice[0,0]:.8f} {lattice[0,1]:.8f} {lattice[0,2]:.8f} \n")
        f.write(f"{lattice[1,0]:.8f} {lattice[1,1]:.8f} {lattice[1,2]:.8f} \n")
        f.write(f"{lattice[2,0]:.8f} {lattice[2,1]:.8f} {lattice[2,2]:.8f} \n")

        atom=""
        for iatom in atom_name:
            atom=atom+iatom+' '
        f.write(atom+"\n")
        atom=""
        for iatom in atom_num:
            atom=atom+str(iatom)+' '
        f.write(atom+"\n")
        f.write(coordinateformat+"\n")

        for x,y,z in coords:
            f.write(f"{x:.16f} {y:.16f} {z:.16f} \n")
    f.close()

def read_POSCAR(filename):
    coords=[]
    atom_name=[]
    atom_num=[]
    lattice=[]
    with open( filename , 'r') as f:
        lines = f.readlines()
        commentiline = lines[0].strip()
        scalingfactor = float(lines[1].strip())
        lattice=[[[0.0] for i in range(3)] for j in range(3)]
        for i in range(3):
            line = lines[2+i].strip()
            for j in range(len(line.split())):
                lattice[j][i]=float(line.split()[j])
        lattice=np.array(lattice)
    
        line = lines[5].strip()
        for i in line.split():
            atom_name.append(i)
        atom_name=np.array(atom_name)
    
        line = lines[6].strip()
        for i in line.split():
            atom_num.append(int(i))
        atom_num=np.array(atom_num)
    
        coordinateformat = lines[7].strip()
    
        for line in lines[8:8+sum(atom_num)]:
            line = line.strip()
            x = line.split()[0]
            y = line.split()[1]
            z = line.split()[2]
            coords.append([float(x),float(y),float(z)])
        coords=np.array(coords)

    return commentiline,scalingfactor,lattice,atom_name,atom_num,coordinateformat,coords

def bit2id(bit,atom_num):
    if len(atom_num)==2:
        atom_id=np.array(bit)
    if len(atom_num)==3:
        bit=np.reshape(bit,[sum(atom_num),len(atom_num)])
        atom_id=[]
        for x1,x2,x3 in bit:
            if x1==1:
                atom_id.append(0)
            elif x2==1:
                atom_id.append(1)
            elif x3==1:
                atom_id.append(2)
            else:
                print('Error')
    return np.array(atom_id)

def get_POSCAR(reference_POSCAR_file,bit,iloop):
    commentiline,scalingfactor,lattice,atom_name,atom_num,coordinateformat,coords=read_POSCAR(reference_POSCAR_file)
    atom_id=bit2id(bit,atom_num)
    id_sort=np.argsort(atom_id)
    write_POSCAR(commentiline,scalingfactor,lattice,atom_name,atom_num,coordinateformat,coords[id_sort],iloop)
