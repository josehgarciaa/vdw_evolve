import ase.db
import json
import numpy as np
import os.path
import spglib


def dbelem2array(x):
    shape, dtype, arr  = x['__ndarray__'];
    return np.array(arr, dtype=dtype).reshape(shape);   

def get_online_data(uid, verbose= True):
    dataurl = 'https://cmrdb.fysik.dtu.dk/c2db/row/'+uid+'/all_data';
    file = "data/"+uid+"_data.json"
    if os.path.isfile(file):
        if verbose:    
            print("file: ",file, "found")
    else:
        if verbose:    
            print("downloading from:",dataurl);
        
    with open('./data/'+uid+'_data.json', 'r') as file:
        syst_data = json.loads(next(file));
    return syst_data

def get_lattice_vectors(syst_data):
    return dbelem2array( syst_data['structure.json']['1']['cell']['array'] );

def get_atom_pos(syst_data):
    species = dbelem2array( syst_data['structure.json']['1']['numbers'] );
    positions= dbelem2array(syst_data['structure.json']['1']['positions'] );
    return list(zip(species,positions));

def get_cell(syst_data):
    lattice= get_lattice_vectors(syst_data);
    tofrac = np.linalg.inv(lattice);
    atoms  = get_atom_pos(syst_data);
    numbers, positions = zip(*atoms);
    positions = np.dot(positions,tofrac); positions-=lattice[2]/2;
    return ( np.round(lattice,5), np.round(positions,5), np.array(numbers,dtype=int) ); 

def get_symmetry(cell, tol=1e-3):
    return spglib.get_symmetry_dataset(cell, symprec=tol) ;

def get_symmetry_matrices(cell, tol=1e-3):
    return spglib.get_symmetry(cell, symprec=tol) ;

def C(a,n):
    ax,ay,az = a/np.linalg.norm(a);
    th = 2*np.pi/n;
    return np.array(
[[ np.cos(th)+ ax**2*( 1 - np.cos(th))    , ax*ay*( 1 - np.cos(th)) - az*np.sin(th), ax*az*( 1 - np.cos(th)) + ay*np.sin(th) ],
 [ ay*ax*( 1 - np.cos(th)) + az*np.sin(th), np.cos(th)+ ay**2*( 1 - np.cos(th))    , ay*az*( 1 - np.cos(th)) - ax*np.sin(th) ],
 [ az*ax*( 1 - np.cos(th)) - ay*np.sin(th), az*ay*( 1 - np.cos(th)) + ax*np.sin(th), np.cos(th)+ az**2*( 1 - np.cos(th))     ]]).T;

def sigma(a):
    return -C(a,2);

def symmetry_directions(cell):
    x,y,z   = np.eye(3);
    lat_vec, pos, numb = cell;
    a2,a3,a1= lat_vec;
    d2,d3,  = [ np.cross(a,z) for a in (a2,a3) ];
    return dict([ ("v1",a1), ("v2",a2), ("v3",a3), ("d2",d2), ("d3",d3)  ]);

def sym_label(Rn, cell):
    directions = symmetry_directions(cell);
    lat_vec, pos, numb = cell;
    toabs   = lat_vec;
    tofrac  = np.linalg.inv(toabs);
    
    R = (tofrac.dot(Rn.T).dot( toabs ));    
    for dir_idx0 in directions:
        if np.linalg.norm( R+ np.eye(3) )< 1e-3:
            return "Inversion";
        if np.linalg.norm( R- np.eye(3) )< 1e-3:
            return "Identity";
        
        u0 = directions[dir_idx0];
        if np.linalg.norm( R-sigma(u0))< 1e-3:
            return "M_"+dir_idx0;
        for n in (6,3,4,2):
            Cu0n = C(u0,n);
            Cu0nl="C"+str(n)+"_"+dir_idx0;
            if np.linalg.norm( R-Cu0n)< 1e-3:
                return Cu0nl;
            if np.linalg.norm( R+Cu0n)< 1e-3:
                return "inv"+Cu0nl;
            
    for dir_idx0 in directions:
        u0 = directions[dir_idx0];
        for n in (6,3,4,2):
            Cu0n = C(u0,n);
            Cu0nl="C"+str(n)+"_"+dir_idx0;
            for dir_idx1 in directions:
                u1 = directions[dir_idx1];
                if np.linalg.norm( R-sigma(u1).dot(Cu0n))< 1e-3:
                    return "M_"+dir_idx1+Cu0nl;
                if np.linalg.norm( R+sigma(u1).dot(Cu0n))< 1e-3:
                    return "invM_"+dir_idx1+Cu0nl;
    return "None";
