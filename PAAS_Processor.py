# this file will compute all the PAAS for all proteins 
# Only the PAAS are computed, meaning the aminoacids that appear together
# No Coordinates whatsoever

import multiprocessing as mp
import numpy as np
import scipy
from glob import glob
import pickle, time
import gzip
import os


def GetAAData():
    lines=open("aas.txt", "rt").readlines()
    aas= [lin.strip().split("\t") for lin in lines[1:]]
    return aas

def extractData(data, D_3to1):
    #creates distance matrix, atom coordinates and atoms -> amino acids
    atoms_aas=[]
    coords=[]
    for line in data:
        #fs=line.strip().split()
        ATOM = line[:4]
        if ATOM=="ATOM":
            #ATOM, aid, atm, aa, lixo, aaid, x,y,z, s1, s2, atm2=fs
            aid = line[6:11].strip()
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            aa = line[17:20].strip()
            aaid = line[22:26].strip()
            aid=int(aid)
            #x,y,z, = float(x), float(y), float(z)
            atoms_aas.append(D_3to1[aa]+aaid)
            coords.append([x,y,z])

    coords=np.array(coords)
    atoms_aas=np.array(atoms_aas)        
    DM = scipy.spatial.distance_matrix(coords, coords)
    return DM, coords, atoms_aas



def compute_fps(DM, coords, atoms_aas, max_dist=5):
    fps=[]
    N=coords.shape[0]
    for i in range(N) :
        fp=set(atoms_aas[DM[i]<max_dist])
        if len(fp)>2:
            L = sorted(list(fp))
            L2=[aa[0] for aa in L]
            fps.append(tuple(L2))
    return fps



def getFPs(fname, D, max_dist=5):
    #data=open(fname, "rt").readlines()
    prot=fname.split("/")[-1].split(".")[0] 
    with gzip.open(fname, "rb") as gz_file:
        d_bytes = gz_file.read()
    d_bytes = d_bytes.decode('utf-8')
    data    = d_bytes.split("\n")

    dm, crds, aaas = extractData(data, D)
    fps=compute_fps(dm, crds, aaas, max_dist)
    unique_fps=set(fps)
    #this is better! - a string for each PAAS
    unique_fps=["".join(i) for i in unique_fps]
    prot=fname.split("/")[-1].split(".")[0]
    #print(prot)
    return {prot: unique_fps}




def process_file(file_path, D, max_dist):
    #print("starting file:", file_path)
    R= getFPs(file_path, D, max_dist = max_dist)
    #print("... done", file_path)
    return R


counter=0
if __name__=="__main__":
    from sys import argv
    if len(argv) == 1: 
        print("(v.2) Need a distance")
        exit()
    max_dist=float(argv[1])

    mypath="paas_%3.1f_pickles" % max_dist
    os.makedirs(mypath, exist_ok=True)
    
    print("[INFO] Computing PAAS with max dist of %3.1f Angstroms" % max_dist)
    print("[INFO] getting all files...", end=" ")
    files=glob("sp_structs/*.gz")
    print(len(files), "files read - done!")
    print("[INFO] Getting AA data and processing...", end="")
    #first get all the individual amino acid data and vectorize them
    aas = GetAAData()
    #this is a useful struct to convert 3lettered AAs into 1 letter AAs
    D_3to1={cod3.upper(): cod1 for cod1, cod3, nam in aas}
    print("done")
    t1=time.time()
    print("[INFO] Processing with #Cores: %d" % (mp.cpu_count()-2))
    chnk=10000
    start=0
    for end in range(chnk,len(files), chnk):
        files_togo=files[start:end]
        pool=mp.Pool(mp.cpu_count()-2)
        res = [pool.apply_async(process_file, args=(f, D_3to1, max_dist)) for f in files_togo]
        pool.close()
        pool.join()
        results = [r.get() for r in res]
        mypickle="/paasD_r%3.1f_%d.pickle" % (max_dist, start)
        pickle.dump(results, open(mypath+mypickle, "wb"))
        print("[INFO] created file", mypickle)
        start=end
    #finalizing last bits
    print("[INFO] Last bits")
    files_togo=files[start:]
    pool=mp.Pool(mp.cpu_count()-2)
    res = [pool.apply_async(process_file, args=(f, D_3to1, max_dist)) for f in files_togo]
    pool.close()
    pool.join()
    results = [r.get() for r in res]
    mypickle="/paasD_r%3.1f_%d.pickle" % (max_dist, start)
    pickle.dump(results, open(mypath+mypickle, "wb"))
    print("[INFO] created file", mypickle)
    
    print("[INFO] Time to finish:", time.time()-t1)
