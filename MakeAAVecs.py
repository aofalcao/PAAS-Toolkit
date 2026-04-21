import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import scipy
import pickle
from skbio import DistanceMatrix
from scipy.spatial import distance_matrix
from skbio.stats.ordination import pcoa

def GetAAData_MDS(N_pcs=None): #use all dimensions by default
    lines=open("aas.txt", "rt").readlines()
    aas= [lin.strip().split("\t") for lin in lines[1:]]
    
    lines=open("blosum62.txt", "rt").readlines()
    aas_blosum=lines[0].strip().split("\t")
        
    N_aa=len(aas_blosum)
    if N_pcs is None: N_pcs=N_aa
    blosum62=np.zeros((N_aa, N_aa), dtype=int)
    for row, lin in enumerate(lines[1:]):
        cols = lin.strip().split("\t")
        blosum62[row]=[int(v) for v in cols[1:]]
    
    for row, lin in enumerate(lines[1:]):
        cols = lin.strip().split("\t")
        blosum62[row]=[int(v) for v in cols[1:]]

    DM1=distance_matrix(blosum62, blosum62)
    dm1 = DistanceMatrix(DM1, ids=aas_blosum)
    res1 = pcoa(dm1, number_of_dimensions=N_pcs)   #true to life PCooA and not that scikit fake stuff

    aas_t=res1.samples.values
    aa_coords={}
    for aa_idx, aa in enumerate(aas_blosum): aa_coords[aa] = aas_t[aa_idx, :N_pcs]
    return aas, aa_coords



if __name__=="__main__":
    print("[INFO] This code will create the full AA Vectors  to be used in all PAAS")
    print("       A file called AAVecs.pickle will be created - do not change it")

    aas, aa_coords=GetAAData_MDS()
    pickle.dump((aas, aa_coords), open("AAVecs.pickle", "wb"))
    print("Done")
    
    