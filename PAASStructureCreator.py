import numpy as np
import scipy
from glob import glob
import pickle, time
import gzip
import os
import sys
from Bio.PDB import PDBParser, PDBIO

# TODO There should be a file with parameters so that we should not have to rely on
#      PAAS scores and PAAS IDX Clusters and Radius from the beginning     

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



def compute_paas(DM, coords, atoms_aas, max_dist=5):
    paas=[]
    paas_raw=[]
    N=coords.shape[0]
    unique_raws=set()
    for i in range(N) :
        fp=set(atoms_aas[DM[i]<max_dist])
        if len(fp)>2:
            L = sorted(list(fp))
            L2=tuple([aa[0] for aa in L])
            L_r=tuple([str(aa) for aa in L]) ## this includes the Amino Acid Position in the sequence!!

            #it is important not to repeat PAAS (even raw) so that atoms will not be counted more than one time
            if L_r not in unique_raws:
                paas.append(L2)
                paas_raw.append(L_r)
                unique_raws.add(L_r)
    #print(fps_raw)
    #paas contains the actual PAAS (in tuple format)
    #paas_raw is also a list of tuples but for each amino acid
    return paas, paas_raw



def get_paas(fname, D, max_dist=5):
    #data=open(fname, "rt").readlines()
    prot=fname.split("/")[-1].split(".")[0] 
    with gzip.open(fname, "rb") as gz_file:
        d_bytes = gz_file.read()
    d_bytes = d_bytes.decode('utf-8')
    data    = d_bytes.split("\n")

    dm, crds, aaas = extractData(data, D)
    paas, paas_raw = compute_paas(dm, crds, aaas, max_dist)
    unique_paas=set(paas)

    # Now this is where things get strange - we will create a dictionary for each unique_paas and will attach to each, each AA
    # from paas_raw and this should be easy because paas and paas raw are actually in the same order
    paas_aas={}
    aas_paas={}
    for i, p in enumerate(paas):
        unique_p = "".join(p)
        for aa_idx in paas_raw[i]:
            paas_aas.setdefault(unique_p,[]).append(aa_idx)
            aas_paas.setdefault(aa_idx,[]).append(unique_p)
    #all paas as strings. 
    unique_paas=["".join(i) for i in unique_paas]    
    prot=fname.split("/")[-1].split(".")[0]
    return unique_paas, paas_aas, aas_paas

def compute_aa_scores(aas_paas, paas2idx):
    #compute the AA scores for a full protein
    #we are getting the MAX score a residue can have
    
    aa_scores=[None]*(len(aas_paas)+1) #we will start at 1 for consistency with PDB
    for aan in aas_paas:
        aa = aan[0]
        n = int(aan[1:])
        max_score=0
        for paas in aas_paas[aan]:
            score=paas_scores[paas2idx[paas]]
            if score > max_score: 
                max_score=score
        aa_scores[n]=(max_score, aa)
    return aa_scores


def write_residue_scores_to_bfactor(input_pdb_gz, output_pdb, residue_scores, chain_id=None):
    """
    input_pdb_gz: path to .pdb.gz file
    output_pdb:   path to output .pdb file
    residue_scores: dict mapping residue number -> score
                    e.g. {1: 0.0, 2: 1.3, 3: 4.8, ...}
    chain_id: optional chain filter, e.g. 'A'
    """
    parser = PDBParser(QUIET=True)

    with gzip.open(input_pdb_gz, "rt") as handle:
        structure = parser.get_structure("prot", handle)

    for model in structure:
        for chain in model:
            if chain_id is not None and chain.id != chain_id:
                continue

            for residue in chain:
                # skip hetero/water residues
                hetflag, resseq, icode = residue.id
                resname = residue.get_resname()
                
                if hetflag != " ":
                    continue
                if residue_scores[resseq] is not None:
                    score = residue_scores[resseq][0]
                else:
                    score = 0
                
                for atom in residue: atom.set_bfactor(score)
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)


if __name__=="__main__":
    if len(sys.argv)!=5:
        print("""
This script requires 4 sequential parameters:
   1. filename with protein structure (.pdb.gz)
   2. filename with output (.pdb)
   3. filename with PAAS IDX scores produced earlier
   4. filename with PAAS clusters""")
        exit()
    fin=sys.argv[1]
    fout=sys.argv[2]
    fscores=sys.argv[3]
    fclusters=sys.argv[4]

    #get IDXScores
    #paas_scores = pickle.load(open("../PAASScores_trimmed.pickle", "rb"))
    paas_scores = pickle.load(open(fscores, "rb"))

    #get the actual clusters for reverse finding the PAAS corresponding to each IDX
    #data_labels, centroids, idx= pickle.load(open("../PAAS_Clusters.pickle", "rb"))
    data_labels, centroids, idx= pickle.load(open(fclusters, "rb"))
    paas2idx=dict(zip(data_labels,idx.tolist()))

    #here we start processing our protein
    max_dist=4.0 #WEIRD Default <- will have to do for now. 
    
    aas=GetAAData()
    D_3to1={cod3.upper(): cod1 for cod1, cod3, nam in aas}
    #prot=fin.split("/")[-1].split(".")[0]

    #get all paas data in different formats
    unique_paas, paas_aas, aas_paas = get_paas(fin, D_3to1, max_dist = max_dist)

    #compute the scores for each residue
    aa_scores = compute_aa_scores(aas_paas, paas2idx)

    #write it into a file
    write_residue_scores_to_bfactor(fin, fout, aa_scores, chain_id=None)
