import pickle
import numpy as np
import glob
import sys

def get_prots_from_file(fname, paas2idx):
    # Get the data as a list of 1 KEY dictionaries
    # and transform it into a proper PAAS Dictionary for each protein in the file
    
    prots = pickle.load(open(fname, "rb"))
    all_prots={}
    for p in prots:
        key=list(p.keys())[0]
        #all_prots[key] = np.array([paas2idx[paas] for paas in p[key]], dtype=np.int32)  
        #only UNIQUES!!! - Think well!!!
        all_prots[key] = np.unique(np.array([paas2idx[paas] for paas in p[key]], dtype = np.int32))

    return all_prots


if __name__=="__main__":
    if len(sys.argv)!=2:
        print("OOps! - Need the folder where the PAAS pickles are!")
        exit()

    path=sys.argv[1]
    #first get the centroid data and PAAS with existing IDs
    all_paas, centroids, clusters =pickle.load(open("PAAS_Clusters.pickle", "rb"))
    paas2idx=dict(zip(all_paas, clusters))
    
    fnames = glob.glob(path+"/*.pickle")

    all_prot_paas={}
    for fname in fnames:
        print("[INFO] - Processing ", fname)
        prot_paas = get_prots_from_file(fname, paas2idx)
        all_prot_paas.update(prot_paas)

    print("[INFO] processed", len(all_prot_paas), "proteins")
    print("[INFO] Saving the Pickle with all the paas")
    pickle.dump(all_prot_paas, open("AllProtPaas.pickle", "wb"))
    print("[INFO] finished")
    