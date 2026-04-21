import numpy as np
#from sklearn.decomposition import PCA
#from sklearn.preprocessing import StandardScaler
import scipy
from glob import glob
import pickle, time
#from skbio import DistanceMatrix
#from scipy.spatial import distance_matrix
#from skbio.stats.ordination import pcoa
import gzip
import sys
from scipy.spatial import cKDTree
from sklearn.cluster import KMeans

PAAS_DIM =  16

#this is a UNIQUE approach <- only the UNIQUE PAAS are preserved irrespective of the quantity
def get_unique_paas(fname):
    data=pickle.load(open(fname, "rb"))
    all_paas=set()
    for prot in data:
        key=list(prot.keys())[0]
        name=key.split("-")[1]
        paas=prot[key]
        #print(name, len(paas), paas[:5])
        all_paas|=set(paas)
    return all_paas

def compute_coords(all_paas, aa_coords):
    L_paas=list(all_paas)
    COORDS_PAAS={}
    for paas in L_paas:
        coords = np.zeros([1,PAAS_DIM])
        for aa in paas: coords+=aa_coords[aa]
        COORDS_PAAS[paas]=coords
    return COORDS_PAAS

def get_all_paas(path):
    files=glob(path+"/*.pickle")

    all_paas = set()
    for file in files:
        print("   processing", file, end=" ")
        file_paas = get_unique_paas(file)
        all_paas |= file_paas
        print("Done:",len(all_paas))
    return all_paas


if __name__=="__main__":
    if len(sys.argv)!=2: 
        print("Ooops. Need a path where the pickles are!")
        exit()
    
    path=sys.argv[1]
    
    #this is now FIXED!!! - All Coordinates come from here
    aas, aa_coords_full=pickle.load(open("AAVecs.pickle", "rb"))
    #but now we trim them accordint to the PAAS_DIM 
    aa_coords={}
    for aa in aa_coords_full: aa_coords[aa]=aa_coords_full[aa][:PAAS_DIM]


    print("[INFO] STEP 1 - get all PAAS for clustering")
    all_paas = get_all_paas(path)
    
    print("[INFO] Step 2 - Compute all PAAS Coordinates")
    paas_coords=compute_coords(all_paas, aa_coords)
    data_coords=[]
    data_labels=[]
    #not taking any chances on weirsd orders of dictionaries
    for paas in paas_coords:
        data_labels.append(paas) 
        data_coords.append(paas_coords[paas][0])
    data_coords=np.array(data_coords)
    print("   N: Of Unique PAAS:", data_coords.shape[0])

    print("[INFO] Step 3 - Do the clustering")
    print("[INFO] Step 3.1 - First, Level 0 Clustering") 
    kmeans_0 = KMeans(n_clusters=144, random_state=0, n_init="auto", max_iter=9999)
    kmeans_0.fit(data_coords)
    
    # now for each PAAS, assign a cluster
    print("[INFO] Step 3.2 - Assign the PAAS to a Centroid of Level 0")
    clusters_0 = kmeans_0.predict(data_coords)

    #here I will separate all the PAAS according to the Cluster_0 each one belongs
    paas_clusters_0={}
    for i, cl in enumerate(clusters_0):
        paas_clusters_0.setdefault(cl, []).append(data_coords[i])

    #############################################
    #THERE IS A MEMORY LEAK IN KMEANS. IF THIS IS NOT DONE, WE GET NASTY WARNINGS
    import warnings
    warnings.filterwarnings("ignore")
    #############################################

    print("[INFO] Step 3.3 - Do the Level 1 Clustering")
    centroids=[]
    for i in range(144):
        kmeans_1 = KMeans(n_clusters=144, random_state=0, n_init="auto", max_iter=9999)
        data_coords_1=np.array(paas_clusters_0[i])
        kmeans_1.fit(data_coords_1)
        centroids.append(kmeans_1.cluster_centers_)
    
    all_centroids=centroids[0]
    for i in range(1,144):
        all_centroids=np.vstack((all_centroids, centroids[i]))

    print("[INFO] Step 4 - Assign each PAAS to the Level 1 Clusters")

    tree = cKDTree(all_centroids)
    dist, idx = tree.query(data_coords, k=1, workers=16)
    
    print("[INFO] Step 5 - Save the PAAS, Centroids and assignments as pickles")
    pickle.dump((data_labels, centroids, idx), open("PAAS_Clusters.pickle", "wb"))
    print("Done")
    
    