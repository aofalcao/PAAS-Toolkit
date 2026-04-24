# this file is just a modification of PAASTrimmer with the exact same scoring function
# It will compute scores for all PAAS IDx with the sole purpose of 


import math
from collections import Counter, defaultdict
 
import numpy as np
import pandas as pd
import pickle
import sys


def compute_scores(protein_to_labels, alpha= 1.0,beta = 1.0,tau = 5.0,ic_cap = 6.0,
                   gamma = 1.0):
    # Compute rarity/reliability scores for labels from a dictionary:

    protein_ids = list(protein_to_labels.keys())
    N = len(protein_ids)
    label_counts = Counter()
    
    print("[INFO] get and count the unique labels from all the data")
    for protein_id, labels in protein_to_labels.items():
        if labels is None: labels = []
        if isinstance(labels, np.ndarray): labels = labels.tolist()

        # Usually repeated labels within one protein should not count multiple times
        labels_for_counting = set(labels)
        label_counts.update(labels_for_counting)

    all_labels = sorted(label_counts.keys())
    rows = []
    print("[INFO] compute IC and scores")
    for label in all_labels:
        n_l = int(label_counts[label])
        p_raw = n_l / N
        p_tilde = (n_l + alpha) / (N + alpha + beta)
        IC = -math.log(p_tilde)
        IC_capped = min(IC, ic_cap) if ic_cap is not None else IC
        R_l = 1.0 - math.exp(-n_l / tau)
        g_l = 1
        D_l = 1.0
        score = IC_capped * R_l * (D_l ** gamma)

        rows.append({ "label": label, "n_l": n_l, "p_raw": p_raw, "p_tilde": p_tilde,
                      "IC": IC, "IC_capped": IC_capped, "R_l": R_l,
                      "g_l": g_l,"D_l": D_l, "score": score})

    result = pd.DataFrame(rows).sort_values("score",
                                            ascending=False).reset_index(drop=True)
    return result

def compute_view_scores(scores, ic_labels, bad_labels):
    #applies a power transform for smoothing the low end and getting a range between [0,1]
    d=np.array(scores['IC'])
    mx = d.max()
    mn = d.min()
    a=2.5
    paas_scores={}
    for ic, label in ic_labels:
        if label not in bad_labels:
            score=float(((ic-mn)/(mx-mn))**a)
        else:
            score=0
        paas_scores[label] = score
    return paas_scores



if __name__ == "__main__":
    if len(sys.argv)!=3:
        print("Ooops! Need the prunning size (#labels to remove)  and the output file name (pickle file)")
        exit()
    prunning_size = int(sys.argv[1])
    outfname=sys.argv[2]
    
    print("[INFO] reading 'AllProtPaas.pickle' file with all the Idxs for the database")
    all_paas=pickle.load(open("AllProtPaas.pickle", "rb"))

    scores = compute_scores(all_paas, gamma=0.0)

    # this next part is actually irrelevant in most cases. The scores are already sorted
    # and are directly related to the IC in the low end (the one we will remove)
    
    print("[INFO] sorting by Information Content and identifying the labels to remove") 
    ic=np.array(scores['IC']).tolist()
    label=np.array(scores['label']).tolist()
    ic_labels=list(zip(ic,label))
    ic_labels.sort()
    ic_labels.reverse()
    
    #this is the low end we will remove
    bad_labels=[ label for ic, label in ic_labels[-prunning_size:]]
    bad_labels=set(bad_labels)
    
    print("[INFO] trimming low end of low IC") 
    good_labels={}
    for prot in all_paas:
        full_paas=set(all_paas[prot].tolist())
        trimmed_paas=np.array(list(full_paas - bad_labels))
        trimmed_paas.sort()
        good_labels[prot]=trimmed_paas

    
    print("[INFO] computing viewing scores (power transformed between 0-1") 
    paas_scores = compute_view_scores(scores, ic_labels, bad_labels)
    
    print("[INFO] Saving %s file with scores for PAAS" %(outfname))
    pickle.dump(paas_scores, open(outfname, "wb"))
    print("Done")
