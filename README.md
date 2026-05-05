# PAAS-Toolkit 

**This is the Instruction Manual of PAAS**

*All the tools needed for creating Proximal Amino Acid Sets. form protein structure to protein fingerprint construction and searching*


## Part 1. Preparation

### 1.1. Create the vectors for each aminoacid  

Requires files aas.txt and blosum62.txt

   `$ MakeAAVecs.py  --> creates a Pickle with ALL the dimensions of the aminoacid`

### 1.2. Create the PAAS for al the proteins

requires aas.txt, the name of a folder where all the ALphaFold Structures are stored (sp_structs/ currently) and an input is the Radius in Angstroms

    `$ PAAS_processor [DIST] -> Creates a set of files paasD_[DIST].pickle inside a folder with all the PAAS are computed`

### 1.3. Create the PAAS Indexes and the centroids for dataset transformation. Also all the known PAAS will be given a index without the need for computation

    `$ MakePAASCentroids  -> creates a pickle with the exiting PAAS and corresponding indexes and the coordinates of the centroids for future assignments`



## PART 2. Attribute PAAS Indexes to proteins

This part has 2 possible components: 
a) creating the PAAS for the proteins in the dataset. That is easy, because we already computed the IDXs for each Exiting PAAS
b) creating the PAAS Indexes for a new file - This is more complex as there may exist some PAAS never computed before, so we need to compute the vectors of each unknown PAAS and assign it to the appropriate cluster

### 2.1. Assign the PAAS IDXs to a file with already defined PAAS (the dataset)

   `$ ProcessPAASidxs [folder with Pickles with existing PAAS] ---> creates AllProtsPaas.pickle`

This will read the PAAS idxs and Centroids defined at 1.3 and go for all file in the pickle folder to process all the PAAS


#### 2.1.1.  It may be important to trim the least significant PAAS - this is especially true for Searching and PC modeling

`$ python PAASTrimmer.py [#bits to remove] [trimmed file name]`

With this tool  it was created a file that for each protein has only the relevant PAAS with the exact same structure of `AllProtPaas.pickle`



## PART 3.  SEARCHING

A draft of the searching algorithm is in Searcher-Draft.ipynb

It has been implemented as a webapp and is VERY fast and is looks very robust
 

## PART 4. Structure inspection

### 4.1. Create PAAS Scores

The first step is to create a score for each PAAS IDX. This is done using a modification of PAASTrimmer [this pipleine will probably change in the future]. It creates a score, and then transforms that score with a scaled power transform for [0, 1] output.

`$ python PAASScorer.py [#bits to remove] [trimmed file name]`

### 4.2 Create PDB files for visualization in PDB visualization software (PyMol, Chimera, ...)

Currently there is just `PAASStructureCreator.py`. This is a file that using a .gz PDB file is able to create a new PDB file with PAAS scores inputed for each amino acid residue. 

This script requires 4 sequential parameters:
   1. filename with protein structure (.pdb.gz)
   2. filename with output (.pdb)
   3. filename with PAAS IDX scores produced in 4.1
   4. filename with PAAS clusters

Example

`$ python PAASStructureCreator.py sample_files/AF-P22374-F1-model_v4.pdb.gz \
          AF-P22374-F1-model_v4.paas.pdb \
          PAASScores_trimmed.pickle \
          PAAS_Clusters.pickle`





