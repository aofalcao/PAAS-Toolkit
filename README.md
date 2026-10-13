# PAAS-Toolkit

**PAAS-Toolkit** contains the reference implementation used to construct and analyse **Proximal Amino Acid Sets (PAAS)**, a protein representation based on local, sequence-order-independent three-dimensional amino-acid neighbourhoods.

The toolkit implements the workflow used in the accompanying PAAS study:

1. construct a substitution-aware numerical representation of amino acids from BLOSUM62;
2. extract local PAAS neighbourhoods from protein structures;
3. construct a discrete PAAS vocabulary by two-level K-means clustering;
4. assign PAAS labels to proteins;
5. compute label information-content scores;
6. project PAAS scores back onto protein residues for structural inspection.

The repository also contains `Sample_Data/`, with 20 protein structures that can be used to verify the PAAS extraction workflow on a small dataset.

> **Important:** the 20-protein sample is intended as a usage example. It is **not large enough to reconstruct the 20,736-label vocabulary used in the paper**, because the production clustering procedure uses 144 first-level clusters followed by 144 second-level clusters within each first-level cluster. Vocabulary construction should therefore be run on a large structural corpus.

---

## Repository contents

| File | Purpose |
|---|---|
| `MakeAAVecs.py` | Builds the BLOSUM62-derived amino-acid coordinate representation and writes `AAVecs.pickle`. |
| `PAAS_Processor.py` | Extracts unique PAAS neighbourhoods from gzipped PDB files at a user-specified distance cutoff. |
| `MakePAASCentroids.py` | Builds the two-level 144 × 144 PAAS vocabulary and writes `PAAS_Clusters.pickle`. |
| `ProcessPAASidxs.py` | Assigns vocabulary labels to the PAAS previously extracted for each protein and writes `AllProtPaas.pickle`. |
| `PAASTrimmer.py` | Removes a specified number of the lowest-information PAAS labels from `AllProtPaas.pickle`. |
| `PAASScorer.py` | Computes label information-content scores and writes a label-to-score pickle for residue visualization. |
| `PAASStructureCreator.py` | Projects PAAS scores onto protein residues and writes the scores into the PDB B-factor field for visualization. |
| `aas.txt` | Amino-acid lookup table used by the toolkit. |
| `blosum62.txt` | BLOSUM62 substitution matrix used to construct the amino-acid representation. |
| `Sample_Data/` | Twenty example protein structures for testing PAAS extraction. |

---

## Software environment

The analyses reported in the study were performed using two Python environments.

### Main PAAS analysis environment

Used for amino-acid vectorization, PAAS extraction, clustering, scoring and downstream processing:

- Python 3.12.9
- NumPy 2.0.1
- SciPy 1.15.1
- pandas 2.2.3
- scikit-bio 0.6.3
- scikit-learn 1.6.1

The corresponding dependencies are listed in:

```text
requirements-main.txt
```

A compatible environment can be created, for example, with:

```bash
python3.12 -m venv paas-main
source paas-main/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements-main.txt
```

### Structure-writing environment

`PAASStructureCreator.py` uses Biopython to read and write PDB structures.  
The study used:

- Python 3.14.4
- Biopython 1.87
- NumPy 2.4.4
- SciPy 1.17.1

Dependencies are listed in:

```text
requirements-pdb.txt
```

---

# Quick start with `Sample_Data`

The sample dataset contains 20 structures and is intended to demonstrate **PAAS extraction**.

The current `PAAS_Processor.py` expects gzipped structures under a directory named `sp_structs/`:

```python
glob("sp_structs/*.gz")
```

Therefore, from the repository root, expose the supplied sample directory under that name. On Linux/macOS the simplest approach is a symbolic link:

```bash
ln -s Sample_Data sp_structs
```

Alternatively, copy the sample structures:

```bash
mkdir -p sp_structs
cp Sample_Data/*.gz sp_structs/
```

The structures supplied to `PAAS_Processor.py` must be gzip-compressed PDB files (`*.pdb.gz` or equivalent `*.gz` PDB files).

The PAAS study used a **4.0 Å neighbourhood radius**, so the sample extraction can be run with:

```bash
python PAAS_Processor.py 4.0
```

The program processes every `*.gz` structure in `sp_structs/` and creates:

```text
paas_4.0_pickles/
```

For the 20 supplied structures, the output will normally contain a single pickle such as:

```text
paas_4.0_pickles/paasD_r4.0_0.pickle
```

Each protein is represented by its set of unique local PAAS strings.

This step is sufficient to verify that protein structures can be converted into PAAS neighbourhoods.

---

# Full PAAS vocabulary construction

The following steps reproduce the **methodological pipeline** used to build a PAAS vocabulary from a large structural corpus.

They should not be interpreted as an instruction to fit the production 20,736-label vocabulary from the 20-protein sample dataset.

## 1. Construct amino-acid vectors

`MakeAAVecs.py` reads `aas.txt` and `blosum62.txt`, converts the BLOSUM62 substitution relationships into a continuous amino-acid coordinate representation using principal coordinates analysis, and writes the result to `AAVecs.pickle`.

Run from the repository root:

```bash
python MakeAAVecs.py
```

Output:

```text
AAVecs.pickle
```

The clustering code subsequently retains the first 16 dimensions of this representation.

---

## 2. Extract PAAS neighbourhoods from structures

Place the structural corpus as gzip-compressed PDB files in:

```text
sp_structs/
```

Then run:

```bash
python PAAS_Processor.py 4.0
```

where `4.0` is the neighbourhood radius in Å.

For radius `R`, output files are written to:

```text
paas_R_pickles/
```

For the study setting:

```text
paas_4.0_pickles/
```

Large collections are written in batches of 10,000 proteins.

A PAAS is stored as an order-independent amino-acid string representing the unique amino-acid composition of a local three-dimensional neighbourhood.

---

## 3. Construct the PAAS vocabulary

After `AAVecs.pickle` and the corpus PAAS pickle files have been generated, build the PAAS vocabulary with:

```bash
python MakePAASCentroids.py paas_4.0_pickles
```

Output:

```text
PAAS_Clusters.pickle
```

The current implementation performs:

1. first-level K-means clustering into 144 clusters;
2. second-level K-means clustering into 144 clusters inside each first-level cluster;
3. assignment of each observed PAAS to the closest final centroid.

This produces a maximum vocabulary of:

```text
144 × 144 = 20,736 labels
```

Both K-means stages use:

```text
random_state = 0
n_init = "auto"
max_iter = 9999
```

### Dataset-size requirement

Each first-level cluster must contain enough PAAS observations to support a 144-cluster second-level K-means fit. A small demonstration set such as `Sample_Data/` is therefore not suitable for rebuilding the production vocabulary.

---

## 4. Assign PAAS labels to proteins

Once `PAAS_Clusters.pickle` exists, assign vocabulary labels to every protein in the PAAS pickle directory:

```bash
python ProcessPAASidxs.py paas_4.0_pickles
```

Output:

```text
AllProtPaas.pickle
```

This file is a dictionary mapping each protein identifier to the set of unique PAAS label indices assigned to that protein.

---

# Information-content scoring

PAAS labels differ in their prevalence across the structural corpus.  
The toolkit derives label scores from their protein-level prevalence, so common labels receive low information content and rarer labels receive higher information content.

## 5. Optional trimming of common labels

To remove the `N` lowest-information labels:

```bash
python PAASTrimmer.py N OUTPUT.pickle
```

For example:

```bash
python PAASTrimmer.py 1000 AllProtPaas_trimmed.pickle
```

Input is read automatically from:

```text
AllProtPaas.pickle
```

The resulting file preserves the protein-to-label structure while excluding the specified number of lowest-information labels.

Trimming is optional and is not required simply to compute residue scores.

---

## 6. Compute PAAS scores for visualization

Run:

```bash
python PAASScorer.py N PAASScores.pickle
```

where `N` is the number of lowest-information labels to suppress.

For example:

```bash
python PAASScorer.py 1000 PAASScores_trimmed.pickle
```

`PAASScorer.py` reads:

```text
AllProtPaas.pickle
```

and writes a dictionary mapping PAAS label identifiers to visualization scores.

The underlying information content is based on label prevalence across proteins. For structural visualization, the scores are transformed to the interval `[0,1]`; labels excluded by the pruning step receive a score of zero.

---

# Residue-level structural visualization

`PAASStructureCreator.py` maps PAAS scores back to individual residues. Each residue receives the **maximum score of the PAAS neighbourhoods in which it participates**.

The resulting score is written to the PDB **B-factor field**, allowing direct visualization in PyMOL, ChimeraX or other molecular-graphics programs.

Usage:

```bash
python PAASStructureCreator.py \
    INPUT.pdb.gz \
    OUTPUT.pdb \
    PAASScores.pickle \
    PAAS_Clusters.pickle
```

Example:

```bash
python PAASStructureCreator.py \
    Sample_Data/AF-P22374-F1-model_v4.pdb.gz \
    AF-P22374-F1-model_v4.paas.pdb \
    PAASScores_trimmed.pickle \
    PAAS_Clusters.pickle
```

The current implementation uses a fixed **4.0 Å** PAAS radius for this residue-mapping step, matching the setting used in the study.

The resulting PDB file can then be coloured by B-factor in a molecular viewer.

---

# Minimal workflow summary

For a **large structural corpus**, the complete workflow is:

```bash
# 1. Amino-acid representation
python MakeAAVecs.py

# 2. Extract 4 Å PAAS
python PAAS_Processor.py 4.0

# 3. Build the 20,736-label vocabulary
python MakePAASCentroids.py paas_4.0_pickles

# 4. Assign labels to proteins
python ProcessPAASidxs.py paas_4.0_pickles

# 5. Compute information-content scores
python PAASScorer.py 1000 PAASScores_trimmed.pickle

# 6. Map scores to a structure
python PAASStructureCreator.py \
    INPUT.pdb.gz \
    OUTPUT.pdb \
    PAASScores_trimmed.pickle \
    PAAS_Clusters.pickle
```

For the supplied **20-protein `Sample_Data/`**, use steps 1–2 to verify PAAS extraction. Reconstructing the full 144 × 144 vocabulary requires a substantially larger structural corpus.

---

# Main output files

| Output | Produced by | Description |
|---|---|---|
| `AAVecs.pickle` | `MakeAAVecs.py` | BLOSUM62-derived amino-acid vectors. |
| `paas_4.0_pickles/*.pickle` | `PAAS_Processor.py` | PAAS extracted from each protein structure. |
| `PAAS_Clusters.pickle` | `MakePAASCentroids.py` | PAAS strings, hierarchical cluster centroids and assigned label indices. |
| `AllProtPaas.pickle` | `ProcessPAASidxs.py` | Protein-to-PAAS-label mapping. |
| user-defined trimmed pickle | `PAASTrimmer.py` | Protein-to-label mapping after low-information labels are removed. |
| user-defined score pickle | `PAASScorer.py` | PAAS label-to-score mapping for structural visualization. |
| user-defined `.pdb` | `PAASStructureCreator.py` | Structure with residue PAAS score stored in the B-factor field. |

---

# Reproducibility notes

The implementation used for the PAAS study used the following fixed methodological settings:

- local neighbourhood radius: **4.0 Å**;
- amino-acid representation derived from **BLOSUM62**;
- retained PAAS vector dimensionality: **16**;
- first clustering level: **144 K-means clusters**;
- second clustering level: **144 K-means clusters per first-level cluster**;
- total vocabulary size: up to **20,736 PAAS labels**;
- K-means random seed: **0**.

The production PAAS vocabulary reported in the paper was constructed from the large reviewed UniProtKB/Swiss-Prot AlphaFold structural corpus, not from `Sample_Data/`.

The sample structures are provided only to make the structure-to-PAAS portion of the workflow readily testable.

---

# Data and code availability

Source code is available from this repository under the MIT license.

The processed datasets associated with the PAAS study are archived on Zenodo:

**DOI: 10.5281/zenodo.20082135**

The software snapshot used for the study is archived on Zenodo:

**DOI: 10.5281/zenodo.21724410**

---
# Hardware and Operating system

Most computations were performed on a **Dell Power Edge R530 with 2x Xeon E5-2630 v4 @ 2.2 GHz and 64 Gb of RAM** running Debian 12.8

### Hardware requirements
There are no specialized hardware requirements, and the code has been run on Windows 11, MacOS (M processors) and Linux

### Installation time

The PAAS-Toolkit itself does not require compilation or installation and can be obtained directly with:

`git clone https://github.com/aofalcao/PAAS-Toolkit`

Python dependencies must be installed separately using the supplied requirements-main.txt and requirements-pdb.txt files. Installation time therefore depends mainly on the Python environment, package manager and network connection.

### Demo Runtime

The principal computational cost in the PAAS-Toolkit is extraction of PAAS neighbourhoods from protein structures; runtime depends on the number and size of the structures and on the neighbourhood radius. Using the supplied 20-structure Sample_Data dataset and a 4.0 Angstrom radius, all 20 structures were processed in 21.0 seconds on the Dell PowerEdge R530 described above, corresponding to approximately **1 second per structure**.

**Note:** This runtime refers only to PAAS extraction from the 20 supplied example structures. The sample dataset is intended to demonstrate and verify the structure-to-PAAS workflow and is not sufficiently large to reconstruct the full 20,736-label PAAS vocabulary used in the study.


# Citation

If you use PAAS or PAAS-Toolkit in published work, please cite the associated PAAS publication and the archived software release.

---

# License

PAAS-Toolkit is distributed under the **MIT License**. See `LICENSE` for details.
