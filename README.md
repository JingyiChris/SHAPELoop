# SHAPELoop

SHAPELoop is an RNA secondary structure prediction tool based on characteristic SHAPE patterns for various loop motifs.

The SHAPELoop framework consists of:

- Identify characteristic SHAPE patterns.
- Generate guidance and candidate structures. 
- Calculate penalties for loops in guidance and candidate structures.
- Evaluate loops in the guidance structure.
- Select candidates.

## Prerequisites

- Python (>=3.6.10)
- [RNAstructure](https://rna.urmc.rochester.edu/RNAstructure.html) (default)
- [MC-Fold](https://major.iric.ca/MajorLabEn/MC-Tools.html) (if non-canonical base pairs are considered)

Please make sure that RNAstructure and its data tables are in your environment variable.
```
export PATH=$PATH:/path/to/RNAstructure/exe/
export DATAPATH=/path/to/RNAstructure/data_tables/
```

## Before running

It is **recommended** to define the environment variables by adding the following to your bash profiles:
```
export PATH=$PATH:/path/to/SHAPELoop/bin/
```
After defining the `$PATH` environment, you can simply run `SHAPELoop` without specifying the entire path.

A helper message is shown:

```
----------------------------------------------------------------------------------------------------
SHAPELoop: version 1.0
This step is the main procedure of SHAPELoop.
----------------------------------------------------------------------------------------------------

Usage:
	 SHAPELoop [Options] -s RNA_sequence -r SHAPE_reactivities -o output_dir

Required Parameters:
	 <RNA_sequence> <SHAPE_reactivities> <output_dir>

Options With Parameters:
	-n	<string>  	The name of this study, default: test

	-s	<filename>	Sequence file, FASTA format

	-r	<filename>	Specify a SHAPE restraints file

	-o	<directory>	Output directory

	-d	<filename>	Dataset for SHAPE pattern identification

	-c	<filename>	Dataset to be combined for SHAPE pattern identification

	-m	<filename>	Candidate structures predicted by MC-Fold, dot file

	-N	<int>     	Specify the size of Boltzmann-weighted candidate ensemble. Default is 1000 structures

	-b	<int>     	Specify the number of flanking base pairs of loops. Default is two base pairs

Options Without Parameters:
	-a	          	Provide all suboptimal structures

	-h	          	Help

----------------------------------------------------------------------------------------------------
```

## Usage of SHAPELoop

### Step 1: Prepare input files
Required input files:
- <RNA_sequence> : FASTA format sequence input.
- <SHAPE_reactivities> : The file format comprises two columns. The first column is the nucleotide number (1-based), and the second is the reactivities. Nucleotides without SHAPE data can be set as less than -500. Columns are separated by a tab. It may look like this:

 Nucleotide | Reactivity 
:--:|:--:
 1 | -999
 2 | -999
 3 | 0.9755 
 4 | 0.2680 
 5 | 0.1520 

Optional input files:
- Additional_dataset_for_SHAPE_pattern (with "-d" or "-c" options). Its format is shown below:
```
>novel1  ##name
AGGGUGAGAGUCCCGAACUGUGAAGGCAGAAGUAACAGUUAGCCUAACGCAAGGGUGUCCGUGGCGACAUGGAAUCUGAAGGAAGCGGACGGCA  ##sequence
.(((.......)))..(((((.............)))))..(((...(((..((((.((((((....)))))))))).......)))...))).  ##structure
0.420042,0.233648,0.095457,0.196786,0.394410,0.327837,0.317320,0.325740,0.199715,0.149509,0.315817,...  ##SHAPE reactivity
>novel2 
AACCUUCGGUCUGAGGAACACGAACUUCAUAUGAGGCUAGGUAUCAAUGGAUGAGUUUGCAUAACAAAACAAAGUCCUUUCUGCCAAAGUUGGUACAGAGUAAAUGAAGCAGAUUGAUGAAGGGA
..(((((((((((.(.....)...(((((......(((....((((((((((..((((........))))...))))...........))))))....)))...))))))))))...))))))..
0.734082,0.588421,0.215076,0.053556,0.215826,0.228872,0.274630,0.158168,0.140292,0.126801,...
...
```
- MC-Fold_candidates (with the "-m" option). Its format is shown below:
```
>test2
GCCGUGAUAGUUUAAUGGUCAGAAUGGGCGCUUGUCGCGUGCCAGAUCGGGGUUCAAUUCCCCGUCGCGGCGCCA
((((((((................(((((((.....)))).)))....((((......)))).))))))))....
((((((((........((((....(((((((.....)))).)))))))((((......)))).))))))))....
...
```
### Step 2: Predict RNA secondary structures

You can use the provided example data to run SHAPELoop:
```
SHAPELoop -s /path/to/SHAPELoop/example/test.seq -r /path/to/SHAPELoop/example/test.shape -o /workspace/test
```

If you want to update the characteristic SHAPE patterns for loop motifs, you can use the '-c' option to combine your own SHAPE data with SHAPELoop provided data:
```
SHAPELoop -s /path/to/SHAPELoop/example/test.seq -r /path/to/SHAPELoop/example/test.shape -o /workspace/test -c /path/to/SHAPELoop/example/novel_data
```

Alternatively, you can use SHAPELoop to identify characteristic SHAPE patterns only based on your data.
```
SHAPELoop -s /path/to/SHAPELoop/example/test.seq -r /path/to/SHAPELoop/example/test.shape -o /workspace/test -d /path/to/SHAPELoop/example/novel_data
```

If you want to use MC-Fold predicted structures as the candidate structures, please use the '-m' option:
```
SHAPELoop -s /path/to/SHAPELoop/example/test.seq -r /path/to/SHAPELoop/example/test.shape -o /workspace/test -m /path/to/SHAPELoop/example/mcfold.dot
```


#### Output files:
```
test/
├── guidance_structure
│   ├── test.guidance.ct
│   └── test.guidance.dot
│
├── candidate_ensemble
│   ├── test.ensemble.ct  
│   ├── test.ensemble.dot
│   └── test.ensemble.dot.tmp
│
├── penalty
│   ├── test.guidance.penalty
│   └── test.candidates.penalty
│
└── test.SHAPELoop.dot
```
> **Note:**
> * The `guidance_structure` folder contains the MFE structure predicted with SHAPE restraints.
> * The `candidate_ensemble` folder contains candidates sampled with Boltzmann conditional probabilities and SHAPE restraints.
> * The `penalty` folder contains loop penalties for the guidance and candidate structures and may look like this:
 
| Column  | Description   |
|---------|---------------|
| Column1 | Candidate ID  |
| Column2 | Loop type     |
| Column3 | Loop length   |
| Column4 | Loop position |
| Column5 | Penalty       |

> * The `test.SHAPELoop.dot` contains the predicted structures.

