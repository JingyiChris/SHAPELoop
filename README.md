# SHAPELoop

SHAPELoop is an RNA secondary structure prediction tool developed based on conserved SHAPE patterns of various loop motifs.

The SHAPELoop framework consists of:

- Generate guidance and candidate structures. 
- Calculate penalties for loops in guidance and candidate structures.
- Classify loops in the guidance structure.
- Select candidates.

## Prerequisites

- Python (version >= 2.7.15)
- RNAstructure (version 6.1)

Please make sure that RNAstructure and its data tables are in your environment variable.
```
export PATH=$PATH:/path/to/RNAstructure/exe/
export DATAPATH=/path/to/RNAstructure/data_tables/
```

## Before runnng

It is **recommended** to define the environment variables by add the following to your bash profiles:
```
export PATH=$PATH:/path/to/SHAPELoop/bin/
```
After defining the `$PATH` environment, you can simply run `SHAPELoop` without specifying the entire path.

A helper message is shown:

```
----------------------------------------------------------------------------------------------------
SHAPELoop: version 1.0
This step is the main precedure of SHAPELoop.
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

	-N	<int>     	Specify the size of Boltzmann-weighted candidate ensemble. Default is 1000 structures

	-b	<int>     	Specify the number of flanking base pairs of loops. Default is 2 base pairs

Options Without Parameters:
	-a	          	Provide all suboptimal structures

	-h	          	Help

----------------------------------------------------------------------------------------------------
```

## Usage of SHAPELoop

### Step 1: Prepare input files

- <RNA_sequence> : FASTA format sequence input.
- <SHAPE_reactivities> : The file format comprises two columns. The first column is the nucleotide number (1-based), and the second is the reactivities. Nucleotides without SHAPE reactivities can be set as values less than -500. Columns are separated by tab. It may look like this:

 Nucleotide | Reactivity 
------------|------------
 1 | -999
 2 | -999
 3 | 0.9755 
 4 | 0.2680 
 5 | 0.1520 

### Step 2: Predict RNA secondary structures

You can use the provided example data to run SHAPELoop:
```
SHAPELoop -s /path/to/SHAPELoop/example/test.seq -r /path/to/SHAPELoop/example/test.shape -o /workspace/test
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
> * The `penalty` folder contains loop penalties for the guidance strcucture and candidate structures, and may look like this:
 
| Column  | Description   |
|---------|---------------|
| Column1 | Candidate ID  |
| Column2 | Loop type     |
| Column3 | Loop length   |
| Column4 | Loop position |
| Column5 | Penalty       |


> * The `test.SHAPELoop.dot` contains the predicted structures.
