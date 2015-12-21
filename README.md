# IsoMut: a robust method for calling unique mutations (SNVs and small indels) from multiple isogenic samples

This is a C implementation of the mutation calling algorith described here (article to come ... ), and python wappers fro parallel exectution, and the BAQ based filtering of SNVs. 


---

### For the impatient:

- compile IsoMut
- modify the pathnames for your data in example_script_isomut_w_pp.py, and run the script
- check out test_isomut_w_pp.ipynb, for necessities, and details

---

### Compile IsoMut:

```
gcc -O3 -c isomut_lib.c fisher.c  -W -Wall
gcc -O3 -o isomut isomut.c isomut_lib.o  fisher.o -lm -W -Wall

```

---

### Files:

For the standalone C application:
- isomut_lib.h isomut_lib.c  the library
- isomut.c the application

For the parallel wrapper:
- isomut_wrappers.py the python functions for parallel execution, and BAQ filtering
- example_script_isomut_w_pp.py, example_script_parallel_isomut.py, example scripts  to run wrappaers

Testing and demonstration:
- test.ipynb jupyter notebook showing a test run for the standalone C program
- test_isomut_parallel.ipynb jupyter notebook showing a parallel test run
- test_isomut_parallel.ipynb jupyter notebook showing a parallel test run, with the BAQ based filtering

---

### Dependencies

- samtools needs to be installed


---

### Understanding the output

- The output is a tab separated file, each line is a detected mutation. The columns are the following:
	- #sample_idx: the index of the sample in the input list of samples given which contain the mutation. 
	- chr: the chromosome/contig name of the mutation
	- pos: the 1-based position of the mutation for SNV, and the previous base before INS/DEL
	- type: the mutation type: values can be: SNV/INS/DEL
	- score: Quality score for the mutation
	- ref: reference base/sequence of the mutation (- for insertions)
	- mut: mutated base/sequence of the mutation (- for deletions)
	- cov: coverage of the mutated position (after base quality filtering!)
	- mut_freq: frequency of reads containing the mutation
	- cleanliness: lowest reference nucleotide frequency among the other samples

---

### Understading the mutation quality score

- We generalized the scoring used by VarScan , the Fisher's exact p-value for multiple samples. The score is calculated taking the counts in the mutated samples, and the counts from the noisiest sample among the other samples. This value is realted to the probability of the mutation being a false positive, but it has not direct phisical interpretation. Therefore, to avoid confusion we use the negative 10 based logarithm of the Fisher's exact p-value calculated. This means, that the higher the score, the better the mutation is.
-  We have found that the score can be used to significantly improve the IsoMuts hard filters, and it allows for a very precise tuning between sensitivity and sepcificity of IsoMut.

---

### Workflow:

- 1, IsoMut calls samtools mpileup wih -B flag for all samples together
- 2, It processes every line and selects the ones which meet the following criteria
	- reliability: coverage is higher than a given threshold
	- enough evidence: mutated base frequency is higher than a given threshold
	- low noise: other samples minimum frequency is higher than a given threshold
	- no gaps are found closer than a given radius 
- 3, For efficient proximal indel induced false positive SNV filtering, IsoMut repeats the the above for the suspected SNVs without the -B flag.


---

### Authors

Main code was written by me (Dezso Ribli), but the research process is done in the collaboration of  
- SE enzimology insitute: David Szuts, Janos Molnar
- ELTE department of complex physics: Istvan Csabai, Orsi Pipek, Dezso Ribli


Fisher's exact test from this github repo, ( thanks Christopher Chang ) :
- https://github.com/chrchang/stats

---
