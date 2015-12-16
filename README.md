# IsoMut: a robust method for calling unique mutations from multiple isogenic samples

This is a C implementation of the mutation calling method, and python wappers fro parallel exectution, and the BAQ based filtering of SNVs

---

### For the impatient:
- modify the pathinames for your data in example_script_isomut_w_pp.py, and run the script
	- check out test_isomut_w_pp.ipynb, for necessities, and details

--

### Files:

For the standalone C application:
- isomut_lib.h unique_mutatiob_lib.c  the library
- isomut.c the application

For the parallel wrapper:
- isomut_wrappers.py the python functions for parallel execution, and BAQ filtering
- example_script_isomut_w_pp.py, example_script_parallel_isomut.py, example scripts  to run wrappaers

Testing and demonstration:
- test.ipynb jupyter notebook showing a test run for the standalone C program
- test_isomut_parallel.ipynb jupyter notebook showing a parallel test run
- test_isomut_parallel.ipynb jupyter notebook showing a parallel test run, with the BAQ based filtering

---

Main code was written by me (Dezso Ribli), but the research process is done in the collaboration of  
- SE enzimology insitute: David Szuts, Janos Molnar
- ELTE department of complex physics: Istvan Csabai, Orsi Pipek, Dezso Ribli


Fisher's exact test from this github repo, ( thanks Christopher Chang ) :
- https://github.com/chrchang/stats

---

