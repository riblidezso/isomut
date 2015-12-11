# Robust method for calling unique mutations from multiple isogenic samples

This is a C implementation of the mutation calling method, and a python wapper fro parallel exectution.

---

### Files:

For the standalone C application:
- unique_mutation_lib.h unique_mutatiob_lib.c  the library
- unique_mutation_app.c the application

For the parallel wrapper:
- isomut_parallel.py the python functions for parallel execution
- isomut_sample_run.py an example script

Testing and demonstration:
- unique_mutation_C_lib.ipynb jupyter notebook showing a test run
- test_Isomut_parallel.ipynb jupyter notebook showing a parallel test run

---

Main code was written by me (Dezso Ribli), but the research process is done in the collaboration of  
- SE enzimology insitute: David Szuts, Janos Molnar
- ELTE department of complex physics: Istvan Csabai, Orsi Pipek, Dezso Ribli


Fisher's exact test from this github repo, ( thanks Christopher Chang ) :
- https://github.com/chrchang/stats

---

