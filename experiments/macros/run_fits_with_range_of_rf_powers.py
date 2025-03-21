from instruments import   station, qdac,  Triton, zurich,exp

drives=[1.5e-3,2e-3,3e-3,4e-3,5e-3]

for drive in drives:
    zurich.output1_amp1(drive) 
   # %run experiments/cs_linesweeps/chargesensing2d_qdac_zurich_fitv9.py