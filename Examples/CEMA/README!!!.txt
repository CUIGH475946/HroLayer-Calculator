# First, cd to the namd directory
# Then, run the .conf files in order

./namd3 +p28 ../CEMA/min.conf > ../CEMA/min.log
./namd3 +p28 ../CEMA/equ.conf > ../CEMA/equ.log
./namd3 +p28 ../CEMA/run.conf > ../CEMA/run.log
./namd3 +p28 ../CEMA/mrt.conf > ../CEMA/mrt.log

# ./namd3 stands for running the program, in some older versions it is ./namd2
# +p28 stands for 28 threads calling the CPU, on some hosts it needs to be adjusted itself
# ../CEMA/mrt.conf stands for Document Catalog

# !!!
# Run in strict order of min, equ, run, mrt !!!

# Commands for extracting data

evaluate project.psf run.dcd -ref "type OG2D2" -first 10000 -ODF 1 -HB 1 -MSD 1
evaluate project.psf run.dcd -ref "type OG2D2" -first 10000 -last 10999 -RDF 1
evaluate project.psf mrt.dcd -ref "type OG2D2" -first 20000 -MRT 1
