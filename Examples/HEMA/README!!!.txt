# First, cd to the namd directory
# Then, run the .conf files in order

./namd3 +p28 ../HEMA/min.conf > ../HEMA/min.log
./namd3 +p28 ../HEMA/equ.conf > ../HEMA/equ.log
./namd3 +p28 ../HEMA/run.conf > ../HEMA/run.log
./namd3 +p28 ../HEMA/mrt.conf > ../HEMA/mrt.log

# ./namd3 stands for running the program, in some older versions it is ./namd2
# +p28 stands for 28 threads calling the CPU, on some hosts it needs to be adjusted itself
# ../HEMA/mrt.conf stands for Document Catalog

# !!!
# Run in strict order of min, equ, run, mrt !!!

# Commands for extracting data

evaluate project.psf run.dcd -ref "type OG311" -first 10000 -ODF 1 -HB 1 -MSD 1
evaluate project.psf run.dcd -ref "type OG311" -first 10000 -last 10999 -RDF 1
evaluate project.psf mrt.dcd -ref "type OG311" -first 20000 -MRT 1
