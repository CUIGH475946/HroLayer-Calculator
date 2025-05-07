# First, cd to the namd directory
# Then, run the .conf files in order

./namd3 +p28 ../Bulkwater/min.conf > ../Bulkwater/min.log
./namd3 +p28 ../Bulkwater/equ.conf > ../Bulkwater/equ.log
./namd3 +p28 ../Bulkwater/run.conf > ../Bulkwater/run.log
./namd3 +p28 ../Bulkwater/mrt.conf > ../Bulkwater/mrt.log

# ./namd3 stands for running the program, in some older versions it is ./namd2
# +p28 stands for 28 threads calling the CPU, on some hosts it needs to be adjusted itself
# ../Bulkwater/mrt.conf stands for Document Catalog

# !!!
# Run in strict order of min, equ, run, mrt !!!

# Commands for extracting data

evaluate project.psf run.dcd -ref "index 22005" -cal "type OT and (not index 22005)" -first 10000 -ODF 1 -HB 1 -MSD 1
evaluate project.psf run.dcd -ref "index 22005" -cal "type OT and (not index 22005)" -first 10000 -last 10999 -RDF 1
evaluate project.psf mrt.dcd -ref "index 22005" -cal "type OT and (not index 22005)" -first 20000 -MRT 1
