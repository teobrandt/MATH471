#!/bin/bash                                                                    
#SBATCH -J homework7           # job name                                          
#SBATCH -o homework7.o%j       # output and error file name (%j expands to jobID)  
#SBATCH -n 4              # total number of mpi tasks requested                
#SBATCH -p development     # queue (partition) -- normal, development, etc.    
#SBATCH -t 01:30:00        # run time (hh:mm:ss) - 1.5 hours                   
#SBATCH --mail-user=username@tacc.utexas.edu                                   
#SBATCH --mail-type=begin  # email me when the job starts                      
#SBATCH --mail-type=end    # email me when the job finishes                    
ibrun ./homework7.x              # run the MPI executable named a.out