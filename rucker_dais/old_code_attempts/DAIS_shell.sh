#PBS -l nodes=4:ppn=16
#PBS -l walltime=144:00:00
#PBS -j oe
#PBS -o dais_error.txt
#PBS -m abe
#PBS -M klr324@psu.edu

cd $PBS_O_WORKDIR

Rscript DAIScali_hetero_model_iid_mcmc.R
