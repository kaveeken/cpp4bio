#!/bin/bash
#SBATCH --job-name=cpp_proj
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2000
#SBATCH --array=1-9

module load GCC/7.3.0-2.30
pwd


# take next configuration filename
fileName=$(ls run* | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
# make a directory to run simulation in and go there
mkdir ${fileName}dir
cd ${fileName}dir
pwd
# convert configuration file to program-readable format 
# by deleting "commented" lines and put it here
sed '/#/d' ../$fileName > config.txt
# copy the program to the simulation directory and run it
cp ../executable .
time ./executable

# normalize & visualize 
# take arguments from parameters.txt
#maxS=$(cat parameters.txt | tail -n 3 | head -n 1 | awk '{print $3}')
#maxR=$(cat parameters.txt | tail -n 2 | head -n 1 | awk '{print $3}')
maxX=$(cat parameters.txt | tail -n 1 | awk '{print $3}')
# put it into a wolfram script because arguments dont work(?)
cp ../maakplaatjes.wls .
sed -i 's/val/$maxS/g' maakplaatjes.wls
# run the script to normalize and make pngs
#wolfram -script maakplaatjes.wls

# add colour
# lists of images
#imageR=$(ls -1 *R.png)
#imageS=$(ls -1 *S.png)
# convert
#for file in $imageR; do
#	convert $file +level-colors blue,white blu${file}
#done

#for file in $imageS; do
#	convert $file +level-colors red,white red${file}
#done

# go back to parent directory
cd ../

