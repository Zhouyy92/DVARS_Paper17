#This script make an appropriate bash files and submit them to the cluster for given time series and variance het. 
#
# SA, UoW, 2017
# srafyouni@gmail.com
#
####REFERENCES
#
#   Afyouni S. & Nichols T.E., Insights and inference for DVARS, 2017
#   http://www.biorxiv.org/content/early/2017/04/06/125021.1

# WallClock=${3}
# T_Idx=${1}
# SD_Idx=${2}
# nRLz=${4}
# Mem=${5}

FileName="tmp_Sim_DVRARS_T${1}_SD${2}.sh"
cat > $FileName << EOF
#!/bin/bash
#$ -cwd
#$ -j y
#$ -o /data/home/mpx143/DVARS/Sim/logs
#$ -l h_rt=${3}:00:00
#$ -l h_vmem=${5}G
#$ -N DVARS_Sim_${1}_${2}
#$ -r y
# #$ -t 1-{4}

cd /data/home/mpx143/DVARS/Sims/

#echo "submit:${hh},${mm} -- From ${Site} ${Atlas} ${GSRStat} -- Spec ${ToNum} ${mat_filename}"

module load matlab
#run a script file
matlab -nodisplay -singleCompThread -nosplash -nojvm -nodesktop -r "T_idx=${1};SD_idx=${2};nRlz=${4};Sim_DVARS_Inference"
EOF

qsub $FileName
