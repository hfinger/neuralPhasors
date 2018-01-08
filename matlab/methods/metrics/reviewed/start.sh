#!/bin/bash
#$ -t 2
#$ -N opt
#$ -l mem=2000M
#$ -pe matlab 4
#$ -m n
#$ -l h=!ramsauer.ikw.uni-osnabrueck.de
#$ -wd /net/store/nbp/projects/phasesim/src_holger/matlab/methods/metrics/reviewed
echo $(hostname)
matlab -nodisplay -r "addpath('/net/store/nbp/projects/phasesim/src_holger/matlab'); addScriptPaths(); taskId=$SGE_TASK_ID; initOptimizeSC"
exit 0