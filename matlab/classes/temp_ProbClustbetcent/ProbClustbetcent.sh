#!/bin/bash
#$ -t 1:660
#$ -N ProbClustbetcent
#$ -l mem=1000M
#$ -pe matlab 3
#$ -m n
#$ -l h=!(*ramsauer*)
#$ -q nbp.q
#$ -wd /net/store/nbp/projects/phasesim/src_Arushi/matlab/classes/temp_ProbClustbetcent
echo $(hostname)
matlab -nodisplay -r "addpath('/net/store/nbp/projects/phasesim/src_Arushi/matlab'); addScriptPaths(); Gridjob.startJobid('/net/store/nbp/projects/phasesim/src_Arushi/matlab/classes/temp_ProbClustbetcent/jobDesc.mat',$SGE_TASK_ID);"
exit 0