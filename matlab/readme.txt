Edit the paths in "include/dataPaths.m" to match your system setup...


start a shell/terminal

%change your directory to our project folder on the IKW server:
cd /net/store/nbp/phasesim/

%checkout the source code into a new folder where you will edit your code. 
%Replace YOURNAME with your name, so that everyone has their own copy of the source code:
git clone file:////net/store/nbp/phasesim/remote/repo.git/ src_YOURNAME

%start matlab (you might also want to create a shortcut on your desktop):
matlab

%In matlab change the directory to your checked out source code folder (replace YOURNAME):
cd /net/store/nbp/phasesim/src_YOURNAME/matlab

%To add all required folders to your matlab path run:
addScriptPaths();

%Change to the directory with the parameters for simulations:
cd parameters

%create a new subfolder where you will store the parameters of your simulations (replace YOURNAME):
mkdir YOURNAME

%copy an example parameter file to the directory you just created (replace YOURNAME):
copyfile('20140501_Connectome/20150217_ExampleSimulation','YOURNAME/20150217_ExampleSimulation')

%change to that directory (replace YOURNAME):
cd YOURNAME/20150217_ExampleSimulation

%open the file which contains all paramters...
edit script_SAR

%run the simulation:
script_SAR

%the results of that simulation can then be found in the directory: 
/net/store/nbp/phasesim/workdir/YOURNAME/20150217_ExampleSimulation/CompareWithEEG

