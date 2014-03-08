Edit the paths in "include/dataPaths.m" to match your system setup...




%change your directory to our project folder on the IKW server:
cd /net/store/nbp/phasesim/

%checkout the source code into a new folder where you will edit your code. 
%Replace yourName with your name, so that everyone has their own copy of the source code:
git clone file:////net/store/nbp/phasesim/remote/repo.git/ src_yourName

%start matlab (you might also want to create a shortcut on your desktop):
matlab

%In matlab change the directory to your checked out source code folder:
cd net/store/nbp/phasesim/src_yourName/matlab

%To add all required folders to your matlab path run:
addScriptPaths();

%Change to the directory with the parameters for simulations:
cd parameters

%create a new subfolder where you will store the parameters of your simulations:
mkdir yourname

%copy an example parameter file to the directory you just created:
copyfile('20140501_Connectome/20150217_ExampleSimulation','yourname/20150217_ExampleSimulation')

%change to that directory:
cd yourname/20150217_ExampleSimulation

%run the simulation:
script_SAR

%edit the parameters in the file... maybe set a breakpoint in the last line of that file where the job is started, so that you can step through the code:
edit script_SAR