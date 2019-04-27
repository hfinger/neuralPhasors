clear paramsAll;

clear params;

params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 5000; % measured 3400 MB
params.Gridjob.wc_host = '!(*saturn*|*c01grid*|*ramsauer*|*kuma*|*taygete*|*helike*|*leda*|*jupiter*)';
params.Gridjob.jobname = 'Connectome';
params.Gridjob.continue = true;
params.Gridjob.initRandStreamWithJobid = true;
params.Gridjob.runOnlyJobIds = [
734;
735;
736;
739;
740;
743;
745;
746;
747;
748;
751;
753;
757;
758;
759;
760;
762;
763;
764;
766;
768;
770;
771;
772;
773;
775;
779;
781;
782;
784;
785;
787;
791;
794;
795;
796;
797;
798;
799;
800;
803;
804;
806;
808;
809;
811;
813;
816;
819;
820;
821;
822;
824;
826;
828;
830;
831;
833;
834;
836;
839;
841;
843;
844;
845;
848;
849;
850;
852;
854;
855;
856;
859;
862;
863;
866;
868;
871;
873;
874;
875;
876;
878;
880;
881;
882;
883;
886;
890;
891;
892;
894;
895;
898;
900;
901;
903;
904;
907;
908;
909;
912;
914;
915;
918;
920;
921;
922;
923;
924;
926;
928;
932;
933;
935;
937;
938;
941;
942;
945;
946;
947;
948;
950;
951;
952;
955;
956;
958;
959;
960;
961;
965;
969;
972;
974;
975;
976;
977;
978;
980;
981;
982;
983;
984;
987;
988;
990;
991;
993;
996;
997;
1000;
1001;
1003;
1004;
1005;
1006;
1009;
1011;
1012;
1015;
1017;
1018;
1020;
1023;
1024;
1026;
1029;
1030;
1031;
1032;
1033;
1036;
1037;
1039;
1040;
1041;
1042;
1044;
1048;
1050;
1052;
1054;
1055;
1058;
1059;
1060;
1061;
1063;
1065;
1066;
1067;
1069;
1070;
1071;
1075;
1077;
1079;
1080;
1081;
1083;
1087;
1088;
1089;
1090;
1092;
1095;
1096;
1098;
1100;
1102;
1103;
1106;
1108;
1109;
1110;
1113;
1115;
1116;
1118;
1119;
1121;
1122;
1124;
1127;
1128;
1130;
1133;
1134;
1136;
1139;
1141;
1142;
1143;
1144;
1145;
1147;
1149;
1152;
1154;    
];
params.Gridjob.combParallel = false;
params.Gridjob.walltime = '00:29:00';
params.Gridjob.requiredThreads = '3';

params.JansenRitConnectomePaper.p = 1; %defines what kind of p norm to use for normalization of structural connectivity
params.JansenRitConnectomePaper.k = 20; %global connection strength scaling
params.JansenRitConnectomePaper.v = 3; % velocity [m/s]
params.JansenRitConnectomePaper.tMax = 605; %max simulation time [seconds]
params.JansenRitConnectomePaper.dt = 0.0005; % simulation step size [seconds]
params.JansenRitConnectomePaper.sampling = 2; % sampling every x steps
params.JansenRitConnectomePaper.noiseVar = 0; % variance of noise used to drive neural masses (22 was mentioned in msc thesis)
params.JansenRitConnectomePaper.noiseMu = [120, 320]; %num2cell(round(20:20:140)); %220; % mean of noise used to drive neural masses (220 was mentioned in msc thesis)
params.JansenRitConnectomePaper.runningAvg = false; % if true, substract running average from input to neural masses
params.JansenRitConnectomePaper.netInp = [1,0,0]; % scales the input to [pyramidal, excitatory, inhibitory] neurons from connectome
params.JansenRitConnectomePaper.subInp = [1,0,0]; % scales the input to [pyramidal, excitatory, inhibitory] neurons from connectome
params.JansenRitConnectomePaper.initSampRem = 124; %initial interval to remove [seconds]
params.JansenRitConnectomePaper.verbose = false; % if we want to print time steps to console

c_tmp = 135;
params.JansenRitConnectomePaper.cs = [c_tmp, c_tmp*0.8, c_tmp*0.25, c_tmp*0.25, 0]; % connectivity strength (can sometimes be interpreted as average synaptic contacts)

params.JansenRitConnectomePaper.fTarget = 'drivFreq'; % [Hz]

params.JansenRitConnectomePaper.FCMeasure = {{'Coherence'}}; % FC measure to use
params.JansenRitConnectomePaper.fullFC = true; %whether to store full coherence matrix or only coherence of driven region
params.JansenRitConnectomePaper.corrSimFC = false; % if true, compare simulated coherence to empirical coherence
params.JansenRitConnectomePaper.nWindows = 1; %number of timewindows over which to calculate mean coherence
params.JansenRitConnectomePaper.nBins = 32; % number of bins to use to discretize phase signal (only for MI and SE)
params.JansenRitConnectomePaper.storeY = false; %if true, store raw signal on simResults
params.JansenRitConnectomePaper.filterSig = true; % if true, raw PSPs will be bandpass-filtered before coherence is calculated

params.JansenRitConnectomePaper.u0 = 6e-3; % membrane voltage for which 50 % of maximum mean firing rate is observed [V].. was 0 in master thesis
params.JansenRitConnectomePaper.He = 3.25e-3; % Average synaptic gain for excitatory synapses [V]
params.JansenRitConnectomePaper.Hi = 22e-3; % Average synaptic gain for inhibitory synapses [V]
params.JansenRitConnectomePaper.Te = 10e-3; % Average time constant for excitatory signal transfer (synaptic delays,..) [s]
params.JansenRitConnectomePaper.Ti = 20e-3; % Average time constant for inhibitory signal transfer (synaptic delays,..) [s]
params.JansenRitConnectomePaper.e0 = 2.5; % determines maximum mean firing rate [1/s]
params.JansenRitConnectomePaper.r = 560; % steepness of the sigmoid transfer function from mean psp to mean firing rate [1/V]

params.JansenRitConnectomePaper.calcFC_nwin1 = false;
params.JansenRitConnectomePaper.subtract_S0 = false;
params.JansenRitConnectomePaper.use_moran = false;
params.JansenRitConnectomePaper.use_out_psp = false;
params.JansenRitConnectomePaper.use_sigm_as_out = false;
params.JansenRitConnectomePaper.use_inpP_as_out = true;
params.JansenRitConnectomePaper.use_sigm_y0_as_out = false;
params.JansenRitConnectomePaper.calcCohWithDriver = true;
params.JansenRitConnectomePaper.saveSpectrum = false;

% connectivity matrix
num_nodes = 4;
C = zeros(num_nodes,num_nodes);
C(1,2) = 0.1;
C(2,1) = 0.1;
C(1,3) = 0.1;
C(3,1) = 0.1;
C(2,4) = 0.1;
C(4,2) = 0.1;
C(3,4) = 0.1;
C(4,3) = 0.1;

% concat many to do multiple simulations simulatanously:
num_repeats = 16;
C_full = zeros(num_repeats*num_nodes, num_repeats*num_nodes);
for k=1:num_repeats
    start_idx=num_nodes*(k-1)+1;
    C_full(start_idx:start_idx+num_nodes-1, start_idx:start_idx+num_nodes-1) = C;
end
C = C_full;

% delay matrix
d213 = 0:10:200; % fix ratio
d243 = 100; % fix total length
d24 = 0:10:100;

Delays = cell(length(d213),length(d24));
for d213_idx=1:length(d213)
    for d24_idx=1:length(d24)
        D = zeros(num_nodes, num_nodes);
        D(1,2) = 0.25 * d213(d213_idx);
        D(2,1) = 0.25 * d213(d213_idx);
        D(1,3) = 0.75 * d213(d213_idx);
        D(3,1) = 0.75 * d213(d213_idx);
        D(2,4) = d24(d24_idx);
        D(4,2) = d24(d24_idx);
        D(3,4) = d243 - d24(d24_idx);
        D(4,3) = d243 - d24(d24_idx);
        D_full = zeros(num_repeats*num_nodes, num_repeats*num_nodes);
        for k=1:num_repeats
            start_idx=num_nodes*(k-1)+1;
            D_full(start_idx:start_idx+num_nodes-1, start_idx:start_idx+num_nodes-1) = D;
        end
        Delays{d213_idx,d24_idx} = D_full;
    end
end
Delays = Delays(:)';

params.JansenRitConnectomePaper.C = C; % connectivity matrix
params.JansenRitConnectomePaper.D = Delays; % distance matrix

drivPos = [2, 3];
drivPosPerSPO = cell(1,num_repeats);
for k=1:num_repeats
    start_idx=num_nodes*(k-1)+1;
    drivPosPerSPO{k} = drivPos+start_idx-1;
end
drivPos = cell2mat(drivPosPerSPO);

%% set driver params
params.JansenRitConnectomePaper.drivPos = drivPos; % indices of network nodes to be driven
params.JansenRitConnectomePaper.drivScale = num2cell(ones(1,5)*0.25); % 0.25, 0.5 driver strength [mV]
phase_offsets = [0:0.0625:(1-0.0625)]*2*pi;
phase_offsets = cat(1, phase_offsets, zeros(size(phase_offsets)));
phase_offsets = phase_offsets(:);
params.JansenRitConnectomePaper.drivPO = phase_offsets'; % phase offset of drivers
params.JansenRitConnectomePaper.drivFreq = 11; % 10 and 11 frequency of driver [Hz]
params.JansenRitConnectomePaper.drivDur = 605; % duration of driver [s]

paramsAll{1} = params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);

%data = dataPaths();
%job = JansenRitConnectomePaper(fullfile(data.workdir, 'Holger/2018_JR/ToyModels/WithDriver/4NodesNetwork/VarRatioAndVarTot/sim4_grid_like2_repeat5', 'temp_Connectome/jobDesc.mat'));
%job.finishJobs()