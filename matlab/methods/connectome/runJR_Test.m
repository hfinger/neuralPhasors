clear all;

x = zeros(1,13);
dt = 0.0001;
Tmax = 20;

% Default Jansen Rit parameters
He = 3.25e-3; % A;       Average synaptic gain for excitatory synapses [V]
Hi = 22e-3; % B;         Average synaptic gain for inhibitory synapses [V]
Te = 10e-3; % 1/a;       Average time constant for excitatory signal transfer (synaptic delays,..) [s]
Ti = 20e-3; % 1/b;       Average time constant for inhibitory signal transfer (synaptic delays,..) [s]
c = 135;
%C = [c,c*0.8,c*0.25,c*0.25]; % connectivity strength (can sometimes be interpreted as average synaptic contacts)
C = [128,128,64,64,16];
Cpe = C(1); % Connection from pyramidal cells to excitatory interneurons
Cpi = C(3); % Connection from pyramidal cells to inhibitory interneurons
Cep = C(2); % Connection from excitatory interneurons to pyramidal cells
Cip = C(4); % Connection from inhibitory interneurons to pyramidal cells
Cii = C(5); % recurrent intrinsic connection from inhibitory interneurons to themselfes
u0 = 6e-3; % membrane voltage for which 50 % of maximum mean firing rate is observed [V]
e0 = 2.5; % determines maximum mean firing rate [1/s]
r = 560; % steepness of the sigmoid transfer function from mean psp to mean firing rate [1/V]
Ka = 1.9531; % adaptation rate constant [1/s]

% alternative parameters leading to gamma-oscillations (23 Hz)
He = 15e-3;
Hi = 22e-3;
Te = 3.5e-3;
Ti = 10e-3;


% pre-calculate JR constants
Ke = He/Te;
Ki = Hi/Ti;
de1 = -2/Te;
de2 = -1/Te^2;
di1 = -2/Ti;
di2 = -1/Ti^2;
S0 = WTP(0, e0, u0, r);

% driver parameters
drivStrength = 0.1;
drivFreq = 30;
drivEnd = Tmax/2;

%% Run simulation
numIters = ceil(Tmax/dt);
x_all = zeros(numIters,13);
driv_all = zeros(numIters,1);
a = 0;
tic()
t = 0;
for i=1:numIters
  t= t + dt;
  
  x_all(i,:) = x;
  Pe = 220 + 22 * randn(1);
  Pp = drivStrength * sin(t*2*pi*drivFreq);
  driv_all(i) = Pp;
  if t > drivEnd
      Pp = 0;
  end
  
  [ dx, da ] = JansenRitUpdate( x, a, Pp, Pe, Ke, Ki, Ka, Cpe, Cpi, Cep, Cip, Cii, de1, de2, di1, di2, e0, u0, r, S0 );
  x = x + dx * dt;
  a = a + da * dt;
  
end
toc()

%% plot stuff

% spectrogram
Fs = 1/dt;
Hs = spectrum.periodogram;
psd(Hs,x_all(:,1),'Fs',Fs)
hold on
xlim([0.01,0.07])
ylim([-90,0])
xvec = [0.01:0.0001:0.07];
y = zeros(size(xvec));
line([0.0295,0.0295],[-90,0],'Color','r','LineStyle','--')
line([0.01,0.07],[-32.1,-32.1],'Color','r','LineStyle','--')
line([0.05905,0.05905],[-90,0],'Color','g','LineStyle','-.')
line([0.01,0.07],[-54.5,-54.5],'Color','g','LineStyle','-.')
hold off

% phase difference between driver and pyramidal cells
env.t_rm = 1; % time to start bandpassing the signal at [s]
env.sigBandpass(1).Fst1 = 25.5; % end stop band [Hz] 
env.sigBandpass(1).Fp1 = 28; % start pass band [Hz]
env.sigBandpass(1).Fp2 = 32; % [Hz] end pass band
env.sigBandpass(1).Fst2 = 34.5; % [Hz] start stop band
env.sigBandpass(1).Ast1 = 30; % frequency attenuation in first stopband
env.sigBandpass(1).Ap = 1; % passband ripples
env.sigBandpass(1).Ast2 = 30; % frequency attenuation in second stopband

bp = env.sigBandpass(1);
dBandPass = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',bp.Fst1,bp.Fp1,bp.Fp2,bp.Fst2,bp.Ast1,bp.Ap,bp.Ast2,Fs);
HdBandPass = design(dBandPass,'butter');

data = vertcat(x_all(:,1)',driv_all');
sourceBP = zeros(size(data));
for k=1:size(data,1)
    sourceBP(k,:) = filter(HdBandPass,squeeze(data(k,:)),2);
end
sigHilbert = zeros(size(sourceBP));
for n=1:size(sigHilbert,1)
    sigHilbert(n,:) = hilbert(squeeze(sourceBP(n,:)));
end
phaseBP = angle(sigHilbert);
phaseDiff = mod(bsxfun(@minus,phaseBP(1,:),phaseBP(2,:)), 2*pi);

figure(2)
xvec = 0:size(phaseDiff, 2) - 1;
xvec = (xvec - drivEnd/dt)*dt;
plot(xvec,phaseDiff)
ylabel('phase')
xlabel('Seconds since driver was turned off')
title('Phase Difference between Driver and NMM')

% PSPs of the 3 different cell populations
% figure(3)
% plot(0:dt:(Tmax-dt),x_all(:,1:3))
% figure(2)
% 
% legend('Pyr','Exc','Inh')