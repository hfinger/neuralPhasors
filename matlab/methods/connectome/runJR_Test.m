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
drivStrength = 100;
drivFreq = 30;

%%
numIters = ceil(Tmax/dt);
x_all = zeros(numIters,13);
a = 0;
tic()
t = 0;
for i=1:numIters
  t= t + dt;
  
  x_all(i,:) = x;
  Pe = 220 + 22 * randn(1);
  Pp = drivStrength * sin(t*2*pi*drivFreq);
  
  [ dx, da ] = JansenRitUpdate( x, a, Pp, Pe, Ke, Ki, Ka, Cpe, Cpi, Cep, Cip, Cii, de1, de2, di1, di2, e0, u0, r, S0 );
  x = x + dx * dt;
  a = a + da * dt;
  
end
toc()
%%
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

% figure(2)
% plot(0:dt:(Tmax-dt),x_all(:,1:3))
% figure(2)
% 
% legend('Pyr','Exc','Inh')