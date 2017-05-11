clear all;
close all;

x = zeros(1,13);
dt = 0.0001;
Tmax = 60;

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
drivStrengths = {0, 0.1, 0.5}; % search in range 0.01 up to 0.1
drivFreqs = {23, 28};%{27, 28, 29};

counter = 1;
for k=1:length(drivStrengths)
  for l=1:length(drivFreqs)
  
    drivStrength = drivStrengths{k};
    drivFreq = drivFreqs{l};

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
    peakFreq = 29.55;
    figure(counter);
    clf;
    [pxx,f] = periodogram(x_all(:,1),[],200000,1/dt);
    plot(f(1:2000),10*log10(pxx(1:2000)),'k')
    xlabel('Frequency (Hz)')
    ylabel('Magnitude (dB)')
    xlim([0 40])
    ylim([-120 -20])
    %title(['drivStrength=' num2str(drivStrength) ' drivFreq=' num2str(drivFreq)])
    hold on;
    hDrivFreq = line([drivFreq,drivFreq],[-200,-90],'Color','r','LineStyle','--','linewidth',2);
    hPeakFreq = line([peakFreq,peakFreq],[-200,-90],'Color','b','LineStyle','--','linewidth',2);
    set(gca,'fontsize',18)
    set(gca,'linewidth',2)
    if drivStrength>0
      legend(hDrivFreq,{['driver at ' num2str(drivFreq) ' Hz']},'Location','northwest')
    else
      legend(hPeakFreq,{['peak at ' num2str(peakFreq) ' Hz']},'Location','northwest')
    end
    saveas(gcf,['runJR_Test_holger_drivStrength' num2str(drivStrength) '_drivFreq' num2str(drivFreq) '.pdf'])
    
    
    counter = counter + 1;

  end
end

% figure(2)
% plot(0:dt:(Tmax-dt),x_all(:,1:3))
% figure(2)
% 
% legend('Pyr','Exc','Inh')