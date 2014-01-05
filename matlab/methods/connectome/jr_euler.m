function jr_euler()
    clc; clear; close all;
    
    global W T1 T2 I K ke ki c1 c2 c3 c4 iu istd he taue hi taui v0 e0 r p  T ...
        cpp cpe cpi cep cee cei cip cie cii h
    
    % Connectivity strength
    c = 135; c1 = c; c2 = 0.8*c; c3 = 0.25*c; c4 = 0.25*c;   
    cpp = 0; cpe = c1; cpi = c3;
    cep = c2; cee = 0; cei = 0;
    cip = c4; cie = 0; cii = 0;
    
    tstart = 0;
    tend = 2.5;
    
    % Parameters of the nonlinear sigmoid function
    v0 = 6E-3; %[v0]=V   
    e0 = 2.5;  %[e0]=1/s 
    r = 560;   %[r]=1/V  
    he = 3.25E-3; %[HeJR]=V 
    hi = -22E-3; %[HiJR]=V 
    taue = 10e-3;
    taui = 20e-3;
    Fs = 1000;
    T = tstart:1/Fs:tend;
    % input
    iu = 220; % input expectation value %220 JR'95
    istd = 22; % input standard deviation %100/sqrt(3) JR'95
    
    p = iu + istd * rand(1,length(T));
    h = 220;
    
    % Connectivity matrix
    ke = he/taue; ki = hi/taui;
    K = [ke; ke; ki];
    W = [0  cep  cip; 
         cpe  0  0; 
         cpi  0  0]; 
    T1 =  [-2/taue 0 0; 0 -2/taue 0; 0 0 -2/taui];
    T2 =  [-1/taue^2 0 0; 0 -1/taue^2 0; 0 0 -1/taui^2];
    I = [0; ke*h; 0];
     
    ystart = zeros(6,1);
    
    % Euler solver
    stepsize = 1e-6;
    y = ystart; ynow = ystart();
    for t = tstart:stepsize:tend
      if mod(t,0.1)==0
        disp(t);
      end
        dy = odesys(t,ynow);
        ynext = ynow + stepsize*dy;
        ynow = ynext;
        y = [y ynow];
    end

    %% Plot
%     t = tstart:stepsize:tend;
%     figure
%     for i = 1:min(size(y))
%        hold on;subplot 311; plot(t,y(i,1:length(t)),'Color',rand(1,3)); 
%     end
%     hold off;
%     tt = t(ceil(end/2):10:end); 
%     u = y(1,:);
%     u = u(ceil(end/2):10:end); u = u - mean(u);
%     subplot 312; plot(tt,u);axis tight
%     %%spectrum
%     Fs = 2/stepsize;
%     T = 1/Fs;                     % Sample time
%     L = 1000;                     % Length of signal
%     NFFT = 2^nextpow2(L);         % Next power of 2 from length of y
%     U = fft(u,NFFT)/L;
%     f = Fs/2*linspace(0,1,NFFT/2);
%     subplot 313; plot(f,2*abs(U(1:NFFT/2))); xlim([0 50]);     
    
    save(['euler_stepsize=' num2str(stepsize,'%.e') '.mat']);
   
end

function dy = odesys(t,y)
    global W T1 T2 I K
    
%     disp(t);       

    S = sigm(y(1:3));
    S = K.*S;
    
    dy(1:3) = y(4:6);
    dy(4:6) = W*S + T1*y(4:6) + T2*y(1:3) + I;

    dy = dy(:);

end

function Q = sigm(y)
    global v0 e0 r
    Q = 2*e0 ./ (1 + exp(r*(v0-(y))));
end