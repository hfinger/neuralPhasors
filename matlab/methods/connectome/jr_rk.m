function jr_rk()
%%
%% Solving JR ODE system using Runge-Kutta 
%% algorithm with adative stepsize
%% @MNT MPI Leipzig 2011
%%
    % Clean
    clc;
    clear;
    close all;
    
    % Global variables
    global W T1 T2 I v0 e0 r p  T h p
    
    % Connectivity strength from JR95
    c = 135; c1 = c; c2 = 0.8*c; 
    c3 = 0.25*c; c4 = 0.25*c;   
    
    cpp = 0; cpe = c1; cpi = c3;
    cep = c2; cee = 0; cei = 0;
    cip = c4; cie = 0; cii = 0;
    
    % Integration time (s)
    tstart = 0;  tend = 5; stepsize = 0.01;
    
    % Parameters of the nonlinear sigmoid function
    v0 = 6E-3;    %[v0]=V   
    e0 = 2.5;     %[e0]=1/s 
    r = 560;      %[r]=1/V  
    
    % Dendrite parameter
    he = 3.25E-3; %[HeJR]=V 
    hi = -22E-3;  %[HiJR]=V 
    taue = 10e-3; % s
    taui = 20e-3; % s
    
    % input
    iu = 220; % input expectation value %220 JR'95
    istd = 22; % input standard deviation %100/sqrt(3) JR'95    
    p = iu + istd * rand(1,length(T));
    h = 220;
    
    % Connectivity matrix
    kp = he/taue; ke = he/taue; ki = hi/taui;
    W = [0  ke*cep  ki*cip; ke*cpe  0  0; ke*cpi  0  0]; 
    
    T1 =  [-2/taue 0 0; 0 -2/taue 0; 0 0 -2/taui];
    T2 =  [-1/taue^2 0 0; 0 -1/taue^2 0; 0 0 -1/taui^2];
    I = [0; ke*h; 0];
    
    % Initial condition
    ystart = zeros(6,1);
    
    % Runge-Kutta Solver
    options = odeset('RelTol',1e-9,'AbsTol',1e-9); 
    [t,y] = ode45(@odesys,tstart:stepsize:tend,ystart,options);   
    
    % Plot
    figure
    for i = 1:min(size(y))
        hold on; plot(t,y(:,i),'Color',rand(1,3)); 
    end
    hold off;
end
  
%% ODE system
function dy = odesys(t,y)
    global W T1 T2 I K    
    disp(t);       
    S = sigm(y(1:3));
    dy(1:3) = y(4:6);
    dy(4:6) = W*S + T1*y(4:6) + T2*y(1:3) + I;
    dy = dy(:);
end

%% Sigmoid function
function Q = sigm(y)
    global v0 e0 r
    Q = 2*e0 ./ (1 + exp(r*(v0-(y))));
end