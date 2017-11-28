function [ dx, da ] = JansenRitUpdate( x, a, Pp, Psc, Inp, Ke, Ki, Ka, Cpe, Cpi, Cep, Cip, Cii, de1, de2, di1, di2, e0, u0, r, S0 )
% JANSENRITUPDATE Function that updates the 13 state variables of the
% Jansen-Rit-Modell as well as a spike rate adaption variable a.
% Follows the architecture as proposed by Moran et al. 2008.
%
%   Input Parameters:
%
%       x   -   2-dim array with 1.dim = Nodes and 2.dim = states
%       a   -   spike rate adaption (scalar)
%       Pp  -   External input to the pyramidal cells. Vector with length =
%               number of nodes [V]
%       Psc -   Input from subcortical regions to target neurons. 
%               Vector with length = number of nodes [mean firing rate]
%       Inp -   Input from network to all cells. Vector with
%               length = number of nodes [mean firing rate]
%       Ke  -   Average synaptic gain for excitatory synapses per second.
%               Scalar [V/s]
%       Ki  -   Average synaptic gain for inhibitory synapses per second.
%               Scalar [V/s]
%       Ka  -   Adaptation rate constant [1/s]
%       Cpe -   Connectivity from pyramidal cells to excitatory interneurons
%       Cpi -   Connectivity from pyramidal cells to inhibitory interneurons
%       Cep -   Connectivity from excitatory interneurons to pyramidal cells
%       Cip -   Connectivity from inhibitory interneurons to pyramidal cells
%       de1 -   Time constant of excitatory synapses for first derivative [s]
%       de2 -   Time constant of excitatory synapses for second derivative [s]
%       di1 -   Time constant of inhibitory synapses for first derivative [s]
%       di2 -   Time constant of inhibitory synapses for second derivative [s]
%       e0  -   Determines maximum mean firing rate [1/s]
%       u0  -   Membrane voltage for which 50 % of maximum mean firing rate is observed [V]
%       r   -   Steepness of the sigmoid transfer function from mean psp to mean firing rate [1/V]
%       S0  -   Mean firing rate constant offset (baseline firing rate)
%
%   Output:
%
%       dx  -   Derivative of the 13 state variables
%       da  -   Derivative of the spike rate adaption

%% pre-calculate constants

dx = zeros(size(x));

% mean spike rates
Se = WTP(x(:,2), e0, u0, r) - S0;
Si = WTP(x(:,3), e0, u0, r) - S0;
Spe = WTP(x(:,1) - a, e0, u0, r) - S0;
Spi = WTP(x(:,1), e0, u0, r) - S0;

% pre-synaptic inputs
Mep = Cep * Ke * Se;
Mip = Cip * Ki * Si;
Mpe = Cpe * Ke * Spe;
Mpi = Cpi * Ke * Spi;
Mii = Cii * Ki * Si;
Mex = Ke * Inp;
Msc = Ke * Psc;

%% voltage/current changes

% excitatory pyramidal cells
dx(:,1) = x(:,4) - x(:,9) + Pp;

    % excitatory input
    dx(:,7) = x(:,4);
    dx(:,4) = Msc(:,1) + Mex(:,1) + Mep + de1*x(:,4) + de2*x(:,7);
    
    % inhibitory input
    dx(:,8) = x(:,9);
    dx(:,9) = Mip + di1*x(:,9) + di2*x(:,8);

% excitatory interneurons
dx(:,2) = x(:,5);
dx(:,5) = Msc(:,2) + Mex(:,2) + Mpe + de1*x(:,5) + de2*x(:,2);

% inhibitory interneurons
dx(:,3) = x(:,11) - x(:,13);

    % excitatory input
    dx(:,10) = x(:,11);
    dx(:,11) = Msc(:,3) + Mex(:,3) + Mpi + de1*x(:,11) + de2*x(:,10);
    
    % recurrent intrinsic inhibitory connection
    dx(:,12) = x(:,13);
    dx(:,13) = Mii + di1*x(:,13) + di2*x(:,12);

%% spike-frequency adaption

da = Ka * (Spe - a);

end

