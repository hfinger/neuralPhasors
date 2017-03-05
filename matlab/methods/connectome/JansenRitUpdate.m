function [ dx, da ] = JansenRitUpdate( x, a, Pp, Pe, Ke, Ki, Ka, Cpe, Cpi, Cep, Cip, Cii, de1, de2, di1, di2, e0, u0, r, S0 )
% JANSENRITUPDATE Function that updates the 6 state variables of the
% standard Jansen-Rit-Modell (Jansen & Rit, 1995)
% Parameters:
%   x   -   2-dim array with 1.dim = Nodes and 2.dim = states
%   Pp  -   external input to the pyramidal cells. Vector with lenght =
%           number of nodes
%   Pe  -   input from network to excitatory interneurons. Vector with
%           length = number of nodes

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

%% voltage/current changes

% excitatory pyramidal cells
dx(:,1) = x(:,4) - x(:,9) + Pp;

    % excitatory input
    dx(:,7) = x(:,4);
    dx(:,4) = Mep + de1*x(:,4) + de2*x(:,7);
    
    % inhibitory input
    dx(:,8) = x(:,9);
    dx(:,9) = Mip + di1*x(:,9) + di2*x(:,8);

% excitatory interneurons
dx(:,2) = x(:,5);
dx(:,5) = Ke*Pe + Mpe + de1*x(:,5) + de2*x(:,2);

% inhibitory interneurons
dx(:,3) = x(:,11) - x(:,13);

    % excitatory input
    dx(:,10) = x(:,11);
    dx(:,11) = Mpi + de1*x(:,11) + de2*x(:,10);
    
    % recurrent intrinsic inhibitory connection
    dx(:,12) = x(:,13);
    dx(:,13) = Mii + di1*x(:,13) + di2*x(:,12);

%% spike-frequency adaption

da = Ka * (Spe - a);

end

