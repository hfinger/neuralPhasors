function b = mBOLDs(r,dt,d)

% The Balloon-Windkessel Hemodynamic Model
% Used in:
%    Cabral J, Hugues E, Sporns O, Deco G.
%    Role of local network oscillations in resting-state functional connectivity.
%    Neuroimage. 2011 Jul 1;57(1):130-9. Epub 2011 Apr 12.
%
%    Cabral J, Kringelbach ML, Deco G
%    Functional Graph Alterations in Schizophrenia: A Result from a Global Anatomic Decoupling?
%    Pharmacopsychiatry 45 (1), 57
%
% Note: If using phase oscillators, r=sin(ths). 
%
% b = BOLDs(r,dt,r)
%
% r       : neural activity
% dt      : time step (s)
% d       : initial seconds to remove (s)
%
% b       : BOLD signal

if nargin<3, d = 0; end
if size(r,2)*dt<=d, error(strcat('Simulation time must be greater than ',num2str(d),' seconds')); end

% BOLD model parameters

taus   = 0.65; % 0.8;    % time unit (s)
tauf   = 0.41; % 0.4;    % time unit (s)
tauo   = 0.98; % 1;      % mean transit time (s)
alpha  = 0.32; % 0.2;    % a stiffness exponent
itaus  = 1/taus;
itauf  = 1/tauf;
itauo  = 1/tauo;
ialpha = 1/alpha;
Eo     = 0.34; % 0.8;    % resting oxygen extraction fraction
vo     = 0.02;
k1     = 7*Eo; 
k2     = 2; 
k3     = 2*Eo-0.2;

% Time integration
R      = size(r,1);
T      = size(r,2)*dt;
n_t    = round(T/dt)+1;
b      = zeros(R,n_t-round(d/dt));

for i = 1:R
    x = zeros(n_t,4);
    x(1,:) = [0 1 1 1];
    for n = 1:n_t-1
        x(n+1,1) = x(n,1) + dt*( r(i,n)-itaus*x(n,1)-itauf*(x(n,2)-1) );
        x(n+1,2) = x(n,2) + dt*x(n,1);
        x(n+1,3) = x(n,3) + dt*itauo*(x(n,2)-x(n,3)^ialpha);
        x(n+1,4) = x(n,4) + dt*itauo*(x(n,2)*(1-(1-Eo)^(1/x(n,2)))/Eo - (x(n,3)^ialpha)*x(n,4)/x(n,3));
        if sum(isfinite(x(n+1,:)))<4 || isreal(x(n+1,:))==0, error('bold calculation failed !!!'); end
    end
    b_tmp = 100/Eo*vo*( k1.*(1-x(:,4)) + k2*(1-x(:,4)./x(:,3)) + k3*(1-x(:,3)) );
    b(i,:) = b_tmp(round(d/dt)+1:end);
end
