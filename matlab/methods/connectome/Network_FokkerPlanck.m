function ys = Network_FokkerPlanck(C,D,k,tau,v,t_max,dt,sampling,sig_n,d,verbose)

% FokkerPlanck network simulation
% 
% C                    = Matrix of coupling weights (NxN) between pairs of regions (can be directed and weighted) 
% D                    = Matrix of distances (in mm) (NxN)
% k                    = Coupling strength
% tau                  = Time scale (ms)
% v                    = Velocity (in m/s)
% t_max                = Total time of simulated activity (in seconds)
% dt                   = integration step (smaller than smaller delays) (in seconds) (eg. 0.0001 s = 0.1 ms)
% sampling             = sampling for saving simulated activity (eg. 10) => 10*dt  = 1ms
% sig_n                = standard deviation of noise in the phases (can be zero, in radians) 
% d                    = initial seconds to remove (transient dynamics) (in seconds)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<11, verbose = 0; end
if nargin<10, d = 0; end
if t_max<=d, error(strcat('Simulation time must be greater than ',num2str(d),' seconds')); end

n_t_total = ceil(t_max/dt);           % total number of time steps e.g 300*0.0001 = 3000000
n_ts      = ceil(n_t_total/sampling); % same with downsampling ex. 300000
n_tcycle  = 10000;                    % number of time steps in each cycle (for RAM memory reasons)
n_tcs     = n_tcycle/sampling;        % same with downsampling ex. 1000
n_cycles  = ceil(n_t_total/n_tcycle); % number of cycles ex. 300

tau        = tau/1e3;
N          = size(C,1);
C          = k*C*dt/tau;            % Scale the coupling strengths per time step
%v          = mean(D(C>0))/tau*1e3;  % conduction velocity (m/s)
stepsDelay = round(D/(v*dt*1e3));   % number of time steps for delays
sig_noise  = sig_n*sqrt(dt/tau);    % Scale noise per time step


if v>0
    n_td = fix(max(stepsDelay(:)))+10; % number of time steps for maximal delays
    n_tp = n_td+n_tcycle;              % number of time steps in one cycle with the time for delays
else
    n_td = 1;
    n_tp = n_tcycle;
end


y   = zeros(N,n_tp,'single');     % initialize phase timeseries for one cycle
ys  = zeros(N,n_ts,'single');     % initialize phase timeseries to save

id = cell(1,N); for i = 1:N, id{i} = find(C(i,:)); end


% Initialization

    y(:,1:n_td) = sig_n*randn(N,n_td);

% Equations integration

if verbose, tic; end

if v>0
    for c = 1:n_cycles
        if verbose, fprintf(['Cycle = ' num2str(c) ' of ' num2str(n_cycles)]); end
        y(:,n_td+1:n_tp) = 0;

        if c < n_cycles          % total number os steps in this cycle
            n_tpc = n_tp;    % normal nr of steps 10000
        else
            n_tpc = n_t_total-(n_cycles-1)*n_tcycle+n_td; % nr of steps to complete total time
        end

        for t = n_td:n_tpc-1
            dy = sig_noise*randn(N,1);
            for n = 1:N
                for p = id{n}
                    dy(n) = dy(n) + C(n,p)*y(p,t-stepsDelay(n,p));
                end
            end
            y(:,t+1) = (1-dt/tau)*y(:,t)+dy;
        end
        ni = (c-1)*n_tcs;
        ns = ceil((n_tpc-n_td)/sampling);
        ys(:,ni+1:ni+ns) = y(:,n_td+1:sampling:n_tpc);

        y(:,1:n_td)      = y(:,n_tp-n_td+1:n_tp);
        if verbose, fprintf(['   (elapsed time = ',num2str(toc),' s)\n']); end
    end
else
    for c = 1:n_cycles
        if verbose, fprintf(['Cycle = ' num2str(c) ' of ' num2str(n_cycles)]); end
        y(:,2:n_tcycle) = 0;

        if c < n_cycles          % total number os steps in this cycle
            n_tpc = n_tcycle;    % normal nr of steps 10000
        else
            n_tpc = n_t_total-(n_cycles-1)*n_tcycle; % nr of steps to complete total time
        end

        for t = 1:n_tpc-1
            y(:,t+1) = (1-dt/tau)*y(:,t)+C*y(:,t) + sig_noise*randn(N,1);
        end
        ni = (c-1)*n_tcs;
        ns = ceil(n_tpc/sampling);
        ys(:,ni+1:ni+ns) = y(:,1:sampling:n_tpc);

        y(:,1)      = y(:,n_tcycle);
        if verbose, fprintf(['   (elapsed time = ',num2str(toc),' s)\n']); end
    end    
end
if verbose, fprintf(['elapsed time = ',num2str(toc),' s\n']); end

% if max(abs(ys(:)))>3*sig_n, ys=ys./max(abs(ys(:))); end
% [B,A] = butter(5,2*dt*sampling*.15,'low'); for n = 1:N, ys(n,:) = filtfilt(B,A,double(ys(n,:))); end

ys(:,1:round(d/(dt*sampling))) = [];  % remove initial d seconds of simulations to exclude transient dynamics
