function Y = ode3new(odefun,tspan,Y,options)
%ODE3  Solve differential equations with a non-adaptive method of order 3.
%   Y = ODE3(ODEFUN,TSPAN,Y) with TSPAN = [T1, T2, T3, ... TN] integrates 
%   the system of differential equations y' = f(t,y) by stepping from T0 to 
%   T1 to TN. Function ODEFUN(T,Y) must return f(t,y) in a column vector.
%   The vector Y is the initial conditions at T0. Each row in the solution 
%   array Y corresponds to a time specified in TSPAN.
%
%   Y = ODE3(ODEFUN,TSPAN,Y,P1,P2...) passes the additional parameters 
%   P1,P2... to the derivative function as ODEFUN(T,Y,P1,P2...). 
%
%   This is a non-adaptive solver. The step sequence is determined by TSPAN
%   but the derivative function ODEFUN is evaluated multiple times per step.
%   The solver implements the Bogacki-Shampine Runge-Kutta method of order 3.
%
%   Example 
%         tspan = 0:0.1:20;
%         y = ode3(@vdp1,tspan,[2 0]);  
%         plot(tspan,y(:,1));
%     solves the system y' = vdp1(t,y) with a constant step size of 0.1, 
%     and plots the first component of the solution.   
%

if nargin < 4
  options = odeset();
end

if ~isnumeric(tspan)
  error('TSPAN should be a vector of integration steps.');
end

if ~isnumeric(Y)
  error('Y should be a vector of initial conditions.');
end

h = diff(tspan);
if any(sign(h(1))*h <= 0)
  error('Entries of TSPAN are not in order.') 
end  

N = length(tspan);

if ~isempty(options.OutputFcn)
  feval(options.OutputFcn,tspan,Y,'init');
end

for i = 2:N
  ti = tspan(i-1);
  hi = h(i-1);
  F1 = feval(odefun,ti,Y);
  F2 = feval(odefun,ti+0.5*hi,Y+0.5*hi*F1);
  F3 = feval(odefun,ti+0.75*hi,Y+0.75*hi*F2);
  Y = Y + (hi/9)*(2*F1 + 3*F2 + 4*F3);
  if ~isempty(options.OutputFcn)
    stop = feval(options.OutputFcn,tspan(i),Y,[]);
    if stop
      break;
    end
  end
end

if ~isempty(options.OutputFcn)
  feval(options.OutputFcn,[],[],'done');
end