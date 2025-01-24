function [eq, xveccell] = EoM_guided(par)
% Equation of motion to be investigated. A first order delay-differential
% equation (DDE) should be provided in a symbolic form, the origin of which
% is a fixed point.
% Input:
%   par - structure of parameters including the delay par.tau (which may
%   also be zero). Note that par.tau will only be used in other functions.
% Outputs:
%   eq - nonlinear equation of motion (symbolic)
%   xveccell - symbolic state vectors of appropriate dimension
%
% The current equation of motion describes the guidance of a human-driven
% vehicle (HV) via an automated vehicle (AV):
%
% dh(t) = dv(t)-dv_{-1}(t)
% dv_{-1}(t) = alpha*(V(h(t-tau))-v_{-1}(t-tau))+beta*(v(t-tau)-v_{-1}(t-tau))
% dv(t) = betahat*(v(t-tau)-v_{ref})+betamin1*(v(t-tau)-v_{-1}(t-tau))
%
% where
% h : distance between the vehicles
% v_{-1} : velocity of the HV
% v : velocity of the AV
% V(h) = vmax*(3*hgo-hst-2*h)*(h-hst)^2/(hgo-hst)^3 : range policy function
% (see Szaksz et. al. (2023): Guided control of a human driver via an 
% automated vehicle, IFAC PapersOnline): 

%% dimension of the system and the corresponding state vectors
dim = 3;    % dimension of the system

x = sym('x',[dim,1]);             % state vector (symbolic)
dx = sym('dx',[dim,1]);           % first derivative of the state vector (symbolic)
xdel = sym('xdel',[dim,1]);       % delayed state vector (symbolic)
xveccell = {x,dx,xdel};

%% Nonlinear equation of motion
% Carry out a Taylor expansion of the range policy function V(h) at the
% steady state
syms h
Vh = par.vmax*(3*par.hgo-par.hst-2*h)*(h-par.hst)^2/(par.hgo-par.hst)^3;
solh = sort(double(solve(Vh==par.vref,h,'MaxDegree',3)));
if abs(solh(2)-real(solh(2)))<1e-10
    hstar = real(solh(2));
else
    error('The steady state is in the saturated domain.');
end
Vhtilde = subs(Vh,h,xdel(1)+hstar)-par.vref;
kappa = double(subs(diff(Vhtilde,xdel(1)),xdel(1),0));
[cVht,tVht] = coeffs(expand(Vhtilde));

% Equation of motion: 
% dx(t) = A0*x(t)+Atau*x(t-tau)+Gnlin(x(t-tau))
A0 = [0,-1,1
    0, 0, 0
    0, 0, 0];

Atau = [0,0,0
    par.alpha*kappa, -(par.alpha+par.beta), par.beta
    0, par.betamin1, -(par.betahat+par.betamin1)];

Gnlin = [0;par.alpha*(Vhtilde-cVht(end)*tVht(end)-cVht(end-1)*tVht(end-1));0];

eq = dx == A0*x+Atau*xdel+Gnlin;              % equation of motion
