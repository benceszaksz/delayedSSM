function [eq, xveccell] = EoM_scalar(par)
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
% The current equation of motion is the scalar DDE: 
% dx(t)= -x(t-tau)+x(t)^3

%% dimension of the system and the corresponding state vectors
dim = 1;

x = sym('x',[dim,1]);             % state vector (symbolic)
dx = sym('dx',[dim,1]);           % first derivative of the state vector (symbolic)
xdel = sym('xdel',[dim,1]);       % delayed state vector (symbolic)
xveccell = {x,dx,xdel};
%% Nonlinear equation of motion
eq = dx== -xdel+x^3;              % equation of motion
