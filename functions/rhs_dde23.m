function dy = rhs_dde23(t,y,Z,L,R,fNL)
% Rewrites the right-hand side of the EoM in a form, which is compatible
% with the dde23 solver
% Inputs:
%   t - time
%   y - actual state vector
%   Z - delayed state vector
%   L - coefficient matrix of the actual state vector y
%   R - coefficient matrix of the delayed state vector Z
%   fNL - structure of nonlinear terms including fNL.fnlin_func
% Output
%   dy - time derivative of the state vector

dy = L*y+R*Z+fNL.fnlin_func(y.',Z.');