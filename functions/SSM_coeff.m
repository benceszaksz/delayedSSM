function [W,beta,p,nW] = SSM_coeff(L,R,fNL,tau,lambda,nSSM,varargin)
% Calculates the coefficients of the SSM and the corresponding reduced
% dynamics. Both the forms 
% [W,beta,p,nW] = SSM_coeff(L,R,fNL,tau,lambda)
% and the form
% [W,beta,p] = SSM_coeff(L,R,fNL,tau,lambda)
% are working. (The latter may be used for example in case of real dominant
% eigenvalues)
% Inputs:
%   L - coefficient matrix of the linear non-delayed state
%   R - coefficient matrix of the linear delayed state
%   fNL - structure containing the nonlinearities
%   tau - time delay of the system
%   lambda - eigenvalue, corresponding to which the projection is carried
%   nSSM - order of the spectral submanifold
%
%   additional inputs (varargin):
%   nlin_type - type of nonlinearity (optional)
%               options: 'non-delayed', 'delayed', 'combined'
%               by default, it is 'combined'
%   red_dyn_style - style of the reduced dynamics (optional)
%               options: 'graph', 'normal_form' (only in case of complex conjugate roots), 'manual'
%               by default, it is 'graph'
% Outputs:
%   W - structure of the coefficients of the SSM
%   beta - coefficient(s) in the reduced dynamics. In case of real lambda,
%          It is a vector containing the coefficients of the 2nd and 3rd
%          order terms in the reduced dynamics. In case of complex lambda, 
%          it contains the coefficient of the rezonant z^2*conj(z) term.
%   p - left eigenvector of the characteristic matrix
%   nW - structure containing the scalar product of the eigenfunction n and
%        the coefficients of the SSM

%% type of nonlinearity and reduced dynamiscs
nlin_type = 'combined';  % by default it is 'combined'
red_dyn_style = 'graph';  % by default it is 'graph'
for i=1:2:length(varargin)
    switch lower(varargin{i})
        case 'nlin_type'
            nlin_type = varargin{i+1};
        case 'red_dyn_style'
            red_dyn_style = varargin{i+1};
     otherwise
      error(message('SSM_coeff:UnknownParameter', varargin{ i }))
    end
end
%% Characteristic matrix and eigenvectors
dim = size(L,1);              % dimension of the system

Delta = @(lambda) lambda*eye(dim)-L-R*exp(-lambda*tau); % Characteristic matrix
DeltaL1 = Delta(lambda);
DeltaL1prime = eye(dim)+tau*R*exp(-lambda*tau);

%Calculation of the right eigenvector
q2end = inv(DeltaL1(2:end,2:end))*(-1*DeltaL1(2:end,1));
q = [1;q2end];    %right eigenvector

%Calculation of the left eigenvector
p2end = -DeltaL1(1,2:end)*inv(DeltaL1(2:end,2:end));
p1 = 1/([1,p2end]*DeltaL1prime*q);  %normalization
p = p1*[1,p2end];    %left eigenvector

%% Calculation of the coefficients of the SSM

if isreal(lambda) % lambda is real
    [W,beta,p,nW] = SSM_balance_real(lambda,p,q,R,Delta,fNL,tau,nSSM,nlin_type,red_dyn_style);
else    % lambda is complex
    [W,beta,p,nW] = SSM_balance_complex(lambda,p,q,R,Delta,fNL,tau,nSSM,nlin_type,red_dyn_style);
end