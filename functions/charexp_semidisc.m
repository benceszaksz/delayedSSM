function charexpout = charexp_semidisc(A,R,par,Neig,varargin)
% Function for the calculation of the characteristic exponents with the
% help of the semidiscretization method, which is followed by a
% Newton-Raphson iteration.
% Inputs:
%   A - coefficient matrix of the linear non-delayed state
%   R - coefficient matrix of the linear delayed state
%   par - structure of parameters including the delay par.tau (which may
%   also be zero)
%   Neig - number of the eigenvalues to be returned
%
%   additional inputs (varargin):
%   disc_num - delay discretization number
%   tol - tolerance for the Newton-Raphson iteration
% Output:
%   charexpout - vector of the dominant characteristic exponents

% Assign optional parameters
if par.tau == 0
    r = 0;      % in the delay-free case there is no sampling
elseif par.tau ~= 0
    r = 100;   % the default number of the sampling number is 100
end
tol = 1e-12;
for i=1:2:length(varargin)
    switch lower(varargin{i})
        case 'disc_num'
            r = varargin{i+1};
        case 'tol'
            tol = varargin{i+1};
     otherwise
      error(message('SSM_coeff:UnknownParameter', varargin{ i }))
    end
end
%% Semidiscretization
dim = size(A,1);
dt = par.tau/(r+1/2);
syms s lambda
if r == 0 && par.tau == 0 % no time delay, just nonlinearity
    A_tau0 = A+R;
    Phi = expm(A_tau0*dt);
elseif r ~= 0 && par.tau ~= 0
    Ad = expm(A*dt);
    % If the condition number of matrix A is small, it is faster to
    % calculate int(expm(-s*A),0,dt) utilizing matrix decomposition
    if cond(A)<1e6 
        [P,J] = eig(A);
        intexpA = real(P*double(int(expm(-s*J),0,dt))/P);
    else
        intexpA = double(int(expm(-s*A),0,dt));
    end

    Bd = Ad*intexpA*R;
    Phi = [Ad,zeros(dim,dim*(r-1)),Bd;
        eye(dim*r,dim*r),zeros(dim*r,dim)];
else
    error('The parameters r and par.tau are not consistent. In case of a delay free system, r should be 0.')
end

%% calculation of the characteristic exponents
eigvec = eig(Phi);      % characteristic multipliers
charexp = (log(eigvec))/dt; % characteristic exponents
[~,ind] = sort(-real(charexp));
charexp = charexp(ind);

charexpout = charexp(1:Neig);

%% Newton-Raphson iteration for more accurate characteristic exponents
% characteristic function and its derivative with respect to lambda
charfunc = det(lambda*eye(size(A))-A-R*exp(-lambda*par.tau));
dcharfunc = diff(charfunc,lambda);
charfunc = matlabFunction(charfunc);
dcharfunc = matlabFunction(dcharfunc);

err = max(abs(charfunc(charexpout)));
iiter = 1;
while err>tol && iiter<100
    charexpout = charexpout-charfunc(charexpout)./dcharfunc(charexpout);
    err = max(abs(charfunc(charexpout)));
    iiter = iiter+1;
end
if err>tol
    warning('The iteration for the tuning of the characteristic exponents did not converge under the tolerance within 100 steps.')
end
