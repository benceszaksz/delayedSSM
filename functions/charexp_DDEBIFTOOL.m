function charexp = charexp_DDEBIFTOOL(L,R,par,minreal)
% Function for the calculation of the characteristic exponents with the
% help of the MATLAB package DDE-BIFTOOL. Requires the installation of the
% DDE-BIFTOOL package.
% Inputs:
%   L - coefficient matrix of the linear non-delayed state
%   R - coefficient matrix of the linear delayed state
%   par - structure of parameters including the delay par.tau (which may
%   also be zero)
%   minreal - minimum value of the real part of the eigenvalues.
% Output:
%   charexp - vector of the dominant characteristic exponents
 
funcs=set_funcs(...
    'sys_rhs',@(xx,par)L*xx(:,1)+R*xx(:,2),...
    'sys_tau',@()[1]);

stst.kind='stst';
stst.parameter=par.tau;
stst.x=zeros(size(L,1),1);

method=df_mthod(funcs,'stst');
method.stability.minimal_real_part=minreal;
[stst,success]=p_correc(funcs,stst,[],[],method.point);
% compute its stability:
stst.stability=p_stabil(funcs,stst,method.stability);

charexp = stst.stability.l1;