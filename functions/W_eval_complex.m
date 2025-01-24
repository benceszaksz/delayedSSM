function [Wsurface,wsurface_xi,Rexi_grid,Imxi_grid] = W_eval_complex(W,Rez,Imz,outputind,nW)
% Evaluation of the SSM over the real and imaginary parts of the 
% parametrization variable in case of complex dominant eigenvalue.
% Inputs:
%   W - coefficients of the SSM
%   Rez - real part of the parametrization variable (may be a matrix -> grid)
%   Imz - imag part of the parametrization variable (may be a matrix -> grid)
%   outputind - index of that component of the state vector, for which the
%               evaluation should be carried out
%   nW - structure containing the scalar product of the eigenfunction n and
%        the coefficients of the SSM
% Output:
%   Wsurface - SSM over z
%   wsurface_xi - nonlinear part of the SSM above xi
%   Rexi_grid - real part of the parametrization variable xi 
%   Imxi_grid - real part of the parametrization variable xi 
% The outputs are matrices if the imputs Rez and Imz are matrices.
% The parametrization variables z and xi are related to each other through
% the near identity transformation (see near_Id_transf.m).
% The function may be used in all the following forms:
% [Wsurface,wsurface_xi,Rexi_grid,Imxi_grid] = Wsurf(W,Rez,Imz,outputind,nW)
% and
% Wsurface = Wsurf(W,Rez,Imz,outputind,nW)
% or
% Wsurface = Wsurf(W,Rez,Imz,outputind)

FieldList = fieldnames(W);
Wsurface = zeros(size(Rez));
for iField = 1 : length(FieldList)
    WFieldi = FieldList{iField};
    zind = get_zind_from_varname(WFieldi);
    Wsurface = Wsurface + 1/factorial(zind(1))/factorial(zind(2))*W.(WFieldi)(outputind)*(Rez+1i*Imz).^zind(1).*(Rez-1i*Imz).^zind(2);
end
Wsurface = real(Wsurface);

if nargin == 5
    xi_grid = near_Id_transf(Rez,Imz,nW);
    Rexi_grid = real(xi_grid);
    Imxi_grid = imag(xi_grid);

    wsurface_xi = Wsurface-W.W_1_0(outputind)*xi_grid-W.W_0_1(outputind)*conj(xi_grid);
end