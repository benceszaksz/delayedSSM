function xi = near_Id_transf(Rez,Imz,nW)
% Near identity transformation xi(z)
% Inputs:
%   Rez - real part of the variable z (may be a matrix)
%   Imz - imaginary part of the variable z (may be a matrix)
%   nW - structure containing the scalar product of the eigenfunction n and
%        the coefficients of the SSM
% Output:
%   xi - parametrization variable after the transformation

FieldList = fieldnames(nW);
xi = 0;
for iField = 1 : length(FieldList)
    nWFieldi = FieldList{iField};
    zind = get_zind_from_varname(nWFieldi);
    xi = xi + nW.(nWFieldi)*(Rez+1i*Imz).^zind(1).*(Rez-1i*Imz).^zind(2);
end
    