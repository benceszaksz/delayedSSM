function varname_out = varname_zind(varname_in,zind)
% This function returns the variable name of a coefficient, which includes 
% the information of the corresponding exponents of the parametrization
% variable(s). For Example:  varname_zind('W',[5]) yields 'W_5' 
% or varname_zind('W',[3,2]) yields 'W_3_2


if length(zind) == 1
    varname_out = [varname_in '_' num2str(zind(1,1))];
elseif length(zind) == 2
    varname_out = [varname_in '_' num2str(zind(1,1)) '_' num2str(zind(1,2))];
else
    error('Error in function varname_zind, the dimension of input zind is wrong.')
end