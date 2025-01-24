function zind = get_zind_from_varname(varname)
% This function returns the exponent of the parametrization variable(s)
% obtained from the variable name. For example: W_5 yields zind = 5, or
% W_3_2 yields zind = [3,2];

varnamesplit = strsplit(varname,'_');

if length(varnamesplit) == 2
    zind = str2double(varnamesplit(2));
elseif length(varnamesplit) == 3
    zind1 = str2double(varnamesplit(2));
    zind2 = str2double(varnamesplit(3));
    zind = [zind1,zind2];
end