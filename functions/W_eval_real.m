function Wz = W_eval_real(W,z,outputind)
% Evaluation of the SSM over the parametrization variable z in case of real
% dominant eigenvalue.
% Inputs:
%   W - coefficients of the SSM
%   z - parametrization variable (may be a vector)
%   outputind - index of that component of the state vector, for which the
%               evaluation should be carried out
% Output:
%   Wz - SSM over z

if nargin == 2
    outputind = 1: length(W.W_1);      % by default, the whole space vector is used
end

if size(z,1)>1
    z = z.';
end

FieldList = fieldnames(W);
Wz = zeros(length(outputind),size(z,2));
for iField = 1 : length(FieldList)
    WFieldi = FieldList{iField};
    zind = get_zind_from_varname(WFieldi);
    Wz = Wz + 1/factorial(zind)*W.(WFieldi)(outputind)*z.^zind;
end