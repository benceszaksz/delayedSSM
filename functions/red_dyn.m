function red_dyn_act = red_dyn(lambda,beta)

syms z
FieldList = fieldnames(beta);
rhs = lambda*z;
for iField = 1 : length(FieldList)
    betaFieldi = FieldList{iField};
    zind = get_zind_from_varname(betaFieldi);
    rhs = rhs + beta.(betaFieldi)*z.^zind(1).*conj(z).^zind(2);
end
red_dyn_act = matlabFunction(rhs);