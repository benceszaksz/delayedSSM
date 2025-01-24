function zout = inv_near_ID_transf(yin,nW)
syms z cz

FieldList = fieldnames(nW);
y = 0;
cy = 0;
for iField = 1 : length(FieldList)
    nWFieldi = FieldList{iField};
    zind = get_zind_from_varname(nWFieldi);
    ycoeffField = varname_zind('ycoeff',zind);
    ycoeff.(ycoeffField) = sym((ycoeffField));
    y = y+nW.(nWFieldi)*z.^zind(1).*cz^zind(2);
    cy = cy + conj(nW.(nWFieldi))*z.^zind(2).*cz^zind(1);
end

zsum = -z;
for iField = 1 : length(FieldList)
    nWFieldi = FieldList{iField};
    zind = get_zind_from_varname(nWFieldi);
    ycoeffField = varname_zind('ycoeff',zind);
    zsum = zsum + ycoeff.(ycoeffField)*y^zind(1)*cy^zind(2);
end

SSM_order = zind(1)+zind(2);

C = coeffs(expand(zsum),[z,cz],"All");

for iField = 1 : length(FieldList)
    nWFieldi = FieldList{iField};
    zind = get_zind_from_varname(nWFieldi);
    ycoeffField = varname_zind('ycoeff',zind);
    ycoeffnum.(ycoeffField) = double(vpasolve(C(end-zind(1),end-zind(2))==0,ycoeff.(ycoeffField)));
    C = subs(C,ycoeff.(ycoeffField),ycoeffnum.(ycoeffField));
end        

zout = zeros(size(yin));
for iField = 1 : length(FieldList)
    nWFieldi = FieldList{iField};
    zind = get_zind_from_varname(nWFieldi);
    ycoeffField = varname_zind('ycoeff',zind);
    zout = zout + ycoeffnum.(ycoeffField)*yin.^zind(1).*conj(yin).^zind(2);
end 



