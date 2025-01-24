function NL = getNL_complex(nSSM,zind,tau,fNL,W,nlin_type)

dim = size(fNL.H2,1);

zind_comb2 = comb(0:nSSM,0:nSSM); % combinations
zind_comb2sum = sum(zind_comb2,2);
zind_comb3 = comb(0:nSSM,0:nSSM,0:nSSM); % combinations
zind_comb3sum = sum(zind_comb3,2);

% indices for the nonlinearity
% b(Wij,Wkl)
NL_b_ik = zind_comb2(zind_comb2sum == zind(1),:);
NL_b_jl = zind_comb2(zind_comb2sum == zind(2),:);
NL_b_ij = comb(NL_b_ik(:,1),NL_b_jl(:,1));
NL_b_kl = comb(NL_b_ik(:,2),NL_b_jl(:,2));
[b_Wind_ij,b_Wind_kl,b_mult] = bind_complex(NL_b_ij,NL_b_kl);

% c(Wij,Wkl,Wmn)
NL_c_ikm = zind_comb3(zind_comb3sum == zind(1),:);
NL_c_jln = zind_comb3(zind_comb3sum == zind(2),:);
NL_c_ij = comb(NL_c_ikm(:,1),NL_c_jln(:,1));
NL_c_kl = comb(NL_c_ikm(:,2),NL_c_jln(:,2));
NL_c_mn = comb(NL_c_ikm(:,3),NL_c_jln(:,3));
[c_Wind_ij,c_Wind_kl,c_Wind_mn,c_mult] = cind_complex(NL_c_ij,NL_c_kl,NL_c_mn);

% Calculation of the nonlinearity term on the right-hand-side
NL = zeros(dim,1);
for ib = 1 : size(b_Wind_ij,1)
    mult = b_mult(ib,1);
    ind1 = b_Wind_ij(ib,:);
    ind2 = b_Wind_kl(ib,:);
    WField1 = varname_zind('W',ind1);
    WField2 = varname_zind('W',ind2);
    Phi1tilde = 1/factorial(ind1(1))/factorial(ind1(2))*state_extension(W.(WField1).sym,tau,nlin_type);
    Phi2tilde = 1/factorial(ind2(1))/factorial(ind2(2))*state_extension(W.(WField2).sym,tau,nlin_type);
    NL = NL + mult*1/2*Bcalc(Phi1tilde,Phi2tilde,fNL.H2);
end
for ic = 1 : size(c_Wind_ij,1)
    mult = c_mult(ic,1);
    ind1 = c_Wind_ij(ic,:);
    ind2 = c_Wind_kl(ic,:);
    ind3 = c_Wind_mn(ic,:);
    WField1 = varname_zind('W',ind1);
    WField2 = varname_zind('W',ind2);
    WField3 = varname_zind('W',ind3);
    Phi1tilde = 1/factorial(ind1(1))/factorial(ind1(2))*state_extension(W.(WField1).sym,tau,nlin_type);
    Phi2tilde = 1/factorial(ind2(1))/factorial(ind2(2))*state_extension(W.(WField2).sym,tau,nlin_type);
    Phi3tilde = 1/factorial(ind3(1))/factorial(ind3(2))*state_extension(W.(WField3).sym,tau,nlin_type);
    NL = NL + mult*1/6*Ccalc(Phi1tilde,Phi2tilde,Phi3tilde,fNL.H3);
end