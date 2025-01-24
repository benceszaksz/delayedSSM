function [indWij,indWkl,indWmn,mult] = cind_complex(NL_c_ij,NL_c_kl,NL_c_mn)

zero_row = ismember(NL_c_ij,[0,0],'rows');
NL_c_ij = NL_c_ij(~zero_row,:);
NL_c_kl = NL_c_kl(~zero_row,:);
NL_c_mn = NL_c_mn(~zero_row,:);
zero_row = ismember(NL_c_kl,[0,0],'rows');
NL_c_ij = NL_c_ij(~zero_row,:);
NL_c_kl = NL_c_kl(~zero_row,:);
NL_c_mn = NL_c_mn(~zero_row,:);
zero_row = ismember(NL_c_mn,[0,0],'rows');
NL_c_ij = NL_c_ij(~zero_row,:);
NL_c_kl = NL_c_kl(~zero_row,:);
NL_c_mn = NL_c_mn(~zero_row,:);

NL_c1 = NL_c_ij(:,1)+1i*NL_c_ij(:,2);
NL_c2 = NL_c_kl(:,1)+1i*NL_c_kl(:,2);
NL_c3 = NL_c_mn(:,1)+1i*NL_c_mn(:,2);
NL_c_Wind = [NL_c1,NL_c2,NL_c3];

NL_c_Wind = sort(NL_c_Wind,2);
[NL_c_Wind_unique,~,ic] = unique(NL_c_Wind,'rows');
mult = accumarray(ic, 1);

indWij = [real(NL_c_Wind_unique(:,1)), imag(NL_c_Wind_unique(:,1))];
indWkl = [real(NL_c_Wind_unique(:,2)), imag(NL_c_Wind_unique(:,2))];
indWmn = [real(NL_c_Wind_unique(:,3)), imag(NL_c_Wind_unique(:,3))];