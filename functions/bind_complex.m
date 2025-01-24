function [indWij,indWkl,mult] = bind_complex(NL_b_ij,NL_b_kl)

zero_row = ismember(NL_b_ij,[0,0],'rows');
NL_b_ij = NL_b_ij(~zero_row,:);
NL_b_kl = NL_b_kl(~zero_row,:);
zero_row = ismember(NL_b_kl,[0,0],'rows');
NL_b_ij = NL_b_ij(~zero_row,:);
NL_b_kl = NL_b_kl(~zero_row,:);

NL_b1 = NL_b_ij(:,1)+1i*NL_b_ij(:,2);
NL_b2 = NL_b_kl(:,1)+1i*NL_b_kl(:,2);
NL_b_Wind = [NL_b1,NL_b2];

 NL_b_Wind = sort(NL_b_Wind,2);
[NL_b_Wind_unique,~,ic] = unique(NL_b_Wind,'rows','stable');
mult = accumarray(ic, 1);

indWij = [real(NL_b_Wind_unique(:,1)), imag(NL_b_Wind_unique(:,1))];
indWkl = [real(NL_b_Wind_unique(:,2)), imag(NL_b_Wind_unique(:,2))];