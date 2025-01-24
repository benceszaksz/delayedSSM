function [solrho,solom] = SSM_per(lambda,beta)

if imag(lambda)~=0
FieldList = fieldnames(beta);
norder = 0;
for iField = 1 : length(FieldList)
    betaFieldi = FieldList{iField};
    zind = get_zind_from_varname(betaFieldi);
    zind1 = zind(1);
    zind2 = zind(2);
    norder = max(norder,zind1+zind2);
    %check whether beta is in normal form
    if beta.(betaFieldi)~=0 && zind1~=(zind2+1)
        error('The reduced dynamics should be in Poincare normal form.')
    end
end
else
    error('Lambda is real')
end

if norder < 5
    if real(lambda)/real(beta.beta_2_1)<0
        solrho = sqrt(-real(lambda)/real(beta.beta_2_1));
        solom = imag(lambda)+imag(beta.beta_2_1)*solrho^2;
    else
        solrho = [];
        solom = [];
        warning('No limit cycle was found.')
    end
elseif norder < 7
    solrho = sqrt((-real(beta.beta_2_1)+sqrt(real(beta.beta_2_1)^2-4*real(lambda)*real(beta.beta_3_2)))/(2*real(beta.beta_3_2)));
    if isreal(solrho) && solrho>0
        solrho2 = sqrt((-real(beta.beta_2_1)-sqrt(real(beta.beta_2_1)^2-4*real(lambda)*real(beta.beta_3_2)))/(2*real(beta.beta_3_2)));
        if isreal(solrho2) && solrho2>0
            solrho = [solrho;solrho2];
        end
        solom = imag(lambda)+imag(beta.beta_2_1)*solrho.^2+imag(beta.beta_3_2)*solrho.^4;
    else
        solrho = [];
        solom = [];
        warning('No limit cycle was found.')
    end
else
    solrho_O5_1 = sqrt((-real(beta.beta_2_1)+sqrt(real(beta.beta_2_1)^2-4*real(lambda)*real(beta.beta_3_2)))/(2*real(beta.beta_3_2)));
    solrho_O5_2 = sqrt((-real(beta.beta_2_1)-sqrt(real(beta.beta_2_1)^2-4*real(lambda)*real(beta.beta_3_2)))/(2*real(beta.beta_3_2)));
    solrho_O5 = [solrho_O5_1;solrho_O5_2];
    solrho_O5 = solrho_O5(imag(solrho_O5)==0 & real(solrho_O5)>0);

    Reeq = real(lambda);
    syms rho2
    assume(rho2,'positive')
    for i = 2 : norder/2+0.5
        betaFieldi = varname_zind('beta',[i,i-1]);
        Reeq = Reeq + real(beta.(betaFieldi))*rho2^(i-1);
    end
    if isempty(solrho_O5)
        solrho_all = sqrt(double(vpasolve(Reeq == 0, rho2)));
    elseif length(solrho_O5) == 1
        solrho_all = sqrt(double(vpasolve(Reeq == 0, rho2,solrho_O5^2)));
    elseif length(solrho_O5) == 2
        solrho_all1 = sqrt(double(vpasolve(Reeq == 0, rho2,solrho_O5(1)^2)));
        solrho_all2 = sqrt(double(vpasolve(Reeq == 0, rho2,solrho_O5(2)^2)));
        solrho_all = unique([solrho_all1;solrho_all2]);
    end
    solrho = solrho_all(imag(solrho_all)==0);
    solrho = solrho(solrho>0);
    if isempty(solrho)
        solrho = [];
        solom = [];
        warning('No limit cycle was found.')
    else
        solom = imag(lambda);
        for i = 2 : norder/2+0.5
            betaFieldi = varname_zind('beta',[i,i-1]);
            solom = solom + imag(beta.(betaFieldi))*solrho.^(2*i-2);
        end
    end
end

