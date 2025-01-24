function [W,beta,p,nW] = SSM_balance_complex(lambda,p,q,R,Delta,fNL,tau,nSSM,nlin_type,red_dyn_style)

%%
syms theta
W.W_1_0.c = q;
W.W_1_0.lambda = lambda;
W.W_1_0.sym = W.W_1_0.c*exp(W.W_1_0.lambda*theta);
W.W_0_1 = Wconj(W.W_1_0);
nW.nW_1_0 = 1;
nW.nW_0_1 = 0;

zind_comb2 = comb(0:nSSM,0:nSSM); % combinations
zind_comb2sum = sum(zind_comb2,2);

zind_all = zind_comb2((zind_comb2sum >zeros(size(zind_comb2sum))) & (zind_comb2sum<=nSSM*ones(size(zind_comb2sum))),:);
zindsort = sortrows([zind_all,sum(zind_all,2)],3);
zind_all = zindsort(:,[2,1]);
zind_basis = zind_all(zind_all(:,1)>=zind_all(:,2),:);

for iField = 1 : length(zind_all)-2
    zind = zind_all(iField+2,:); % the first element is known (W10 = s) 
    betaField_act = varname_zind('beta',zind);
    beta.(betaField_act) = 0;
    betalist(iField) = string(betaField_act);
end

switch red_dyn_style
    case 'manual'
        [indx,tf] = listdlg('PromptString',{'Select the coefficients in the reduced dynamics.',...
            'Multiple selection is possible with the Shift/Ctrl buttons or the with the select all option','','',''},'ListString',betalist,'InitialValue',5);
        if tf == 1
            for i = 1 : length(indx)
                beta.(betalist(indx(i))) = [];
            end
        else
            error('The user terminated the calculation.')
        end
    case 'normal_form'
        for i = 2 : nSSM/2+0.5
            zind = [i,i-1];
            betaFieldi = varname_zind('beta',zind);
            beta.(betaFieldi) = [];
        end
    case 'graph'
        for iField = 1 : length(zind_all)-2
            zind = zind_all(iField+2,:); % the first element is known (W10 = s) 
            betaField_act = varname_zind('beta',zind);
            beta.(betaField_act) = [];
        end
    otherwise
        error('This reduced dynamics style is not supported in case of projection to complex conjugate eigenvalues. Supported options: ''graph'', ''normal_form'', or ''manual''')
end


for iField = 1 : length(zind_basis)-1
    zind =  zind_basis(iField+1,:);
    betaField_act = varname_zind('beta',zind);
    WField_act = varname_zind('W',zind);
    
    [dWdz_Wind,dWdz_betaind,dWdcz_Wind,dWdcz_cbetaind] = get_Wbeta_ind_complex(nSSM,zind);
    
    NL = getNL_complex(nSSM,zind,tau,fNL,W,nlin_type);
    
    % beta
    if isempty(beta.(betaField_act))
        % dWdz*(beta20*z^2+beta11*z*cz+beta02*cz^2+...)
        ommit_row = ismember(dWdz_Wind,[1,0],'rows'); % this row corresponds to that beta parameter which is to be obtained 
        dWdz_Wind_prev = dWdz_Wind(~ommit_row,:);
        dWdz_betaind_prev = dWdz_betaind(~ommit_row,:);

        n_betaprev = size(dWdz_betaind_prev,1);
        betanPhi = 0;
        for ibeta = 1 : n_betaprev
            WFieldi = varname_zind('W',dWdz_Wind_prev(ibeta,:));
            betaFieldi = varname_zind('beta',dWdz_betaind_prev(ibeta,:));
            Phi = W.(WFieldi);
            betanPhi = betanPhi+1/factorial(dWdz_Wind_prev(ibeta,1)-1)/factorial(dWdz_Wind_prev(ibeta,2))*beta.(betaFieldi)*pairing_nPhi(p,R,tau,lambda,Phi);
        end

        % dWdcz*conj(beta20*z^2+beta11*z*cz+beta02*cz^2+...)
        ommit_row = ismember(dWdcz_Wind,[0,1],'rows'); % <n,W01>=0
        dWdcz_Wind_prev = dWdcz_Wind(~ommit_row,:);
        dWdcz_cbetaind_prev = dWdcz_cbetaind(~ommit_row,:);

        n_cbetaprev = size(dWdcz_cbetaind_prev,1);
        cbetanPhi = 0;
        for ibeta = 1 : n_cbetaprev
            WFieldi = varname_zind('W',dWdcz_Wind_prev(ibeta,:));
            betaFieldi = varname_zind('beta',dWdcz_cbetaind_prev(ibeta,:));
            Phi = W.(WFieldi);
            cbetanPhi = cbetanPhi+1/factorial(dWdcz_Wind_prev(ibeta,1))/factorial(dWdcz_Wind_prev(ibeta,2)-1)*conj(beta.(betaFieldi))*pairing_nPhi(p,R,tau,lambda,Phi);
        end
        beta.(betaField_act) = p*NL-betanPhi-cbetanPhi;
    end

    % beta_ij is already obtained, but beta_ji is also needed for the
    % determination of W_ij
    betaField_mirror = varname_zind('beta',flip(zind));
    if  isempty(beta.(betaField_mirror))
        zind_mirror = flip(zind);
        NL_mirror = getNL_complex(nSSM,zind_mirror,tau,fNL,W,nlin_type);
        [dWdz_Wind_mirror,dWdz_betaind_mirror,dWdcz_Wind_mirror,dWdcz_cbetaind_mirror] = get_Wbeta_ind_complex(nSSM,zind_mirror);
        
        % dWdz*(beta20*z^2+beta11*z*cz+beta02*cz^2+...)
        ommit_row = ismember(dWdz_Wind_mirror,[1,0],'rows'); % this rho corresponds to that beta parameter which is to be obtained 
        dWdz_Wind_prev = dWdz_Wind_mirror(~ommit_row,:);
        dWdz_betaind_prev = dWdz_betaind_mirror(~ommit_row,:);

        n_betaprev = size(dWdz_betaind_prev,1);
        betanPhi = 0;
        for ibeta = 1 : n_betaprev
            WFieldi = varname_zind('W',dWdz_Wind_prev(ibeta,:));
            betaFieldi = varname_zind('beta',dWdz_betaind_prev(ibeta,:));
            Phi = W.(WFieldi);
            betanPhi = betanPhi+1/factorial(dWdz_Wind_prev(ibeta,1)-1)/factorial(dWdz_Wind_prev(ibeta,2))*beta.(betaFieldi)*pairing_nPhi(p,R,tau,lambda,Phi);
        end

        % dWdcz*conj(beta20*z^2+beta11*z*cz+beta02*cz^2+...)
        ommit_row = ismember(dWdcz_Wind_mirror,[0,1],'rows'); % <n,W01>=0
        dWdcz_Wind_prev = dWdcz_Wind_mirror(~ommit_row,:);
        dWdcz_cbetaind_prev = dWdcz_cbetaind_mirror(~ommit_row,:);

        n_cbetaprev = size(dWdcz_cbetaind_prev,1);
        cbetanPhi = 0;
        for ibeta = 1 : n_cbetaprev
            WFieldi = varname_zind('W',dWdcz_Wind_prev(ibeta,:));
            betaFieldi = varname_zind('beta',dWdcz_cbetaind_prev(ibeta,:));
            Phi = W.(WFieldi);
            cbetanPhi = cbetanPhi+1/factorial(dWdcz_Wind_prev(ibeta,1))/factorial(dWdcz_Wind_prev(ibeta,2)-1)*conj(beta.(betaFieldi))*pairing_nPhi(p,R,tau,lambda,Phi);
        end
        beta.(betaField_mirror) = p*NL_mirror-betanPhi-cbetanPhi;
    end

    % W
    Psi_betaW.c = [];
    Psi_betaW.lambda = [];
    for ibetaW = 1 : size(dWdz_betaind,1)
        betaFieldi = varname_zind('beta',dWdz_betaind(ibetaW,:));
        WFieldi = varname_zind('W',dWdz_Wind(ibetaW,:));
        Psi_betaW.c = [Psi_betaW.c, -1/factorial(dWdz_Wind(ibetaW,1)-1)/factorial(dWdz_Wind(ibetaW,2))*beta.(betaFieldi)*W.(WFieldi).c];
        Psi_betaW.lambda = [Psi_betaW.lambda, W.(WFieldi).lambda];
    end
    for ibetaW = 1 : size(dWdcz_cbetaind,1)
        betaFieldi = varname_zind('beta',dWdcz_cbetaind(ibetaW,:));
        WFieldi = varname_zind('W',dWdcz_Wind(ibetaW,:));
        Psi_betaW.c = [Psi_betaW.c, -1/factorial(dWdcz_Wind(ibetaW,1))/factorial(dWdcz_Wind(ibetaW,2)-1)*conj(beta.(betaFieldi))*W.(WFieldi).c];
        Psi_betaW.lambda = [Psi_betaW.lambda, W.(WFieldi).lambda];
    end
    lambdatilde = zind(1)*lambda+zind(2)*conj(lambda);
    Wsol = solBVP(lambdatilde,Delta,NL,Psi_betaW);
    Wact.c = factorial(zind(1))*factorial(zind(2))*Wsol.c;
    Wact.sym = factorial(zind(1))*factorial(zind(2))*Wsol.sym;
    Wact.lambda = Wsol.lambda;
    W.(WField_act) = Wact;
    if zind(1)~=zind(2)
        WField_mirror = varname_zind('W',flip(zind));
        W.(WField_mirror) = Wconj(Wact);
    end
end

%nW
for iField = 1 : length(zind_all)-2
    zind =  zind_all(iField+2,:);
    betaField_act = varname_zind('beta',zind);
    nWField_act = varname_zind('nW',zind);
    if beta.(betaField_act) == 0
        WField_act = varname_zind('W',zind);
        Phi = W.(WField_act);
        nW.(nWField_act) = 1/factorial(zind(1))/factorial(zind(2))*pairing_nPhi(p,R,tau,lambda,Phi);    
    else
      nW.(nWField_act) = 0;
    end
end


