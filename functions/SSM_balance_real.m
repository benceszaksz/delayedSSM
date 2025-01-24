function [W,beta,p,nW] = SSM_balance_real(lambda,p,q,R,Delta,fNL,tau,nSSM,nlin_type,red_dyn_style)

%%
dim = size(q,1);

syms theta
W.W_1.c = q;
W.W_1.lambda = lambda;
W.W_1.sym = W.W_1.c*exp(W.W_1.lambda*theta);
nW.nW_1 = 1;

for iField = 1 : nSSM-1
    zind = iField+1;
    betaFieldi =  varname_zind('beta',zind);
    betalist(iField) = string(betaFieldi);
    beta.(betaFieldi) = 0;
end

switch red_dyn_style
    case 'manual'
        [indx,tf] = listdlg('PromptString',{'Select the coefficients in the reduced dynamics.',...
            'Multiple selection is possible with the Shift/Ctrl buttons or the with the select all option','','',''},'ListString',betalist,'InitialValue',1:nSSM-1);
        if tf == 1
            for i = 1 : length(indx)
                beta.(betalist(indx(i))) = [];
            end
        else
            error('The user terminated the calculation.')
        end
    case 'graph'
        for i = 1 : nSSM-1
            betaFieldi = varname_zind('beta',i+1);
            beta.(betaFieldi) = [];
        end
    otherwise
    error('This reduced dynamics style is not supported in case of projection to a real eigenvalue. Supported options: ''graph'' or ''manual''')
end

indcombinations2 = comb(1:nSSM,1:nSSM); % combinations
sumindcombinations2 = sum(indcombinations2,2);

indcombinations3 = comb(1:nSSM,1:nSSM,1:nSSM); % combinations
sumindcombinations3 = sum(indcombinations3,2);

for iField = 1 : nSSM-1
    zind =  iField+1;
    betaField_act = varname_zind('beta',zind);
    WField_act = varname_zind('W',zind);

    % calculate the indeces necessary for the following calculations
    betaW_betaind = (2:zind).'; % because beta_1=lambda is known
    betaW_Wind = zind-betaW_betaind+1; % +1: because of the derivation
    
    NL_b_Wind = indcombinations2(sumindcombinations2 == zind,:);
    NL_b_Wind = sort(NL_b_Wind,2);
    [NL_b_Wind_unique,~,ic] = unique(NL_b_Wind,'rows','stable'); % index/indices
    b_mult = accumarray(ic, 1);  % multiplier(s)
    
    NL_c_Wind = indcombinations3(sumindcombinations3 == zind,:); 
    NL_c_Wind = sort(NL_c_Wind,2);
    [NL_c_Wind_unique,~,ic] = unique(NL_c_Wind,'rows'); % index/indices
    c_mult = accumarray(ic, 1); % multiplier(s)

    %NL
    NL = zeros(dim,1);
    for ib = 1 : size(NL_b_Wind_unique,1)
        mult = b_mult(ib,1);
        ind1 = NL_b_Wind_unique(ib,1);
        ind2 = NL_b_Wind_unique(ib,2);
        WField1 = varname_zind('W',ind1);
        WField2 = varname_zind('W',ind2);
        Phi1tilde = state_extension(W.(WField1).sym,tau,nlin_type);
        Phi2tilde = state_extension(W.(WField2).sym,tau,nlin_type);
        NL = NL + mult*1/2*1/factorial(ind1)*1/factorial(ind2)*Bcalc(Phi1tilde,Phi2tilde,fNL.H2);
    end
    for ic = 1 : size(NL_c_Wind_unique,1)
        mult = c_mult(ic,1);
        ind1 = NL_c_Wind_unique(ic,1);
        ind2 = NL_c_Wind_unique(ic,2);
        ind3 = NL_c_Wind_unique(ic,3);
        WField1 = varname_zind('W',ind1);
        WField2 = varname_zind('W',ind2);
        WField3 = varname_zind('W',ind3);
        Phi1tilde = 1/factorial(ind1)*state_extension(W.(WField1).sym,tau,nlin_type);
        Phi2tilde = 1/factorial(ind2)*state_extension(W.(WField2).sym,tau,nlin_type);
        Phi3tilde = 1/factorial(ind3)*state_extension(W.(WField3).sym,tau,nlin_type);
        NL = NL + mult*1/6*Ccalc(Phi1tilde,Phi2tilde,Phi3tilde,fNL.H3);
    end
    
    %beta
    if isempty(beta.(betaField_act))
        betaprev_loc = find(betaW_betaind<zind); %get the previous beta-s
        betaprev_ind = betaW_betaind(betaprev_loc);
        betaW_Wind_prev = betaW_Wind(betaprev_loc);
        n_betaprev = length(betaprev_loc);
        betanPhi = 0;
        for ibeta = 1 : n_betaprev
            WFieldi = varname_zind('W',betaW_Wind_prev(ibeta));
            betaFieldi = varname_zind('beta',betaprev_ind(ibeta));
            Phi = W.(WFieldi);
            betanPhi = betanPhi+1/factorial(betaW_Wind_prev(ibeta)-1)*beta.(betaFieldi)*pairing_nPhi(p,R,tau,lambda,Phi);
        end
        beta.(betaField_act) = p*NL-betanPhi;
    end

    %W
    Psi_betaW.c = [];
    Psi_betaW.lambda = [];
    for ibetaW = 1 : length(betaW_betaind)
        betaFieldi = varname_zind('beta',betaW_betaind(ibetaW));
        WFieldi = varname_zind('W',betaW_Wind(ibetaW));
        Psi_betaW.c = [Psi_betaW.c, -1/factorial(betaW_Wind(ibetaW)-1)*beta.(betaFieldi)*W.(WFieldi).c];
        Psi_betaW.lambda = [Psi_betaW.lambda, W.(WFieldi).lambda];
    end
    lambdatilde = zind*lambda;
    Wact = solBVP(lambdatilde,Delta,NL,Psi_betaW);
    Wact.c = factorial(zind)*Wact.c;
    Wact.sym = factorial(zind)*Wact.sym;
    W.(WField_act) = Wact;    
end

%nW
for iField = 1 : nSSM-1
    zind =  iField+1;
    betaField_act = varname_zind('beta',zind);
    nWField_act = varname_zind('nW',zind);
    if beta.(betaField_act) == 0
        WField_act = varname_zind('W',zind);
        Phi = W.(WField_act);
        nW.(nWField_act) = 1/factorial(zind)*pairing_nPhi(p,R,tau,lambda,Phi);
    else
       nW.(nWField_act) = 0;
    end
end


