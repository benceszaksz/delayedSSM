function nPhi = pairing_nPhi(p,R,tau,lambda,Phi)
nPhi = 0;
nterms = length(Phi.lambda);
for iterm = 1 : nterms
    if lambda == Phi.lambda(iterm)
        nPhi = nPhi + p*(eye(size(R))+R*exp(-lambda*tau)*tau)*Phi.c(:,iterm);
    else
        nPhi = nPhi + p*(eye(size(R))+R*exp(-lambda*tau)*(exp((lambda-Phi.lambda(iterm))*tau)-1)/(lambda-Phi.lambda(iterm)))*Phi.c(:,iterm);
    end
end
