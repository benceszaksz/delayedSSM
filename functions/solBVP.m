function Phiout = solBVP(lambdatilde,Delta,NL,Psi)

syms theta
ncterms = length(Psi.lambda);
dim = size(Psi.c,1);

Phi.c(:,1) = Delta(lambdatilde)\NL;
Phi.lambda(1,1) = lambdatilde;
for iterm = 1 : ncterms
    Phi.c(:,1) = Phi.c(:,1) + Delta(lambdatilde)\(-Delta(Psi.lambda(iterm))*Psi.c(:,iterm)/(lambdatilde-Psi.lambda(iterm)));
    Phi.c(:,iterm+1) = Psi.c(:,iterm)/(lambdatilde-Psi.lambda(iterm));
    Phi.lambda(1,iterm+1) = Psi.lambda(iterm);
end

[lambda_u,~,i_u] = unique(Phi.lambda,'stable');
Phi2.lambda = lambda_u;
Phi2.c = zeros(dim,length(lambda_u));
for i = 1 : length(i_u)
    Phi2.c (:,i_u(i)) = Phi2.c (:,i_u(i))+Phi.c(:,i);
end

Phiout.c = [];
Phiout.lambda = [];
for i = 1 : size(Phi2.c,2)
    if Phi2.c(:,i) ~= zeros(dim,1)
        Phiout.c = [Phiout.c, Phi2.c(:,i)];
        Phiout.lambda = [Phiout.lambda, Phi2.lambda(i)];
    end
end

Phiout.sym = zeros(dim,1);
for iterm = 1 : length(Phiout.lambda)
    Phiout.sym = Phiout.sym+Phiout.c(:,iterm)*exp(Phiout.lambda(1,iterm)*theta);
end