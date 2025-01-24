function [dWdz_Wind_ij,dWdz_betaind_kl,dWdcz_Wind_ij,dWdcz_cbetaind_kl] = get_Wbeta_ind_complex(nSSM,zind)

zind_comb2 = comb(0:nSSM,0:nSSM); % combinations
zind_comb2sum = sum(zind_comb2,2);

% indices for dWdz*(beta20*z^2+beta11*z*cz+beta02*cz^2+...)
% W_ij and beta_kl
dWdzbeta_ik = zind_comb2(zind_comb2sum == (zind(1)+1),:);
dWdzbeta_jl = zind_comb2(zind_comb2sum == (zind(2)),:);
dWdzb_ij = comb(dWdzbeta_ik(:,1),dWdzbeta_jl(:,1));
dWdzb_kl = comb(dWdzbeta_ik(:,2),dWdzbeta_jl(:,2));
[dWdz_Wind_ij,dWdz_betaind_kl] = Wbetaind_complex(dWdzb_ij,dWdzb_kl,'dWdz');
% indices for dWdcz*conj(beta20*z^2+beta11*z*cz+beta02*cz^2+...)
% W_ij and conj(beta_kl)
dWdczcbeta_il = zind_comb2(zind_comb2sum == (zind(1)),:);
dWdczcbeta_jk = zind_comb2(zind_comb2sum == (zind(2)+1),:);
dWdczcb_ij = comb(dWdczcbeta_il(:,1),dWdczcbeta_jk(:,1));
dWdczcb_lk = comb(dWdczcbeta_il(:,2),dWdczcbeta_jk(:,2));
dWdczcb_kl = [dWdczcb_lk(:,2),dWdczcb_lk(:,1)];
[dWdcz_Wind_ij,dWdcz_cbetaind_kl] = Wbetaind_complex(dWdczcb_ij,dWdczcb_kl,'dWdconjz');

    function [Wind_ij,betaind_kl] = Wbetaind_complex(indij,indkl,difftype)
        % W_ij and beta_kl
        % there is no beta_00, beta_10, beta_01
        ommit_row = ismember(indkl,[0,0],'rows');
        indij = indij(~ommit_row,:);
        indkl = indkl(~ommit_row,:);
        ommit_row = ismember(indkl,[1,0],'rows');
        indij = indij(~ommit_row,:);
        indkl = indkl(~ommit_row,:);
        ommit_row = ismember(indkl,[0,1],'rows');
        indij = indij(~ommit_row,:);
        indkl = indkl(~ommit_row,:);

        switch difftype
            case 'dWdz'
                ommit_row = ismember(indij(:,1),0,'rows');   %d(W_0j*conj(z)^j)/dz = 0
                indij = indij(~ommit_row,:);
                indkl = indkl(~ommit_row,:);
            case 'dWdconjz'
                ommit_row = ismember(indij(:,2),0,'rows');   %d(W_i0*z^i)/d(conj(z)) = 0
                indij = indij(~ommit_row,:);
                indkl = indkl(~ommit_row,:);
            otherwise
                error('The variable difftype is wrong within the arguments of the function Wbetaind_complex(indij,indkl,difftype).')
        end

        Wind_ij = indij;
        betaind_kl = indkl;
    end
end