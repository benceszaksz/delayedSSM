function xi = proj_to_manif(xt,lambda,p,R,r,dt)
% Projection of a sampled trajectory to the SSM.
% Inputs:
%   xt - sampled trajectory
%   lambda - dominant eigenvalue
%   p - left eigenvector of the characteristic matrix
%   R - coefficient of the delayed terms in the linearized dynamics
%   r - sampling delay number
%   dt - sampling time
% Output:
%   xi - coordinate along the corresponding eigenvector

dim = size(R,1);
xt_reverse = flip(xt,2);
xt_reverse_vec = reshape(xt_reverse,[],1);
Tinv = zeros(1,dim*(r+1));
Tinv(1,1:dim) = p;
for i = 1 : r
    Tinv(1,i*dim+1:(i+1)*dim) = p*R*exp(lambda*(i-r)*dt)*dt;
end
xi_reverse = zeros(1,length(xt)-r);
for i = 1 : length(xt)-r
xi_reverse(:,i)= (Tinv*xt_reverse_vec((i-1)*dim+1:(i-1)*dim+dim*(r+1),1));
end
xi = flip(xi_reverse,2);