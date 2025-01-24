function Si = Bcalc(phi1,phi2,H2)
    dim = size(H2,1);
    dim_ext = size(H2,2);
    Sj = zeros(dim,dim_ext);
    Si = zeros(dim,1);
    for i = 1 : dim_ext
        for j = 1:dim_ext
            Sj(:,i) = Sj(:,i)+H2(:,i,j)*phi2(j);
        end
        Si(:,1) = Si(:)+Sj(:,i)*phi1(i);
    end
end