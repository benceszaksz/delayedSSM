function Si = Ccalc(phi1,phi2,phi3,H3)
    dim = size(H3,1);
    dim_ext = size(H3,2);
    Sk = zeros(dim,dim_ext,dim_ext);
    Sj = zeros(dim,dim_ext);
    Si = zeros(dim,1);
    for i = 1 : dim_ext
        for j = 1:dim_ext
            for k = 1:dim_ext
                Sk(:,i,j) = Sk(:,i,j)+H3(:,i,j,k)*phi3(k);
            end
            Sj(:,i) = Sj(:,i)+Sk(:,i,j)*phi2(j);
        end  
        Si(:,1) = Si(:)+Sj(:,i)*phi1(i);
    end
end