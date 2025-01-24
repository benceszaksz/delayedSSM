function out = comb(vec1,vec2,vec3)

if nargin == 2
    vec3 = [];
end

if size(vec1,2)~=1
    vec1 = vec1.';
end
if size(vec2,2)~=1
    vec2 = vec2.';
end
if size(vec3,2)~=1
    vec3 = vec3.';
end

n1 = size(vec1,1);
n2 = size(vec2,1);
n3 = size(vec3,1);

if n3 ~= 0
    out = nan(n1*n2*n3,3);
    ii = 0;
    for i1 = 1 : n1
        for i2 = 1 : n2
            for i3 = 1 : n3
                ii = ii+1;
                out(ii,:) = [vec1(i1),vec2(i2),vec3(i3)];
            end
        end
    end
else
    out = nan(n1*n2,2);
    ii = 0;
    for i1 = 1 : n1
        for i2 = 1 : n2
            ii = ii+1;
            out(ii,:) = [vec1(i1),vec2(i2)];
        end
    end
end