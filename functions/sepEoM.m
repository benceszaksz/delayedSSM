function  [A,R,fNL] = sepEoM(par,EoMfname,nlin_type)
% Separates the equation of motion (EoM) into linear (non-delayed, and 
% delayed) and nonlinear parts (2nd and 3rd order nonlinearities).
% Inputs:
%   par - structure of parameters including the delay par.tau (which may
%   also be zero). Note that par.tau will only be used in other functions.
%   EoMfname - the name of the function, which contains the EoM
%   nlin_type - type of nonlinearity (optional)
%               options: 'non-delayed', 'delayed', 'combined'
%               by default, it is 'combined'
% Outputs:
%   A - coefficient matrix of x(t)
%   R - coefficient matrix of x(t-tau)
%   fNL - structure including the nonlinear terms 
%         fNL.H2 and fNL.H3 are the coefficients of the 2nd and 3rd order
%         terms
%         fNL.fnlin_func is the nonlinearity as a function (for faster
%         computations)
%% type of nonlinearity
if nargin == 2
    nlin_type = 'combined';     % by default it is 'combined'
end
%% Call the EoM
callEoM = [EoMfname, '(par)'];
[eq, xveccell] = eval(callEoM);

% symbolic state vectors
x = xveccell{1};
dx = xveccell{2};
xdel = xveccell{3};
dim = size(x,1);

switch nlin_type
    case 'non-delayed'
        xtilde = x;  % extended state vector
    case 'delayed'
        xtilde = xdel;  % extended state vector;
    case 'combined'
        xtilde = [x;xdel];  % extended state vector;
    otherwise
        warning('The type of nonlinearity is wrong, the code calculates with the combined option. The calculation may be longer but it is still accurate.')
        nlin_type = 'combined';
        xtilde = [x;xdel];  % extended state vector;        
end


%% Separate the EoM
% Express the delayed state to the left-hand side
if dim==1    
    dxsoltemp = solve(eq,dx);
    dxsol.dx1 = dxsoltemp;
else
    dxsol = solve(eq,dx);
end

graddxsol_x = cellfun(@(z)gradient(z,x),struct2cell(dxsol),'UniformOutput',false);
graddxsol_xdel = cellfun(@(z)gradient(z,xdel),struct2cell(dxsol),'UniformOutput',false);

% coefficient matrix of the linear non-delayed state vector x(t)
A = subs([graddxsol_x{:}].',[x;xdel],zeros(size([x;xdel])));

% coefficient matrix of the linear delayed state vector x(t-tau)
R = subs([graddxsol_xdel{:}].',[x;xdel],zeros(size([x;xdel])));

% nonlinear part of the right-hand side
fnlin = simplify(struct2cell(dxsol)-A*x-R*xdel);
A = double(A);
R = double(R);
%% Calculate the gradients corresponding the 2nd and 3rd order nonlinearities
% get the dimension of the extended state
if strcmp(nlin_type,'combined')
    dim_ext = 2*dim;
else
    dim_ext = dim;
end
fNL.H2 =sym(zeros(dim,dim_ext,dim_ext));
fNL.H3 = sym(zeros(dim,dim_ext,dim_ext,dim_ext));
for j = 1 : dim_ext
    for i = 1 : dim_ext
        fNL.H2(:,i,j) = diff(diff(fnlin,xtilde(i)),xtilde(j));
    end
end
for k = 1 : dim_ext
    for j = 1 : dim_ext
        for i = 1 : dim_ext
            fNL.H3(:,i,j,k) = diff(diff(diff(fnlin,xtilde(i)),xtilde(j)),xtilde(k));
        end
    end
end

%% Create the fNL structure
fNL.H2 = double(subs(fNL.H2,xtilde,zeros(dim_ext,1)));
fNL.H3 = double(subs(fNL.H3,xtilde,zeros(dim_ext,1)));
fNL.fnlin_func = matlabFunction(fnlin,'Vars',{x.',xdel.'});
end