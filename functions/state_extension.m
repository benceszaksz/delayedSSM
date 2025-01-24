function Phitilde = state_extension(Phi,tau,nlin_type)

if Phi == zeros(size(Phi))
    Phitilde = [Phi;Phi];    
else
    Phi = matlabFunction(Phi);
    switch nlin_type
        case 'non-delayed'
            Phitilde = Phi(0);  % extended state vector
        case 'delayed'
            Phitilde = Phi(-tau);  % extended state vector;
        case 'combined'
            Phitilde = [Phi(0);Phi(-tau)];  % extended state vector;
        otherwise
            warning('The type of nonlinearity is wrong, the code calculates with the combined option. The calculation may be longer but it is still accurate.')
            Phitilde = [Phi(0);Phi(-tau)];  % extended state vector;        
    end
end