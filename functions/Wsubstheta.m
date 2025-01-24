function W_theta = Wsubstheta(W,theta_act, FieldList)
% Evaluates the SSM coefficients at a given value of theta
% Inputs:
%   W - structure of the SSM coefficients
%   theta_act - value of theta
%   FieldList - list of fields, in which the substitution should be carried
%               out (optional)
% Output:
%   W_theta - structure of the SSM coefficients after the substitution

syms theta
if nargin == 2
   FieldList = fieldnames(W);
end 
for iField = 1:numel(FieldList)
   Field    = FieldList{iField};
   W_theta.(Field) = double(subs(W.(Field).sym,theta,theta_act));
end
end