% I am glad you are using my code. I will be happy to receive your 
% feedback and help you if needed. If necessary, please contact me at
% szaksz@mm.bme.hu.
%
% This code is an example for the calculation of the spectral submanifold 
% (SSM) in a scalar nonlinear time delay system. More specifically, it is
% foxusing to the case of complex conjugate dominant eigenvalues. 
%
% The current equation of motion is the scalar DDE: 
% dx(t)= -x(t-tau)+x^2(t)-0.5*x^3(t)
%%
clear;
%close all;

set(0,'defaulttextinterpreter','latex')
set(0, 'DefaultLegendInterpreter', 'latex')
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex')

mycolor.green = [20 150 0]/255;
mycolor.brown = [210,105,30]/255;
mycolor.orange = [255,165,0]/255;
%% parameters
% time delay
par.tau = 1.55; % time delay

% type of nonlinearity (may be 'non-delayed', 'delayed' or 'combined')
nlin_type = 'non-delayed';

% style of reduced dynamics 'graph', 'normal_form' (only in case of complex conjugate roots), 'manual'
red_dyn_style = 'normal_form';

% order of the Spectral submanifold and the corresponding reduced dynamics
SSMorder = 7;

% the index of the state vector for visualization
outputind = 1;

% initial functions for two particular trajectories (determined with dde23)
initf1 = @(t) 0.6;
initf2 = @(t) 0.95;

% length of the dde23 simulation
T = 100;

% surface coloring needed for 3D visualization
map = interp1([0,0.5,1],[[0.6275    0.3216    0.1765]/2;[205,133,63]/255;1 1 1],[0:0.01:1]);

%% EoM

[L,R,fNL] = sepEoM(par,'EoM_scalar',nlin_type);

%% Spectrum of the linearized dynamics
% determine the dominant characteristic root with semi-discretization method
Neig = 6;
Ndisc = 30;
charexp = charexp_semidisc(L,R,par,Neig,'disc_num',Ndisc,'tol',1e-10);

% Dominant eigenvalue
lambda = charexp(1);

% plot the right-most eigenvalues
figure(1);
plot(real(charexp),imag(charexp),'x')
xlim([-2.5,1])
xlabel('Re$\lambda$')
ylabel('Im$\lambda$')
set(gca,"FontSize",12)

if isreal(lambda)
    warning('The dominant root is real. Only, the coefficients of the SSM will be calculated. For the visualization, see the file: SSM_main_scalar_real.m')
end

%% Calculation of the SSM corresponding to lambda
[W,beta,p,nW] = SSM_coeff(L,R,fNL,par.tau,lambda,SSMorder,'nlin_type',nlin_type,'red_dyn_style',red_dyn_style);

% Take the coefficients of the SSM at theta=0 and at theta=-tau
W0 = Wsubstheta(W,0);
Wmintau = Wsubstheta(W,-par.tau);

%% Obtain the trajectoris initiated from 'initf1/2' with the help of the built in dde23 solver
options = ddeset('RelTol',1e-4,'AbsTol',1e-8);
sol1 = dde23(@(t,y,Z) rhs_dde23(t,y,Z,L,R,fNL),par.tau,initf1,[0, T],options);
sol2 = dde23(@(t,y,Z) rhs_dde23(t,y,Z,L,R,fNL),par.tau,initf2,[0, T],options);

% discretization of the trajectory (dt=tau/r)
r = 1000;
dt = par.tau/r;

% time vectors and the corresponding trajectories
tvec1 = 0:dt:sol1.x(end);
tvec2 = 0:dt:sol2.x(end);
xt1 = deval(sol1,tvec1);
xt2 = deval(sol2,tvec2);

% coordinate of the trajectory along the eigenvector
xi1 = proj_to_manif(xt1,lambda,p,R,r,dt);
xi2 = proj_to_manif(xt2,lambda,p,R,r,dt);

% Plot the trajectories
figure(2);
plot(tvec1,xt1,'Color',mycolor.green,'LineWidth',1)
hold on;
plot(tvec2,xt2,'r','LineWidth',1)
title('Time evolution of the trajectories')
legend('Traj. \#1','Traj. \#2')
xlabel('$t$')
ylabel('$x$')
ylim([-4,4])
set(gca,"FontSize",12)

%% visualization of the SSM together with the obtained trajectory
if imag(lambda)~=0 % the dominant eigenvalues form a complex conjugate pair
    %% Check for the existence of a limit cycle and visualize it
    [rho,om] = SSM_per(lambda,beta); % rho: amplitude in the reduced dynamics, om: angular frequency of the limit cycle
    tvec1 = linspace(0,2*pi/om(1),100);         % time vector for the limit cycle
    tvec2 = linspace(0,2*pi/om(2),100);         % time vector for the limit cycle

    % limit cycle on the SSM
    zpervec_t1 = rho(1).*exp(1i*om(1)*tvec1);                    % actual state (t) 
    zpervec_t2 = rho(2).*exp(1i*om(2)*tvec2);                    % actual state (t) 

    % limit cycle in the original coordinates 
    Wper_t1 = W_eval_complex(W0,real(zpervec_t1),imag(zpervec_t1),outputind);               % actual state (t) 
    Wper_tmintau1 = W_eval_complex(Wmintau,real(zpervec_t1),imag(zpervec_t1),outputind);    % delayed state (t-tau) 
    Wper_t2 = W_eval_complex(W0,real(zpervec_t2),imag(zpervec_t2),outputind);               % actual state (t) 
    Wper_tmintau2 = W_eval_complex(Wmintau,real(zpervec_t2),imag(zpervec_t2),outputind);    % delayed state (t-tau) 
    
    % Add the approximated limit cycle to the figure
    figure(4)
    hold on
    plot(Wper_tmintau1,Wper_t1,'r--','LineWidth',2)
    hold on
    plot(Wper_tmintau2,Wper_t2,'--','Color',mycolor.green,'LineWidth',2)
    legend('unstable limit cycle','stable limit cycle')
    xlabel('$x$')
    ylabel('$\dot{x}$')
    
end
