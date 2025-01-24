% I am glad you are using my code. I will be happy to receive your 
% feedback and help you if needed. If necessary, please contact me at
% szaksz@mm.bme.hu.
%
% This code is an example for the calculation of the spectral submanifold 
% (SSM) in a 3-dimensional nonlinear time delay system. More specifically, 
% it is foxusing to the case of a real dominant eigenvalue. For the case of
% complex conjugate eigenvalues see the code: SSM_main_guided_complex.m
%
% The current equation of motion describes the guidance of a human-driven
% vehicle (HV) via an automated vehicle (AV):
%
% dh(t) = dv(t)-dv_{-1}(t)
% dv_{-1}(t) = alpha*(V(h(t-tau))-v_{-1}(t-tau))+beta*(v(t-tau)-v_{-1}(t-tau))
% dv(t) = betahat*(v(t-tau)-v_{ref})+betamin1*(v(t-tau)-v_{-1}(t-tau))
%
% where
% h : distance between the vehicles
% v_{-1} : velocity of the HV
% v : velocity of the AV
% V(h) = vmax*(3*hgo-hst-2*h)*(h-hst)^2/(hgo-hst)^3 : range policy function
% (see Szaksz et. al. (2023): Guided control of a human driver via an 
% automated vehicle, IFAC PapersOnline): 
%%
clear;
close all;

set(0,'defaulttextinterpreter','latex')
set(0, 'DefaultLegendInterpreter', 'latex')
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex')

mycolor.green = [20 150 0]/255;
mycolor.brown = [210,105,30]/255;
%% parameters
% parameters for the human driver
par.alpha = 0.3;
par.beta = 0.4;
par.vmax = 30;
par.vref = 26.55;
par.hgo = 55;
par.hst = 5; 

% parameters for the automated vehicle leading to real dominant root
par.betamin1 = 1;
par.betahat = 0.2;

% time delay
par.tau = 0.4;

% type of nonlinearity (may be 'non-delayed', 'delayed' or 'combined')
nlin_type = 'delayed';

% style of reduced dynamics 'graph', 'manual'
red_dyn_style = 'graph';

% order of the Spectral submanifold and the corresponding reduced dynamics
SSMorder = 3;

% the index of the state vector for visualization
outputind = 2;

% initial function for a particular trajectory (determined with dde23)
initf = @(t) [2;0;0];

% length of the dde23 simulation
T = 100; 

%% Equation of motion

[L,R,fNL] = sepEoM(par,'EoM_guided',nlin_type);

%% Spectrum of the linearized dynamics
% determine the dominant characteristic root with semi-discretization method 
Neig = 6;
Ndisc = 30;
charexp = charexp_semidisc(L,R,par,Neig,'disc_num',Ndisc,'tol',1e-10);

% Dominant eigenvalue
lambda = charexp(1);

% plot the right-most eigenvalues
figure;
plot(real(charexp),imag(charexp),'x')
xlim([-10,1])
xlabel('Re$\lambda$')
ylabel('Im$\lambda$')
set(gca,'FontSize',12)

if imag(lambda)~=0
    warning('The dominant root is complex with nonzero imaginary part. Only, the coefficients of the SSM will be calculated. For the visualization, see the file: SSM_main_guided_complex.m')
end

%% Calculation of the SSM corresponding to lambda
[W,beta,p] = SSM_coeff(L,R,fNL,par.tau,lambda,SSMorder,'nlin_type',nlin_type,'red_dyn_style',red_dyn_style);

% Take the coefficients of the SSM at theta=0 and at theta=-tau 
W0 = Wsubstheta(W,0);
Wmintau = Wsubstheta(W,-par.tau);

%% Obtain the trajectory initiated from 'initf' with the help of the built in dde23 solver
options = ddeset('RelTol',1e-4,'AbsTol',1e-8);
sol = dde23(@(t,y,Z) rhs_dde23(t,y,Z,L,R,fNL),par.tau,initf,[0, T],options);

% discretization of the trajectory (dt=tau/r)
r = 1000; 
dt = par.tau/r;

% time vector and the corresponding trajectory
tvec_dde23 = 0:dt:sol.x(end);         
xt = deval(sol,tvec_dde23);

% coordinate of the trajectory along the eigenvector
xi = proj_to_manif(xt,lambda,p,R,r,dt);

% Plot the trajectory
figure;
subplot(2,1,1)
plot(tvec_dde23,xt(1,:),'LineWidth',1)
xlabel('$t$ [s]')
ylabel('$h$ [m]')
set(gca,'FontSize',12)
subplot(2,1,2)
plot(tvec_dde23,xt([2,3],:),'LineWidth',1)
xlabel('$t$ [s]')
ylabel('$v$ [m/s]')
legend('$v_{-1}$','$v$')
set(gca,'FontSize',12)

%% visualization of the SSM together with the obtained trajectory 
if isreal(lambda) % if the dominant eigenvalue is real
    % horizontal axis
    zvec = linspace(-2,2,1001);

    % SSM at theta=0 and at theta=-tau
    W0vec = W_eval_real(W0,zvec);
    Wmintauvec = W_eval_real(Wmintau,zvec);

    % orthogonal part of the SSM
    w0vec = W0vec-W0.W_1*zvec;

    % transformed trajectory
    xt_transf = xt(:,r+1:end)-W0.W_1*xi;

    %% visualization
    figure;
    % SSM and the transformed trajectory as a function of the
    % parametrization variable xi(z)=z
    subplot(1,2,1) 
    plot(zvec,w0vec(outputind,:),'Color',mycolor.brown,'LineWidth',2)
    hold on; 
    plot(xi,xt_transf(outputind,:),'Color',mycolor.green,'LineWidth',2)
    plot(0,0,'g.','MarkerSize',20)
    legend('SSM','Traj.')
    xlabel('$y(z)$')
    ylabel('$w_2(z;0)$')
    xlim([-2,2])
    ylim([-1,1])
    set(gca,'FontSize',12)
    grid on
    
    % SSM and trajectory in the plane of the delayed and actual states
    subplot(1,2,2) 
    plot(Wmintauvec(outputind,:),W0vec(outputind,:),'Color',mycolor.brown,'LineWidth',2)
    hold on;
    plot(xt(outputind,1:end-r),xt(outputind,r+1:end),'Color',mycolor.green,'LineWidth',2)
    plot(0,0,'g.','MarkerSize',20)
    legend('SSM','Traj.')
    xlabel('$v_{-1}(t-\tau)$ [m/s]')
    ylabel('$v_{-1}(t)$ [m/s]')
    xlim([-0.6,0.6])
    ylim([-0.6,0.6])
    set(gca,'FontSize',12)
    grid on
end