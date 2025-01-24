% I am glad you are using my code. I will be happy to receive your 
% feedback and help you if needed. If necessary, please contact me at
% szaksz@mm.bme.hu.
%
% This code is an example for the calculation of the spectral submanifold 
% (SSM) in a scalar nonlinear time delay system. More specifically, it is
% foxusing to the case of a real dominant eigenvalue. For the case of
% complex conjugate eigenvalues see the code: SSM_main_scalar_complex.m
%
% The current equation of motion is the scalar DDE: 
% dx(t)= -x(t-tau)+x(t)^3
%%
clear;
close all;

set(0,'defaulttextinterpreter','latex')
set(0, 'DefaultLegendInterpreter', 'latex')
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex')

mycolor.blue = [0 102 234]/255;
mycolor.green = [20 150 0]/255;
mycolor.brown = [210,105,30]/255;
%% parameters
% time delay
par.tau = 0.25; % time delay

% type of nonlinearity (may be 'non-delayed', 'delayed' or 'combined')
nlin_type = 'non-delayed';

% style of reduced dynamics 'graph', 'normal_form' (only in case of complex conjugate roots), 'manual'
red_dyn_style = 'graph';

% order of the Spectral submanifold and the corresponding reduced dynamics
SSMorder = 3;

% the index of the state vector for visualization
outputind = 1;

% initial functions for two particular trajectories (determined with dde23)
initf1 = @(t) 4*t;
initf2 = @(t) -4*t;

% length of the dde23 simulation
T = 50;

%% Equation of motion

[L,R,fNL] = sepEoM(par,'EoM_scalar',nlin_type);

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
xlim([-20,1])
xlabel('Re$\lambda$')
ylabel('Im$\lambda$')
set(gca,"FontSize",12)

if imag(lambda)~=0
    warning('The dominant root is complex with nonzero imaginary part. Only, the coefficients of the SSM will be calculated. For the visualization, see the file: SSM_main_scalar_complex.m')
end

%% Calculation of the SSM corresponding to lambda
[W,beta,p] = SSM_coeff(L,R,fNL,par.tau,lambda,SSMorder,'nlin_type',nlin_type,'red_dyn_style',red_dyn_style);

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
figure;
plot(tvec1,xt1,'Color',mycolor.blue,'LineWidth',1)
hold on;
plot(tvec2,xt2,'Color',mycolor.green,'LineWidth',1)
title('Time evolution of the trajectories')
legend('Traj. \#1','Traj. \#2')
xlabel('$t$')
ylabel('$x$')
set(gca,"FontSize",12)

%% visualization of the SSM together with the obtained trajectory 
if isreal(lambda) % if the dominant eigenvalue is real
    % horizontal axis
    zvec = linspace(-1,1,1000);
    
    % SSM at theta=0 and at theta=-tau
    W0vec = W_eval_real(W0,zvec);
    Wmintauvec = W_eval_real(Wmintau,zvec);

    % orthogonal part of the SSM
    w0vec = W0vec-W0.W_1*zvec;

    % transformed trajectories
    xt1_transf = xt1(:,r+1:end)-W0.W_1*xi1;
    xt2_transf = xt2(:,r+1:end)-W0.W_1*xi2;

    %% visualization
    figure;
    % SSM and the transformed trajectory as a function of the
    % parametrization variable xi(z)=z
    subplot(1,2,1)
    plot(zvec,w0vec(outputind,:),'Color',mycolor.brown,'LineWidth',2)
    hold on;
    plot(xi1,xt1_transf(outputind,:),'Color',mycolor.blue,'LineWidth',2)
    plot(xi2,xt2_transf(outputind,:),'Color',mycolor.green,'LineWidth',2)
    plot(0,0,'g.','MarkerSize',20)
    legend('SSM','Traj. \#1','Traj. \#2')
    xlabel('$y(z)$')
    ylabel('$\mathbf{w}(z;0)$')
    title(['$\tau=$ ' num2str(par.tau)])
    xlim([-0.2,0.2])
    ylim([-0.03,0.03])
    set(gca,'FontSize',12)
    grid on

    % SSM and trajectory in the plane of the delayed and actual states
    subplot(1,2,2)
    plot(Wmintauvec(outputind,:),W0vec(outputind,:),'Color',mycolor.brown,'LineWidth',2)
    hold on;
    plot(xt1(outputind,1:end-r),xt1(outputind,r+1:end),'Color',mycolor.blue,'LineWidth',2)
    plot(xt2(outputind,1:end-r),xt2(outputind,r+1:end),'Color',mycolor.green,'LineWidth',2)
    plot(0,0, 'g.','MarkerSize',20)
    legend('SSM','Traj. \#1','Traj. \#2')
    xlabel('$x(t-\tau)$')
    ylabel('$x(t)$')
    title(['$\tau=$ ' num2str(par.tau)])
    xlim([-0.2,0.2])
    ylim([-0.2,0.2])
    set(gca,'FontSize',12)
    grid on
end

