% I am glad you are using my code. I will be happy to receive your 
% feedback and help you if needed. If necessary, please contact me at
% szaksz@mm.bme.hu.
%
% This code is an example for the calculation of the spectral submanifold 
% (SSM) in a scalar nonlinear time delay system. More specifically, it is
% foxusing to the case of complex conjugate dominant eigenvalues. 
% For the case of a real eigenvalue see the code: SSM_main_scalar_real.m
%
% The current equation of motion is the scalar DDE: 
% dx(t)= -x(t-tau)+x(t)^3
%%
clear;
close all;

set(0,'defaulttextinterpreter','latex')
set(0, 'DefaultLegendInterpreter', 'latex')
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex')

mycolor.green = [20 150 0]/255;
mycolor.brown = [210,105,30]/255;
mycolor.orange = [255,165,0]/255;
%% parameters
% time delay
par.tau = 1.3; % time delay

% type of nonlinearity (may be 'non-delayed', 'delayed' or 'combined')
nlin_type = 'non-delayed';

% style of reduced dynamics 'graph', 'normal_form' (only in case of complex conjugate roots), 'manual'
red_dyn_style = 'normal_form';

% order of the Spectral submanifold and the corresponding reduced dynamics
SSMorder = 3;

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
    % horizontal grid
    [Rez,Imz] = meshgrid(-0.3:0.001:0.3,-0.3:0.001:0.3);
    
    % SSM at the desired coorinate (outputind)
    [Wsurface,wsurface_xi,Rexi_grid,Imxi_grid] = W_eval_complex(W0,Rez,Imz,outputind,nW);

    % transformed trajectories
    xt1_transf = xt1(:,r+1:end)-W0.W_1_0*xi1-W0.W_0_1*conj(xi1);
    xt2_transf = xt2(:,r+1:end)-W0.W_1_0*xi2-W0.W_0_1*conj(xi2);

    %% One may either project the trajectory to the manifold and visualize
    % that, or project only its initial function to the manifold, and
    % obtain a trajectory on the SSM with the corresponding reduced
    % dynamics. In the original version of the code, the latter one is
    % used.
    
    % projection to manifold:
%     W0vec1 = W_eval_complex(W0,real(xi1),imag(xi1),outputind); % projection to the SSM at theta=0
%     Wmintauvec1 = W_eval_complex(Wmintau,real(xi1),imag(xi1),outputind); %  projection to the SSM at theta=-tau
%     W0vec2 = W_eval_complex(W0,real(xi2),imag(xi2),outputind); % projection to the SSM at theta=0
%     Wmintauvec2 = W_eval_complex(Wmintau,real(xi2),imag(xi2),outputind); %  projection to the SSM at theta=-tau

    % reduced dynamics on the manifold:
    % Traj. #1
    z0_1 = xi1(1);
    red_dyn_act = red_dyn(lambda,beta);
    solSSM1 = ode45(@(t,z) red_dyn_act(z), [0,T],z0_1);
    tvec_SSM1 = 0:dt:solSSM1.x(end);
    zsol1 = deval(solSSM1,tvec_SSM1);
    W0vec1 = W_eval_complex(W0,real(zsol1),imag(zsol1),outputind); % reduced dynamics on the SSM at theta=0
    Wmintauvec1 = W_eval_complex(Wmintau,real(zsol1),imag(zsol1),outputind); % reduced dynamics on the SSM at theta=-tau
    % Traj. #2
    z0_2 = xi2(1);
    solSSM2 = ode45(@(t,z) red_dyn_act(z), [0,T],z0_2);
    tvec_SSM2 = 0:dt:solSSM2.x(end);
    zsol2 = deval(solSSM2,tvec_SSM2);
    W0vec2 = W_eval_complex(W0,real(zsol2),imag(zsol2),outputind); % reduced dynamics on the SSM at theta=0
    Wmintauvec2 = W_eval_complex(Wmintau,real(zsol2),imag(zsol2),outputind); % reduced dynamics on the SSM at theta=-tau

    %% Generate the figure
    figure(3)
    % SSM above the plane of Re(xi) and Im(xi), and the transformed
    % trajectory
    subplot(1,2,1)
    hold off
    sl = surfl(Rexi_grid,Imxi_grid,wsurface_xi,[-20 20], [.7 .8 0.1 10]);
    sl.FaceAlpha = 0.8;
    colormap(map)
    shading interp
    hold on
    plot3(real(xi1),imag(xi1),xt1_transf(outputind,:),'Color',mycolor.green,'LineWidth',2)
    plot3(0,0,0,'.', 'color', 'g','MarkerSize',20)
    xlabel('Re($y(z,\overline{z})$)')
    ylabel('Im($y(z,\overline{z})$)')
    zlabel('$\mathbf{w}(z,\overline{z};0)$')
    title(['$\tau=$' num2str(par.tau)])
    xlim([-0.3,0.3])
    ylim([-0.3,0.3])
    zlim([-0.015,0.015])
    view(10,15)
    set(gca,'FontSize',12)
    
    % SSM and trajectory in the plane of the delayed and actual states
    subplot(1,2,2)
    hold off
    plot(xt1(outputind,1:end-r),xt1(outputind,r+1:end),'Color',mycolor.green,'LineWidth',2)
    hold on
    plot(xt2(outputind,1:end-r),xt2(outputind,r+1:end),'r','LineWidth',2)
    plot(Wmintauvec1,W0vec1,'--','Color',mycolor.brown,'LineWidth',2)
    plot(Wmintauvec2,W0vec2,'--','Color',mycolor.brown,'LineWidth',2)
    plot(0,0,'.', 'color', 'g','MarkerSize',20)
    plot([-1,1],[-1,1],'rx','MarkerSize',14,'LineWidth',2)
    legend('Traj. \#1', 'Traj. \#2', 'SSM','Location', 'northwest')
    xlabel('$x(t-\tau)$')
    ylabel('$x(t)$')
    title(['$\tau=$' num2str(par.tau)])
    xlim([-1.1,1.1])
    ylim([-1.1,1.1])
    set(gca,'FontSize',12)
    grid on

    %% Check for the existence of a limit cycle and visualize it
    [rho,om] = SSM_per(lambda,beta); % rho: amplitude in the reduced dynamics, om: angular frequency of the limit cycle
    tvec = linspace(0,2*pi/om,100);         % time vector for the limit cycle
    
    % limit cycle on the SSM
    zpervec_t = rho*exp(1i*om*tvec);                    % actual state (t) 

    % limit cycle in the original coordinates 
    Wper_t = W_eval_complex(W0,real(zpervec_t),imag(zpervec_t),outputind);               % actual state (t) 
    Wper_tmintau = W_eval_complex(Wmintau,real(zpervec_t),imag(zpervec_t),outputind);    % delayed state (t-tau) 
    
    % Add the approximated limit cycle to the figure
    figure(3)
    subplot(1,2,2)
    plot(Wper_tmintau,Wper_t,'--','Color',mycolor.orange,'LineWidth',2,'DisplayName','SSM limit cycle')
    
end
