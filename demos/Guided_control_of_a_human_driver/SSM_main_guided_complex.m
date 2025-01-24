% I am glad you are using my code. I will be happy to receive your 
% feedback and help you if needed. If necessary, please contact me at
% szaksz@mm.bme.hu.
%
% This code is an example for the calculation of the spectral submanifold 
% (SSM) in a 3-dimensional nonlinear time delay system. More specifically, 
% it is foxusing to the case of complex conjugate dominant eigenvalues. 
% For the case of a real eigenvalue see the code: SSM_main_guided_real.m
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

% parameters for the automated vehicle leading to a pair of complex conjugate roots
par.betamin1 = -0.8;
par.betahat = 0.8;

% time delay
par.tau = 0.4;

% type of nonlinearity (may be 'non-delayed', 'delayed' or 'combined')
nlin_type = 'delayed';

% style of reduced dynamics 'graph', 'normal_form' (only in case of complex conjugate roots), 'manual'
red_dyn_style = 'normal_form';

% order of the Spectral submanifold and the corresponding reduced dynamics
SSMorder = 3;

% the index of the state vector for visualization
outputind = 2;

% initial function for a particular trajectory (determined with dde23)
initf = @(t) [2;0;0];

% length of the dde23 simulation
T = 100; 

% surface coloring needed for 3D visualization
map = interp1([0,0.5,1],[[0.6275    0.3216    0.1765]/2;[205,133,63]/255;1 1 1],[0:0.01:1]);
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

if isreal(lambda)
    warning('The dominant root is real. Only, the coefficients of the SSM will be calculated. For the visualization, see the file: SSM_main_guided_real.m')
end

%% Calculation of the SSM corresponding to lambda
[W,beta,p,nW] = SSM_coeff(L,R,fNL,par.tau,lambda,SSMorder,'nlin_type',nlin_type,'red_dyn_style',red_dyn_style);

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
if imag(lambda)~=0 % the dominant eigenvalues form a complex conjugate pair
    % horizontal grid
    [Rez,Imz] = meshgrid(-0.6:0.01:0.6,-0.6:0.01:0.6);

    % SSM at the desired coorinate (outputind)
    [Wsurface,wsurface_xi,Rexi_grid,Imxi_grid] = W_eval_complex(W0,Rez,Imz,outputind,nW);

    % transformed trajectory
    xt_transf = xt(:,r+1:end)-W0.W_1_0*xi-W0.W_0_1*conj(xi);

    %% One may either project the trajectory to the manifold and visualize
    % that, or project only its initial function to the manifold, and
    % obtain a trajectory on the SSM with the corresponding reduced
    % dynamics. In the original version of the code, the latter one is
    % used.
    
    % projection to manifold:
%     W0vec = W_eval_complex(W0,real(xi),imag(xi),outputind); % projection to the SSM at theta=0
%     Wmintauvec = W_eval_complex(Wmintau,real(xi),imag(xi),outputind); %  projection to the SSM at theta=-tau

    % reduced dynamics on the manifold:
    z0 = xi(1);
    red_dyn_act = red_dyn(lambda,beta);
    solSSM = ode45(@(t,z) red_dyn_act(z), [0,T],z0);
    tvec_SSM = 0:dt:solSSM.x(end);
    zsol = deval(solSSM,tvec_SSM);
    W0vec = W_eval_complex(W0,real(zsol),imag(zsol),outputind); % reduced dynamics on the SSM at theta=0
    Wmintauvec = W_eval_complex(Wmintau,real(zsol),imag(zsol),outputind); % reduced dynamics on the SSM at theta=-tau

           
    %% Generate the figure
    figure;
    % SSM above the plane of Re(xi) and Im(xi), and the transformed
    % trajectory
    subplot(1,3,1)
    hold off
    sl = surfl(Rexi_grid,Imxi_grid,wsurface_xi,[-20 20], [.7 .8 0.1 10]);
    sl.FaceAlpha = 0.8;
    view(-60,18)
    colormap(map)
    shading interp
    xlabel('Re($y(z,\overline{z})$)')
    ylabel('Im($y(z,\overline{z})$)')
    zlabel('$w_2(z,\overline{z};0)$ [m/s]') % the 2nd component of w if outputind=2
    hold on
    plot3(real(xi),imag(xi),xt_transf(outputind,:),'Color',mycolor.green,'LineWidth',2)
    plot3(0,0,0,'.', 'color', 'g','MarkerSize',20)
    set(gca,'FontSize',14)
    xlim([-0.6,0.6])
    ylim([-0.6,0.6])
    zlim([-0.3,0.05])
    
    % the same plot as the previous one, just enlarged
    subplot(1,3,2)
    hold off
    sl = surfl(Rexi_grid,Imxi_grid,wsurface_xi,[-20 20], [.7 .8 0.1 10]);
    sl.FaceAlpha = 0.8;
    view(-60,18)
    colormap(map)
    shading interp
    xlabel('Re($y(z,\overline{z})$)')
    ylabel('Im($y(z,\overline{z})$)')
    zlabel('$w_2(z,\overline{z};0)$ [m/s]') % the 2nd component of w if outputind=2
    hold on
    plot3(real(xi),imag(xi),xt_transf(outputind,:),'Color',mycolor.green,'LineWidth',2)
    plot3(0,0,0,'.', 'color', 'g','MarkerSize',20)
    set(gca,'FontSize',14)
    xlim([-0.29,0.29])
    ylim([-0.29,0.29])
    zlim([-0.004,0.002])
    
    % SSM and trajectory in the plane of the delayed and actual states
    subplot(1,3,3)
    hold off
    plot(xt(outputind,1:end-r),xt(outputind,r+1:end),'Color',mycolor.green,'LineWidth',2)
    hold on; 
    plot(Wmintauvec,W0vec,'--','Color',mycolor.brown,'LineWidth',2)
    plot(0,0,'.', 'color', 'g','MarkerSize',20)
    legend('SSM','Traj.','Location','NorthWest')
    xlabel('$v_{-1}(t-\tau)$ [m/s]')
    ylabel('$v_{-1}(t)$ [m/s]')
    set(gca,'FontSize',14)
    xlim([-0.6,0.6])
    ylim([-0.6,0.6])
    grid on
end