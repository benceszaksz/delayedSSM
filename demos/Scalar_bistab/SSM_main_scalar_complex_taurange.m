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
tauvec = linspace(1.53,1.6,10);
par.tau = tauvec(1);

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

for itau = 1 : length(tauvec)
    disp([num2str(itau) '/' num2str(length(tauvec))])
    par.tau = tauvec(itau); % time delay

%% Spectrum of the linearized dynamics
% determine the dominant characteristic root with the semi-discretization method
Neig = 6;
Ndisc = 30;
charexp = charexp_semidisc(L,R,par,Neig,'disc_num',Ndisc,'tol',1e-10);

% Dominant eigenvalue
lambda = charexp(1);

% plot the right-most eigenvalues
% figure(1);
% plot(real(charexp),imag(charexp),'x')
% xlim([-2.5,1])
% xlabel('Re$\lambda$')
% ylabel('Im$\lambda$')
% set(gca,"FontSize",12)

if isreal(lambda)
    warning('The dominant root is real. Only, the coefficients of the SSM will be calculated. For the visualization, see the file: SSM_main_scalar_real.m')
end

%% Calculation of the SSM corresponding to lambda
[W,beta,p,nW] = SSM_coeff(L,R,fNL,par.tau,lambda,SSMorder,'nlin_type',nlin_type,'red_dyn_style',red_dyn_style);

% Take the coefficients of the SSM at theta=0 and at theta=-tau
W0 = Wsubstheta(W,0);
Wmintau = Wsubstheta(W,-par.tau);


    %% Check for the existence of a limit cycle and visualize it
    [rho,om] = SSM_per(lambda,beta); % rho: amplitude in the reduced dynamics, om: angular frequency of the limit cycle
    if ~isempty(rho)
        rho = flip(rho);
        om = flip(rho);
        Nrepmat = 10; % to obtain more periods of the limit cycle
        tvec1 = linspace(0,Nrepmat*2*pi/om(1),Nrepmat*100);         % time vector for the limit cycle
        % limit cycle on the SSM
        zpervec_t1 = rho(1).*exp(1i*om(1)*tvec1);                    % actual state (t) 
        % limit cycle in the original coordinates 
        Wper_t1 = W_eval_complex(W0,real(zpervec_t1),imag(zpervec_t1),outputind);               % actual state (t) 
        Wper_tmintau1 = W_eval_complex(Wmintau,real(zpervec_t1),imag(zpervec_t1),outputind);    % delayed state (t-tau) 
    
    
        Y_SSM = fft(Wper_t1);
        L_FFT = length(Y_SSM);
        f_SSM = (0:(L_FFT/2));
        P2 = abs(Y_SSM/L_FFT);
        P1_SSM = P2(1:L_FFT/2+1);
        P1_SSM(2:end-1) = 2*P1_SSM(2:end-1);
        FFT_om_SSM1(itau) = P1_SSM(Nrepmat+1);
        FFT_3om_SSM1(itau) = P1_SSM(3*Nrepmat+1);
        FFT_5om_SSM1(itau) = P1_SSM(5*Nrepmat+1);
    
        if size(rho,1)>1
            tvec2 = linspace(0,Nrepmat*2*pi/om(2),Nrepmat*100);         % time vector for the limit cycle
        
            zpervec_t2 = rho(2).*exp(1i*om(2)*tvec2);                    % actual state (t) 
        
            Wper_t2 = W_eval_complex(W0,real(zpervec_t2),imag(zpervec_t2),outputind);               % actual state (t) 
            Wper_tmintau2 = W_eval_complex(Wmintau,real(zpervec_t2),imag(zpervec_t2),outputind);    % delayed state (t-tau) 
    
    
            Y_SSM = fft(Wper_t2);
            L_FFT = length(Y_SSM);
            f_SSM = (0:(L_FFT/2));
            P2 = abs(Y_SSM/L_FFT);
            P1_SSM = P2(1:L_FFT/2+1);
            P1_SSM(2:end-1) = 2*P1_SSM(2:end-1);
            FFT_om_SSM2(itau) = P1_SSM(Nrepmat+1);
            FFT_3om_SSM2(itau) = P1_SSM(3*Nrepmat+1);
            FFT_5om_SSM2(itau) = P1_SSM(5*Nrepmat+1);
        else
            FFT_om_SSM2(itau) = nan;
            FFT_3om_SSM2(itau) = nan;
            FFT_5om_SSM2(itau) = nan;
        end         
    else
            FFT_om_SSM1(itau) = nan;
            FFT_3om_SSM1(itau) = nan;
            FFT_5om_SSM1(itau) = nan;
            FFT_om_SSM2(itau) = nan;
            FFT_3om_SSM2(itau) = nan;
            FFT_5om_SSM2(itau) = nan;
    end

end

figure(1);
subplot(1,3,1)
plot(tauvec,FFT_om_SSM1,'--','Color',mycolor.green,'LineWidth',2)
hold on
plot(tauvec,FFT_om_SSM2,'r--','LineWidth',2)
xlabel('$\tau$')
ylabel('amplitude of base harmonic ($\omega t$)')
%xticks([1.2,1.3,1.4,1.5,1.6])
set(gca,"FontSize",14)
grid on
subplot(1,3,2)
plot(tauvec,FFT_3om_SSM1,'--','Color',mycolor.green,'LineWidth',2)
hold on
plot(tauvec,FFT_3om_SSM2,'r--','LineWidth',2)
xlabel('$\tau$')
ylabel('amplitude of third harmonic ($3\omega t$)')
%xticks([1.2,1.3,1.4,1.5,1.6])
set(gca,"FontSize",14)
grid on
subplot(1,3,3)
plot(tauvec,FFT_5om_SSM1,'--','Color',mycolor.green,'LineWidth',2)
hold on
plot(tauvec,FFT_5om_SSM2,'r--','LineWidth',2)
xlabel('$\tau$')
ylabel('amplitude of fifth harmonic ($5\omega t$)')
%xticks([1.2,1.3,1.4,1.5,1.6])
set(gca,"FontSize",14)
grid on