% Different initial values
% biting rate changed from 0.5 to highest 1.5
% prob of transmission changed from 0.4 to highest 0.75
% initial values
x0 = [8000000,... % Sh[0]
        1000000,... % Eh[0]
        900000,... % Ih[0]
        100000,... % Rh[0]
        2500000,... % Awf[0]
        250000,... % Swf[0]
        50000,... % Iwf[0]
        250000,... % Mwf[0]
        500000,... % Awi[0]
        250000,... % Swi[0]
        0,... % Iwi[0]
        250000]; % Mwi[0] 
  
tspan = [0, 600];
ops = odeset('OutputFcn',@odetpbar);
[t, x] = ode45(@SEIR, tspan, x0, ops);
I_h = x(length(t),3); % Ih[final];
I_wf = x(length(t),7); % Iwf[final]
I_wi = x(length(t),11); % Iwi[final]
S_wf = x(length(t),6); %Swf[final]
S_wi = x(length(t),10); %Swi[final] 
[dx, params, R, R_w, M, R_0z] = SEIR(t,x);
this_file = mfilename;
%save_folder = 'Z:\Matlab\graphics\';
save_folder = pwd;


figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,3,1)
plot(t,x(:,1:4),'LineWidth',2)
legend('Sh', 'Eh', 'Ih', 'Rh')
title('Humans')
subplot(1,3,2)
plot(t,x(:,5:7),'LineWidth',2)
legend('Awf', 'Swf', 'Iwf')
title('Wolbachia free mosq.')
subplot(1,3,3)
plot(t,x(:,9:11),'LineWidth',2)
legend('Awi', 'Swi', 'Iwi')
title('Wolbachia infected mosq.')
xlabel('Time in days')
ylabel('Population') 
cleanfigure()
matlab2tikz(strcat(save_folder,'\',this_file, '_main.tex')); 

figure('DefaultAxesFontSize',26)
plot(t,x(:,3),'Color',[0.9290, 0.6940, 0.1250],'LineWidth',2)
legend('Ih')
title('Inf humans')
xlim([0 600])
ylim([0 1*10^5])
cleanfigure()
matlab2tikz(strcat(save_folder,'\',this_file, '_ih.tex')); 

figure('DefaultAxesFontSize',26)
plot(t,x(:,7),'Color',[0.9290, 0.6940, 0.1250],'LineWidth',2)
legend('Iwf')
title('Inf Wolb. free')
xlim([0 600])
ylim([0 5*10^6])    
cleanfigure()
matlab2tikz(strcat(save_folder,'\',this_file, '_iwf.tex')); 

figure('DefaultAxesFontSize',26)
plot(t,x(:,11),'Color',[0.9290, 0.6940, 0.1250],'LineWidth',2)
legend('Iwf')
title('Inf Wolb. inf')
xlim([0 600])
ylim([0 6*10^6])    
cleanfigure()
matlab2tikz(strcat(save_folder,'\',this_file, '_iwi.tex')); 

function [dx, params, R, R_w, M, R_0z] = SEIR(~,x)
Sh = x(1);
Eh = x(2);
Ih = x(3);
Rh = x(4);
Awf = x(5);
Swf = x(6);
Iwf = x(7);
Mwf = x(8);
Awi = x(9);
Swi = x(10);
Iwi = x(11);
Mwi = x(12);

% parameters
Lambda = 20000000/(78.8*365);
mu_h = 1/(78.8*365);
nu_h = 1/10;
gamma_h = 1/5;%changed from 
beta_vh = 1.25*0.75; %biting rate changed from 0.5 to highest 1.5
%beta_hh = 0.05;
beta_hh = 0;
beta_hv = 1.25*0.75; %prob of transmission changed from 0.4 to highest 0.75
beta_hvw = 0.042*1.25*0.75;
beta_vhw = 0.042*1.25*0.75;
k = 10^9;
mu_v = 0.061;
mu_vi = 0.068;
mu_Ai = 0.2;
gamma_wf = 0.11;
gamma_wi = 0.11;
mu_A = 0.02;
alpha = 0.5;
eta = 13;
q1 = 11;
q2 = 11;

params = [Lambda, mu_h, nu_h, gamma_h, beta_vh,...
 beta_hh, beta_hv, beta_hvw, beta_vhw,...
 k, mu_v, mu_vi, mu_Ai, gamma_wf, gamma_wi,...
 mu_A, alpha, eta, q1, q2];

R = (eta*alpha*(1-alpha)*gamma_wf)/(mu_v*(gamma_wf + mu_A));
oldC = ((gamma_wi + mu_Ai)/(gamma_wf + mu_A) - (q2*gamma_wi*mu_v)/(eta*gamma_wf*mu_vi)) / (q1*(gamma_wi/gamma_wf)^2 * (mu_v/mu_vi)^2);
A_wfstar = (1 -  (mu_v*(gamma_wf + mu_A)*(gamma_wf * mu_vi + oldC*gamma_wi*mu_v)) / (eta*alpha*(1-alpha)*gamma_wf^2*mu_vi))  *  (k/(1 + oldC));
M_wfstar = ((1-alpha)*gamma_wf/mu_v)*A_wfstar;
R_wi = (alpha*gamma_wi*q2*M_wfstar*(1 - A_wfstar/k))/(gamma_wf*A_wfstar*(gamma_wi + mu_Ai));
R_w = (q2*alpha*(1-alpha)*gamma_wi)/(mu_vi*(gamma_wi + mu_Ai));
R_d = (nu_h*beta_hh)/((mu_h + nu_h)*(gamma_h + mu_h));
R_z_wf = (alpha*beta_vh*beta_hv*gamma_wf*mu_h*nu_h*k*(1-1/R))/(Lambda*mu_vi*(mu_h + nu_h)*(gamma_h + mu_h));
R_z = R_z_wf + R_d;
M = (q1*alpha*(1-alpha)*gamma_wi)/(mu_vi*(gamma_wi + mu_Ai));
R_z_wi = (alpha*beta_vhw*beta_hvw*gamma_wi*mu_h*nu_h*k*(1-1/M))/(Lambda*mu_vi^2*(mu_h + nu_h)*(gamma_h + mu_h));
R_zi = R_z_wi + R_d;
C0 = 2;
C=(1-R_w/R)*(gamma_wf*mu_vi*R)/(gamma_wi*mu_v*M);
q = 1/(1+C)*(1-1/R-1/M*(1-R_w/R))*(R_z_wf/(1-1/R) + C*R_z_wi/(1 - 1/M));
R_0z = (R_d+sqrt(R_d*R_d+4*q))/2;

dSh = Lambda - (beta_vh*Iwf + beta_vhw*Iwi +beta_hh*Ih)*Sh*(mu_h/Lambda) -  mu_h*Sh;
dEh = (beta_vh*Iwf + beta_vhw*Iwi + beta_hh*Ih)*Sh*(mu_h/Lambda) - (nu_h+mu_h)*Eh;
dIh = nu_h*Eh - (gamma_h + mu_h)*Ih;
dRh = gamma_h*Ih - mu_h*Rh;
dAwf = eta*((Swf*Mwf)/(Swf + Iwf + Mwf + Swi + Iwi + Mwi))*(1 - (Awi + Awf)/k) - (gamma_wf + mu_A)*Awf;
dSwf = alpha*gamma_wf*Awf - beta_hv*Ih*(mu_h/Lambda)*Swf - mu_v*Swf;
dIwf = beta_hv*Ih*(mu_h/Lambda)*Swf - mu_v*Iwf;
dMwf = (1-alpha)*gamma_wf*Awf - mu_v*Mwf;
dAwi = Swi*(q1*Mwi + q2*Mwf)/(Swf + Iwf + Mwf + Swi + Iwi + Mwi)*(1-(Awf + Awi)/k) - (gamma_wi + mu_Ai)*Awi;
dSwi = alpha*gamma_wi*Awi - beta_hvw*Ih*(mu_h/Lambda)*Swi- mu_v*Swi;
dIwi = beta_hvw*Ih*(mu_h/Lambda)*Swi - mu_vi*Iwi;
dMwi = (1-alpha)*gamma_wi*Awi - mu_vi*Mwi;

dx = [dSh; dEh; dIh; dRh; dAwf; dSwf; dIwf; dMwf; dAwi; dSwi; dIwi; dMwi];
end



