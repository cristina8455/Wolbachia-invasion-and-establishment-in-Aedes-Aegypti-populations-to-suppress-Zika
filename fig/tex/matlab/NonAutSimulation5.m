% Awi changed from 500,000 to 1,000,000 and Mwi changed from 250,000 to
% 1,100,000
%initial values
x0 = [10000000,... % Sh[0]
        3000000,... % Eh[0]
        2500000,... % Ih[0]
        10000,... % Rh[0]
        500000,... % Awf[0]
        250000,... % Swf[0]
        50000,... % Iwf[0]
        250000,... % Mwf[0]
        1000000,... % Awi[0]
        250000,... % Swi[0]
        50000,... % Iwi[0]
        1100000]; % Mwi[0]
  
tspan = [0, 1095];
ops = odeset('OutputFcn',@odetpbar);
[t, x] = ode45(@SEIR, tspan, x0, ops);
I_h = x(length(t),3); % Ih[final]
I_wf = x(length(t),7); % Iwf[final]
I_wi = x(length(t),11); % Iwi[final]
[dx, params, R, R_w, R_z] = SEIR(t,x);


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

str = {...
    %['$\Lambda = $' num2str(params(1))],...
    %['$\mu_h = $' num2str(params(2))],...
    %['$\nu_h = $' num2str(params(3))],...
    %['$\gamma_h = $' num2str(params(4))],...
    %['$\beta_{vh} = $' num2str(params(5))],...
    %['$\beta_{hh} = $' num2str(params(6))],...
    %['$\beta_{hv} = $' num2str(params(7))],...
    %['$\beta_{hvw} = $' num2str(params(8))],...
    %['$\beta_{vhw} = $' num2str(params(9))],...
    %['$k = $' num2str(params(10))],...
    %['$\mu_v = $' num2str(params(11))],...
    %['$\mu_{vi} = $' num2str(params(12))],...
    %['$\mu_{Ai} = $' num2str(params(13))],...
    %['$\gamma_{wf} = $' num2str(params(14))],...
    %['$\gamma_{wi} = $' num2str(params(15))],...
    %['$\mu_A = $' num2str(params(16))],...
    %['$\alpha = $' num2str(params(17))],...
    %['$\eta = $' num2str(params(18))],...
    %['$q_1 = $' num2str(params(19))],...
    %['$q_2 = $' num2str(params(20))],...
    ['$R = $' num2str(R)],...
    ['$R_{w}/R = $' num2str(R_w/R)],...
    ['$R_{z} = $' num2str(R_z)]};
    %['$I_h = $' num2str(I_h)],...
    %['$I_{wf} = $' num2str(I_wf)],...
    %['$I_{wi} = $' num2str(I_wi)]};
annotation('textbox',[0.8, 0.7, 0.1, 0.1],'interpreter','latex','String',str,'FitBoxToText','on');
    
function [dx, params, R, R_w, R_z] = SEIR(t,x)
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
Lambda = 21500000/(78.8*365);
mu_h = 1/(78.8*365);
nu_h = 1/10;
gamma_h = 1/5;
beta_vh = 0.5*0.4;
beta_hh = 0.05;
beta_hv = 0.5*0.4;
beta_hvw = 0.042*0.5*0.4;
beta_vhw = 0.042*0.5*0.4;
k = 10^6;
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
C = ((gamma_wi + mu_Ai)/(gamma_wf + mu_A) - (q2*gamma_wi*mu_v)/(eta*gamma_wf*mu_vi)) / (q1*(gamma_wi/gamma_wf)^2 * (mu_v/mu_vi)^2);
A_wfstar = (1 -  (mu_v*(gamma_wf + mu_A)*(gamma_wf * mu_vi + C*gamma_wi*mu_v)) / (eta*alpha*(1-alpha)*gamma_wf^2*mu_vi))  *  (k/(1 + C));
M_wfstar = ((1-alpha)*gamma_wf/mu_v)*A_wfstar;
R_wi = (alpha*gamma_wi*q2*M_wfstar*(1 - A_wfstar/k))/(gamma_wf*A_wfstar*(gamma_wi + mu_Ai));
R_w = (q2*alpha*(1-alpha)*gamma_wi)/(mu_vi*(gamma_wi + mu_Ai));
R_z = (alpha*beta_vh*beta_hv*gamma_wf*mu_h*nu_h*k*(1-1/R))/(Lambda*mu_vi*(mu_h + nu_h)*(gamma_h + mu_h)) + (nu_h*beta_hh)/((mu_h + nu_h)*(gamma_h + mu_h));
M = (q1*alpha*(1-alpha)*gamma_wi)/(mu_vi*(gamma_wi + mu_Ai))

dSh = Lambda - (beta_vh*Iwf + beta_vhw*Iwi +beta_hh*Ih)*Sh*(mu_h/Lambda) -  mu_h*Sh;
dEh = (beta_vh*Iwf + beta_vhw*Iwi + beta_hh*Ih)*Sh*(mu_h/Lambda) - (nu_h+mu_h)*Eh;
dIh = nu_h*Eh - (gamma_h + mu_h)*Ih;
dRh = gamma_h*Ih - mu_h*Rh;
dAwf = (11*sin(0.0172*t-91*0.0172)+13)*((Swf*Mwf)/(Swf + Iwf + Mwf + Swi + Iwi + Mwi))*(1 - (Awi + Awf)/k) - (gamma_wf + mu_A)*Awf;
dSwf = alpha*gamma_wf*Awf - beta_hv*Ih*(mu_h/Lambda)*Swf - (0.035*sin(0.0172*t-91*1.0172)+0.061)*Swf;
dIwf = beta_hv*Ih*(mu_h/Lambda)*Swf - (0.035*sin(0.0172*t-91*1.0172)+0.061)*Iwf;
dMwf = (1-alpha)*gamma_wf*Awf - (0.035*sin(0.0172*t-91*1.0172)+0.061)*Mwf;
dAwi = Swi*((9*sin(0.0172*t-91*0.0172)+11)*Mwi + (9*sin(0.0172*t-91*0.0172)+11)*Mwf)/(Swf + Iwf + Mwf + Swi + Iwi + Mwi)*(1-(Awf + Awi)/k) - (gamma_wi + mu_Ai)*Awi;
dSwi = alpha*gamma_wi*Awi - beta_hvw*Ih*(mu_h/Lambda)*Swi- (0.055*sin(0.0172*t-91*1.0172)+0.068)*Swi;
dIwi = beta_hvw*Ih*(mu_h/Lambda)*Swi - (0.055*sin(0.0172*t-91*1.0172)+0.068)*Iwi;
dMwi = (1-alpha)*gamma_wi*Awi - (0.055*sin(0.0172*t-91*1.0172)+0.068)*Mwi;

dx = [dSh; dEh; dIh; dRh; dAwf; dSwf; dIwf; dMwf; dAwi; dSwi; dIwi; dMwi];
end



