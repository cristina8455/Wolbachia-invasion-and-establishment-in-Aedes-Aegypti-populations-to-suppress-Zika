% same as 6b but with 0.42 instead of 0.042 as decrease of vector
% compentence
% Different initial values
% biting rate changed from 0.5 to highest 1.5
% prob of transmission changed from 0.4 to highest 0.75
% From 6a to 6b: changing Awf from 2,500,000 to 500,000 and Awi from
% 500,000 to 1,500,000 we get Wolbachia inf persists and disease dies
% initial values
%initial values
x0 = [8000000,... % Sh[0]
        1000000,... % Eh[0]
        900000,... % Ih[0]
        100000,... % Rh[0]
        500000,... % Awf[0]
        250000,... % Swf[0]
        50000,... % Iwf[0]
        250000,... % Mwf[0]
        1500000,... % Awi[0]
        250000,... % Swi[0]
        0,... % Iwi[0]
        250000]; % Mwi[0] 
  
tspan = [0, 1095];
ops = odeset('OutputFcn',@odetpbar);
[t, x] = ode45(@SEIR, tspan, x0, ops);
I_h = x(length(t),3); % Ih[final];
I_wf = x(length(t),7); % Iwf[final]
I_wi = x(length(t),11); % Iwi[final]
S_wf = x(length(t),6); %Swf[final]
S_wi = x(length(t),10); %Swi[final] 
[dx, params, R, R_w, M, R_0z] = SEIR(t,x);

% plot settings
legend_sz = 54;
title_sz = 54;
axis_sz = 40;
axis_tick_sz = 24;
line_wdth = 2;
this_file = mfilename;
%save_folder = 'Z:\Matlab\graphics\';
save_folder = pwd;
colors = ["#0072BD" "#D95319" "#7E2F8E" "#EDB120"];

% plotting fcn args - see explanations at bottom
% plot_subplot_panel(t, x, idx, x_vars, x_lim, x_cols, ...,
%    line_wdth, legend_sz, axis_tick_sz, title_str, title_sz, title_pos, legend_text)
% plot_inset_panel(t, x, panel_pos, x_vars, x_lim, x_cols, ...,
%   y_lim, line_wdth, legend_sz, axis_tick_sz, title_str, title_sz, title_pos, legend_text)


figure('units','normalized','position',[0 0 1 0.6])
setappdata(gcf, 'SubplotDefaultAxesLocation', [0, 0, 1, 1]);
set(gcf, 'PaperUnits', 'Inches', 'PaperSize', [6.29, 1.65])

plot_subplot_panel(t, x, [1,3,1], [1,2,4], [0 1095], colors, line_wdth, legend_sz, axis_tick_sz,...
    'Humans', title_sz, [550, 1.1*10^7, 10], ["$S_h$", "$E_h$", "$R_h$"])

plot_inset_panel(t, x, [0.18 0.25 0.13 0.25], 3, [0 1095], colors(4), [0 1*10^5], line_wdth, legend_sz, axis_tick_sz,...
    '', title_sz, [550, 8.5*10^4, 10], "$I_h$")

plot_subplot_panel(t, x, [1,3,2], [5, 6], [0 1095], colors(1:2), line_wdth, legend_sz, axis_tick_sz,...
    '\emph{Wolb}-free mosq.', title_sz, [550, 8.4*10^6, 10], ["$A_{wf}$", "$S_{wf}$"])

plot_inset_panel(t, x, [0.49 0.25 0.13 0.25], 7, [0 1095], colors(4), [0 1.5*10^6], line_wdth, legend_sz, axis_tick_sz,...
    '', title_sz, [550, 1*10^6, 10], "$I_{wf}$")

plot_subplot_panel(t, x, [1,3,3], [9, 10], [0 1095], colors(1:2), line_wdth, legend_sz, axis_tick_sz,...
    '\emph{Wolb}-inf. mosq.', title_sz, [550, 2.8*10^9, 10], ["$A_{wi}$", "$S_{wi}$"])

plot_inset_panel(t, x, [0.82 0.25 0.13 0.25], 11, [0 1095], colors(4), [0 1.75*10^8], line_wdth, legend_sz, axis_tick_sz,...
    '', title_sz, [550, 1*10^6, 10], "$I_{wi}$")


set(gcf,'color','w');
%export_fig(strcat(save_folder,'\',this_file, '_main.pdf'))%, '-q101', '-m2')
%exportgraphics(gcf,strcat(save_folder,'\',this_file, '_main_ML.pdf'),'ContentType','vector')
cleanfigure()
matlab2tikz('myfile.tex'); 

function [dx, params, R, R_w, M, R_0z] = SEIR(t,x)
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
    gamma_h = 1/5;%changed from 1/5
    beta_vh = 1.25*0.75; %biting rate changed from 0.5 to highest 1.25=.9375
    beta_hh = 0.05;
    beta_hv = 1.25*0.75; %prob of transmission changed from 0.4 to highest 0.75
    beta_hvw = 0.42*1.25*0.75; %increase it to .39375 from .039375;
    beta_vhw = 0.42*1.25*0.75;
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
    C = 2;
    q = (1/(1+C))*(1 - (eta*M + (R - R_w))/(eta*R*M))*(R_z_wf/(1-1/R) + C*R_z_wi/(1 - 1/M));
    R_0z = R_d*(2/(3*(9*q + sqrt(3)*sqrt(27*q^2 - 4*R_d^3))))^(1/3) + ((9*q+sqrt(3)*sqrt(27*q^2 - 4*R_d^3))/18)^(1/3);
    
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

% Plotting notes:
% For all functionality, need to install two addons:
% export_fig from: https://www.mathworks.com/matlabcentral/fileexchange/23629-export_fig
% subplot_tight from: https://www.mathworks.com/matlabcentral/fileexchange/30884-controllable-tight-subplot


% In the functions plot_inset_panel() and plot_subplot_panel(),
% t = ind var, x = dep var
% idx = location parameters for subplot/subplot_tight
% panel_pos = location for inset panel
%     essentially it is [horizontal loc, vertical loc, width, height],
%     expressed as percentages of the overall plot/figure size
% x_vars = which components of x to plot
% x_lim = obvious, x_cols = colors to use for the xs
% line_wdth, legend_sz, axis_tick_sz = size params for these
% title_str = string for the title
% legend_text = string array for legend entries


