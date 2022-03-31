x=[0:0.1: 365];

% plot settings
legend_sz = 26;
title_sz = 32;
axis_sz = 32;
axis_tick_sz = 18;
line_wdth = 3;
this_file = mfilename;
%save_folder = 'Z:\Matlab\graphics\';
save_folder = pwd;

%seasonlity of eta - egg laying rate of wild mosquitoes
y1=11*sin(0.0172*x-91*0.0172)+13;


%seasonlity of q1 and q2 - egg laying rate of Wolbachia mosquitoes
y2=9*sin(0.0172*x-91*0.0172)+11;
figure()

plot(x,y1,x,y2,'LineWidth',line_wdth);
xlim([0 365])
ax = gca;
ax.YAxis.FontSize = axis_tick_sz;
ax.XAxis.FontSize = axis_tick_sz;
legend('$\eta(t)$','$q_1(t),q_2(t)$')
set(legend,'FontSize',legend_sz,'interpreter','latex')
title('Birth Rates', 'interpreter', 'latex','FontSize', title_sz);
%ylabel({'Rate','of mosquitoes'},'interpreter', 'latex','FontSize',
%axis_sz) % multiline example
ylabel('Rate','interpreter', 'latex','FontSize',axis_sz);
xlabel('Time (days)','interpreter', 'latex','FontSize', axis_sz);
set(gcf,'color','w');
export_fig(strcat(save_folder,'\',this_file, 'seasonality_birth.png'), '-m0.83')
cleanfigure()
matlab2tikz(strcat(save_folder,'\',this_file, 'seasonality_birth.tex')); 

%seasonlity of muvi - death rate of Wolbachia mosquitoes
y4=0.055*sin(0.0172*x-91*1.0172)+0.068;

figure()
%seasonlity of muv - death rate of wild mosquitoes
y3=0.035*sin(0.0172*x-91*1.0172)+0.061;
plot(x,y3,x,y4,'LineWidth',line_wdth);
xlim([0 365])
ax = gca;
ax.YAxis.FontSize = axis_tick_sz;
ax.XAxis.FontSize = axis_tick_sz;
legend('$\mu_v(t)$','$\mu_{vi}(t)$')
set(legend,'FontSize',legend_sz,'interpreter','latex')
xlabel('Time (days)','interpreter', 'latex','FontSize', axis_sz);
ylabel('Rate', 'interpreter', 'latex', 'FontSize', axis_sz);
title('Death Rates', 'interpreter', 'latex','FontSize', title_sz);
set(gcf,'color','w');
export_fig(strcat(save_folder,'\',this_file, 'seasonality_death.png'), '-m0.83')
cleanfigure()
matlab2tikz(strcat(save_folder,'\',this_file, 'seasonality_death.tex')); 


figure()
%seasonlity of gammawf and wi - transition rates of mosquitoes
y5=0.1*sin(2*pi/365*(x-91))+0.12;
y6=0.1*sin(2*pi/365*(x-91))+0.12;
plot(x,y5,'LineWidth',line_wdth);
xlim([0 365])
ax = gca;
ax.YAxis.FontSize = axis_tick_sz;
ax.XAxis.FontSize = axis_tick_sz;
legend('$\gamma_{wf}(t),\gamma_{wi}(t)$')
set(legend,'FontSize',legend_sz,'interpreter','latex')
xlabel('Time (days)','interpreter', 'latex','FontSize', axis_sz);
ylabel('Rate','interpreter', 'latex','FontSize', axis_sz);
title('Transition Rates', 'interpreter', 'latex','FontSize', title_sz);
set(gcf,'color','w');
export_fig(strcat(save_folder, '\', this_file, 'seasonality_trans.png'), '-m0.83')
cleanfigure()
matlab2tikz(strcat(save_folder,'\',this_file, 'seasonality_trans.tex')); 