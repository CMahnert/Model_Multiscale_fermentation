clear all, clc

load('Result_adjust_W3110.mat');
load('Result_adjust_VAL22.mat');
load('Result_adjust_VAL23.mat');
load('Result_adjust_VAL24.mat');
load('Result_adjut_monod_W3110.mat');
load('Result_adjut_monod_VAL22.mat');
load('Result_adjut_monod_VAL23.mat');
load('Result_adjut_monod_VAL24.mat');

%% Figure
%Biomasa vs tiempo
figure(1),clf
subplot (2,1,1)
%Strain W3110
plot(time_fmincon_W3110/60,Biomasa_obs_fmincon_W3110,'o','LineWidth',2,'Color',[0.4660 0.6740 0.1880])
hold on 
plot(t_fmincon_W3110/60,biomasa_model_fmincon_W3110,'LineWidth',2,'Color',[0.4660 0.6740 0.1880])
xlim([0,16]),ylim([0,4])
%Strain V22
plot(time_fmincon_V22/60,Biomasa_obs_fmincon_V22,'o','LineWidth',2,'Color',[0.9290 0.6940 0.1250])
plot(t_fmincon_V22/60,biomasa_model_fmincon_V22,'LineWidth',2,'Color',[0.9290 0.6940 0.1250])
%Strain V23
plot(time_fmincon_V23/60,Biomasa_obs_fmincon_V23,'o','LineWidth',2,'Color',[0.8500 0.3250 0.0980])
plot(t_fmincon_V23/60,biomasa_model_fmincon_V23,'LineWidth',2,'Color',[0.8500 0.3250 0.0980])
%Strain V24
plot(time_fmincon_V24/60,Biomasa_obs_fmincon_V24,'o','LineWidth',2,'Color',[0 0.4470 0.7410])
plot(t_fmincon_V24/60,biomasa_model_fmincon_V24,'LineWidth',2,'Color',[0 0.4470 0.7410])
hold off
ylabel('Biomass [g/L]'),xlabel('time [h]')
legend('','W3110','','Val22','','Val23','','Val24')
set(gca,'fontsize',30),ylim([0,4])
% 
%Substrate
subplot(2,1,2)
%Strain W3110
plot(time_fmincon_W3110/60,Sustrate_obs_fmincon_W3110,'o','LineWidth',2,'Color',[0.4660 0.6740 0.1880])
hold on 
plot(t_fmincon_W3110/60,sustrate_model_fmincon_W3110,'LineWidth',2,'Color',[0.4660 0.6740 0.1880])
xlim([0,16]),ylim([0,5])
%Strain V22
plot(time_fmincon_V22/60,Sustrate_obs_fmincon_V22,'o','LineWidth',2,'Color',[0.9290 0.6940 0.1250])
plot(t_fmincon_V22/60,sustrate_model_fmincon_V22,'LineWidth',2,'Color',[0.9290 0.6940 0.1250])
%Strain V23
plot(time_fmincon_V23/60,Sustrate_obs_fmincon_V23,'o','LineWidth',2,'Color',[0.8500 0.3250 0.0980])
plot(t_fmincon_V23/60,sustrate_model_fmincon_V23,'LineWidth',2,'Color',[0.8500 0.3250 0.0980])
%Strain V24
plot(time_fmincon_V24/60,Sustrate_obs_fmincon_V24,'o','LineWidth',2,'Color',[0 0.4470 0.7410])
plot(t_fmincon_V24/60,sustrate_model_fmincon_V24,'LineWidth',2,'Color',[0 0.4470 0.7410])
ylabel('Substrate [g/L]'),xlabel('time [h]')
set(gca,'fontsize',30)
ylim([0,6])
hold off
% 
%Producto
figure(2),clf
subplot(2,1,1)
%Strain W3110
plot(time_fmincon_W3110/60,Producto_obs_fmincon_W3110,'o','LineWidth',2,'Color',[0.4660 0.6740 0.1880])
hold on 
plot(t_fmincon_W3110/60,producto_model_fmincon_W3110,'LineWidth',2,'Color',[0.4660 0.6740 0.1880])
xlim([0,16]),ylim([0,0.4])
%Strain V22
plot(time_fmincon_V22/60,Producto_obs_fmincon_V22,'o','LineWidth',2,'Color',[0.9290 0.6940 0.1250])
plot(t_fmincon_V22/60,producto_model_fmincon_V22,'LineWidth',2,'Color',[0.9290 0.6940 0.1250])
%Strain V23
plot(time_fmincon_V23/60,Producto_obs_fmincon_V23,'o','LineWidth',2,'Color',[0.8500 0.3250 0.0980])
plot(t_fmincon_V23/60,producto_model_fmincon_V23,'LineWidth',2,'Color',[0.8500 0.3250 0.0980])
%Strain V24
plot(time_fmincon_V24/60,Producto_obs_fmincon_V24,'o','LineWidth',2,'Color',[0 0.4470 0.7410])
plot(t_fmincon_V24/60,producto_model_fmincon_V24,'LineWidth',2,'Color',[0 0.4470 0.7410])
ylabel('GFP [g/L]'),xlabel('time [h]'),ylim([0,0.3])
set(gca,'fontsize',30),legend('','W3110','','Val22','','Val23','','Val24')
hold off
% 
%Acetate
subplot(2,1,2)
%Strain W3110
plot(time_fmincon_W3110/60,acetato_obs_fmincon_W3110,'o','LineWidth',2,'Color',[0.4660 0.6740 0.1880])
hold on 
plot(t_fmincon_W3110/60,acetato_model_fmincon_W3110,'LineWidth',2,'Color',[0.4660 0.6740 0.1880])
xlim([0,16]),ylim([0,1])
%Strain V22
plot(time_fmincon_V22/60,acetato_obs_fmincon_V22,'o','LineWidth',2,'Color',[0.9290 0.6940 0.1250])
plot(t_fmincon_V22/60,acetato_model_fmincon_V22,'LineWidth',2,'Color',[0.9290 0.6940 0.1250])
%Strain V23
plot(time_fmincon_V23/60,acetato_obs_fmincon_V23,'o','LineWidth',2,'Color',[0.8500 0.3250 0.0980])
plot(t_fmincon_V23/60,acetato_model_fmincon_V23,'LineWidth',2,'Color',[0.8500 0.3250 0.0980])
%Strain V24
plot(time_fmincon_V24/60,acetato_obs_fmincon_V24,'o','LineWidth',2,'Color',[0 0.4470 0.7410])
plot(t_fmincon_V24/60,acetato_model_fmincon_V24,'LineWidth',2,'Color',[0 0.4470 0.7410])
ylabel('Acetate [g/L]'),xlabel('time [h]')
set(gca,'fontsize',30)
hold off
%%

figure(3),clf 
colororder({'k','k'})
subplot(1,4,1)
%strain W3110
plot(time_fmincon_W3110/60,acetato_obs_fmincon_W3110,'o','LineWidth',2,'Color',[0.8500 0.3250 0.0980])
hold on 
plot(t_fmincon_W3110/60,acetato_model_fmincon_W3110,'-','LineWidth',2,'Color',[0.8500 0.3250 0.0980])
xlim([0,16])
ylabel('Acetate[g/L]')
yyaxis right
plot(time_fmincon_W3110/60,Producto_obs_fmincon_W3110,'o','LineWidth',2,'Color',[0 0.4470 0.7410])
plot(t_fmincon_W3110/60,producto_model_fmincon_W3110,'-','LineWidth',2,'Color',[0 0.4470 0.7410])
ylabel('gfp[g/L]')
set(gca,'fontsize',16),legend('','Acetate','','gfp')
%strain V22
subplot(1,4,2)
plot(time_fmincon_V22/60,acetato_obs_fmincon_V22,'o','LineWidth',2,'Color',[0.8500 0.3250 0.0980])
hold on 
plot(t_fmincon_V22/60,acetato_model_fmincon_V22,'-','LineWidth',2,'Color',[0.8500 0.3250 0.0980])
xlim([0,12])
ylabel('Acetate[g/L]')
yyaxis right
plot(time_fmincon_V22/60,Producto_obs_fmincon_V22,'o','LineWidth',2,'Color',[0 0.4470 0.7410])
plot(t_fmincon_V22/60,producto_model_fmincon_V22,'-','LineWidth',2,'Color',[0 0.4470 0.7410])
ylabel('gfp[g/L]'),ylim([0,0.055])
set(gca,'fontsize',30)
%strain V23
subplot(1,4,3)
plot(time_fmincon_V23/60,acetato_obs_fmincon_V23,'o','LineWidth',2,'Color',[0.8500 0.3250 0.0980])
hold on 
plot(t_fmincon_V23/60,acetato_model_fmincon_V23,'-','LineWidth',2,'Color',[0.8500 0.3250 0.0980])
xlim([0,13])
ylabel('Acetate[g/L]')
yyaxis right
plot(time_fmincon_V23/60,Producto_obs_fmincon_V23,'o','LineWidth',2,'Color',[0 0.4470 0.7410])
plot(t_fmincon_V23/60,producto_model_fmincon_V23,'-','LineWidth',2,'Color',[0 0.4470 0.7410])
ylabel('gfp[g/L]')
set(gca,'fontsize',30)
%strain V24
subplot(1,4,4)
plot(time_fmincon_V24/60,acetato_obs_fmincon_V24,'o','LineWidth',2,'Color',[0.8500 0.3250 0.0980])
hold on 
plot(t_fmincon_V24/60,acetato_model_fmincon_V24,'-','LineWidth',2,'Color',[0.8500 0.3250 0.0980])
xlim([0,10])
ylabel('Acetate[g/L]')
yyaxis right
plot(time_fmincon_V24/60,Producto_obs_fmincon_V24,'o','LineWidth',2,'Color',[0 0.4470 0.7410])
plot(t_fmincon_V24/60,producto_model_fmincon_V24,'-','LineWidth',2,'Color',[0 0.4470 0.7410])
ylabel('gfp[g/L]')
set(gca,'fontsize',30)
%% Análsis de J

figure(4),clf 
parnames={'W3110','Val22','Val23','Val24'
    };
Y_bar_Model=[ss0_fmincon_W3110,ss0_Monod_W3110
    ss0_fmincon_V22,ss0_Monod_V22
    ss0_fmincon_V23,ss0_Monod_V23
    ss0_fmincon_V24,ss0_Monod_V24];
b=bar(Y_bar_Model,'FaceColor','flat');
b(1).CData=[0.9290 0.6940 0.1250];
b(2).CData=[0.47,0.67,0.19];
set(gca,'xticklabel',parnames,'FontSize',20)
set(gca,'fontsize',30)
ylabel('Good of fit [J]')
legend('Integrated model','Monod Model')

%% Análisis J stack 
figure(5),clf 
% Y_bar_biomasa=[J_suma_biomasa_fmincon_W3110,J_suma_biomasa_Monod_W3110
%                J_suma_biomasa_fmincon_V22,J_suma_biomasa_Monod_V22
%                J_suma_biomasa_fmincon_V23,J_suma_biomasa_Monod_V23
%                J_suma_biomasa_fmincon_V24,J_suma_biomasa_fmincon_V24];
% 
% Y_bar_sustrato=[J_suma_sustrato_fmincon_W3110,J_suma_sustrato_Monod_W3110
%                J_suma_sustrato_fmincon_V22,J_suma_sustrato_Monod_V22
%                J_suma_sustrato_fmincon_V23,J_suma_sustrato_Monod_V23
%                J_suma_sustrato_fmincon_V24,J_suma_sustrato_fmincon_V24];

Y_bar=[J_suma_biomasa_fmincon_W3110,J_suma_sustrato_fmincon_W3110
       J_suma_biomasa_Monod_W3110,J_suma_sustrato_Monod_W3110;
        J_suma_biomasa_fmincon_V22,J_suma_sustrato_fmincon_V22
        J_suma_biomasa_Monod_V22,J_suma_sustrato_Monod_V22;
        J_suma_biomasa_fmincon_V23,J_suma_sustrato_fmincon_V23
        J_suma_biomasa_Monod_V23,J_suma_sustrato_Monod_V23;
        J_suma_biomasa_fmincon_V24,J_suma_sustrato_fmincon_V24
        J_suma_biomasa_Monod_V24,J_suma_sustrato_Monod_V24];

% Y_bar=rand([J_suma_biomasa_fmincon_W3110,J_suma_sustrato_fmincon_W3110,J_suma_biomasa_Monod_W3110,J_suma_sustrato_Monod_W3110;
%         J_suma_biomasa_fmincon_V22,J_suma_sustrato_fmincon_V22,J_suma_biomasa_Monod_V22,J_suma_sustrato_Monod_V22;
%         J_suma_biomasa_fmincon_V23,J_suma_sustrato_fmincon_V23,J_suma_biomasa_Monod_V23,J_suma_sustrato_Monod_V23;
%         J_suma_biomasa_fmincon_V24,J_suma_sustrato_fmincon_V24, J_suma_biomasa_Monod_V24,J_suma_sustrato_Monod_V24]);
% 
% 
% groupLabels = {1,2};
% plotBarStackGroups(Y_bar, groupLabels)
Yb=bar(Y_bar,'stacked','FaceColor','flat')
legend('Biomass','Substrate')




%% Análisis de sensibilidad
load('Model_efast_W3110.mat')
load('Model_efast_V22.mat')
load('Model_efast_V23.mat')
load('Model_efast_V24.mat')

figure(5),clf 
title('Sensibity')
parnames={'K_{\gamma}','w_p','w_{fer}','w_{res}','n_f','w_t'};
%W3110
subplot(4,1,1)
ba=bar([Si_W3110(:,2), Sti_W3110(:,2)],'FaceColor','flat');
ba(1).CData=[0.9290 0.6940 0.1250];
ba(2).CData=[0.47,0.67,0.19];
set(gca,'xticklabel',parnames,'FontSize',30)
set(gca,'fontsize',30)
ylabel('W3110'),ylim([0,1])
legend('first-order S_i','total-order S_{Ti}')

%V22
subplot(4,1,2)
ba2=bar([Si_V22(:,2),Sti_V22(:,2)],'FaceColor','flat');
ba2(1).CData=[0.9290 0.6940 0.1250];
ba2(2).CData=[0.47,0.67,0.19];
set(gca,'xticklabel',parnames,'FontSize',30)
set(gca,'fontsize',30)
ylabel('Val22'),ylim([0,1])

%V23
subplot(4,1,3)
ba3=bar([Si_V23(:,2),Sti_V23(:,2)],'FaceColor','flat');
ba3(1).CData=[0.9290 0.6940 0.1250];
ba3(2).CData=[0.47,0.67,0.19];
set(gca,'xticklabel',parnames,'FontSize',30)
set(gca,'fontsize',30)
ylabel('Val23'),ylim([0,1])

%V24
subplot(4,1,4)
ba4=bar([Si_V24(:,2),Sti_W3110(:,2)],'FaceColor','flat');
ba4(1).CData=[0.9290 0.6940 0.1250];
ba4(2).CData=[0.47,0.67,0.19];
set(gca,'xticklabel',parnames,'FontSize',30)
set(gca,'fontsize',30)
ylabel('Val24'),ylim([0,1])

%% Spider plot 
figure(6),clf 
%W3110
D1_W3110=[Sti_W3110(1,2) Sti_W3110(2,2) Sti_W3110(3,2) Sti_W3110(4,2) Sti_W3110(5,2) Sti_W3110(6,2) ];
D1_V22=[Sti_V22(1,2) Sti_V22(2,2) Sti_V22(3,2) Sti_V22(4,2) Sti_V22(5,2) Sti_V22(6,2)];
D1_V23=[Sti_V23(1,2) Sti_V23(2,2) Sti_V23(3,2) Sti_V23(4,2) Sti_V23(5,2) Sti_V23(6,2) ];
D1_V24=[Sti_V24(1,2) Sti_V24(2,2) Sti_V24(3,2) Sti_V24(4,2) Sti_V24(5,2) Sti_V24(6,2)];

P_strain=[D1_W3110; D1_V22; D1_V23; D1_V24];
%spider plot 
spider_plot(P_strain,'AxesLimits',[0, 0.6, 0, 0.45, 0, 0.05; 0.2, 0.9, 0.25, 0.8, 0.4, 0.55],'AxesPrecision', [1, 1, 1, 1, 1, 1],'AxesInterval',2,'AxesLabels',{'K_{\gamma}','w_p','w_{ef}','w_{er}','n_f','w_t'},'Color', [0.4660 0.6740 0.1880; 0.9290 0.6940 0.1250;0.8500 0.3250 0.0980;0 0.4470 0.7410] ); 
legend('W3110','Val22','Val23','Val24')

figure(5),clf 
spider_plot(P_strain,'AxesLimits',[0, 0, 0, 0, 0, 0; 1,1, 1, 1, 1, 1],'AxesPrecision', [1, 1, 1, 1, 1, 1],'AxesInterval',3,'AxesLabels',{'K_{\gamma}','w_p','w_{ef}','w_{er}','n_f','w_t'},'Color', [0.4660 0.6740 0.1880; 0.9290 0.6940 0.1250;0.8500 0.3250 0.0980;0 0.4470 0.7410],'AxesDisplay', 'one' ); 
legend('W3110','Val22','Val23','Val24')


