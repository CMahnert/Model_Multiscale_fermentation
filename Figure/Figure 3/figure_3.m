clear all, clc 

load('Result_wp_W3110')

%% Figure
figure(1),clf 
colororder({'k','k'})
yyaxis right
semilogx(wp,biomasa_model_Wp_W3110,'LineWidth',3,'Color',[0 0.4470 0.7410])
xlabel('Wp [molecs/cell/min]'),ylabel('Biomass')
legend('Biomass')
hold on 
yyaxis left
semilogx(wp,producto_model_Wp_W3110,'LineWidth',3,'Color',[0.8500 0.3250 0.0980])
semilogx(wp,acetato_model_max_Wp_W3110,'lineStyle','-','LineWidth',3,'Color',[0.9290 0.6940 0.1250])
xlabel('Wp [molecs/cell/min]'),ylabel('Acetate, GFP [g/L] and growth rate [1/h]')
semilogx(wp,max_lam_Wp_W3110,'lineStyle','-','LineWidth',3,'Color',[0.4660 0.6740 0.1880])
legend('GFP','Acetate','Growth rate','Biomass [g/L]')
hold off 
set(gca,'fontsize',30)

figure(2),clf 
plot(productividad_producto_Wp_W3110*60,yield_growth_rate,'LineWidth',3),xlabel('Protein [g/L]'),ylabel('Growth yield')
set(gca,'fontsize',30),xlabel('gfp Productivity [g/L/h]'),ylabel('growth yield')

% subplot(1,4,3)
% plot(producto_model_Wp_W3110,max_lam_Wp_W3110,'LineWidth',3),xlabel('Protein [g/L]'),ylabel('growth rate [1/h]')
% subplot(1,4,4)
% scatter(max_lam_Wp_W3110,yield_growth_rate,'filled'),xlabel('growth rate [1/h]'),ylabel('Growth Yield')

% hold on 
% max_lam_Wp_W3110_varios_puntos=[max_lam_Wp_W3110(1,1), max_lam_Wp_W3110(1,670) max_lam_Wp_W3110(1,end)];
% yield_growth_rate_varios_puntos=[yield_growth_rate(1,1), yield_growth_rate(1,670),yield_growth_rate(1,end)];
% scatter(max_lam_Wp_W3110_varios_puntos,yield_growth_rate_varios_puntos,'filled')
% hold off


%% 
load('Result_wp_0_W3110.mat')
load('Result_wp_1565_W3110.mat')
frac=[frac_R_0 frac_R_1565; frac_p_0 frac_p_1565;frac_q_0 frac_q_1565 ;frac_et_0 frac_et_1565;frac_er_0 frac_er_1565;frac_eAc_0 frac_eAc_1565]; 
parnames={'\phi_{R}', '\Phi_P', '\phi_q','\phi_{e_t}','\phi_{e_{res}}','\phi_{e_{fer}}'};

figure(3),clf
bar(frac)
set(gca,'fontsize',30)
set(gca,'xticklabel',parnames,'FontSize',30)


figure(4),clf 
frac_molec=[frac_R_0 frac_R_1565; frac_p_0 frac_p_1565;frac_q_0 frac_q_1565];
frac_enzimes=[frac_er_0 frac_er_1565;frac_et_0 frac_et_1565;frac_eAc_0 frac_eAc_1565];
% subplot(1,2,1)
% parnames_frac={'\phi_{R}','\phi_p','\phi_q'};
% bar(frac_molec)
% set(gca,'xticklabel',parnames_frac,'FontSize',20)
% set(gca,'fontsize',18)
% subplot(1,2,2)
parnames_enz={'\phi_{e_{res}}','\phi_{e_{t}}','\phi_{e_{fer}}'};
bar(frac_enzimes)
set(gca,'xticklabel',parnames_enz,'FontSize',30)
set(gca,'fontsize',30)
%%
load('Result_wp_RBS_W3110')

figure(5),clf 

blue_end = [0, 0.4470, 0.7410];

% Crear un mapa de colores personalizado desde blanco hasta azul claro
nSteps = 100;  % Número de pasos de transición
white_color = [0.95,1,1];
miMapaColores = [linspace(white_color(1), blue_end(1), nSteps)', ...
                 linspace(white_color(2), blue_end(2), nSteps)', ...
                 linspace(white_color(3), blue_end(3), nSteps)'];


subplot(2,2,1)
contourf(A_RBS_Wp_diff_W3110,B_RBS_Wp_diff_W3110,max_lam_RBS_Wp_W3110*60,'LineStyle','none')
hcb=colorbar;
hcb.Label.String='Growth rate [1/h]';
set(gca,'XScale','log','Yscale','log'),xlabel('wp [molecs/cell/min]'),ylabel('RBS')
set(gca,'FontSize',30),ylim([10^0 10^3]);
subplot(2,2,2)
contourf(A_RBS_Wp_diff_W3110,B_RBS_Wp_diff_W3110,productividad_biomasa_RBS_Wp_Diff_W3110,'LineStyle','none')
hcb=colorbar;
hcb.Label.String='Biomass Productivity [g/L/h]';
set(gca,'XScale','log','Yscale','log'),xlabel('wp [molecs/cell/min]'),ylabel('RBS')
set(gca,'FontSize',30),ylim([10^0 10^3])
 colormap (miMapaColores)
subplot(2,2,3)
contourf(A_RBS_Wp_diff_W3110,B_RBS_Wp_diff_W3110,productividad_producto_RBS_Wp_Diff_W3110,'LineStyle','none')
hcb=colorbar;
hcb.Label.String='Protein productivity [g/L/h]';
set(gca,'XScale','log','Yscale','log'),xlabel('wp [molecs/cell/min]'),ylabel('RBS')
set(gca,'FontSize',30),ylim([10^0 10^3])

subplot(2,2,4)
levels = linspace(0, 40, 5);
Time_Wp_RBS_diff_W3110(Time_Wp_RBS_diff_W3110>=40)=40;
Time_Wp_RBS_diff_W3110(Time_Wp_RBS_diff_W3110==0)=40;
contourf(A_RBS_Wp_diff_W3110,B_RBS_Wp_diff_W3110,Time_Wp_RBS_diff_W3110,levels,'LineStyle','none')
hcb=colorbar;ylim([10^0 10^3])
hcb.Label.String='time[h]';
set(gca,'XScale','log','Yscale','log'),xlabel('wp [molecs/cell/min]'),ylabel('RBS')
set(gca,'FontSize',30,'Colormap',flipud(colormap))
%%
clear all 

load('Result_optimization_W3110.mat')
figure(6),clf 
% subplot(1,3,1)
% loglog(pp_optimizacion_W3110(:,1),RBS_W3110,'LineWidth',3),xlabel('W_p'),ylabel('RBS')
% hold on
% Wp_W3110_puntos=[pp_optimizacion_W3110(1,1),pp_optimizacion_W3110(25,1),pp_optimizacion_W3110(50,1)];
% RBS_W3110_puntos=[RBS_W3110(1,1),RBS_W3110(25,1),RBS_W3110(50,1)];
% alfa_puntos=[alfa_1(1,1),alfa_1(1,25),alfa_1(1,50)];
% scatter(Wp_W3110_puntos,RBS_W3110_puntos,70,"filled")
% set(gca,'FontSize',16),ylim([1e-6,1e1])
% hold off

figure(7),clf 
semilogy(pp_optimizacion_W3110(:,1),RBS_W3110,'LineWidth',3),xlabel('W_p'),ylabel('RBS')
hold on
Wp_W3110_puntos=[pp_optimizacion_W3110(1,1),pp_optimizacion_W3110(25,1),pp_optimizacion_W3110(50,1)];
RBS_W3110_puntos=[RBS_W3110(1,1),RBS_W3110(25,1),RBS_W3110(50,1)];
alfa_puntos=[alfa_1(1,1),alfa_1(1,25),alfa_1(1,50)];
scatter(Wp_W3110_puntos,RBS_W3110_puntos,80,"filled")
set(gca,'FontSize',30),xlim([88.2349901092848,1732.38823940775])
hold off

figure(8),clf
biomasa_model_max_RBS_Wp_W3110=round(biomasa_model_max_RBS_Wp_W3110,4);
biomasa_model_max_RBS_Wp_W3110(biomasa_model_max_RBS_Wp_W3110== 2.35230000000000)=2.5;
biomasa_model_max_RBS_Wp_W3110(biomasa_model_max_RBS_Wp_W3110== 2.49200000000000)=2.5;
biomasa_model_max_RBS_Wp_W3110(biomasa_model_max_RBS_Wp_W3110>=2.5)=2.5;
plot(producto_model_max_RBS_Wp_W3110,biomasa_model_max_RBS_Wp_W3110,'LineWidth',3),ylabel('Biomass [g/L]'),xlabel('Protein [g/L]')
hold on
biomasa_W3110_puntos=[biomasa_model_max_RBS_Wp_W3110(1,1),biomasa_model_max_RBS_Wp_W3110(25,1),biomasa_model_max_RBS_Wp_W3110(50,1)];
proteina_W3110_puntos=[producto_model_max_RBS_Wp_W3110(1,1),producto_model_max_RBS_Wp_W3110(25,1),producto_model_max_RBS_Wp_W3110(50,1)];
scatter(proteina_W3110_puntos,biomasa_W3110_puntos,80,'filled')
hold off
set(gca,'FontSize',30)
%%
figure(9),clf 
semilogy(pp_optimizacion_W3110(:,1),RBS_W3110,'LineWidth',3),xlabel('W_p'),ylabel('RBS')
hold on
scatter(pp_optimizacion_W3110(:,1),RBS_W3110,100,alfa_1,'filled'),xlabel('W_p'),ylabel('RBS')

% Wp_W3110_puntos=[pp_optimizacion_W3110(1,1),pp_optimizacion_W3110(25,1),pp_optimizacion_W3110(50,1)];
% RBS_W3110_puntos=[RBS_W3110(1,1),RBS_W3110(25,1),RBS_W3110(50,1)];
% alfa_puntos=[alfa_1(1,1),alfa_1(1,25),alfa_1(1,50)];
% scatter(Wp_W3110_puntos,RBS_W3110_puntos,80,"filled")
set(gca,'FontSize',30),xlim([88.2349901092848,1732.38823940775])
hold off
title('\alpha');

colorbar 
%%
%eSTA ES LA FIGURA USADA EN EL PAPER
load('Result_wp_diff_W3110.mat')
figure(10),clf 

colororder({'k','k'})
yyaxis right
semilogx(wp,productividad_proteina_Wp_diff_W3110,'LineWidth',3,'Color',[0.8500 0.3250 0.0980])
legend('gfp')
hold on 
yyaxis left
semilogx(wp,productividad_biomasa_Wp_diff_W3110,'lineStyle','-','LineWidth',3,'Color',[0 0.4470 0.7410])
semilogx(wp,productividad_acetato_Wp_diff_W3110,'lineStyle','-','LineWidth',3,'Color',[0.9290 0.6940 0.1250])
semilogx(wp,max_lam_Wp_W3110_diff,'lineStyle','-','LineWidth',3,'Color',[0.4660 0.6740 0.1880])
hold off 
legend('Biomass','Acetate','growth rate specific','gfp')
set(gca,'FontSize',30)
ylabel('Productivity [g/L/h] and growth rate specific [1/h]'),xlabel('W_p [molecs/cell/min]')

figure(11),clf 
semilogx(wp,productividad_proteina_Wp_diff_W3110,'LineWidth',3,'Color',[0.8500 0.3250 0.0980])
legend('gfp')
hold on 
semilogx(wp,productividad_biomasa_Wp_diff_W3110,'lineStyle','-','LineWidth',3,'Color',[0 0.4470 0.7410])
semilogx(wp,productividad_acetato_Wp_diff_W3110,'lineStyle','-','LineWidth',3,'Color',[0.9290 0.6940 0.1250])
semilogx(wp,max_lam_Wp_W3110_diff,'lineStyle','-','LineWidth',3,'Color',[0.4660 0.6740 0.1880])
hold off 
legend('gfp','Biomass','Acetate','growth rate specific')
set(gca,'FontSize',30)
ylabel('Productivity [g/L/h] and growth rate specific [1/h]'),xlabel('W_p [molecs/cell/min]')

figure(12),clf 
% scatter(productividad_proteina_Wp_diff_W3110,yield_growth,[100],wp,'filled')
hold on 
plot(productividad_proteina_Wp_diff_W3110,yield_growth,'Color',[0 0.4470 0.7410],'LineWidth',3)
hold off
xlabel('gfp productivity [g/L/h]'),ylabel('Biomass yield')
set(gca,'FontSize',30),ylim([0.47,0.011])
%%
clear all 
load('Result_productivity_wp_kb_ku.mat')
figure(13),clf 
plot(productividad_PR_Wp_W3110,productividad_biomasa_Wp_W3100,'LineWidth',3)
hold on 
scatter(productividad_PR_Wp_W3110,productividad_biomasa_Wp_W3100,[120],alfa_1,'filled')

hold off 
xlabel('gfp productivity'), ylabel('Biomass productivity')
set(gca,'FontSize',30)