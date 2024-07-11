clear all, clc
load('Result_wp_diff_W3110.mat')
figure(1),clf 

x_wp_productividad_proteina=[productividad_proteina_Wp_diff_W3110(1,1),productividad_proteina_Wp_diff_W3110(1,19),productividad_proteina_Wp_diff_W3110(1,30)]
y_wp_productividad_acetato=[productividad_acetato_Wp_diff_W3110(1,1),productividad_acetato_Wp_diff_W3110(1,19),productividad_acetato_Wp_diff_W3110(1,30)]
plot(productividad_proteina_Wp_diff_W3110,productividad_acetato_Wp_diff_W3110,'LineWidth',3,'Color','black')
hold on 
scatter(x_wp_productividad_proteina,y_wp_productividad_acetato,200,"red","filled")
hold off
xlabel('GFP productivity [g/L/h]'),ylabel('Acetate productivity [g/L/h]')
set(gca,'fontsize',30)

figure(2),clf
plot(max_lam_Wp_W3110_diff,productividad_acetato_Wp_diff_W3110,'LineWidth',3)
hold on 
x_wp_lam=[max_lam_Wp_W3110_diff(1,1),max_lam_Wp_W3110_diff(1,19),max_lam_Wp_W3110_diff(1,30)]
y_wp_productividad_acetato=[productividad_acetato_Wp_diff_W3110(1,1),productividad_acetato_Wp_diff_W3110(1,19),productividad_acetato_Wp_diff_W3110(1,30)]
scatter(x_wp_lam,y_wp_productividad_acetato,200,'filled')
hold off
xlabel('Specific growth rate [1/h]'),ylabel('Acetate productivity [g/L/h]')
set(gca,'fontsize',30)

%%
clear all, clc
load('Result_we_wp_0.mat')
max_lam_We_W3110_0=max_lam_We_W3110;
Productividad_acetato_We_W3110_diff_0=Productividad_acetato_We_W3110_diff;
v_produccion_Ac_We_W3110_0=v_produccion_Ac_We_W3110;
load('Result_we_wp_500.mat')
max_lam_We_W3110_500=max_lam_We_W3110;
Productividad_acetato_We_W3110_diff_500=Productividad_acetato_We_W3110_diff;
v_produccion_Ac_We_W3110_500=v_produccion_Ac_We_W3110;
load('Result_we_wp_1000.mat')
max_lam_We_W3110_1000=max_lam_We_W3110;
Productividad_acetato_We_W3110_diff_1000=Productividad_acetato_We_W3110_diff;
v_produccion_Ac_We_W3110_1000=v_produccion_Ac_We_W3110;
load('Result_we_wp_1500.mat')
max_lam_We_W3110_1500=max_lam_We_W3110;
Productividad_acetato_We_W3110_diff_1500=Productividad_acetato_We_W3110_diff;
v_produccion_Ac_We_W3110_1500=v_produccion_Ac_We_W3110;
load('Result_we_wp_2000.mat')
max_lam_We_W3110_2000=max_lam_We_W3110;
Productividad_acetato_We_W3110_diff_2000=Productividad_acetato_We_W3110_diff;
v_produccion_Ac_We_W3110_2000=v_produccion_Ac_We_W3110;
load('Result_we_wp_3000.mat')
max_lam_We_W3110_3000=max_lam_We_W3110;
Productividad_acetato_We_W3110_diff_3000=Productividad_acetato_We_W3110_diff;
v_produccion_Ac_We_W3110_3000=v_produccion_Ac_We_W3110;

figure(3),clf 
plot(max_lam_We_W3110_0(1:22),Productividad_acetato_We_W3110_diff_0(1:22),'LineWidth',4,'Color',[0.4660 0.6740 0.1880])
hold on 
%plot(max_lam_We_W3110_500(1:23),Productividad_acetato_We_W3110_diff_500(1:23),'LineWidth',2)

plot(max_lam_We_W3110_1000(1:22),Productividad_acetato_We_W3110_diff_1000(1:22),'LineWidth',4,'Color',[0.8500 0.3250 0.0980])
%plot(max_lam_We_W3110_1500(1:23),Productividad_acetato_We_W3110_diff_1500(1:23),'LineWidth',2)
plot(max_lam_We_W3110_2000(1:22),Productividad_acetato_We_W3110_diff_2000(1:22),'LineWidth',4,'Color',[0.9290 0.6940 0.1250])
%plot(max_lam_We_W3110_3000(1:22),Productividad_acetato_We_W3110_diff_3000(1:22),'LineWidth',4,'Color',[0.9290 0.6940 0.1250])

hold off
xlabel('Specific growth rate [1/h]'),ylabel('Acetate productivity')
set(gca,'fontsize',30)
xlim([0.13,0.5])
figure(4),clf 
plot(max_lam_We_W3110_0(1:24),v_produccion_Ac_We_W3110_0(1:24)*1/(mpp*Na*V)*60,'LineWidth',3)
hold on 
plot(max_lam_We_W3110_500(1:23),v_produccion_Ac_We_W3110_500(1:23)*1/(mpp*Na*V)*60,'LineWidth',2)

plot(max_lam_We_W3110_1000(1:23),v_produccion_Ac_We_W3110_1000(1:23)*1/(mpp*Na*V)*60,'LineWidth',2)
plot(max_lam_We_W3110_1500(1:23),v_produccion_Ac_We_W3110_1500(1:23)*1/(mpp*Na*V)*60,'LineWidth',2)
plot(max_lam_We_W3110_2000(1:23),v_produccion_Ac_We_W3110_2000(1:23)*1/(mpp*Na*V)*60,'LineWidth',2)
hold off
xlabel('Specific growth rate [1/h]'),ylabel('V_{prod_{Ac}}')
set(gca,'fontsize',30)

figure(5),clf 
semilogy(we,v_importacion_S_We_cell_W3110,'LineWidth',5,'Color','black')
set(gca,'fontsize',30)
xlabel('w_t'), ylabel('Substrate import rate')
%%
load('Result_wet_wp.mat')

figure(5),clf
blue_end = [0, 0.4470, 0.7410];

% Crear un mapa de colores personalizado desde blanco hasta azul claro
nSteps = 100;  % Número de pasos de transición
white_color = [0.95,1,1];
miMapaColores = [linspace(white_color(1), blue_end(1), nSteps)', ...
                 linspace(white_color(2), blue_end(2), nSteps)', ...
                 linspace(white_color(3), blue_end(3), nSteps)'];

subplot(1,4,1)
contourf(A_We_Wp_W3110,B_We_Wp_W3110,max_lam_We_Wp_W3110,'LineStyle','none')
hcb=colorbar;
hcb.Label.String='Specific growth rate [1/h]';
set(gca,'XScale','log'),xlabel('wp [molecs/cell/min]'),ylabel('We [molecs/cell/min]')
 set(gca,'FontSize',20)

subplot(1,4,2)
contourf(A_We_Wp_W3110,B_We_Wp_W3110,productividad_biomasa_We_Wp_Diff_W3110,'LineStyle','none')
hcb=colorbar;
hcb.Label.String='Biomass productivity [g/L/h]';
set(gca,'XScale','log'),xlabel('wp [molecs/cell/min]'),ylabel('we [molecs/cell/min]')
 set(gca,'FontSize',20)

subplot(1,4,3)
contourf(A_We_Wp_W3110,B_We_Wp_W3110,productividad_producto_We_Wp_Diff_W3110,'LineStyle','none')
xlabel('wp'),ylabel('we')
hcb=colorbar;
hcb.Label.String='gfp propductivity [g/L/h]';
set(gca,'XScale','log')
 set(gca,'FontSize',20)

subplot(1,4,4)
contourf(A_We_Wp_W3110,B_We_Wp_W3110,productividad_acetato_We_Wp_Diff_W3110,'LineStyle','none')
xlabel('wp [molecs/cell/min]'),ylabel('we[molecs/cell/min]')
hcb=colorbar;
hcb.Label.String='Acetate productivity [g/L/h]';
set(gca,'XScale','log')
set(gca,'FontSize',20)
colormap (miMapaColores)


