clear all, clc
%kb=1;
global columna alfa productividad_biomasa col Time_diff productividad_PR N_T_diff time_ast productividad_Ac concentracion_proteina concentracion_biomasa
alfa_1=linspace(0,1,30);%[0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85  0.9 0.95 1];
%alfa_1=0.2;%[0 0.2 0.6 0.8 1];
 %    wp   ku   kb
lb_0=[1e0 1e-6 0.8];
ub_0=[1e5 1 1];
wp=1e2;
ku=1;
kb=1;
k00=[wp ku kb];
for i=1:length(alfa_1)
    
    alfa=alfa_1(i);

parametro_estimar=3; %Corresponde a los tres parámetros a estimar para optimizar (Wp Ku Kb)
myfun=@Adjust_optimization_W3110;
%options = gaoptimset('PopulationSize', 50, 'Generations', 100);
%[pp_optimizacion_W3110(i,:),ss0(i,:)]=ga(myfun,parametro_estimar,[],[],[],[],lb_0,ub_0,[],[]); %Entrega muy buenos resultados. Mejor que el resto debido a los pesos que considera
[pp_optimizacion_W3110(i,:),ss0(i,:)]=gamultiobj(myfun,parametro_estimar,[],[],[],[],lb_0,ub_0,[],[]); %Entrega muy buenos resultados. Mejor que el resto debido a los pesos que considera
%[pp_optimizacion_W3110(i,:),ss0(i,:)]=fmincon(myfun,k00,[],[],[],[],lb_0,ub_0,[]);

lb_0=[1e0 1e-6 1e-5]; 
%lb_0=[1e0 pp_optimizacion_W3110(i,2) 1e-6]; 

%ub_0=[pp_optimizacion_W3110(i,1) 1 pp_optimizacion_W3110(i,3)];
ub_0=[pp_optimizacion_W3110(i,1) 1 pp_optimizacion_W3110(i,3)];

%ub_0=[pp_optimizacion_W3110(i,1) 1 1];

    productividad_biomasa_Wp_W3100(i,:)=productividad_biomasa;
    productividad_PR_Wp_W3110(i,:)=productividad_PR;
    concentracion_biomasa_W3110(i,:)=concentracion_biomasa;
    concentracion_proteina_W3110(i,:)=concentracion_proteina;
    Time_diff_W3110(i,:)=Time_diff;
end

%% Figure 
%Alfa=0
% wp_W3110_alfa_0=pp_optimizacion_W3110(1,1);
% ku_W3110_alfa_0=pp_optimizacion_W3110(1,2);
% kb_W3110_alfa_0=pp_optimizacion_W3110(1,3);
% 
% %Alfa=0.5
% wp_W3110_alfa_05=pp_optimizacion_W3110(6,1);
% ku_W3110_alfa_05=pp_optimizacion_W3110(6,2);
% kb_W3110_alfa_05=pp_optimizacion_W3110(6,3);
% 
% %Alfa=1
% wp_W3110_alfa_1=pp_optimizacion_W3110(11,1);
% ku_W3110_alfa_1=pp_optimizacion_W3110(11,2);
% kb_W3110_alfa_1=pp_optimizacion_W3110(11,3);

RBS_W3110=pp_optimizacion_W3110(:,3)./pp_optimizacion_W3110(:,2);
figure(1),clf 
semilogx(pp_optimizacion_W3110(:,1),RBS_W3110,'LineWidth',2),xlabel('W_p'),ylabel('RBS')
%%
figure(2),clf 
plot3(pp_optimizacion_W3110(:,1),pp_optimizacion_W3110(:,2),pp_optimizacion_W3110(:,3),'LineWidth',2,'LineStyle','-')
% hold on 
% scatter3(wp_W3110_alfa_0,ku_W3110_alfa_0,kb_W3110_alfa_0,'filled','MarkerFaceColor',[1 0 0])
% scatter3(wp_W3110_alfa_05,ku_W3110_alfa_05,kb_W3110_alfa_05,'filled','LineWidth',4,'MarkerFaceColor',[1 0 0])
% scatter3(wp_W3110_alfa_1,ku_W3110_alfa_1,kb_W3110_alfa_1,'filled','LineWidth',4,'MarkerFaceColor',[1 0 0])

xlabel('Wp'),ylabel('Ku'),zlabel('Kb'),grid on
hold off
%%
figure(3),clf
scatter3(pp_optimizacion_W3110(:,1),pp_optimizacion_W3110(:,2),pp_optimizacion_W3110(:,3),[],alfa_1,'filled')
s.SizeData = 11;
colorbar

figure(4),clf 
subplot(1,2,1)
plot(productividad_PR_Wp_W3110,productividad_biomasa_Wp_W3100),xlabel('gfp productivity'),ylabel('biomass productivity')
subplot(1,2,2)
plot(concentracion_proteina_W3110,concentracion_biomasa_W3110)