clear all, clc

wp=logspace(0,5,30);

for z=1:length(wp)
%%
b= 0;
dm= 0.1;
kb=1;
ku=1;
f= 0;
ds= 0;
dn=0; 
rates= [b dm kb ku f ds dn];

V=0.4;
vol_cell=10^-15;
Na=6.022e+23; 
PMgluc=180; 
mpp=3e-13;
thetar=427;
s0= 4.8;
gmax= 1260;
thetax= 4.38;
Kt=0.003 ;
M= 1.0e8;
Km= 96000; 
vm= 4300; 
nx= 300;
Kq= 152219;
vt= 10800;
wr= 930;
wq= 949;
nq= 4; 
nr= 7459;
yy=0.45;
np=238; 
ns=0.2;
n_xAc=nx;
theta_xAc=thetax;
Kcat_Ac=93000;
Km_Ac=140000;
Kcat_Ac_in=85200;
Km_Ac_in=14408000;
 Kgamma=172968.9731;
 wAc=0.53652;
 wer=4.9671;
 nf=0.092207;
 we=4.3247;


parametro= [ Na PMgluc mpp thetar s0 gmax thetax Kt M Km vm nx Kq vt wr wq nq nr V yy np ns n_xAc theta_xAc Kcat_Ac Km_Ac Kcat_Ac_in Km_Ac_in Kgamma wp(z) wAc wer nf we] ;

%Condición inicial with R-prot
   rmr_0= 0;
    em_0= 0;
    rmp_0= 0;
    rmAc_0=0;
    rmq_0= 0;
    rmt_0= 0;
    et_0= 0;
    rmm_0= 0;
    mt_0= 0;
    mm_0= 0;
    q_0= 0;
    p_0= 0;
    Ac_0=0;
    si_0= 0;
    mq_0= 0;
    mp_0= 0;
    mAc_0=0;
    mr_0= 0;
    r_0= 10;
    a_0= 1000;
  
inic= [rmr_0 em_0 rmp_0 rmAc_0 rmq_0 rmt_0 et_0 rmm_0 mt_0 mm_0 q_0 p_0 Ac_0 si_0 mq_0 mp_0 mAc_0 mr_0 r_0 a_0];

timee=linspace(0,1e6);%(0,1e5,1e4);
options = odeset('RelTol',1e-4,'AbsTol',1e-4);%,'NonNegative',[1:20],'Stats','off'); %RelTol=error relativo--- AbsTol=error absoluto

%Condicion con strain wit R-prot
[tt,xx]= ode15s(@(t,x)Model_multiscale_initial_condition_analysis(t,x,rates,parametro), timee, inic,options);
rmr_0=xx(end,1);   
em_0= xx(end,2);
    rmp_0=0;% xx(end,3);
    rmAc_0=xx(end,4);
    rmq_0= xx(end,5);
    rmt_0= xx(end,6);
    et_0= xx(end,7);
    rmm_0= xx(end,8);
    mt_0= xx(end,9);
    mm_0= xx(end,10);
    q_0= xx(end,11);
    p_0= 0;%xx(end,12);
    Ac_0=xx(end,13);
    si_0= xx(end,14);
    Ac_int_0=0;
    mq_0= xx(end,15);
    mp_0= 0;%xx(end,16);
    mAc_0=xx(end,17);
    mr_0= xx(end,18);
    r_0= xx(end,19);
    a_0= xx(end,20);
   N_0= 1.8667e+11; %[cell]
biomasa_0=N_0*mpp/V;
s_0=4.8;
PR_0=0.006/27000;%5.9259e-06;%0.16/27000;%5.9259e-06;
Ac_ext_0=0;


inic= [rmr_0 em_0 rmp_0 rmAc_0 rmq_0 rmt_0 et_0 rmm_0 mt_0 mm_0 q_0 p_0 Ac_0 si_0 Ac_int_0 mq_0 mp_0 mAc_0 mr_0 r_0 a_0 N_0 biomasa_0 s_0 PR_0 Ac_ext_0];


times=linspace(0,5e4,5e3);%(0,1e5,1e4);

[t,x]= ode15s(@(t,x)Model_multiscale_analysis(t,x,rates,parametro), times, inic,options);

[dxdt,lam]=cellfun(@(t,x) Model_multiscale_analysis(t,x,rates,parametro),num2cell(t),num2cell(x,2),'uni',0);%función para poder ver variables que no forman parte del ode y permite que se guarden en el ODE 
lam=cell2mat(lam);
max_lam_Wp_W3110_diff(z)=max(lam).*60;
 N_Wp_diff_W3110(:,z)=diff(x(:,23));
 N_Wp_diff_W3110_round(:,z)=round(N_Wp_diff_W3110(:,z),12);%rnd 3
time_ast_Wp_diff_W3110(:,z)=t(1:length(N_Wp_diff_W3110_round(:,z))); %Function para invertir columnas por filas

col{z}=find(N_Wp_diff_W3110_round(:,z)<=0.001);
idx = ~cellfun('isempty',col);
out = zeros(size(col));
out(idx) = cellfun(@(v)v(1),col(idx)); 
col_time=[out].';
%col_time(col_time==0)=1;
Time_Wp_diff_W3110=time_ast_Wp_diff_W3110(col_time,z)./60;
Time_Wp_diff_W3110(Time_Wp_diff_W3110==0)=23;

Time_Wp_diff_W3110(23)=32;
Time_Wp_diff_W3110(24)=42;
Time_Wp_diff_W3110(25)=59;
Time_Wp_diff_W3110(26)=82;
Time_Wp_diff_W3110(27)=117;
Time_Wp_diff_W3110(28)=176;
Time_Wp_diff_W3110(29)=252;
Time_Wp_diff_W3110(30)=360;
Time_Wp_diff_W3110=round(Time_Wp_diff_W3110,2);

biomasa_model_Wp_diff_W3110(z)=x(end,23); %g/L
sustrate_model_Wp_diff_W3110(z)=x(end,24); %g/L
producto_model_Wp_diff_W3110(z)=max(x(:,25)).*27000; %g/L
acetato_model_Wp_diff_W3110=x(:,26).*60/(Na*V);
acetato_model_max_Wp_diff_W3110(z)=max(acetato_model_Wp_diff_W3110);

concentracion_biomasa_Wp_W3110_diff(z)=x(end,23);
concentracion_proteina_Wp_W3110_diff(z)=max(x(:,25)*27000);
concentracion_acetato_Wp_W3110_diff(z)=max(x(:,26).*60/(Na*V));

productividad_proteina_Wp_diff_W3110(z)=producto_model_Wp_diff_W3110(z)./(Time_Wp_diff_W3110(z));
productividad_biomasa_Wp_diff_W3110(z)=biomasa_model_Wp_diff_W3110(z)./(Time_Wp_diff_W3110(z));
%productividad_biomasa_Wp_diff_W3110(z)=round(productividad_biomasa_Wp_diff_W3110(z),3);
productividad_acetato_Wp_diff_W3110(z)=acetato_model_max_Wp_diff_W3110(z)./(Time_Wp_diff_W3110(z));
v_produccion_Ac_Wp_cell_W3110(z)=max(x(:,13).*Kcat_Ac.*x(:,14)./(Km_Ac+x(:,14)));
v_produccion_Ac_Wp_W3110(z)=max(x(:,13).*Kcat_Ac.*x(:,14)./(Km_Ac+x(:,14))).*x(end,21);

v_importacion_S_Wp_cell_W3110(z)=max(x(:,7).*vt.*x(:,24)./(Kt+x(:,24)));
v_importacion_S_Wp_W3110(z)=max(x(:,7).*vt.*x(:,24)./(Kt+x(:,24))).*x(end,21);

v_consumo_Ac_Wp_cell_W3110(z)=max(x(:,13).*Kcat_Ac_in.*(x(:,26)./x(:,22))./(Km_Ac_in+x(:,26)./x(:,22)));%velocidad
v_consumo_Ac_Wp_W3110(z)=max(x(:,13).*Kcat_Ac_in.*(x(:,26)./x(:,22))./(Km_Ac_in+x(:,26)./x(:,22))).*x(end,22);%velocidad

yield_growth(z)=abs(x(1,23)-x(end,23))./abs(x(1,24)-x(end,24)); 
yield_gfp(z)=abs(x(1,25)-x(end,25))/abs(x(1,24)-x(end,24));
yield_acetate(z)=abs(x(1,26)-max(x(:,26)))/abs(x(1,24)-x(end,24));
productividad_lamda_biomasa(z)=concentracion_biomasa_Wp_W3110_diff(z)./max_lam_Wp_W3110_diff(z);
productividad_lamda_acetato(z)=concentracion_acetato_Wp_W3110_diff(z)./max_lam_Wp_W3110_diff(z);

cell(z)=x(end,21);
energia(z)=x(end,20);
end 

%% FIGURES
figure(1),clf 
subplot(1,3,1)
semilogx(wp,productividad_proteina_Wp_diff_W3110,'LineWidth',3),xlabel('wp'),ylabel('gfp productivity')
subplot(1,3,2)
semilogx(wp,productividad_biomasa_Wp_diff_W3110,'LineWidth',3),xlabel('wp'),ylabel('gfp productivity')
subplot(1,3,3)
semilogx(wp,productividad_acetato_Wp_diff_W3110,'LineWidth',3),xlabel('wp'),ylabel('gfp productivity')

figure(2),clf 
plot(productividad_proteina_Wp_diff_W3110,productividad_acetato_Wp_diff_W3110),xlabel('gfp productivity [g/L/h]'),ylabel('Acetate productivity [g/L/h]')

figure(3),clf 
semilogx(wp,Time_Wp_diff_W3110),xlabel('Wp'),ylabel('Time [h]')

figure(4),clf 
plot(productividad_proteina_Wp_diff_W3110,yield_growth,'LineWidth',7),xlabel('gfp productivity [g/L/h]'),ylabel('Growth yield'),set(gca, 'FontSize', 30)

figure(5),clf 
plot(max_lam_Wp_W3110_diff,acetato_model_max_Wp_diff_W3110)

figure(6),clf 
subplot(1,2,1)
plot(productividad_proteina_Wp_diff_W3110,productividad_acetato_Wp_diff_W3110)
subplot(1,2,2)
plot(concentracion_proteina_Wp_W3110_diff,concentracion_acetato_Wp_W3110_diff)
