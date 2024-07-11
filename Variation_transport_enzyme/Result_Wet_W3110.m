clear all, clc

data.ydata=[
  %time     Biomasa [g/L]     Sustrate [g/L]     Producto [g/L] Ace [g/L]
  
 0*60          0.14              4.8                 0.006      0
 1*60          0.19              4.43                0.007      0.05
 2*60          0.22              4.36                0.01       0.13
 3*60          0.27              4.255               0.011      0.15
 4*60          0.33              3.9                 0.016      0.244
 5*60          0.39              3.47                0.02       0.28
 6*60          0.461             3.19                0.022      0.376
 7*60          0.557             3.04                0.03       0.51
 8*60          0.71              2.51                0.038      0.53
 9*60          0.89              2.09                0.048      0.6
 10*60         1                 1.45                0.06       0.639
 11*60         1.24              0.46                0.069      0.724
  12*60         1.35              0.106               0.079      0.68
   13*60         1.544             0                   0.095     0.536
  14*60         1.75              0                   0.11      0.45
 15*60        1.86              0                   0.12       0.29
 16*60         1.91              0                   0.123     0.27
 ];

time=data.ydata(:,1);
Biomasa_obs=data.ydata(:,2);
Sustrate_obs=data.ydata(:,3);
Producto_obs=data.ydata(:,4);
acetato_obs=data.ydata(:,5);
we=linspace(0.1,10);
for z=1:length(we)
%%
b= 0;
dm= 0.1;
kb=1;
ku=1;
f= 0;
ds= 0;
dn=0; %muerte celular
rates= [b dm kb ku f ds dn];

V=0.4;
Na=6.022e+23; %molec/mol
PMgluc=180; %g/mol
mpp=3e-13;%gdw/cell promedio 
thetar=427;
s0=4.8;
gmax= 1260;
thetax= 4.38;
Kt=0.003 ;%el valor esta en uM%6640*(5e+10);
M= 1.0e8;
Km= 96000; %Brenda, 0.24[mM] con ECMDB de 0.4 molecs/cell y sustrato de 3-phospho-D-glycerate
vm= 4300; %1/min, se considero la enzima más lenta que en este caso corresponde a PGK o 3-PGK
nx= 300;
Kq= 152219;
vt= 10800;%1/min..obtenido de: A steady-state of microbial acclimation to substrate limitation
wr= 930;
wq= 949;
nq= 2;%coeficiente de Hill 
nr= 7459;
yy=0.45;
np=238; %GFP
ns=0.2;
n_xAc=nx;
theta_xAc=thetax;
Kcat_Ac=93000;%53400;
Km_Ac=140000;%64000;%140000;
Kcat_Ac_in=85200;%33600;%69900;%85200;
Km_Ac_in=14408000;%11720000;%2800000;
 Kgamma=172968.9731;
 wp=3000;%1565.5665;
 wAc=0.53652;
 wer=4.9671;
 nf=0.092207;
 %we=4.3247;

parametro= [ Na PMgluc mpp thetar s0 gmax thetax Kt M Km vm nx Kq vt wr wq nq nr V yy np ns n_xAc theta_xAc Kcat_Ac Km_Ac Kcat_Ac_in Km_Ac_in Kgamma wp wAc wer nf we(z)] ;

%Condición inicial with R-prot
   rmr_0= 0;
    em_0= 0;
    rmp_0= 0;
    rmAc_0=0;
    rmq_0= 0;
    rmt_0= 0;
    et_0= 1;
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

timee=linspace(0,1e6);
options = odeset('RelTol',1e-4,'AbsTol',1e-4);%,'NonNegative',[1:20],'Stats','off'); %RelTol=error relativo--- AbsTol=error absoluto

%Condicion con strain wit R-prot
[tt,xx]= ode15s(@(t,x)Model_multiscale_initial_condition_analysis(t,x,rates,parametro), timee, inic,options);
rnr_0=xx(end,1);   
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
s_0=s0;
PR_0=0.006/27000;%5.9259e-06;%0.16/27000;%5.9259e-06;
Ac_ext_0=0;

inic= [rmr_0 em_0 rmp_0 rmAc_0 rmq_0 rmt_0 et_0 rmm_0 mt_0 mm_0 q_0 p_0 Ac_0 si_0 Ac_int_0 mq_0 mp_0 mAc_0 mr_0 r_0 a_0 N_0 biomasa_0 s_0 PR_0 Ac_ext_0];

times=linspace(0,1e4,1e3);
[t,x_we]= ode15s(@(t,x)Model_multiscale_analysis(t,x,rates,parametro), times, inic,options);

[dxdt,lam]=cellfun(@(t,x_we) Model_multiscale_analysis(t,x_we,rates,parametro),num2cell(t),num2cell(x_we,2),'uni',0);%función para poder ver variables que no forman parte del ode y permite que se guarden en el ODE 
lam=cell2mat(lam);
max_lam_We_W3110(z)=max(lam).*60;
lam_We_W3110(z)=lam(1,1).*60;
%max_lam_We_W3110_3000(z)=max(lam).*60;

concentracion_biomasa_We_W3110(z)=max(x_we(:,23));
%concentracion_biomasa_We_W3110_3000(z)=max(x_we(:,22));
concentracion_proteina_We_W3110(z)=max(x_we(:,25)).*27000;
%concentracion_proteina_We_W3110_3000(z)=max(x_we(:,24).*27000);
concentracion_acetato_We_W3110(z)=max(x_we(:,26)).*60/(Na*V);
%concentracion_acetato_We_W3110_3000(z)=max(x_we(:,25).*60/(Na*V));

molec_e_Ac_We_W3110(z)=x_we(end,13);
molec_e_Ac_We_W3110_inic(z)=x_we(1,13);
Substrate_internal_We_W3110(z)=max(x_we(:,14));
molec_Ac_We_S_W3110(z)=max(x_we(:,24));
molec_et_Ac_We_W3110(z)=max(x_we(:,7));
molec_et_final_We_W3110(z)=x_we(1,7);
Substrate_internal_We_W3110_inic(z)=x_we(1,14);
Ac_ext_We_W3110_max(z)=max(x_we(end,26));
N_cell_We_W3110(z)=max(x_we(:,22));

v_produccion_Ac_We_cell_W3110(z)=max(x_we(:,13).*Kcat_Ac.*x_we(:,14)./(Km_Ac+x_we(:,14)));
%v_produccion_Ac_We_cell_W3110_3000(z)=max(x_we(:,13).*Kcat_Ac.*x_we(:,14)./(Km_Ac+x_we(:,14)));
v_produccion_Ac_We_W3110(z)=max(x_we(:,13).*Kcat_Ac.*x_we(:,14)./(Km_Ac+x_we(:,14))).*N_cell_We_W3110(z);
%v_produccion_Ac_We_W3110_3000(z)=max(x_we(:,13).*Kcat_Ac.*x_we(:,14)./(Km_Ac+x_we(:,14))).*N_cell_We_W3110(z);

v_importacion_S_We_cell_W3110(z)=max(x_we(:,7).*vt.*x_we(:,24)./(Kt+x_we(:,24)));
%v_importacion_S_We_cell_W3110_3000(z)=max(x_we(:,7).*vt.*x_we(:,23)./(Kt+x_we(:,23)));

v_importacion_S_We_W3110(z)=max(x_we(:,7).*vt.*x_we(:,24)./(Kt+x_we(:,24))).*N_cell_We_W3110(z);
%v_importacion_S_We_W3110_3000(z)=max(x_we(:,7).*vt.*x_we(:,23)./(Kt+x_we(:,23))).*N_cell_We_W3110(z);

v_importacion_S_We_cell_W3110_final(z)=(x_we(50,7).*vt.*x_we(50,24)./(Kt+x_we(50,24)));
%v_importacion_S_We_cell_W3110_final_3000(z)=(x_we(50,7).*vt.*x_we(50,23)./(Kt+x_we(50,23)));

v_importacion_S_We_cell_W3110_inicial(z)=(x_we(1,7).*vt.*x_we(1,24)./(Kt+x_we(1,24)));


v_consumo_Ac_We_cell_W3110(z)=max(x_we(:,13).*Kcat_Ac_in.*(x_we(:,26)./x_we(:,22))./(Km_Ac_in+x_we(:,26)./x_we(:,22)));%velocidad
%v_consumo_Ac_We_cell_W3110_3000(z)=max(x_we(:,13).*Kcat_Ac_in.*(x_we(:,25)./x_we(:,21))./(Km_Ac_in+x_we(:,25)./x_we(:,21)));%velocidad

v_consumo_Ac_We_W3110(z)=max(x_we(:,13).*Kcat_Ac_in.*(x_we(:,26)./x_we(:,22))./(Km_Ac_in+x_we(:,26)./x_we(:,22))).*N_cell_We_W3110(z);%velocidad
%v_consumo_Ac_We_W3110_3000(z)=max(x_we(:,13).*Kcat_Ac_in.*(x_we(:,25)./x_we(:,21))./(Km_Ac_in+x_we(:,25)./x_we(:,21))).*N_cell_We_W3110(z);%velocidad

%Diff sustrato 
 N_We_diff_W3110(:,z)=diff(x_we(:,23));
 N_We_diff_W3110_round(:,z)=round(N_We_diff_W3110(:,z),12);
time_ast_We_diff_W3110(:,z)=t(1:length(N_We_diff_W3110_round(:,z))); %Function para invertir columnas por filas

col{z}=find(N_We_diff_W3110_round(:,z)<=0.001);
idx = ~cellfun('isempty',col);
out = zeros(size(col));
out(idx) = cellfun(@(v)v(1),col(idx)); %Se guardan los valores para cada iteración
col_time=[out].';
col_time(col_time==0)=1;
Time_We_diff_W3110=time_ast_We_diff_W3110(col_time,z)./60;
% Time_We_diff_W3110(1)=t(999);
%  Time_We_diff_W3110(2)=t(930);
% Time_We_diff_W3110(3)=t(605);
% Time_We_diff_W3110(4)=t(458);
% Time_We_diff_W3110(5)=t(370);
% Time_We_diff_W3110(6)=t(309);
%  Time_We_diff_W3110(7)=t(267);
%   Time_We_diff_W3110(8)=t(234);
%   Time_We_diff_W3110(9)=t(210);
%   Time_We_diff_W3110(10)=t(190);
%   Time_We_diff_W3110(11)=t(173);
%   Time_We_diff_W3110(12)=t(159);

 % Rendimiento
 Rendimiento_biomasa_We_W3110(z)=abs(x_we(1,23)-x_we(end,23))/(abs(x_we(1,24)-x_we(end,24)));
 Rendimiento_proteina_We_W3110(z)=abs(x_we(1,25).*27000-x_we(end,25).*27000)/(abs(x_we(1,24)-x_we(end,24)));
 Rendimiento_acetato_We_W3110(z)=abs(x_we(1,26)-max(x_we(:,26))).*60/(Na*V)/(abs(x_we(1,24)-x_we(end,24)));

% Productividad
Productividad_biomasa_We_W3110_diff(z)=concentracion_biomasa_We_W3110(z)/Time_We_diff_W3110(z);
Productividad_proteina_We_W3110_diff(z)=concentracion_proteina_We_W3110(z)/Time_We_diff_W3110(z);
Productividad_acetato_We_W3110_diff(z)=concentracion_acetato_We_W3110(z)/Time_We_diff_W3110(z);  

end 
%%
figure(1),clf 
subplot(1,2,1)
plot(Rendimiento_biomasa_We_W3110,Productividad_acetato_We_W3110_diff); 
subplot(1,2,2)
plot(Rendimiento_biomasa_We_W3110,Productividad_proteina_We_W3110_diff);

figure(2),clf 
plot(max_lam_We_W3110,Productividad_acetato_We_W3110_diff)

figure(3),clf
plot(max_lam_We_W3110,v_produccion_Ac_We_W3110)

figure(4),clf 
subplot(3,1,1)
plot(we,v_importacion_S_We_cell_W3110),xlabel('we'),ylabel('v_{imp}')
subplot(3,1,2)
plot(v_importacion_S_We_cell_W3110,Productividad_acetato_We_W3110_diff),xlabel('v_{imp}'),ylabel('Acetato productivity')
subplot(3,1,3)
plot(we,Productividad_biomasa_We_W3110_diff),xlabel('we'),ylabel('acetate productivity')

figure(4),clf
subplot(1,2,1)
semilogy(we,v_importacion_S_We_cell_W3110_inicial)
subplot(1,2,2)
plot(we,max_lam_We_W3110)

figure(5),clf
plot(max_lam_We_W3110,Productividad_acetato_We_W3110_diff)