clear all, clc

data.ydata=[
  %time     Biomass [g/L]     Substrate [g/L]     Product [g/L] Ace [g/L]
  
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

time_fmincon_W3110=data.ydata(:,1);
Biomasa_obs_fmincon_W3110=data.ydata(:,2);
Sustrate_obs_fmincon_W3110=data.ydata(:,3);
Producto_obs_fmincon_W3110=data.ydata(:,4);
acetato_obs_fmincon_W3110=data.ydata(:,5);

%%
global J_suma_substrato_W3110 J_suma_biomass_W3110

 Kgamma=1e3;%1e3;%172968.9731;
 wp=1565.5665;
 wAc=0.53652;
 wer=4.9671;
 nf=0.092207;
 we=4.3247;
k00=[Kgamma wp wAc wer nf we];
     %Ky    wp    wf  wr  nf     we 
lb_0=[1e2    1000     0.1  1 0.08    1];
ub_0=[4e5   4000     1   6  0.1    5];
parametro_estimar=6;
options = optimoptions(@fmincon,'Algorithm','interior-point');
%[pp,ss0]=fminsearch(@Adjust_batch_W3110,k00,[],data);
%[pp,ss0]=fmincon(@Adjust_batch_W3110,k00,[],[],[],[],lb_0,ub_0,data);
[pp,ss0_fmincon_W3110]=fmincon(@Adjust_batch_W3110,k00,[],[],[],[],lb_0,ub_0,[],options,data);
%[pp,ss0_fmincon_W3110]=gamultiobj(@Adjust_batch_W3110,parametro_estimar,[],[],[],[],lb_0,ub_0,[],[]);

mse=ss0_fmincon_W3110/(length(data.ydata)-6); %el valor de 3 corresponde a los valores estimados

Kgamma=pp(1);
wp=pp(2);
wAc=pp(3);
wer=pp(4);
nf=pp(5);
we=pp(6);

disp(['Kgamma:' num2str(Kgamma)])
disp(['wp:' num2str(wp)])
disp(['wAc:' num2str(wAc)])
disp(['wer:' num2str(wer)])
disp(['nf:' num2str(nf)])
disp(['we:' num2str(we)])

pp=[Kgamma wp wAc wer nf we];

%%
b= 0;
dm= 0.1;
kb= 1;%1;
ku= 1;
f= 0;
dn=0; %muerte celular
rates= [b dm kb ku f dn];

V=0.4;
Na=6.022e+23; %molec/mol
PMgluc=180; %g/mol
mpp=3e-13;%gdw/cell promedio 
thetar=427;
s0= 4.8;%1.2044e+22;
gmax= 1260;
thetax= 4.38;
Kt=0.003 ;%el valor esta en uM%6640*(5e+10);
M= 1.0e8;
Km= 96000; %Brenda, 0.008[mM] con ECMDB de 0.4 molecs/cell y sustrato de 3-phospho-D-glycerate
vm= 4300; %1/min, se considero la enzima más lenta que en este caso corresponde a PGK o 3-PGK
nx= 300;
Kq= 152219;
vt= 10800;%8000;%10800;%1/min..obtenido de: A steady-state of microbial acclimation to substrate limitation
wr= 930;
wq= 949;
nq= 4;%coeficiente de Hill 
nr= 7459;
yy=0.45;
np=238; %GFP
ns=0.15;%0.15;%0.2;%0.3;
n_xAc=nx;
theta_xAc=thetax;
Kcat_Ac=93000;%53400;
Km_Ac=140000;
Kcat_Ac_in=85200;%69900;%85200;
Km_Ac_in=14408000;%11720000;%2800000;
parametro= [Na PMgluc mpp thetar s0 gmax thetax Kt M Km vm nx Kq vt wr wq nq nr V yy np ns n_xAc theta_xAc Kcat_Ac Km_Ac Kcat_Ac_in Km_Ac_in] ;


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

timee=linspace(0,1000000);
options = odeset('RelTol',1e-4,'AbsTol',1e-4);%,'NonNegative',[1:20],'Stats','off'); %RelTol=error relativo--- AbsTol=error absoluto


%Condicion con strain wit R-prot
[tt,xx]= ode15s(@(t,x)model_multiscale_condition_initial(t,x,rates,parametro,pp), timee, inic,options);
    rmr_0= xx(end,1);
	em_0= xx(end,2);
    rmp_0= xx(end,3);
    rmAc_0=xx(end,4);
    rmq_0= xx(end,5);
    rmt_0= xx(end,6);
    et_0= xx(end,7);
    rmm_0= xx(end,8);
    mt_0= xx(end,9);
    mm_0= xx(end,10);
    q_0= xx(end,11);
    p_0= xx(end,12);
    Ac_0=xx(end,13);
    si_0= xx(end,14);
    Ac_int_0=0;
    mq_0= xx(end,15);
    mp_0= xx(end,16);
    mAc_0=xx(end,17);
    mr_0= xx(end,18);
    r_0= xx(end,19);
    a_0= xx(end,20);
   N_0= 1.8667e+11; %[cell]
biomass_0=N_0*mpp/V;
s_0=4.8; %[molec] 
PR_0=0.006/27000;%5.9259e-06;%0.16/27000;%5.9259e-06;
Ac_ext_0=0;
PR_con_0=0.006;
Ac_ext_con_0=0;
inic= [rmr_0 em_0 rmp_0 rmAc_0 rmq_0 rmt_0 et_0 rmm_0 mt_0 mm_0 q_0 p_0 Ac_0 si_0 Ac_int_0 mq_0 mp_0 mAc_0 mr_0 r_0 a_0 N_0 biomass_0 s_0 PR_0 Ac_ext_0 PR_con_0 Ac_ext_con_0];

times=linspace(0,10000);
[t_fmincon_W3110,x_fmincon_W3110]= ode15s(@(t,x)model_multiscale_batch(t,x,rates,parametro,pp), [0 100*60], inic,options);

gamma=gmax*x_fmincon_W3110(:,21)/(pp(1)+x_fmincon_W3110(:,21));
lamda3=(x_fmincon_W3110(:,1)+x_fmincon_W3110(:,3)+x_fmincon_W3110(:,4)+x_fmincon_W3110(:,5)+x_fmincon_W3110(:,6)+x_fmincon_W3110(:,8)).*gamma/M;
lamda=lamda3.*60;

biomass_model_fmincon_W3110=x_fmincon_W3110(:,23); %g/L
substrate_model_fmincon_W3110=x_fmincon_W3110(:,24); %g/L
product_model_fmincon_W3110=x_fmincon_W3110(:,25).*27000; %g/L
acetate_model_fmincon_W3110=x_fmincon_W3110(:,26).*60/(Na*V);

%%
figure(1),clf
subplot(3,2,1)
plot(time_fmincon_W3110/60,Biomasa_obs_fmincon_W3110,'o',t_fmincon_W3110/60,biomass_model_fmincon_W3110,'r','LineWidth',2);xlabel('time [h]'),ylabel('Biomasa [g/L]');legend('Lara et al 2006','Simulation');grid on
xlim([0,16])
subplot(3,2,2)
plot(time_fmincon_W3110/60,Sustrate_obs_fmincon_W3110,'o',t_fmincon_W3110/60,substrate_model_fmincon_W3110,'r','LineWidth',2);xlabel('time [h]'),ylabel('Sustrate [g/L]');legend('Lara et al 2006','Simulation');  grid on
xlim([0,16])
subplot(3,2,3)
plot(time_fmincon_W3110/60,Producto_obs_fmincon_W3110,'o',t_fmincon_W3110/60,product_model_fmincon_W3110,'r','LineWidth',2);xlabel('time [h]'),ylabel('Producte [g/L]');legend('Lara et al 2006','Simulation');  grid on
xlim([0,16])
subplot(3,2,4)
plot(time_fmincon_W3110/60,acetato_obs_fmincon_W3110,'o',t_fmincon_W3110/60,acetate_model_fmincon_W3110,'r','LineWidth',2);xlabel('time [h]'),ylabel('Acetate [g/L]');legend('Lara et al 2006','Simulation');  grid on
xlim([0,16])

figure(2),clf 
plot(t_fmincon_W3110/60,lamda)
