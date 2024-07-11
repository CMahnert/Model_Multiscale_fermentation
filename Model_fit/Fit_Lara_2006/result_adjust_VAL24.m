clear all, clc

data.ydata=[
  %time     Biomass [g/L]     Substrate [g/L]     Product [g/L] Ace [g/L]
 0        0.121             5                   0.016           0
 1*60      0.228             4.89                0.018           0.068
 2*60      0.304             4.5                0.027            0.08
 3*60      0.42              4.2               0.038            0.25
 4*60      0.55              3.7                0.042           0.32
 5*60      0.94              3.19                0.061           0.4
 6*60      1.53              1.1                 0.09            0.42
 7*60      2.2               0.104               0.13            0.47
 8*60      2.4              0.104               0.16            0.45
 9*60      2.52              0                   0.18            0.30
 10*60     2.52              0                   0.18            0.27
 ];

time_fmincon_W3110=data.ydata(:,1);
Biomasa_obs_fmincon_W3110=data.ydata(:,2);
Sustrate_obs_fmincon_W3110=data.ydata(:,3);
Producto_obs_fmincon_W3110=data.ydata(:,4);
acetato_obs_fmincon_W3110=data.ydata(:,5);

%%
global J_suma_substrato_W3110 J_suma_biomass_W3110

Kgamma=1000;
wp=1000.0048;%1267.7;
wAc=0.250781;%0.25001;%34773;%0.3593;
wer=3.7143 ;
nf=0.080969;%0.09;
we=3.0208;
options = optimoptions(@fmincon,'Algorithm','interior-point','Display','off');
k00=[Kgamma wp wAc wer nf we];

parametro_estimar=6;

lb_0=[10    100  0.2   1   0.08  1];
ub_0=[1000  10000 1    5   0.2   4];
%[pp,ss0]=fminsearch(@Adjust_batch_VAL24,k00,[],data);
%[pp,ss0]=fmincon(@Adjust_batch_VAL24,k00,[],[],[],[],lb_0,ub_0,data);
%[pp,ss0_fmincon_W3110]=fmincon(@Adjust_batch_VAL24,k00,[],[],[],[],lb_0,ub_0,[],options,data);
[pp,ss0_fmincon_W3110]=gamultiobj(@Adjust_batch_VAL24,parametro_estimar,[],[],[],[],lb_0,ub_0,[],[]);

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
ns=0.2;%0.15;%0.2;%0.3;
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
PR_0=0.004/27000;%5.9259e-06;%0.16/27000;%5.9259e-06;
Ac_ext_0=0;
PR_con_0=0.004;
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
