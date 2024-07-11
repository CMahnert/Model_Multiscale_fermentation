function ss=Adjust_optimization_W3110(pp,data)

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

global alfa biomasa_model_max producto_model_max acetato_model_max
time=data.ydata(:,1);
Biomasa_obs=data.ydata(:,2);
Sustrate_obs=data.ydata(:,3);
Producto_obs=data.ydata(:,4);
acetato_obs=data.ydata(:,5);

b= 0;
dm= 0.1;
%kb= 1;%1;
f= 0;
ds= 0;
dn=0; %muerte celular

 rates= [b dm f ds dn];                                   

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
Km= 3200; %Brenda, 0.008[mM] con ECMDB de 0.4 molecs/cell y sustrato de 3-phospho-D-glycerate
vm= 4300; %1/min, se considero la enzima más lenta que en este caso corresponde a PGK o 3-PGK
nx= 300;
Kq= 152219;
vt= 10800;%1/min..obtenido de: A steady-state of microbial acclimation to substrate limitation
wr= 930;
wq= 949;
nq= 4;%coeficiente de Hill 
nr= 7459;
yy=0.45;
np=238; %GFP
ns=0.2;%0.3;
n_xAc=nx;
theta_xAc=thetax;
Kcat_Ac=93000;%53400;
Km_Ac=140000;
Kcat_Ac_in=85200;%69900;%85200;
Km_Ac_in=14408000;%11720000;%2800000;
Kgamma=399999.01;
 wAc=0.54856;
 wer=4.5895;
 nf=0.14696;
 we=4.311;

parametro= [Na PMgluc mpp thetar s0 gmax thetax Kt M Km vm nx Kq vt wr wq nq nr V yy np ns n_xAc theta_xAc Kcat_Ac Km_Ac Kcat_Ac_in Km_Ac_in Kgamma wAc wer nf we alfa] ;

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

timee=linspace(0,1e6);
options = odeset('RelTol',1e-6,'AbsTol',1e-9);%,'NonNegative',[1:20],'Stats','off'); %RelTol=error relativo--- AbsTol=error absoluto

%Condicion con strain wit R-prot
[tt,xx]= ode15s(@(t,x)Model_multiscale_initial_condition_optimization(t,x,rates,parametro,pp), timee, inic,options);

 	rmr_0= xx(end,1);
	em_0= xx(end,2);
    rmp_0= 0;%xx(end,3);
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
    mq_0= xx(end,15);
    mp_0= 0;%xx(end,16);
    mAc_0=xx(end,17);
    mr_0= xx(end,18);
    r_0= xx(end,19);
    a_0= xx(end,20);
   N_0= 1.8667e+11; %[cell]
biomasa_0=N_0*mpp/V;
s_0=4.8; %[molec] 
PR_0=0.006/27000;%5.9259e-06;%0.16/27000;%5.9259e-06;
Ac_ext_0=0;
PR_con_0=0.006;
Ac_ext_con_0=0;
inic= [rmr_0 em_0 rmp_0 rmAc_0 rmq_0 rmt_0 et_0 rmm_0 mt_0 mm_0 q_0 p_0 Ac_0 si_0 mq_0 mp_0 mAc_0 mr_0 r_0 a_0 N_0 biomasa_0 s_0 PR_0 Ac_ext_0 PR_con_0 Ac_ext_con_0];

times=linspace(0,1e4);
[t,x]= ode15s(@(t,x)Model_multiscale_optimization(t,x,rates,parametro,pp), times, inic,options);

Biomasa_max=2.5;
Proteina_max=0.1;
biomasa_model_max=max(x(:,22));
producto_model_max=max(x(:,24)).*27000;
acetato_model_max=max(x(:,25).*60/(Na*V));
A=alfa.*((biomasa_model_max-Biomasa_max)).^2+10*(1-alfa).*(producto_model_max-Proteina_max).^2;

%A=alfa.*((biomasa_model_max-Biomasa_max)).^2+10*(1-alfa).*((producto_model_max-Proteina_max)).^2;
ss=sum(A);
end