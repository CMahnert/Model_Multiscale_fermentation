function ss=Adjust_batch_VAL23(pp,data)


data.ydata=[
  %time     Biomass [g/L]     Substrate [g/L]     Product [g/L] Ace [g/L]
  0         0.157             4.8                0.01           0
  1*60      0.157             4.44               0.014          0.009
  2*60      0.188             4.26               0.024          0.009
  3*60      0.25              4.16               0.027          0.01
  4*60      0.34              3.96               0.03           0.115
  5*60      0.44              3.85               0.035          0.125
  6*60      0.518             3.36               0.042          0.178
  7*60      0.7               3.09               0.051          0.196
  8*60      0.88              2.35               0.06           0.257
  9*60     1.16               1.69               0.083          0.48
  10*60     1.5               0.68               0.1            0.55
  11*60     1.8               0.12               0.12           0.66
  12*60     2.13              0                  0.153          0.566
  13*60     2.13              0                  0.153          0.38
 ];


global J_suma_substrato_W3110 J_suma_biomass_W3110

time=data.ydata(:,1);
biomass_obs=data.ydata(:,2);
substrate_obs=data.ydata(:,3);
product_obs=data.ydata(:,4);
acetate_obs=data.ydata(:,5);

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


%Condition initial
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
    Ac_int_0=0;
    mq_0= xx(end,15);
    mp_0= 0;%xx(end,16);
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

[t,x]= ode15s(@(t,x)model_multiscale_batch(t,x,rates,parametro,pp), time, inic,options);

gamma4=gmax*x(1,21)/(pp(1)+x(1,21));
lamda3=(x(1,1)+x(1,3)+x(1,4)+x(1,5)+x(1,6)+x(1,8)).*gamma4/M;
lamda=lamda3*60;

biomass_model=x(:,23); %g/L
substrate_model=x(:,24); %g/L
product_model=x(:,25).*27000; %g/L
acetate_model=x(:,26).*60/(Na*V);

J_biomass_W3110=(biomass_model-biomass_obs).^2;
J_substrate_W3110=(substrate_model-substrate_obs).^2;
J_suma_biomass_W3110=sum(J_biomass_W3110);
J_suma_substrato_W3110=sum(J_substrate_W3110);

A=1*(biomass_model-biomass_obs).^2+1*(substrate_model-substrate_obs).^2+100*(product_model-product_obs).^2+100*(acetate_model-acetate_obs).^2; 
ss=sum(A);
end