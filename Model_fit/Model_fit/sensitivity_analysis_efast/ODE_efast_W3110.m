%% This ODE represents the HIV model in Section 4.2
function dydt=ODEmodel(t,y,rates,parametro,X,run_num)

%% PARAMETERS %%
%parameter_setting_w3110;

Kgamma=X(run_num,1);
wp=X(run_num,2);
wAc=X(run_num,3);
wer=X(run_num,4);
nf=X(run_num,5);
we=X(run_num,6);

b= 0;
dm= 0.1;
kb= 1;%1;
ku= 1;
f= 0;
ds= 0;
dn=0; %muerte celular
rates= [b dm kb ku f ds dn];

V=0.4;
Na=6.022e+23; %molec/mol
PMgluc=180; %g/mol
mpp=3e-13;%gdw/cell promedio 
thetar=427;
s0= 4.8;%1.2044e+22;
gmax= 1260;
cl= 0;
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
PMgly=92.09382;
PMgal=180.156;
PMlac=342.3;
PMsucc=118.1;
ns=0.2;%0.3;
n_xAc=nx;
theta_xAc=thetax;
Kcat_Ac=93000;%53400;
Km_Ac=140000;
Kcat_Ac_in=85200;%69900;%85200;
Km_Ac_in=14408000;%11720000;%2800000;
parametro= [Na PMgluc mpp thetar s0 gmax thetax Kt M Km vm nx Kq vt wr wq nq nr V yy np ns n_xAc theta_xAc Kcat_Ac Km_Ac Kcat_Ac_in Km_Ac_in] ;
options = odeset('RelTol',1e-3,'AbsTol',1e-3);%,'NonNegative',[1:20],'Stats','off'); %RelTol=error relativo--- AbsTol=error absoluto

%% SUSTRATO VARIABLE

 rmr_0= 1.7e3;
	em_0= 464.55;
    rmp_0= 1.8e3;
    rmAc_0=0.73;
    rmq_0= 1202.53;
    rmt_0= 6.10488269087667;
    et_0= 431.750960544174;
    rmm_0= 6.56865908664122;
    mt_0= 41.8461443297373;
    mm_0= 45.0251167976109;
    q_0= 85046.0100063734;
    p_0= 164050.223565854;
    Ac_0=52.0579540854592;
    si_0= 569574196.140578;
    Ac_int_0=0;
    mq_0= 8242.82511238665;
    mp_0= 13375.2810368323;
    mAc_0=5.04555834091491;
    mr_0= 8833.07248158489;
    r_0= 0.190687439260528;
    a_0= 222688.784529954;
   N_0= 1.8667e+11; %[cell]
biomass_0=N_0*mpp/V;
s_0=4.8; %[molec] 
PR_0=0.006/27000;%5.9259e-06;%0.16/27000;%5.9259e-06;
Ac_ext_0=0;
inic= [rmr_0 em_0 rmp_0 rmAc_0 rmq_0 rmt_0 et_0 rmm_0 mt_0 mm_0 q_0 p_0 Ac_0 si_0 Ac_int_0 mq_0 mp_0 mAc_0 mr_0 r_0 a_0 N_0 biomass_0 s_0 PR_0 Ac_ext_0];

 %Variables del modelo
 	rmr= y(1);
	em= y(2);
	rmp= y(3);
    rmAc=y(4); 	
    rmq= y(5);
	rmt= y(6);
	et= y(7);
	rmm= y(8);
	mt= y(9);
	mm= y(10);
	q= y(11);
	p= y(12);
    Ac=y(13); 
    si= y(14);
    Ac_int=y(15);
	mq= y(16);
	mp= y(17);
    mAc=y(18); 
    mr= y(19);
	r= y(20);
	a= y(21);
    N=y(22);
    biomass=y(23);
    s=y(24);
    PR=y(25);
    Ac_ext=y(26);
   

	gamma= gmax*a/(Kgamma + a);
	ttrate= (rmq + rmr + rmt + rmm)*gamma;
	lam= (ttrate+rmp*gamma+rmAc*gamma)/M;

    nucat= em*vm*si/(Km + si);
 nucat_x=p*vm*si/(Km+si);%*a/(Km+a);
    nuimp=et*vt*s/(Kt+s); 
    nucat_Ac=Ac*vm*si/(Km+si);%*a/(Km*a);
    v_prod_Ac=Ac*Kcat_Ac*si/(Km_Ac+si);
   v_in_Ac=0;%Ac*Kcat_Ac_in*(Ac_ext/N)/(Km_Ac_in+Ac_ext/N);%velocidad de entrada del acetato por la ruta AckA
    nucat_in_Ac=em*121*Ac_int/(220000+Ac_int); 

mc=1e-13*exp(36.904*lam);
   feed=2/100*s0;
   if s>feed
  v_in_Ac=0;
   else 
    v_in_Ac=Ac*Kcat_Ac_in*(Ac_ext/N)/(Km_Ac_in+Ac_ext/N);
   end
   

dydt(size(y,1),1)= 0;

	dydt(1)= +kb*r*mr-ku*rmr-gamma/nr*rmr-f*rmr-lam*rmr;
	dydt(2)= +gamma/nx*rmm-lam*em;
	dydt(3)= +kb*r*mp-ku*rmp-gamma/np*rmp-f*rmp-lam*rmp;
	
    dydt(4)=+kb*r*mAc-ku*rmAc-gamma/n_xAc*rmAc-f*rmAc-lam*rmAc;
    
    dydt(5)= +kb*r*mq-ku*rmq-gamma/nx*rmq-f*rmq-lam*rmq;
	dydt(6)= +kb*r*mt-ku*rmt-gamma/nx*rmt-f*rmt-lam*rmt;
	dydt(7)= +gamma/nx*rmt-lam*et;
	dydt(8)= +kb*r*mm-ku*rmm-gamma/nx*rmm-f*rmm-lam*rmm;
	dydt(9)= +(we*a/(thetax + a))+ku*rmt+gamma/nx*rmt-kb*r*mt-dm*mt-lam*mt;
	dydt(10)= +(wer*a/(thetax + a))+ku*rmm+gamma/nx*rmm-kb*r*mm-dm*mm-lam*mm;
	dydt(11)= +gamma/nx*rmq-lam*q;
	dydt(12)= +gamma/np*rmp-lam*p;

    dydt(13)=+gamma/n_xAc*rmAc-lam*Ac;

	dydt(14)= +nuimp-nucat-nucat_Ac-lam*si;
	dydt(15)=+v_in_Ac-nucat_in_Ac-lam*Ac_int;
    dydt(16)= +(wq*a/(thetax + a)/(1 + (q/Kq)^nq))+ku*rmq+gamma/nx*rmq-kb*r*mq-dm*mq-lam*mq;
	dydt(17)= +(wp*a/(thetax + a))+ku*rmp+gamma/np*rmp-kb*r*mp-dm*mp-lam*mp;
	
    dydt(18)=+(wAc*a/(theta_xAc + a))+ku*rmAc+gamma/n_xAc*rmAc-kb*r*mAc-dm*mAc-lam*mAc;
    
    dydt(19)= +(wr*a/(thetar + a))+ku*rmr+gamma/nr*rmr-kb*r*mr-dm*mr-lam*mr;
	dydt(20)= +ku*rmr+ku*rmt+ku*rmm+ku*rmp+ku*rmAc+ku*rmq+gamma/nr*rmr+gamma/nr*rmr+gamma/nx*rmt+gamma/nx*rmm+gamma/np*rmp+gamma/n_xAc*rmAc+gamma/nx*rmq-kb*r*mr-kb*r*mt-kb*r*mm-kb*r*mp-kb*r*mAc-kb*r*mq-lam*r; %ribosoma

    dydt(21)= +(ns*nucat+nf*(nucat_Ac)+ns*0.26*nucat_in_Ac-ttrate-np*rmp/np*gamma-n_xAc*rmAc/n_xAc*gamma-lam*a);

    %Bioreactor
    dydt(22)=+(lam*N-dn*N);
    dydt(23)=+(dydt(22)*mpp/V);
    dydt(24)=-nuimp*N/yy/Na*PMgluc/V;
    dydt(25)=+(lam*p*N-dn*p*N)/(Na*V);
  
   if s>feed
   dydt(26)=+(v_prod_Ac*N );
   elseif s<feed
           v_in_Ac=Ac*Kcat_Ac_in*(Ac_ext/N)/(Km_Ac_in+Ac_ext/N);
       dydt(26)=(v_prod_Ac*N - v_in_Ac*N);
   end 

    

