clear all, clc


wp=logspace(0,5,30);
we=linspace(0.001,10,30);
for h=1:length(we)

for z=1:length(wp)
%%

b= 0;
dm= 0.1;
ku=1;
kb=1;
f= 0;
ds= 0;
dn=0; %muerte celular
rates= [b dm kb ku f ds dn];

V=0.4;
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

parametro= [Na PMgluc mpp thetar s0 gmax thetax Kt M Km vm nx Kq vt wr wq nq nr V yy np ns n_xAc theta_xAc Kcat_Ac Km_Ac Kcat_Ac_in Km_Ac_in Kgamma wp(z) wAc wer nf we(h)] ;

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
options = odeset('RelTol',1e-4,'AbsTol',1e-4);%,'NonNegative',[1:20],'Stats','off'); %RelTol=error relativo--- AbsTol=error absoluto

[tt,xx]= ode15s(@(t,x)Model_multiscale_initial_condition_analysis(t,x,rates,parametro), timee, inic,options);
rmr_0=xx(end,1);   
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
   N_0= 1.8667e+11; 
biomasa_0=N_0*mpp/V;
s_0=4.8;  
PR_0=0.006/27000;
Ac_ext_0=0;
 inic= [rmr_0 em_0 rmp_0 rmAc_0 rmq_0 rmt_0 et_0 rmm_0 mt_0 mm_0 q_0 p_0 Ac_0 si_0 Ac_int_0 mq_0 mp_0 mAc_0 mr_0 r_0 a_0 N_0 biomasa_0 s_0 PR_0 Ac_ext_0];

times=linspace(0,1e4);
[t,x]= ode15s(@(t,x)Model_multiscale_analysis(t,x,rates,parametro), times, inic,options);

[dxdt,lam]=cellfun(@(t,x) Model_multiscale_analysis(t,x,rates,parametro),num2cell(t),num2cell(x,2),'uni',0);%función para poder ver variables que no forman parte del ode y permite que se guarden en el ODE 
lam=cell2mat(lam);
max_lam_We_Wp_W3110(h,z)=max(lam).*60;
biomasa_model_We_Wp_W3110(z)=x(end,23); %g/L
sustrate_model_We_Wp_W3110(z)=x(end,24); %g/L
producto_model_We_Wp_W3110(z)=max(x(:,25)).*27000; %g/L
acetato_model_We_Wp_W3110=x(:,26).*60/(Na*V);
acetato_model_max_We_Wp_W3110(z)=max(x(:,26).*60/(Na*V));

concentracion_biomasa_We_Wp_W3110(h,z)=biomasa_model_We_Wp_W3110(z);
concentracion_proteina_We_Wp_W3110(h,z)=producto_model_We_Wp_W3110(z);
concentracion_acetato_We_Wp_W3110(h,z)=(acetato_model_max_We_Wp_W3110(z));

N_Wp_We_diff_W3110{h,z}=diff(x(:,23));
 N_Wp_We_diff_W3110_round{h,z}=round(N_Wp_We_diff_W3110{h,z},1);
 time_ast_Wp_We_diff_W3110{h,z}=t(1:length(N_Wp_We_diff_W3110{h,z}));%.'/60; %Function para invertir columnas por filas
 col{h,z}=find(N_Wp_We_diff_W3110_round{h,z}==0);

 idx = ~cellfun('isempty',col);
 out = zeros(size(col));
 out(idx) = cellfun(@(v)v(1),col(idx)); 
 col_time=[out].';
 col_time(col_time==0)=1;
Time_Wp_We_diff_W3110=time_ast_Wp_We_diff_W3110{h,z}(col_time)/60;
Time_Wp_We_diff_W3110(Time_Wp_We_diff_W3110==0)=100000000000000/60;
 

end 
end 
%%
for j=1:length(we)

    for d=1:length(wp)

productividad_biomasa_We_Wp_Diff_W3110(d,j)=concentracion_biomasa_We_Wp_W3110(d,j)/Time_Wp_We_diff_W3110(j,d);
productividad_biomasa_We_Wp_Diff_W3110(productividad_biomasa_We_Wp_Diff_W3110==inf)=0;
productividad_producto_We_Wp_Diff_W3110(d,j)=concentracion_proteina_We_Wp_W3110(d,j)/Time_Wp_We_diff_W3110(j,d);
productividad_producto_We_Wp_Diff_W3110(productividad_producto_We_Wp_Diff_W3110==inf)=0;
productividad_acetato_We_Wp_Diff_W3110(d,j)=concentracion_acetato_We_Wp_W3110(d,j)/Time_Wp_We_diff_W3110(j,d);
productividad_acetato_We_Wp_Diff_W3110(productividad_acetato_We_Wp_Diff_W3110==inf)=0;

    end 
end



%%
[A_We_Wp_W3110,B_We_Wp_W3110]=meshgrid(wp,we);


figure(1),clf
subplot(1,4,1)
contourf(A_We_Wp_W3110,B_We_Wp_W3110,max_lam_We_Wp_W3110,'LineStyle','none')
hcb=colorbar;
hcb.Label.String='Specific growth rate [1/h]';
set(gca,'XScale','log'),xlabel('wp [molecs/cell/min]'),ylabel('We')
set(gca,'FontSize',20)

subplot(1,4,2)
contourf(A_We_Wp_W3110,B_We_Wp_W3110,concentracion_biomasa_We_Wp_W3110,'LineStyle','none')
xlabel('wp'),ylabel('we')
hcb=colorbar;
hcb.Label.String='Biomasa [g/L]';
set(gca,'XScale','log')
set(gca,'FontSize',20)

subplot(1,4,3)
contourf(A_We_Wp_W3110,B_We_Wp_W3110,concentracion_proteina_We_Wp_W3110,'LineStyle','none')
xlabel('wp'),ylabel('we')
hcb=colorbar;
hcb.Label.String='gfp [g/L]';
set(gca,'XScale','log')
set(gca,'FontSize',20)

subplot(1,4,4)
contourf(A_We_Wp_W3110,B_We_Wp_W3110,concentracion_acetato_We_Wp_W3110,'LineStyle','none')
xlabel('wp'),ylabel('we')
hcb=colorbar;
hcb.Label.String='Acetate [g/L]';
set(gca,'XScale','log')
set(gca,'FontSize',20)

figure(2),clf 
subplot(1,3,1)
contourf(A_We_Wp_W3110,B_We_Wp_W3110,productividad_biomasa_We_Wp_Diff_W3110,'LineStyle','none')
hcb=colorbar;
hcb.Label.String='Biomass productivity [g/L/h]';
set(gca,'XScale','log'),xlabel('wp [molecs/cell/min]'),ylabel('we [molecs/cell/min]')
set(gca,'FontSize',20)

subplot(1,3,2)
contourf(A_We_Wp_W3110,B_We_Wp_W3110,productividad_producto_We_Wp_Diff_W3110,'LineStyle','none')
xlabel('wp'),ylabel('we')
hcb=colorbar;
hcb.Label.String='gfp propductivity [g/L/h]';
set(gca,'XScale','log')
set(gca,'FontSize',20)

subplot(1,3,3)
contourf(A_We_Wp_W3110,B_We_Wp_W3110,productividad_acetato_We_Wp_Diff_W3110,'LineStyle','none')
xlabel('wp'),ylabel('we')
hcb=colorbar;
hcb.Label.String='Acetate productivity [g/L/h]';
set(gca,'XScale','log')
set(gca,'FontSize',20)

