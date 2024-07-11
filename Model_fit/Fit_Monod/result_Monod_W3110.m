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

time_Monod_W3110=data.ydata(:,1);
Biomass_obs_Monod_W3110=data.ydata(:,2);
Substrate_obs_Monod_W3110=data.ydata(:,3);
Product_obs_Monod_W3110=data.ydata(:,4);
Acetate_obs_Monod_W3110=data.ydata(:,5);
%%
global J_suma_sustrato_Monod_W3110 J_suma_biomasa_Monod_W3110

yxs=0.4;
umax=0.2/60;
Km_S=0.0002;
k00=[yxs umax Km_S];
lb_0=[0.3 0.2/60 7e-5];
ub_0=[0.5 0.3/60 0.002 ];
options = optimoptions(@fmincon,'Algorithm','interior-point');
[pp,ss0_Monod_W3110]=fmincon(@Adjust_Monod_W3110,k00,[],[],[],[],lb_0,ub_0,[],options,data);
mse=ss0_Monod_W3110/(length(data.ydata)-2); %el valor de 3 corresponde a los valores estimados
yxs=pp(1);
umax=pp(2);
Km_S=pp(3);

disp(['yxs:' num2str(yxs)])
disp(['umax:' num2str(umax)])
disp(['Km_S:' num2str(Km_S)])


pp=[yxs umax Km_S];

dn=0.01;
S0=4.8;
yps=0.5;
yas=0.5;
parametro=[dn S0 yps yas];

S_0=4.8;
N_0=0.14;
P_0=0.006;
Ac_0=0;
inic=[S_0 N_0 P_0 Ac_0];

options = odeset('RelTol',1e-4,'AbsTol',1e-4);%,'NonNegative',[1:20],'Stats','off'); %RelTol=error relativo--- AbsTol=error absoluto

[t_Monod_W3110,x_Monod_W3110]= ode15s(@(t,x)model_Monod(t,x,parametro,pp), time_Monod_W3110, inic,options);

Biomass_model_Monod_W3110=x_Monod_W3110(:,2);
Substrate_model_Monod_W3110=x_Monod_W3110(:,1);
Product_model_Monod_W3110=x_Monod_W3110(:,3);
Acetate_model_Monod_V22=x_Monod_W3110(:,4);

%% figure 
figure(1),clf 
subplot(1,2,1)
plot(time_Monod_W3110/60,Biomass_obs_Monod_W3110,'o',time_Monod_W3110/60,Biomass_model_Monod_W3110,'r','LineWidth',2)
subplot(1,2,2)
plot(time_Monod_W3110/60,Substrate_obs_Monod_W3110,'o',time_Monod_W3110/60,Substrate_model_Monod_W3110,'r','LineWidth',2)
