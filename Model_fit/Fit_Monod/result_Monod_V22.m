clear all, clc

data.ydata=[
  %time     Biomass [g/L]     Suvstrate [g/L]     Product [g/L] Ace [g/L]
 0          0.14              4.8                0.004          0
 1*60       0.2               4.2                0.004          0.14
 2*60       0.25              3.9                0.0071         0.21
 3*60       0.35              3.57               0.0085         0.23
 4*60       0.53              3.22               0.0114         0.3
 5*60       1.11              2.36               0.0214         0.35
 6*60       1.71              1.15               0.0328         0.4
 7*60       2                 0.034              0.0399         0.39
 8*60       1.98              0.0347             0.0514         0.254
 9*60       1.82              0                  0.052          0.096
 10*60      1.73              0                  0.053          0
 11*60      1.64              0                  0.053          0
 12*60      1.66              0                  0.053          0
 ];

time_Monod_V22=data.ydata(:,1);
Biomass_obs_Monod_V22=data.ydata(:,2);
Substrate_obs_Monod_V22=data.ydata(:,3);
Product_obs_Monod_V22=data.ydata(:,4);
Acetate_obs_Monod_V22=data.ydata(:,5);
%%
global J_suma_sustrato_Monod_V22 J_suma_biomasa_Monod_V22

yxs=0.4;
umax=0.4/60;
Km_S=7e-5;
k00=[yxs umax Km_S];

lb_0=[0.2 0.4/60 7e-5];
ub_0=[0.5 0.6/60 0.002 ];
options = optimoptions(@fmincon,'Algorithm','interior-point');
[pp,ss0_Monod_V22]=fmincon(@Adjust_Monod_V22,k00,[],[],[],[],lb_0,ub_0,[],options,data);
mse=ss0_Monod_V22/(length(data.ydata)-3); %el valor de 3 corresponde a los valores estimados

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
P_0=0.004;
Ac_0=0;
inic=[S_0 N_0 P_0 Ac_0];

options = odeset('RelTol',1e-4,'AbsTol',1e-4);%,'NonNegative',[1:20],'Stats','off'); %RelTol=error relativo--- AbsTol=error absoluto

[t_Monod_V22,x_Monod_V22]= ode15s(@(t,x)model_Monod(t,x,parametro,pp), time_Monod_V22, inic,options);

Biomass_model_Monod_V22=x_Monod_V22(:,2);
Substrate_model_Monod_V22=x_Monod_V22(:,1);
Product_model_Monod_V22=x_Monod_V22(:,3);
Acetate_model_Monod_V22=x_Monod_V22(:,4);

%% figure 
figure(1),clf 
subplot(1,2,1)
plot(time_Monod_V22/60,Biomass_obs_Monod_V22,'o',time_Monod_V22/60,Biomass_model_Monod_V22,'r','LineWidth',2)
subplot(1,2,2)
plot(time_Monod_V22/60,Substrate_obs_Monod_V22,'o',time_Monod_V22/60,Substrate_model_Monod_V22,'r','LineWidth',2)




