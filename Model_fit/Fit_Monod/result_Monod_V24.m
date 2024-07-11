clear all, clc

data.ydata=[
  %time     Biomasa [g/L]     Sustrate [g/L]     Producto [g/L] Ace [g/L]
  0        0.121             5                   0.016           0
 1*60      0.228             4.89                0.021           0.068
 2*60      0.304             4.68                0.03            0.14
 3*60      0.42              4.479               0.04            0.25
 4*60      0.55              3.99                0.044           0.32
 5*60      0.94              3.19                0.064           0.41
 6*60      1.53              1.8                 0.11            0.44
 7*60      2.2               0.104               0.15            0.51
 8*60      2.34              0.104               0.19            0.47
 9*60      2.41              0                   0.23            0.32
 10*60     2.44              0                   0.23            0.27
 ];

time_Monod_V22=data.ydata(:,1);
Biomass_obs_Monod_V22=data.ydata(:,2);
Substrate_obs_Monod_V22=data.ydata(:,3);
Product_obs_Monod_V22=data.ydata(:,4);
Acetate_obs_Monod_V22=data.ydata(:,5);
%%
global J_suma_sustrato_Monod_V24 J_suma_biomasa_Monod_V24

yxs=0.2;
umax=0.4/60;
Km_S=7e-5;
k00=[yxs umax Km_S];

lb_0=[0.2 0.4/60 7e-5];
ub_0=[0.5 0.5/60 0.002 ];
options = optimoptions(@fmincon,'Algorithm','interior-point');
[pp,ss0_Monod_V24]=fmincon(@Adjust_Monod_V24,k00,[],[],[],[],lb_0,ub_0,[],options,data);
mse=ss0_Monod_V24/(length(data.ydata)-3); %el valor de 3 corresponde a los valores estimados

yxs=pp(1);
umax=pp(2);
Km_S=pp(3);
disp(['yxs:' num2str(yxs)])
disp(['umax:' num2str(umax)])
disp(['Km_S:' num2str(Km_S)])

pp=[yxs umax Km_S];

dn=0.01;
S0=5;
yps=0.5;
yas=0.5;
parametro=[dn S0 yps yas];

S_0=5;
N_0=0.121;
P_0=0.006;
Ac_0=0;
inic=[S_0 N_0 P_0 Ac_0];

options = odeset('RelTol',1e-4,'AbsTol',1e-4);%,'NonNegative',[1:20],'Stats','off'); %RelTol=error relativo--- AbsTol=error absoluto

[t_Monod_V24,x_Monod_V24]= ode15s(@(t,x)model_Monod(t,x,parametro,pp), time_Monod_V24, inic,options);

Biomass_model_Monod_V24=x_Monod_V24(:,2);
Substrate_model_Monod_V24=x_Monod_V24(:,1);
Product_model_Monod_V24=x_Monod_V24(:,3);
Acetate_model_Monod_V24=x_Monod_V24(:,4);

%% figure 
figure(1),clf 
subplot(1,2,1)
plot(time_Monod_V24/60,Biomass_obs_Monod_V24,'o',time_Monod_V24/60,Biomass_model_Monod_V24,'r','LineWidth',2)
subplot(1,2,2)
plot(time_Monod_V24/60,Substrate_obs_Monod_V24,'o',time_Monod_V24/60,Substrate_model_Monod_V24,'r','LineWidth',2)
