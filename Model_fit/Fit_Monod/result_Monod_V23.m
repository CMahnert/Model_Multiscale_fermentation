clear all, clc

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

time_Monod_V23=data.ydata(:,1);
Biomass_obs_Monod_V23=data.ydata(:,2);
Substrate_obs_Monod_V23=data.ydata(:,3);
Product_obs_Monod_V23=data.ydata(:,4);
acetate_obs_Monod_V23=data.ydata(:,5);
%%
global J_suma_sustrato_Monod_V23 J_suma_biomasa_Monod_V23


yxs=0.3;
umax=0.25/60;
Km_S=0.001;
k00=[yxs umax Km_S];

lb_0=[0.2 0.2/60 7e-5];
ub_0=[0.5 0.3/60 0.002 ];
options = optimoptions(@fmincon,'Algorithm','interior-point');
[pp,ss0_Monod_V23]=fmincon(@Adjust_Monod_V23,k00,[],[],[],[],lb_0,ub_0,[],options,data);
mse=ss0_Monod_V23/(length(data.ydata)-3); %el valor de 3 corresponde a los valores estimados

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
N_0=0.157;
P_0=0.006;
Ac_0=0;
inic=[S_0 N_0 P_0 Ac_0];

options = odeset('RelTol',1e-4,'AbsTol',1e-4);%,'NonNegative',[1:20],'Stats','off'); %RelTol=error relativo--- AbsTol=error absoluto

[t_Monod_V23,x_Monod_V23]= ode15s(@(t,x)model_Monod(t,x,parametro,pp), time_Monod_V23, inic,options);

Biomass_model_Monod_V23=x_Monod_V23(:,2);
Substrate_model_Monod_V23=x_Monod_V23(:,1);
Product_model_Monod_V23=x_Monod_V23(:,3);
Acetate_model_Monod_V23=x_Monod_V23(:,4);

%% figure 
figure(1),clf 
subplot(1,2,1)
plot(time_Monod_V23/60,Biomass_obs_Monod_V23,'o',time_Monod_V23/60,Biomass_model_Monod_V23,'r','LineWidth',2)
subplot(1,2,2)
plot(time_Monod_V23/60,Substrate_obs_Monod_V22,'o',time_Monod_V23/60,Substrate_model_Monod_V23,'r','LineWidth',2)



