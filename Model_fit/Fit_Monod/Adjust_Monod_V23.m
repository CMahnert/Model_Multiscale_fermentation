function ss=Ajuste_Monod_V23(pp,data)

global J_suma_sustrato_Monod_V23 J_suma_biomasa_Monod_V23


time=data.ydata(:,1);
biomass_obs=data.ydata(:,2);
substrate_obs=data.ydata(:,3);
producto_obs=data.ydata(:,4);
acetato_obs=data.ydata(:,5);

dn=0.01;
S0=4.8;
yps=0.5;
yas=0.5;

parametro=[dn S0 yps yas];

%Condición inicial 

S_0=4.8;
N_0=0.157;
P_0=0.006;
Ac_0=0;

inic=[S_0 N_0 P_0 Ac_0];

options = odeset('RelTol',1e-4,'AbsTol',1e-4);%,'NonNegative',[1:20],'Stats','off'); %RelTol=error relativo--- AbsTol=error absoluto

[t_Monod_V23,x_Monod_V23]= ode15s(@(t,x)model_Monod(t,x,parametro,pp), time, inic,options);

Biomass_model_Monod_V23=x_Monod_V23(:,2);
Substrate_model_Monod_V23=x_Monod_V23(:,1);
Product_model_Monod_V23=x_Monod_V23(:,3);
Acetate_model_Monod_V23=x_Monod_V23(:,4);


J_biomass_Monod_V23=(Biomass_model_Monod_V23-biomass_obs).^2;
J_substrate_Monod_V23=(Substrate_model_Monod_V23-substrate_obs).^2;
J_suma_biomass_Monod_V23=sum(J_biomass_Monod_V23);
J_suma_substrate_Monod_V23=sum(J_substrate_Monod_V23);

A=2*(Biomass_model_Monod_V23-biomass_obs).^2+1*(Substrate_model_Monod_V23-substrate_obs).^2;

ss=sum(A);
end


