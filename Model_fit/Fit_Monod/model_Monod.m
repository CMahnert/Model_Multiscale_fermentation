function [dxdt]=model_Monod(t,x,parametro,pp)

dn=parametro(1);
S0=parametro(2);
yps=parametro(3);
yas=parametro(4);

yxs=pp(1);
umax=pp(2);
Km_S=pp(3);

%Variable del modelo 
S=x(1);
N=x(2);
P=x(3);
Ac=x(4);

u=umax*S/(Km_S+S);

dxdt(size(x,1),1)= 0;

dxdt(1)=-u*N/yxs;
dxdt(2)=u*N;
dxdt(3)=u*N/yps;
dxdt(4)=u*N/yas;
