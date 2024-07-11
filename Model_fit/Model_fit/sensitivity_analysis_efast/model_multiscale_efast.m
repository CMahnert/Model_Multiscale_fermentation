function [dydt]= model_multiscale_efast(t, y, rates, parametro,X,run_num)

		b= rates(1);
	dm= rates(2);
	kb= rates(3);
	ku= rates(4);
	f= rates(5);
    dn= rates(6);
   

      %Initial parameter
    Na=parametro (1); %molec/mol
    PMgluc=parametro (2); %g/mol
    mpp=parametro(3);%gdw/cell
	thetar= parametro(4);
	s0= parametro(5);
	gmax= parametro(6);
	thetax= parametro(7);
	Kt= parametro(8);
	M= parametro(9);
	Km= parametro(10);
	vm= parametro(11);
	nx= parametro(12);
	Kq= parametro(13);
	vt= parametro(14);
	wr= parametro(15);
	wq= parametro(16);
	nq= parametro(17);
	nr= parametro(18);
     V=parametro (19);
    yy=parametro (20);
    np=parametro(21);
 ns=parametro (22);
 n_xAc=parametro(23); 
theta_xAc=parametro(24);
Kcat_Ac=parametro(25);
Km_Ac=parametro(26);
Kcat_Ac_in=parametro(27);
Km_Ac_in=parametro(28);

%Parámetros a estimar
 Kgamma=X(1);
 wp=X(2);
 wAc=X(3);
 wer=X(4);
 nf=X(5);
 we=X(6);
    

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

mc=1e-13*exp(36.904*lam);



dydt(size(x,1),1)= 0;

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

    

