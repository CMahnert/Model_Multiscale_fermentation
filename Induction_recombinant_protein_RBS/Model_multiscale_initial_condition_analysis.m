function [dxdt,lam]= Model_multiscale_initial_condition_analysis(t, x, rates, parametro)

	b= rates(1);
	dm= rates(2);
	kb= rates(3);
	ku= rates(4);
	f= rates(5);
    ds= rates(6);
    dn= rates(7);

  Na=parametro (1); 
    PMgluc=parametro (2);
    mpp=parametro(3);
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

Kgamma=parametro(29);
 wp=parametro(30);
 wAc=parametro(31);
 wer=parametro(32);
nf=parametro(33);
we=parametro(34);

 	rmr= x(1);
	em= x(2);
	rmp= x(3);
    rmAc=x(4); 	
    rmq= x(5);
	rmt= x(6);
	et= x(7);
	rmm= x(8);
	mt= x(9);
	mm= x(10);
	q= x(11);
	p= x(12);
    Ac=x(13); 
    si= x(14);
	mq= x(15);
	mp= x(16);
    mAc=x(17); 
    mr= x(18);
	r= x(19);
	a= x(20);
    
   	gamma= gmax*a/(Kgamma + a);
	ttrate= (rmq + rmr + rmt + rmm)*gamma;
	lam= (ttrate+rmp*gamma+rmAc*gamma)/M;
    nucat= em*vm*si/(Km + si);
    nuimp=et*vt*s0/(Kt+s0); 
    nucat_Ac=Ac*vm*si/(Km+si);
  
mc=1e-13*exp(36.904*lam);

dxdt(size(x,1),1)= 0;

	dxdt(1)= +kb*r*mr-ku*rmr-gamma/nr*rmr-f*rmr-lam*rmr;
	dxdt(2)= +gamma/nx*rmm-lam*em;
	dxdt(3)= +kb*r*mp-ku*rmp-gamma/np*rmp-f*rmp-lam*rmp;
	
    dxdt(4)=+kb*r*mAc-ku*rmAc-gamma/n_xAc*rmAc-f*rmAc-lam*rmAc;
    
    dxdt(5)= +kb*r*mq-ku*rmq-gamma/nx*rmq-f*rmq-lam*rmq;
	dxdt(6)= +kb*r*mt-ku*rmt-gamma/nx*rmt-f*rmt-lam*rmt;
	dxdt(7)= +gamma/nx*rmt-lam*et;
	dxdt(8)= +kb*r*mm-ku*rmm-gamma/nx*rmm-f*rmm-lam*rmm;
	dxdt(9)= +(we*a/(thetax + a))+ku*rmt+gamma/nx*rmt-kb*r*mt-dm*mt-lam*mt;
	dxdt(10)= +(wer*a/(thetax + a))+ku*rmm+gamma/nx*rmm-kb*r*mm-dm*mm-lam*mm;
	dxdt(11)= +gamma/nx*rmq-lam*q;
	dxdt(12)= +gamma/np*rmp-lam*p;

    dxdt(13)=+gamma/n_xAc*rmAc-lam*Ac;

	dxdt(14)= +nuimp-nucat-nucat_Ac-lam*si;
	dxdt(15)= +(wq*a/(thetax + a)/(1 + (q/Kq)^nq))+ku*rmq+gamma/nx*rmq-kb*r*mq-dm*mq-lam*mq;
	dxdt(16)= +(wp*a/(thetax + a))+ku*rmp+gamma/np*rmp-kb*r*mp-dm*mp-lam*mp;
	
    dxdt(17)=+(wAc*a/(theta_xAc + a))+ku*rmAc+gamma/n_xAc*rmAc-kb*r*mAc-dm*mAc-lam*mAc;
    
    dxdt(18)= +(wr*a/(thetar + a))+ku*rmr+gamma/nr*rmr-kb*r*mr-dm*mr-lam*mr;
	dxdt(19)= +ku*rmr+ku*rmt+ku*rmm+ku*rmp+ku*rmAc+ku*rmq+gamma/nr*rmr+gamma/nr*rmr+gamma/nx*rmt+gamma/nx*rmm+gamma/np*rmp+gamma/n_xAc*rmAc+gamma/nx*rmq-kb*r*mr-kb*r*mt-kb*r*mm-kb*r*mp-kb*r*mAc-kb*r*mq-lam*r; %ribosoma

    dxdt(20)= +(ns*nucat+nf*(nucat_Ac)-ttrate-np*rmp/np*gamma-n_xAc*rmAc/n_xAc*gamma-lam*a);
  