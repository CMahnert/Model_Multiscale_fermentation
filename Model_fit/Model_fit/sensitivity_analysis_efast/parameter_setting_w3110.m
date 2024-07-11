%% Parameters %%

pmin=[1000, % Kgamma 
100, %wp
0.1, % wAc
1, % wer
0.1, % nf
2 %we
];%we

pmax=[3e5, % Kgamma 
1500, %wp
1, % wAc
4, % wer
0.2, % nf
5
]; %we

% Parameter Labels 
efast_var={'Kgamma', '\wAc', 'wer', ...
    'nf','we'};

% PARAMETER BASELINE VALUES
Kgamma=2e2; 
wp=1300;
wAc=0.51;
wer=4.3;
nf=0.1;
we=4.46;
%% TIME SPAN OF THE SIMULATION
t_end=16*60; % length of the simulations
tspan=(0:1:t_end);   % time points where the output is calculated
time_points=[100 16*60]; % time points of interest for the US analysis

% % INITIAL CONDITION FOR THE ODE MODEL
% T0=1e3;
% T1=0;
% T2=0;
% V=1e-3;
% 
% y0=[T0,T1,T2,V];
% 
% % Variables Labels
 y_var_label={'rmr','em','rmp','rmAc','rmq','rmt','et','rmm','mt','mm','q','p','Ac','si','mq','mp','mAc','mr','r','a','N','biomass','S','Protein','Acetate'};

%  Kgamma=X(run_num,1);
% wp=X(run_num,2);
% wAc=X(run_num,3);
% wer=X(run_num,4);
% nf=X(run_num,5);
% we=X(run_num,6);

