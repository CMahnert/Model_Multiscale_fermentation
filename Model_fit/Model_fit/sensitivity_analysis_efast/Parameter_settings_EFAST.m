%% Parameters %%

pmin=[100, % Kgamma 
0.01, %wp
0.01, % wAc
0.01, % wer
0.01, % nf
0.01]; %we

pmax=[2000, % Kgamma 
3000, %wp
1, % wAc
4, % wer
0.1, % nf
4]; %we

% Parameter Labels 
efast_var={'Kgamma', '\wAc', 'wer', ...
    'nf','we'};

% PARAMETER BASELINE VALUES
Kgamma=2.32E+05; 
wp=1300;
wAc=0.51;
wer=4.3;
nf=0.1;
we=4.46;
%% TIME SPAN OF THE SIMULATION
t_end=4000; % length of the simulations
tspan=(0:1:t_end);   % time points where the output is calculated
time_points=[2000 4000]; % time points of interest for the US analysis

y_var_label={'rmr','em','rmp','rmAc','rmq','rmt','et','rmm','mt','mm','q','p','Ac','si','mq','mp','mAc','mr','r','a','N','biomasa','s','PR','Ac_ext'};

