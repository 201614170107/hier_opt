% main
% Sovle the traking problem of the aggregated profile of a prosumer
% community, in a distributed way through proximal algorithm.

clear 
close all
rng(0)


%% |--------------------------- Set parameter ---------------------------|
    
H = 24*2;               % Timesteps per horizon, considering dayly horizons
t_prox = 1.5;             % proxy constant and inverse of lambda step
selfish = false;        % agents are cooperative
rand_range = 10;        % maximum number of agents per branch
n = 30;                 % number of days in generated data
N = n*H;                % total number of timesteps
DSO_ref = zeros(N,1);   % tracking profile - 0s = quadratic peak shaving 
pb = 0.2*ones(N,1);     % buying price profile 
ps = 0.05*ones(N,1);    % selling price profile

% battery parameters
eta_in = 0.9;               % charging efficiency
eta_out = 0.9;              % discharging efficiency
ts = 24*3600/H;             % time sampling of simulation [s]
tau_sd = 1/(365*24*3600);   % tau of self discharge
lifetime = 20;              % battery lifetime 
DOD = 0.9;                  % depth of discharge
Nc = 3000;                  % number of cicle before end of life
E_0 = 10;                   % nominal capacity [kWh]
c = 1;                      % C factor [kW/kWh]


%% |------------------- Prepare hierarchical structure -------------------| 

grid_scheleton = [];

% top level
grid_scheleton.V_ref = 1; 
grid_scheleton.P_ref = 1; 

% second level
grid_scheleton.q1.V_ref = 0.5; 
grid_scheleton.q1.P_ref = 0.5; 


% third level
grid_scheleton.q1.q1.V_ref = 0.25; 
grid_scheleton.q1.q1.P_ref = 0.25; 


grid_topology = []; % do not specify grid topology 
slack_b = false; % add a slack battery in each branch. This ensure constraints
                 % are always satisfied. If true, Increase t_prox appropriately
                 % to remove bumpy behavior in the convergence.

defauLt_worker_pars = struct('eta_in',eta_in,'eta_out',eta_out,'ts',ts,...
    'tau_sd',tau_sd,'lifetime',lifetime,'DOD',DOD,'Nc',Nc,'E_0',E_0,...
    'c',c,'h',H,'t',t_prox,'f_method','perfect','selfish',selfish);

% create and print grid struct
gridstr = rand_grid_struct(grid_scheleton,rand_range,grid_topology,...
          defauLt_worker_pars,H,N,slack_b);

% create hierarchical aggregator
tol = 1e-3;
timestep = 1;
opts.selfish = false;
opts.fb_split = true;

hi_aggr = hieragg(gridstr,pb,ps,H,t_prox,DSO_ref,tol);
results = hi_aggr.mpc_step(timestep,opts);


save('one_case_results.mat','results')