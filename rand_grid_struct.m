classdef rand_grid_struct < handle
    % Define a grid_struct object. It is typically defined through the scheleton structure, and then randomly populated
    % with default parameters. The active and reactive sensitivity matrices are thought have nonzero elements only in the Queen-Workers connecions
    % if not differently specified.
    properties
        grid_scheleton   % structure describing the hierarchical relation between queens. This is used to build grid_hierarchy
        grid_hierarchy   % structure containing the hierarchy of queens and workers, including topology if available
        grid_topology    % boolean, if true, a grid topology in terms of sensitivity coefficients matrix must be provided for each Queen
        rand_range % max number of workers per queen (in case the grid struct is randomly populated)
        default_worker_pars % structure containing parameters for the default worker (from which other workers are obtained with random_populate_levels method)
        H % horizon of the simulation
        N % number of horizons in the simulation
        f_method %  forecast method, inherited by workers f_method at initialization
        n_levels %  number of levels
        slack_b % if true, add a slack battery at each node
    end
    
    methods
        function obj = rand_grid_struct(grid_scheleton,rand_range,grid_topology,default_worker_pars,H,N,slack_b)
            % Init function
            if nargin<7
                slack_b = false;
            end
            fprintf('\nBe sure each Queen in the structure has a name in the form qx where x is a number'); 
            obj.slack_b = slack_b;
            % check the scheleton            
            obj.grid_scheleton = grid_scheleton;
            check_scheleton(obj,grid_scheleton)
            
            obj.grid_topology = grid_topology;
            obj.rand_range = rand_range;
            obj.default_worker_pars = default_worker_pars;
            obj.H = H; % horizon of the simulation in timesteps
            obj.N = N; %number of horizons in the simulation, needed to randonmly create Worker's profiles
           
            % check if default_worker_pars contains all the required fields
            dwp_fields = fieldnames(default_worker_pars);
            dwp_list = {'eta_in','eta_out','ts','tau_sd','lifetime','DOD','Nc','E_0','c','h','t','f_method','selfish'};
            
            if (length(dwp_list) ~= length(dwp_fields)) || ~all(ismember(dwp_list,dwp_fields))
                error('Some fields are missing in the default workers pars list');
                fprintf('\ndefault_worker_pars must contain:');
                fprintf('%s',dwp_list{:});
            end
            
            obj.f_method = default_worker_pars.f_method;
            n_level = 0; % starting level = root node level
            % random populate the scheleton
            if rand_range > 0
                [obj.grid_hierarchy,obj.n_levels] = random_populate_levels(obj,grid_scheleton,n_level);
            end
            
            % print scheleton and structure
            fprintf('\n%s',repmat('#',10))
            fprintf('\n \t GRID SCHELETON')
            fprintf('\n%s',repmat('#',10))
            disp_struct(obj,grid_scheleton)
            fprintf('\n%s',repmat('#',10))
            fprintf('\n \t GRID HIERARCHY')
            fprintf('\n%s',repmat('#',10))
            disp_struct(obj,obj.grid_hierarchy)
        end
        
        
        function [level,n_level_max] = random_populate_levels(obj,level,n_level)
            % randomly populate the grid based on the description in
            % grid_hierarchy. Use recursion.
            % Input: -level: the structure for current level
            %        -n_level: current level number
%             rng(0);
            level.n_level = n_level;
            n_level = n_level+1;
            labels = fieldnames(level); 
            q_labels = labels(~cellfun(@isempty,(regexp(labels, 'q\d')))); % match any field with name as q* where * is a number
            if ~isempty(q_labels) % if there are Queens, go on recursively
                level = random_populate_workers(obj,level,length(q_labels)+1); % Consider n+(1+nq) connections
                level.Branch = Branch(obj.H,level.V_ref,level.P_ref,obj.f_method); % create a Queen object containing forecasters for voltages and powers
                % go on recursively
                for i=1: length(q_labels)
                    [level.(q_labels{i}),n_level_max(i)] = random_populate_levels(obj,level.(q_labels{i}),n_level);
                end
                n_level_max = max(n_level_max);
            else    % if we are in a terminal queen, populate with random workers
                level = random_populate_workers(obj,level,1);
                level.Branch = Branch(obj.H,level.V_ref,level.P_ref,obj.f_method); % create a Queen object containing forecasters for voltages and powers
                n_level_max = 0;
            end
           n_level_max = max(n_level,n_level_max); 
        end
        
        function level = random_populate_workers(obj,level,n_extra)
            % Randomly populate the level with workers.
            % Inputs: -level: the level structure to populate
            %         -n_extra: number of extra entities to consider when buildin the sensitivity matrices.
            %                   If we are in a terminal Queen, n_extra = 1 (the Queen)
            
            dwp = obj.default_worker_pars;
            Ntot = obj.N;
            n_days = floor(Ntot/obj.H);
            Hor = obj.H;
            num_w = randperm(obj.rand_range,1);
            norm_factor = level.P_ref/(num_w); % each agent is dimensioned to have a reasonable power profile and storage w.r.t. the nominal power of the level
            
            use_real_data = 1;
            if use_real_data
                load('load_pv_data')
                for i=1:num_w
                    Pm(:,i) = mean(reshape(generate_load_and_pv(P_Load,P_PV, 10, 1, 1, [2 10],90),2,[]));
                end
                violPerc=0.1;
                gain = fmincon(@(x) (violPerc-get_violation_perc(x*Pm,level.P_ref))^2,0.001,[],[],[],[],0,10);
                
                for i=1:num_w
                    level.workers(i) = battery(dwp.eta_in,dwp.eta_out,dwp.ts,dwp.tau_sd,dwp.lifetime,dwp.DOD,dwp.Nc,level.P_ref/(num_w),dwp.c,dwp.h,dwp.t,gain*Pm(:,i),dwp.f_method,dwp.selfish); % 3 hours autonomy if the worker had average power consumption
                end
            else
                for i=1:num_w
                    Pm(:,i) =  (sin((2*pi*(1:Ntot)')/obj.H +randn(1)*pi/2)+0.2*randn(Ntot,1)+0.2)*norm_factor;  % pseudo random power profile
                    level.workers(i) = battery(dwp.eta_in,dwp.eta_out,dwp.ts,dwp.tau_sd,dwp.lifetime,dwp.DOD,dwp.Nc,norm_factor,dwp.c,dwp.h,dwp.t,Pm(:,i),dwp.f_method,dwp.selfish); % 3 hours autonomy if the worker had average power consumption
                end
            end
            labels = fieldnames(level); 
            q_labels = labels(~cellfun(@isempty,(regexp(labels, 'q\d'))));
            if obj.slack_b % && isempty(q_labels)
                num_w = num_w+1;
                level.workers(end+1) =  battery(dwp.eta_in,dwp.eta_out,dwp.ts,dwp.tau_sd,dwp.lifetime,dwp.DOD,dwp.Nc,3*norm_factor,dwp.c,dwp.h,dwp.t,zeros(Ntot,1),dwp.f_method,dwp.selfish,true);
            end
            % give random values to sensitivity coefficients. [k_auto,k_subqueens,k_workers]
            k_v_max = 5e-2*level.V_ref/level.P_ref/num_w; % if P = Pmax 
            level.Kp = -k_v_max*0.5*(1+0.2*rand(num_w+n_extra,1)); % active power sensitivity matrix. 
            level.Kq = -k_v_max*0.5*(1+0.2*rand(num_w+n_extra,1)); % reactive power sensitivity matrix. 
            
            
            % Set the per unit reference power and reference voltage
%             level.P_ref = 20e3; % reference power [W]
%             level.V_ref = 230; % reference voltage [V]
            
            % Set grid constraints
            level.Pmax = 1.0*ones(Hor,1)*level.P_ref;
            level.Pmin = -1.0*ones(Hor,1)*level.P_ref;
            level.Vmax = 1.01*ones(Hor,1)*level.V_ref;
            level.Vmin = 0.99*ones(Hor,1)*level.V_ref;
            
            % Set uncontrolled loads in current Queen
            level.P_unc = zeros(Hor,1);
            level.active_constraints = ones(obj.H,4);
            
            % normalized weights for objective function reallocation
            rand_vect = rand(num_w,1)+0.5;
            level.alphas = rand_vect./sum(rand_vect); 
            
            
        end
        function disp_struct(obj,x,indent)
            % Custom function for displaying a structure. Use recursion.
            % Input: -x: a nested structure
            %        -indent: do not specify this input. For recursion
            %        purposes only.
            if nargin<3
                indent='';
            end
            if isstruct(x)
                labels = fieldnames(x);
                for i=1:length(labels)
                    if isstruct(x.(labels{i}))
                        fprintf('\n%s%s',[indent,labels{i}])
                        disp_struct(obj,x.(labels{i}),strcat(indent,'---'))
                    else
                        fprintf('\n%s%s%s',[indent,labels{i},':'])
                        fprintf('\t[%ix%i %s] ',[size(x.(labels{i}),1),size((x.(labels{i})),2),class(x.(labels{i}))])
                    end
                end
            end
            fprintf('\n')
        end
        
        function check_scheleton(obj,scheleton)
            if isstruct(scheleton)
                labels = fieldnames(scheleton);
                if ~all(ismember({'V_ref','P_ref'},labels))
                    error('Each Queen must have the fields {V_ref,P_ref}')
                end
                q_labels = labels(~cellfun(@isempty,(regexp(labels, 'q\d')))); % match any field with name as q* where * is a number
                if ~isempty(q_labels) % if we are not in a terminal queen
                    for i=1:length(q_labels)
                        check_scheleton(obj,scheleton.(q_labels{i}))
                    end
                end
            end
        end 
    end
end

