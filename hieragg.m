classdef hieragg < handle
    properties
        % Parameters
        grid_struct % grid structue, obtained with the grid_struct class
        H % steps ahead for the receding horizon optimization
        t % prox constant
        DSO_ref % DSO reference profile 
        iter % current iteration 
        % Profiles
        pbg % buying price from the grid (N x 1)
        psg % selling price to the grid (N x 1)
        pbg_H % buying price from the grid in the current horizon (H x 1)
        psg_H % selling price from the grid in the current horizon (H x 1)
        DSO_ref_H %reference power profile in the current horizon (H x 1)
        tol % primal dual residual tolerance for convergence
        result_scheleton
        
    end
    
    methods
        function obj = hieragg(grid_struct,pbg,psg,H,t,DSO_ref,tol)
            obj.grid_struct = grid_struct.grid_hierarchy;
            obj.pbg = pbg;
            obj.psg = psg;
            obj.H = H;
            obj.t = t;
            obj.DSO_ref = DSO_ref;
            obj.tol = tol;
            obj.iter = 1;
            obj.result_scheleton = grid_struct.grid_scheleton;
        end
        
        function [results,history] = mpc_step(obj,timestep,opts)
            % perform an MPC step with the hierarchical distributed implicit optimization
            % Inputs: -timestep: current time step of simulation
            %         -opts: option structure containing {'fb_split','selfish'} 
            %                fb_split: if true project onto power
            %                selfish: if true perform selfish iterations
            %                constraints of first level
            
            % Inizialize hierarcy
            forecast_only = true;
            [~,obj.grid_struct] = initialize_slot(obj,obj.grid_struct,timestep,forecast_only);
            forecast_only = false;
            [~,obj.grid_struct] = initialize_slot(obj,obj.grid_struct,timestep,forecast_only);
            obj.pbg_H = obj.pbg(timestep:timestep+obj.H-1);
            obj.psg_H = obj.psg(timestep:timestep+obj.H-1);
            obj.DSO_ref_H = obj.DSO_ref(timestep:timestep+obj.H-1);
            
            % Iterate until convergence
            plot_opts.do_plots = false;
            plot_opts.only_top = false;
            plot_opts.variables = {'P','V','E','D'};
            init_vals  = plot_results(obj,plot_opts,obj.DSO_ref_H);
            plot_opts.h = figure;
            history = [];
            history.primal_res = 2*obj.tol;
            history.conv_factor = 2*obj.tol;
            % Get initial results for dysplaying purposes
%             retrieve_power_and_voltage(obj,obj.grid_struct);
%             initial_values = plot_results(obj,plot_opts,obj.DSO_ref_H);
            t0 = tic;
            while history.primal_res(obj.iter) > obj.tol && history.conv_factor(obj.iter)>obj.tol && obj.iter<200

                % Do the backward pass
                obj.grid_struct = backward(obj,obj.grid_struct,[],[],[],[]);
                
                % Propagate objective in the hierarchy, do the forward pass
                forward(obj,obj.grid_struct,timestep,opts);
                
                % Plot something from  the first level
%                 plot_iter(obj,h1,y0,obj.DSO_ref_H)
                if obj.iter>1 
                    plot_opts.do_plots = true;
                end
                results = plot_results(obj,plot_opts,obj.DSO_ref_H,init_vals,history);
                obj.iter = obj.iter+1;
                history.primal_res(obj.iter) = mean([results.tot_primal_res_s(:);results.tot_primal_res_v(:)]);
                history.violations(obj.iter) = mean([results.tot_violations_p(:);results.tot_violations_v(:)]);
                history.conv_factor(obj.iter) = (history.primal_res(obj.iter)-history.violations(obj.iter));
                %                 obj.t = obj.t/1.01;
            end
            % Group optimization output in a convenient structure
            history.primal_res = history.primal_res(2:end);
            history.violations = history.violations(2:end);
            history.conv_factor = history.conv_factor(2:end);
            results.n_subworkers = obj.grid_struct.n_subworkers;
            results.time = toc(t0);
            results.time_per_agent = results.time/results.n_subworkers;
            results.init_vals = init_vals;
        end
        
        function [n_subworkers,level] = initialize_slot(obj,level,timestep,forecasts_only)
            % Initialize solutions for the current time slot. All agents
            % perform their forecasts for the uncontrolled loads
            % Each agent retrieves its own power forecast, at current
            % timestep. Find number of subworkers for each queen. Use recursion.
            labels = fieldnames(level);
            n_subworkers = 0; % descendants of current branching node
            if ismember('workers',labels) % If thre are Workers, retrieve solution
                for i=1:length(level.workers)
                    level.workers(i).do_forecasts(timestep); % Update Pm_forecast object in each Worker
                    n_subworkers = length(level.workers);
                end
                if forecasts_only
                    level.Branch.do_forecasts(timestep); % Update P_unc_forecast and V_forecast in the Branch
                else
                    level.Branch.do_forecasts(timestep); % Update P_unc_forecast and V_forecast in the Branch
                    [Pinit,Vinit] = retrieve_power_and_voltage(obj,level); % initialize global variable representing sum of P workers
                    level.Branch.y_s = Pinit;
                    level.Branch.y_v = Vinit;
                    level.Branch.lambda_s = zeros(size(level.Branch.y_s)); % initialize lambda of sum of P workers
                    level.Branch.lambda_v = zeros(size(level.Branch.y_s)); % initialize lambda of sum of P workers
                end
            end
            q_labels = labels(~cellfun(@isempty,(regexp(labels, 'q\d')))); % match any field with name as q* where * is a number
            if ~isempty(q_labels) % if there are Branchs, go on recursively
                for i=1:length(q_labels)
                    [n_subworkers_q,level.(q_labels{i})] = initialize_slot(obj,level.(q_labels{i}),timestep,forecasts_only);
                    n_subworkers = n_subworkers+ n_subworkers_q;
                end
                n_branches = length(q_labels);
            else
                n_branches = 0;
            end
            % write the number of sub-Workers in the current level (workers in the level + number of workers in the sub-Branchs)
            level.n_subworkers = n_subworkers;
            % number of branches in the level
            level.n_branches = n_branches;
        end
        
        function forward(obj,level,timestep,opts)
            % Perform the forward passage of the hierarchical optimization.
            % Decompose the overall objective among the different hierarchy
            % levels and Workers that populate these levels. Use recursion.
            % Inputs: -timestep: current timestep of the simulation
            %         -lambda_p,lambda_v: Lagrange multipliers from the
            %         upper level, associated to power and voltage
            %         constraints
            
            
            labels = fieldnames(level);
            
            % Perform a minimization step for the Workers in this sub-Branch
            if ismember('workers',labels) % If thre are Workers, retrieve solution
                % compute the Workers forward pass
                for i=1:length(level.workers)
                    if obj.iter > 0
                        x = level.workers(i).Pm_opt;                                    % agent's actions at previous iteration
                        Ref_s = -(level.R_s -repmat(x,1,size(level.R_s,2)));            % reference for sum
                        K_coeffs = [level.Kp_l,level.Kp(1+level.n_branches + i)];       % sensitivity coefficients of current profile on all constrained points in its genealogy
                        Ref_v = -(level.R_v -repmat(x,1,size(level.R_v,2)).*K_coeffs);  % reference for voltage
                        Ref = [Ref_s,Ref_v];
                        Coeffs =  [ones(1,size(Ref_s,2)),K_coeffs];                     % multiplicative coefficients of the control variable 
                    else
                        Ref = zeros(obj.H,1);
                        Coeffs = 1;
                    end
                    P_hat = level.workers(i).Pm_forecast;
                    [~,~,~,noneconomic_cost,cost] = level.workers(i).solve_cooperative_admm_yalmip(P_hat,obj.pbg_H,obj.psg_H,Ref,Coeffs,opts.selfish); % U Pm_opt and e_opt are automatically saved in the Worker at each iteration
                    level.workers(i).noneconomic_cost = noneconomic_cost;
                    level.workers(i).cost = cost;
                end
            end
            
            % compute the Branch forward pass with projections
            [P,V] = retrieve_power_and_voltage(obj,level);
            nsw = level.n_subworkers;
            if level.n_level == 0
                % proxy of system-level objective
                level.Branch.y_s = prox_agg(obj,P,level.Branch.lambda_s,nsw);
                if opts.fb_split
                    level.Branch.y_s = min(max(level.Branch.y_s, level.Pmin), level.Pmax);
                end
            else
                % project onto constraint sets
                level.Branch.y_s = min(max(P +level.Branch.lambda_s*nsw, level.Pmin), level.Pmax);
            end
            
            level.Branch.y_v = min(max(V+level.Branch.lambda_v*nsw, level.Vmin), level.Vmax);
            % Detect if there are sub-Branchs in this level
            q_labels = labels(~cellfun(@isempty,(regexp(labels, 'q\d')))); % match any field with name as q* where * is a number
            if ~isempty(q_labels) % if there are Branchs, go on recursively
                for i=1:length(q_labels)
                    forward(obj,level.(q_labels{i}),timestep,opts);
                end
            end
        end
        
        function level = backward(obj,level,Rs_upper,Rv_upper,Kp_upper,nsw_norm_upper)
            % Do the backward step, in which the Lagrange multipliers are
            % updated bottom-up. This does not include the Lagrangian of
            % the upper level consensus constraint, which is updated
            % separately. Use recursion
            
            labels = fieldnames(level);
            nsw = level.n_subworkers;
            % Update power and voltage constraint lambdas for the current level
            [P,V] = retrieve_power_and_voltage(obj,level); % sum of total power at the current Branch
            
            % Dual variables update
            level.Branch.r_s = P-level.Branch.y_s;
            level.Branch.r_v = V-level.Branch.y_v;
            level.Branch.lambda_s = level.Branch.lambda_s + level.Branch.r_s/nsw/obj.t;
            level.Branch.lambda_v = level.Branch.lambda_v + level.Branch.r_v/nsw/obj.t;
            
            % assign reference profiles
            rsi = level.Branch.r_s/nsw+level.Branch.lambda_s;
            rvi = level.Branch.r_v/nsw+level.Branch.lambda_v;
            level.R_s = [Rs_upper,rsi];
            level.R_v = [Rv_upper,rvi];
            level.Kp_l = Kp_upper;
            level.nsw_norm = [nsw_norm_upper,nsw];
            
            q_labels = labels(~cellfun(@isempty,(regexp(labels, 'q\d')))); % match any field with name as q* where * is a number
            if ~isempty(q_labels) % if there are Branchs, go on recursively
                for i=1:length(q_labels)
                    % normalize the voltage reference with the coefficient
                    % of the ith branching node (the first is self-influence)
%                     R_level = [level.R,level.Rv/level.Kp(i+1)];
%                     R_level = [level.R,level.Rv];
%                     R_level = level.R;
                    level.(q_labels{i}) = backward(obj,level.(q_labels{i}),level.R_s,level.R_v,[Kp_upper,level.Kp(1+i)],nsw);
                end
            end
        end
        
        function [P,V,P_subqueens_mat,P_workers_mat] = retrieve_power_and_voltage(obj,level)
            
            % Sum all the power in the current and lower levels in the
            % hierarchy. Use recursion
            labels = fieldnames(level);
            
            % Sum total power of the workers
            if ismember('workers',labels)
                if ~isempty([level.workers.Pm_opt])
                    P_workers_mat = [level.workers.Pm_opt]; % Optimized total power profile of workers at current iteration
                else
                    for i=1:length(level.workers)
                        level.workers(i).Pm_opt = level.workers(i).Pm_forecast;
                        P_workers_mat = [level.workers.Pm_forecast];% Forecasted total power profile of workers at current iteration 
                    end
                end
                P_workers = sum(P_workers_mat,2);
                if isempty(P_workers)
                    P_workers = zeros(obj.H,1);
                end
                
            else
                P_workers = 0;
                P_workers_mat = [];
            end
            
            % Sum total power from the Branchs, recursively
            q_labels = labels(~cellfun(@isempty,(regexp(labels, 'q\d')))); % match any field with name as q* where * is a number
            if ~isempty(q_labels) % if there are Branchs, go on recursively
                P_subqueens_mat = zeros(obj.H,length(q_labels));
                for i=1:length(q_labels)
                    P_subqueens_mat(:,i) = retrieve_power_and_voltage(obj,level.(q_labels{i}));
                end
                P_subqueens = sum(P_subqueens_mat,2);
            else
                P_subqueens = 0;
                P_subqueens_mat = [];
            end
            
            % sum of power in the Branch = P workers + P from subqueens + P from non-workers
            P = P_workers + P_subqueens + level.Branch.P_unc_forecast;
            
            % compute V approximation from active power only
            V = level.Branch.V_forecast + [P,P_subqueens_mat,P_workers_mat]*level.Kp;
        end
        
        function y = prox_agg(obj,sum_x,lambda,nsw)
            y = (obj.t*obj.DSO_ref_H+lambda*nsw+sum_x)/(obj.t +1);
        end
        
        function selfish_step(obj,timestep)
            [~,obj.grid_struct] = initialize_slot(obj,obj.grid_struct,timestep);
            lambda = zeros(obj.H,1);
            obj.pbg_H = obj.pbg(timestep:timestep+obj.H-1);
            obj.psg_H = obj.psg(timestep:timestep+obj.H-1);
            selfish = true;
            forward(obj,obj.grid_struct,timestep,lambda,selfish);
        end
        
        function plot_iter(obj,h1,y0,P_target)
            rmse = @(x,y) mean((x-y).^2)^0.5;
            [P,V] = retrieve_power_and_voltage(obj,obj.grid_struct);
            figure(h1);
            primal_norm = mean((obj.grid_struct.Branch.y_s-P).^2)^0.5;
            fprintf('\nDual resid norm = %0.3e',primal_norm);
            subplot 221
            plot(P);hold on; plot(obj.grid_struct.Branch.y_s); plot(obj.grid_struct.Pmax);plot(y0)
            hold off;
            legend('P','y','P max','y0')
            subplot 222
            plot([obj.grid_struct.Branch.lambda_s,obj.grid_struct.Branch.lambda_v]);
            legend('Lambda S','Lambda V')
            title('Lambda')
            subplot 223
            plot(obj.grid_struct.Vmax); hold on; plot(obj.grid_struct.Vmin);
            plot(obj.grid_struct.Branch.V_forecast); plot(V,'--'); hold off;
            legend('V max','V min','V forecast','V forecast + actions')
            subplot 224
            plot([obj.grid_struct.workers.e_opt]);
            suptitle(sprintf('RMSE = %0.2e',rmse(P,P_target)))
            drawnow;
        end
    end
end