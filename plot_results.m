function result_scheleton = plot_results(hi_aggr,plot_opts,P_target,initial_vars,history)
% Get results stored in the hierarchical structure and plot them if
% requested.
% Inputs: hi_aggr: hierarchical aggregator object, see the related class
%         plot_opts: struct with the following fields: {'do_plots','variables'}
%                    do_plots: boolean, if true do plots
%                    variables: string containing a combination of the following
%                    letters: {'P','V','E','D'}. P: plot aggregated power, reference
%                    and limits. V: plot voltage and limits. E: plot workers and subqueens battery energies.
%                    D: plot the Lagrangian dual variables.
%         P_target: power target profile
%
% Outputs: result_scheleton: a copy of the grid_scheleton contained in the hi_aggr, populated with results


if nargin < 2
    plot_opts.do_plots = false;
end

if nargin < 3
    P_target = [];
end

if nargin < 4
    initial_vars = [];
end

if nargin < 5
    history = [];
end

% Populate the grid_scheleton with results
H = hi_aggr.H;
result_scheleton = hi_aggr.result_scheleton;
[~,~,result_scheleton] = retrieve_solutions(hi_aggr.grid_struct,result_scheleton,H,plot_opts.do_plots);

fprintf('\n Iter %i Mean primal residual tot= %0.2e, Mean violations tot = %0.2e',...
    [hi_aggr.iter,mean([result_scheleton.tot_primal_res_v(:);result_scheleton.tot_primal_res_s(:)]),mean([result_scheleton.tot_violations_p(:);result_scheleton.tot_violations_v(:)])]);
if plot_opts.do_plots
    plot_result(result_scheleton,plot_opts,P_target,initial_vars,history)
end

end

function [P,E,result_scheleton,P_all_workers,P_all_workers_forecasts] = retrieve_solutions(level,result_scheleton,H,do_plots)
% Sum all the power in the current and lower levels in the
% hierarchy. Use recursion

labels = fieldnames(level);

% Sum total power of the workers
if ismember('workers',labels)
    if ~isempty([level.workers.Pm_opt])
        P_workers = sum([level.workers.Pm_opt],2); % Optimized total power profile of workers at current iteration
        P_workers_mat = [level.workers.Pm_opt];
    else
        P_workers = sum([level.workers.Pm_forecast],2);% Forecasted total power profile of workers at current iteration 
        P_workers_mat = [level.workers.Pm_forecast];
    end
%     P_workers = sum([current_struct.workers.Pm_opt],2); % Optimized total power profile of workers at current iteration
    if isempty(P_workers)
        P_workers = zeros(H,1);
    end
    P_workers_mat_forecasts = [level.workers.Pm_forecast];
    e_workers_mat = [level.workers.e_opt];
else
    P_workers = 0;
    P_workers_mat = [];
    e_workers_mat = [];
end


% Sum total power from the Branchs, recursively
q_labels = labels(~cellfun(@isempty,(regexp(labels, 'q\d')))); % match any field with name as q* where * is a number
if ~isempty(q_labels) % if there are Branchs, go on recursively
    P_subqueens_mat = zeros(H,length(q_labels));
    e_subqueens_mat = P_subqueens_mat;
    tot_primal_res_v = 0;
    tot_primal_res_s = 0;
    tot_violations_v = 0;
    tot_violations_p = 0;
    for i=1:length(q_labels)
        [P_subqueens_mat(:,i),e_subqueens_mat(:,i),result_scheleton.(q_labels{i}),P_all_workers_i{i},P_workers_mat_forecasts_i{i}] = ...
            retrieve_solutions(level.(q_labels{i}),result_scheleton.(q_labels{i}),H,do_plots);
        tot_primal_res_v = tot_primal_res_v + abs(result_scheleton.(q_labels{i}).tot_primal_res_v);
        tot_primal_res_s = tot_primal_res_s + abs(result_scheleton.(q_labels{i}).tot_primal_res_s);
        tot_violations_v =  tot_violations_v + abs(result_scheleton.(q_labels{i}).tot_violations_v);
        tot_violations_p = tot_violations_p +abs(result_scheleton.(q_labels{i}).tot_violations_p);
    end
    P_subqueens = sum(P_subqueens_mat,2);
else
    P_subqueens = 0;
    P_subqueens_mat = [];
    e_subqueens_mat = [];
    P_all_workers_i{1} = [];
    P_workers_mat_forecasts_i{1} = [];
    tot_primal_res_s = 0;
    tot_primal_res_v = 0;
    tot_violations_p = 0;
    tot_violations_v = 0;
end

% Power and energy
P = P_workers + P_subqueens + level.Branch.P_unc_forecast;
if ~isempty(e_subqueens_mat) && ~isempty(e_workers_mat)
    E = sum(e_workers_mat,2)+sum(e_subqueens_mat,2);
elseif ~isempty(e_subqueens_mat)
    E = sum(e_subqueens_mat,2);
elseif ~isempty(e_workers_mat)
    E = sum(e_workers_mat,2);
else
    E = zeros(size(P));
end



result_scheleton.P = P;
result_scheleton.P_workers = P_workers_mat;
result_scheleton.P_subqueens = P_subqueens_mat;
result_scheleton.p_max = level.Pmax;
result_scheleton.p_min = level.Pmin;
result_scheleton.e_workers = e_workers_mat;
result_scheleton.E = E;
result_scheleton.e_subqueens = e_subqueens_mat;
P_all_workers = [P_workers_mat,[P_all_workers_i{:}]];
result_scheleton.P_all_workers = P_all_workers; 
P_all_workers_forecasts = [P_workers_mat_forecasts,[P_workers_mat_forecasts_i{:}]];
result_scheleton.P_all_workers_forecasts = P_all_workers_forecasts;
result_scheleton.Branch = level.Branch;

% Voltage
Pmat = [P,P_subqueens_mat,P_workers_mat];
V = level.Branch.V_forecast +Pmat*level.Kp;
result_scheleton.V = V;
result_scheleton.v_max = level.Vmax;
result_scheleton.v_min = level.Vmin;

% Primal residual
result_scheleton.primal_residual_s = abs(P-level.Branch.y_s);
result_scheleton.primal_residual_v = abs(result_scheleton.V-level.Branch.y_v);
result_scheleton.tot_primal_res_s = tot_primal_res_s+result_scheleton.primal_residual_s;
result_scheleton.tot_primal_res_v = tot_primal_res_v+result_scheleton.primal_residual_v;

% Violations
result_scheleton.violations_p = abs(P - min(max(P, level.Pmin), level.Pmax));
result_scheleton.violations_v = abs(V - min(max(V, level.Vmin), level.Vmax));
result_scheleton.tot_violations_p = tot_violations_p + result_scheleton.violations_p;
result_scheleton.tot_violations_v = tot_violations_v + result_scheleton.violations_v;


% Global variables and Lagrangian multipliers
result_scheleton.y_s = level.Branch.y_s;
result_scheleton.y_v = level.Branch.y_v;
result_scheleton.lambdas =  [level.Branch.lambda_s,level.Branch.lambda_v];
% if do_plots
% result_scheleton.fig_h = figure;
% end
end


function plot_result(result_scheleton,plot_opt,P_target,initial_vars,history)
% plot results stored in the grid_scheleton
if plot_opt.only_top
    plot_level(result_scheleton,plot_opt,P_target);
else
    plot_level(result_scheleton,plot_opt,P_target,initial_vars);
    labels = fieldnames(result_scheleton);
    q_labels = labels(~cellfun(@isempty,(regexp(labels, 'q\d')))); % match any field with name as q* where * is a number
    for i=1:length(q_labels)
        plot_result(result_scheleton.(q_labels{i}),plot_opt,P_target,initial_vars.(q_labels{i}),history);
    end
end
if ~isempty(history)
figure(plot_opt.h);
yyaxis('left')
plot([history.primal_res(2:end)',history.violations(2:end)'],'--');
ylabel('primal and violations')
yyaxis('right')
plot((history.primal_res(2:end)-history.violations(2:end)));
ylabel('normalized primal err-violations err')
end
end

function plot_level(result_scheleton,plot_opt,P_target,initial_vars)
if nargin < 3
    P_traget = zeros(result_scheleton.H,1);
end
% plot on the result_scheleton default figure (initialize during get_results call)
% figure(result_scheleton.fig_h);

figure(result_scheleton.Branch.h_f);
n_subplot = length(plot_opt.variables);
[r,c] = subplot_topology(n_subplot);
for i=1:n_subplot
    subplot(r,c,i)
    variable_i = plot_opt.variables(i);
    if strcmp(variable_i,'P')
        plot(initial_vars.P);hold on;
        plot([result_scheleton.P,result_scheleton.Branch.y_s,P_target]);
        plot([result_scheleton.p_min,result_scheleton.p_max],'--');
        legend({'P start','P aggregate','y s','P target','P min','P max'});
        hold off;
    elseif strcmp(variable_i,'V')
        plot(initial_vars.V);hold on;
        plot([result_scheleton.V,result_scheleton.Branch.y_v]);
        plot([result_scheleton.v_min,result_scheleton.v_max],'--');
        legend({'V start','V','y v','V min','V max'});
        hold off;
    elseif strcmp(variable_i,'E')
            plot(result_scheleton.e_workers);
            legend('e workers')
    elseif strcmp(variable_i,'D')
        plot(result_scheleton.lambdas);
        legend({'lambda s','lambda v'});
    end
    
    drawnow
end

end

function [r,c] = subplot_topology(n_subplot)
r = floor(sqrt(n_subplot));
c = ceil(n_subplot/r);
end