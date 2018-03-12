classdef battery < handle
    % Agent rapresenting an electric battery
    properties
        eta_in                               % charge efficiency
        eta_out                              % discharge efficiency
        ts                                   % simulation sampling time, [seconds]
        tau_sd                               % self discharge of the battery [1/seconds]
        lifetime                             % battery lifetime [years]
        DOD                                  % Depth of Discharge [-]
        Nc                                   % Number of cycles befor end of life [-]
        E_0                                  % Nominal capacity kWh
        c                                    % Maximum relative charge/discharge [kW/kWh]
        E                                    % Battery state [kWh]
        objfun                               % Objective function of the agent
        Pm_perfect                           % Main power profile, ground truth [kWh]
        f_method                             % Forecast method, in {'perfect','noise',..}
        x_u                                  % Fixed upper state constraints
        x_l                                  % Fixed lower state constraints
        u_u                                  % Fixed upper operational constraints
        u_l                                  % Fixed lower operational constraints
        Ac
        Bc
        Ec
        Ad
        Bd
        Ed
        C
        h                                   % Lenght of the horizon
        k_t
        k_q
        H
        Hu
        H0
        Hd
        D
        t
        Q
        S
        selfish                             % If true, the agent solve selfish problem even in cooperative setting
        U                                   % Optimized agent's action; this is updated at each iteration
        Pm_opt                              % Optimized agent's net power profile, based on forecasts; this is updated at each iteration
        e_opt                               % Optimized energy profile; this is updated at each iteration
        Pm_forecast                         % Forecasted uncontrolled power of workers. This is updated at the beginning of the time slot.
        noneconomic_cost
        cost
        slack % if true battery is a slack battery
    end
    
    
    methods
        
        function obj = battery(eta_in,eta_out,ts,tau_sd,lifetime,DOD,Nc,E_0,c,h,t,Pm_perfect,f_method,selfish,slack)
            if nargin < 15
                slack = false;
            end
            obj.slack = slack;    
            % Constructor
            obj.eta_in = eta_in;                          % charge efficiency
            obj.eta_out = eta_out;                        % discharge efficiency
            obj.ts = ts;                                  % simulation sampling time, [seconds]
            obj.tau_sd = tau_sd;                          % self discharge of the battery [1/seconds]
            obj.lifetime = lifetime;                      % battery lifetime [years]
            obj.DOD = DOD;                                % Depth of Discharge [-]
            obj.Nc = Nc;                                  % Number of cycles befor end of life [-]
            obj.E_0 = E_0;                                % Nominal capacity kWh
            obj.E = E_0*0.3+E_0*0.6*rand(1);                                % Battery state, initialized as half of the nominal charge [kWh]
            obj.c = c;                                    % Maximum relative charge/discharge [kW/kWh]
            obj.Ac = -tau_sd;                             % Continuous state matrix
            obj.Bc = [eta_in/3600 -1/(eta_out*3600)];     % Continuous input matrix considering that the state is in kWh
            % and the inputs in kW
            obj.Ec = 1;                                   % Continuous process noise matrix
            obj.C = 1;                                    % Output matrix
            obj.h = h;                                    % Lenght of the horizon
            obj.t = t;                                    % Proximal coefficient
            obj.Pm_perfect = Pm_perfect;                  % Main power profile, ground truth [kWh]
            obj.f_method = f_method;                      % Forecast method, in {'perfect','noise',..}
            obj.selfish = selfish;                        % If true, the agent solve selfish problem even in cooperative setting
            
            % Aging effect parameters
            lifetime_IS = lifetime*3600*24*365;           % convert lifetime in international system: from years to seconds
            obj.k_t = -0.2/lifetime_IS;                   % calendar aging coefficient [1/s]
            obj.k_q = -0.2/(2*Nc*E_0*DOD);                % coefficient of degradation
            
            % build discrete matrices under constant input assumption
            [obj.Ad,obj.Bd,obj.Ed] = c2d(obj.Ac,obj.Bc,obj.Ec,obj.ts);
            
            % build impulse response matrices for the length of the horizon
            [obj.Hu,obj.H0,obj.Hd] = get_batched(obj.h,obj.Ad,obj.Bd,obj.C,obj.Ed);
            
            % build difference matrix for total power from/to the battery
            D = zeros(h,2*h);
            for i=1:obj.h
                D(i,1+(i-1)*2:2*i) = [1,-1];
            end
            obj.D = D;
            
            % build quadratic cost matrix
            obj.H = [D,zeros(h,h)];
            obj.Q = 0.5*(obj.H'*obj.H);
            
            % Build summation cost matrix
            obj.S = [zeros(1,h*2),ones(1,h)];
            
            
            % build fixed constraints
            obj.x_u = obj.E_0*0.9;
            obj.x_l = obj.E_0*0.1;
            obj.u_l = [0,0];
            obj.u_u = ones(1,2)*obj.E_0*obj.c;            
        end
        
        function [Ain,bin] = power_and_soc_constraints(obj,x0,x_u,x_l,u_u,u_l)
            % Build power and SOC constraints. The decision variable is
            %[Pin,1 Pout,1 Pin,2 Pout,2...Pin,t Pout,t]'.
            hor = obj.h;
            Eye = eye(size(u_l(:),1)*hor,size(u_l(:),1)*hor);
            Ain = [Eye;-Eye;obj.Hu;-obj.Hu];
            bin = [repmat(u_u',hor,1);
                repmat(-u_l',hor,1);
                -obj.H0*x0+repmat(x_u,hor,1);
                obj.H0*x0-repmat(x_l,hor,1)];
        end
        
        function [Ain_aug,bin_aug] = augment_Ain_bin(obj,Ain,bin,Pm,pb,ps)
            % Augment the inequality matrices in the case the objective
            % function includes an economic cost with an ansymmetric
            % buying/selling price: p*(Pm+Pin-Pout) where p=pbuy if
            % Pm+Pin-Pout>0, p=psell if Pm+Pin-Pout<0
            % Augment the decision variable vector as
            % [y,1 y,2...y,t Pin,1 Pout,1 Pin,2 Pout,2...Pin,t Pout,t]'
            
            %Duplicate the buying and selling price
            pbb = reshape(repmat(pb',2,1),1,[]);
            pss = reshape(repmat(ps',2,1),1,[]);
            Ain = [Ain,zeros(size(Ain,1),obj.h)];
            Add = [obj.D.*pbb,-eye(size(Pm,1));obj.D.*pss,-eye(size(Pm,1))];
            Ain_aug = [Add;Ain];
            bin_aug = [-pb.*Pm;-ps.*Pm;bin];
        end
        
        function [Ain,bin] = build_constraints(obj,Pm,pb,ps,x0)
            [Ain,bin] = power_and_soc_constraints(obj,x0,obj.x_u,obj.x_l,obj.u_u,obj.u_l);
            [Ain,bin] = augment_Ain_bin(obj,Ain,bin,Pm,pb,ps);
        end
        
        function [u,Pm_opt,e_opt] = solve_selfish(obj,Pm,pb,ps)
            [Ain,bin] = build_constraints(obj,Pm,pb,ps,obj.E);
            opts = optimoptions('linprog','algorithm','dual-simplex');
            x = linprog(obj.S,Ain,bin,[],[],[],[],opts);
            u = reshape(x(1:2*obj.h),2,[])';
            Pm_opt = Pm +u(:,1) -u(:,2); % Pm +Pin -Pout
            e_opt = obj.H0*obj.E + obj.Hu*x(1:2*obj.h); % impulse response state equation
        end
        
        function [u,Pm_opt,e_opt,x] = solve_cooperative(obj,Pm,Pref,pb,ps,x_ws)
            % x_ws: warm start for the solution
            if nargin <11
                x_ws = [];
            end
            [Ain,bin] = build_constraints(obj,Pm,pb,ps,obj.E);
            Qm = obj.Q/obj.t/2;
            add_regularization = false;
            if add_regularization
                Qm = Qm + 1e-10*eye(size(Qm))+ 1e1*blkdiag( eye(obj.h*2),zeros(obj.h)) ;
            end
            qoptions = optimoptions('quadprog','display','none');
            L = obj.S-(Pref-Pm)'*obj.H/obj.t;
            %             L = obj.S-(Pref)'*obj.H/obj.t;
            x = quadprog(Qm,L,Ain,bin,[],[],[],[],x_ws,qoptions);
            u = reshape(x(1:2*obj.h),2,[])';
            Pm_opt = Pm +u(:,1) -u(:,2); % Pm +Pin -Pout
            e_opt = obj.H0*obj.E + obj.Hu*x(1:2*obj.h); % impulse response state equation
            obj.U = u(:,1) -u(:,2); % store the optimal solution at current iteration
            obj.Pm_opt = Pm_opt;% store optimal Pm at current iteration
            obj.e_opt = e_opt;% store optimal e at current iteration
        end
        
        function [u,Pm_opt,e_opt,f,g] = solve_cooperative_yalmip(obj,Pm,Pref,pb,ps,iter0,K)
            if nargin<7
                K = 1;
            end
            
            yalmip clear
            
            % check exhistence of gurobi
            if exist('gurobi','file') == 3
                sets = sdpsettings('solver','gurobi','verbose',0);
            else
                sets = sdpsettings('solver','quadprog','verbose',0);
            end
            
            x = sdpvar(obj.h*2,1);
            cost =sdpvar(obj.h,1);
            % p in and p out indexes in the solution vector
            pindex = 1:2:(2*obj.h-1);
            poutdex = 2:2:2*obj.h;
            %
            constraints = [x(pindex)<=repmat(obj.u_u(1),obj.h,1)];
            constraints = [constraints, x(poutdex)<=repmat(obj.u_u(2),obj.h,1)];
            constraints = [constraints, x(pindex)>=repmat(obj.u_l(1),obj.h,1)];
            constraints = [constraints, x(poutdex)>=repmat(obj.u_l(2),obj.h,1)];
            
            constraints = [constraints, obj.H0*obj.E + obj.Hu*x <= repmat(obj.x_u,obj.h,1)];
            constraints = [constraints, obj.H0*obj.E + obj.Hu*x >= repmat(obj.x_l,obj.h,1)];
            constraints = [constraints, obj.H0(end,:)*obj.E + obj.Hu(end,:)*x >= obj.E_0/2];
            constraints = [constraints, cost(1:obj.h) >= (Pm(1:obj.h) +x(pindex)-x(poutdex)).*pb];
            constraints = [constraints, cost(1:obj.h) >= (Pm(1:obj.h)+x(pindex)-x(poutdex)).*ps];
            
            
            % objective function
            of = sum(cost)+100*(Pm<Pref)'*x(poutdex);
            if ~iter0
                of = of +K*(Pm-Pref+x(pindex)-x(poutdex))'*(Pm-Pref+x(pindex)-x(poutdex))/obj.t/2;
            end
            
            sol = optimize(constraints,of,sets);
            x=value(x);
            f = value((Pm-Pref+x(pindex)-x(poutdex))'*(Pm-Pref+x(pindex)-x(poutdex))/obj.t/2);
            g = value(sum(cost));
            u = reshape(x(1:2*obj.h),2,[])';
            Pm_opt = Pm +u(:,1) -u(:,2); % Pm +Pin -Pout
            e_opt = obj.H0*obj.E + obj.Hu*x(1:2*obj.h); % impulse response state equation
            % save current solution in the agent
            obj.Pm_opt = Pm_opt; % save optimized net power profile
            obj.e_opt = e_opt; % save optimized energy profile
            obj.U = u; % save optimized battery profiles
            
        end
        
        function [u,Pm_opt,e_opt,f,g] = solve_cooperative_admm_yalmip(obj,Pm,pb,ps,R,K_coeff,iter0,multipliers_ref,multipliers_coeffs,alpha)
            if nargin<7
                iter0=false;
            end
            if nargin<8
                multipliers_ref=[];
            end
            
            if nargin<9
                multipliers_coeffs=[];
            end
            if nargin<10
                alpha=1;
            end
            
            yalmip clear

            % check exhistence of gurobi
            if exist('gurobi','file') == 3
                sets = sdpsettings('solver','gurobi','verbose',0,'gurobi.TimeLimit',1);
            else
                sets = sdpsettings('solver','quadprog','verbose',0);
            end
            
            if ~obj.slack
                
                x = sdpvar(obj.h*2,1);
                cost =sdpvar(obj.h,1);
                % p in and p out indexes in the solution vector
                pindex = 1:2:(2*obj.h-1);
                poutdex = 2:2:2*obj.h;
                %
                constraints = [x(pindex)<=repmat(obj.u_u(1),obj.h,1)];
                constraints = [constraints, x(poutdex)<=repmat(obj.u_u(2),obj.h,1)];
                constraints = [constraints, x(pindex)>=repmat(obj.u_l(1),obj.h,1)];
                constraints = [constraints, x(poutdex)>=repmat(obj.u_l(2),obj.h,1)];
                
                constraints = [constraints, obj.H0*obj.E + obj.Hu*x <= repmat(obj.x_u,obj.h+1,1)];
                constraints = [constraints, obj.H0*obj.E + obj.Hu*x >= repmat(obj.x_l,obj.h+1,1)];
%                 constraints = [constraints,x(end-1:end)== zeros(2,1)];
%                             constraints = [constraints, obj.H0(end,:)*obj.E + obj.Hu(end,:)*x >= obj.E];
                constraints = [constraints, cost(1:obj.h) >= (Pm(1:obj.h) +x(pindex)-x(poutdex)).*pb];
                constraints = [constraints, cost(1:obj.h) >= (Pm(1:obj.h)+x(pindex)-x(poutdex)).*ps];
                
                
                % objective function
                
                of =  sum(cost);
                if ~iter0
                    for i=1:size(R,2)
                        %             of =  of + (W(:,i).*(Pm+x(pindex)-x(poutdex))-R(:,i))'*(W(:,i).*(Pm+x(pindex)-x(poutdex))-R(:,i))/obj.t/2;
                        of = of +100*double(Pm<R(:,i))'*x(poutdex);
                        of =  of + alpha*((Pm+x(pindex)-x(poutdex))*K_coeff(i)-R(:,i))'*((Pm+x(pindex)-x(poutdex))*K_coeff(i)-R(:,i))/obj.t/2;
                    end
                end
                if ~isempty(multipliers_ref)
                    of = of + (multipliers_ref-multipliers_coeffs'*(Pm+x(pindex)-x(poutdex)))*(multipliers_ref-multipliers_coeffs'*(Pm+x(pindex)-x(poutdex)))/obj.t/2;
                end
                sol = optimize(constraints,of,sets);
                x = value(x);
                g = value(sum(cost));
                f = value(of)-g;
                u = reshape(x(1:2*obj.h),2,[])';
                Pm_opt = Pm +u(:,1) -u(:,2); % Pm +Pin -Pout
                e_opt = obj.H0*obj.E + obj.Hu*x(1:2*obj.h); % impulse response state equation
                % save current solution in the agent
                obj.Pm_opt = Pm_opt; % save optimized net power profile
                obj.e_opt = e_opt(1:obj.h); % save optimized energy profile
                obj.U = u; % save optimized battery profiles
            else
                x = sdpvar(obj.h,1);
                of =x'*x*1e-3;
                for i=1:size(R,2)
                    of =  of + ((x)*K_coeff(i)-R(:,i))'*((x)*K_coeff(i)-R(:,i))/obj.t/2;
                end
                sol = optimize([],of,sets);
                u = value(x);
                Pm_opt = value(x);
                e_opt = value(x);
                f = 0;
                g = 0;
                
                obj.Pm_opt = Pm_opt; % save optimized net power profile
                obj.e_opt = e_opt; % save optimized energy profile
                obj.U = u; % save optimized battery profiles
            end
            
            
        end
        
        function Phat = do_forecasts(obj,n)
            % perform forecasts of the power profile at main, at the
            % current timestep n,using method specified in f_method
            
            if strcmp(obj.f_method, 'perfect')
                % Perfect forecasts
                Phat = obj.Pm_perfect(n:n+obj.h-1);
            elseif strcmp(obj.f_method, 'noise')
                % Forecast is the perfect forecast, corrupted with normal
                % noise with a standard deviation of 5% of the perfect
                % forecast
                Phat = obj.Pm_perfect(n:n+obj.h-1)+randn(obj.h,1)*(std(obj.Pm_perfect(n:n+obj.h-1)))*0.05;
            end
            obj.Pm_forecast = Phat;
        end
        function test_selfish(obj,n)
            
            N = n*obj.h;
            t = [1:N]';
            pb = 0.2*sin(n*2*pi*t/N)+0.25; % buying price profile
            ps = 0.05*sin(n*2*pi*t/N)+0.05; % selling price profile
            Pm =  5*sin(n*1.4*pi*t/N)+2.5+0.2*randn(N,1); % main power profile
            u_opt_t = zeros(N,2);
            r_opt_t = u_opt_t;
            e_opt_t = u_opt_t;
            figure
            for i=1:N-obj.h
                Pm_h = Pm(i:i + obj.h-1);
                pb_h = pb(i:i + obj.h-1);
                ps_h = ps(i:i + obj.h-1);
                [u_opt_h,Pm_opt_h,e_opt_h] = solve_selfish(obj,Pm_h,pb_h,ps_h);
                u_opt_t(i,:) = u_opt_h(1,:);
                r_opt_t(i) = Pm_opt_h(1);
                e_opt_t(i) = e_opt_h(2);
                obj.E = e_opt_t(i);
                fprintf('%i\n',i)
                subplot 311
                plot([e_opt_h,u_opt_h])
                legend({'E','Pin','Pout'})
                subplot 312
                plot([Pm_h,Pm_opt_h])
                legend({'Pm','Pm opt'})
                subplot 313
                plot([pb_h,ps_h])
                legend({'pb','ps'})
                drawnow
            end
            
        end
        
        function test_cooperative(obj,n,type)
            iter0 = false;
            N = n*obj.h;
            t = [1:N]';
            %             pb = 0.2*sin(n*2*pi*t/N)+0.25; % buying price profile
            %             ps = 0.05*sin(n*2*pi*t/N)+0.05; % selling price profile
            pb = 0.2*ones(N,1);
            ps = 0.07*ones(N,1);
            
            Pm =  5*sin(n*1.4*pi*t/N)+2.5+0.2*randn(N,1); % main power profile
            Pref =  0*(5*sin(n*1.4*pi*t/N+pi/3)+2.5+0.2*randn(N,1)); % main power profile
            u_opt_t = zeros(N,2);
            r_opt_t = u_opt_t;
            e_opt_t = u_opt_t;
            figure
            for i=1:N
                Pm_h = Pm(i:i + obj.h-1);
                pb_h = pb(i:i + obj.h-1);
                ps_h = ps(i:i + obj.h-1);
                Pref_h = Pref(i:i + obj.h-1);
                if strcmp(type,'qp')
                    [u_opt_h,Pm_opt_h,e_opt_h] = solve_cooperative(obj,Pm_h,Pref_h,pb_h,ps_h,iter0);
                elseif strcmp(type,'miqp')
                    [u_opt_h,Pm_opt_h,e_opt_h] = solve_cooperative_MIQP(obj,Pm_h,Pref_h,pb_h,ps_h,iter0);
                elseif strcmp(type,'yalmip')
                    [u_opt_h,Pm_opt_h,e_opt_h] = solve_cooperative_yalmip(obj,Pm_h,Pref_h,pb_h,ps_h,iter0);
                end
                %                 [u_opt_h,Pm_opt_h,e_opt_h] = solve_selfish(obj,Pm_h,pb_h,ps_h,x0,x_u,x_l,u_u,u_l);
                u_opt_t(i,:) = u_opt_h(1,:);
                r_opt_t(i) = Pm_opt_h(1);
                e_opt_t(i) = e_opt_h(2);
                obj.E = e_opt_t(i);
                fprintf('%i\n',i)
                subplot 311
                plot([e_opt_h,u_opt_h])
                legend({'E','Pin','Pout'})
                subplot 312
                plot([Pm_h,Pm_opt_h,Pref_h])
                legend({'Pm','Pm opt','Pref'})
                subplot 313
                plot([pb_h,ps_h])
                legend({'pb','ps'})
                drawnow
                
            end
        end
    end
end
