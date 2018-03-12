classdef Branch < handle
    % Queen class, contains forecasting methods for voltage and power
    % forecasting and ground truth power and voltage profiles
    properties
        P_unc_perfect   % ground truth for the uncontrolled power from non-Workers
        V_perfect       % ground truth for the voltage at Queen (due to upper level and non-Workers)
        P_unc_forecast  % forecasted P_unc. This is updated at the beginning of the time slot.
        V_forecast      % forecasted V. This is updated at the beginning of the time slot.
        f_method        % Forecast method, in {'perfect','noise',..}
        H               % Length of the horizon
        P_ref           % reference power
        V_ref           % reference voltage
        y_s             % auxiliary variable for sum
        y_v             % auxiliary variable for generalized sum
        lambda_s        % dual variable sum constraints
        lambda_v        % dual variable generalized sum constraints
        gamma           % dual variable budget constraint
        r_s             % residual variable sum
        r_v             % residual variable generalized sum
        r_b             % residual of budget constraint
        
        h_f             % figure handle
    end
    methods
        function obj = Branch(H,V_ref,P_ref,f_method,P_unc_perfect,V_perfect)
            % Inputs:
            %        - workers: vector of Workers in the Queen. Needed for
            %        retrieving their
            
            if nargin<5
                obj.P_unc_perfect = zeros(H,1); % Base case: no non-Workers in the Queen
            else
                obj.P_unc_perfect = P_unc_perfect;
            end
            
            if nargin<6
                obj.V_perfect = V_ref*ones(H,1); % We forecast a value equal to V_ref
            else
                obj.V_perfect = V_perfect;
            end
            
            obj.f_method = f_method;
            obj.V_ref = V_ref;
            obj.P_ref = P_ref;
            obj.H = H;
            obj.h_f = figure('visible','off');
        end
        
        function [Phat,Vhat] = do_forecasts(obj,n)
            % perform forecasts of the power profile at main, at the
            % current timestep n,using method specified in f_method
            
            if strcmp(obj.f_method, 'perfect')
                % Perfect forecasts
                Phat = obj.P_unc_perfect(n:n+obj.H-1);
                Vhat = obj.V_perfect(n:n+obj.H-1);
            elseif strcmp(obj.f_method, 'noise')
                % Forecast is the perfect forecast, corrupted with normal
                % noise with a standard deviation of 5% of the perfect
                % forecast
                Phat = obj.P_unc_perfect(n:n+obj.H-1)+randn(obj.H,1)*(std(obj.P_unc_perfect(n:n+obj.H-1)))*0.05;
                Vhat = obj.V_perfect(n:n+obj.H-1)+randn(obj.H,1)*(std(obj.V_perfect(n:n+obj.H-1)))*0.05;
            end
            obj.P_unc_forecast = Phat;
            obj.V_forecast = Vhat;
        end
        
    end
end
