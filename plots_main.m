clear 
close all
load('one_case_results.mat')


figure;
% Power -----------------------------------------------------------------

subplot(3,3,1)
P0 = results.init_vals.P;
p_max = results.p_max;
p_min = results.p_min;
P = results.P;
plot([P0,P]); hold on; plot([p_min,p_max],'--');
xlim([0,47]);
xlabel('time [30 min]')
ylabel('Pn [-]')


subplot(3,3,4)
P0 = results.init_vals.q1.P;
p_max = results.q1.p_max;
p_min = results.q1.p_min;
P = results.q1.P;
plot([P0,P]); hold on; plot([p_min,p_max],'--');
xlim([0,47]);
xlabel('time [30 min]');
ylabel('Pn [-]')


subplot(3,3,7)
P0 = results.init_vals.q1.q1.P;
p_max = results.q1.q1.p_max;
p_min = results.q1.q1.p_min;
P = results.q1.q1.P;
plot([P0,P]); hold on; plot([p_min,p_max],'--');
xlim([0,47]);
xlabel('time [30 min]') 
ylabel('Pn [-]')


% voltage-----------------------------------------------------------------
subplot(3,3,2)
P0 = results.init_vals.V;
p_max = results.v_max;
p_min = results.v_min;
P = results.V;
plot([P0,P]/results.V_ref); hold on; plot([p_min,p_max]/results.V_ref,'--');
xlim([0,47]);
ylim([0.85,1.08])
xlabel('time [30 min]') 
ylabel('Voltage [p.u.]')

subplot(3,3,5)
P0 = results.init_vals.q1.V;
p_max = results.q1.v_max;
p_min = results.q1.v_min;
P = results.q1.V;
plot([P0,P]/results.q1.V_ref); hold on; plot([p_min,p_max]/results.q1.V_ref,'--');
xlim([0,47]);
ylim([0.85,1.08])
xlabel('time [30 min]') 
ylabel('Voltage [p.u.]')

subplot(3,3,8)
P0 = results.init_vals.q1.q1.V;
p_max = results.q1.q1.v_max;
p_min = results.q1.q1.v_min;
P = results.q1.q1.V;
plot([P0,P]/results.q1.q1.V_ref); hold on; plot([p_min,p_max]/results.q1.q1.V_ref,'--');
xlim([0,47]);
ylim([0.85,1.08])
xlabel('time [30 min]') 
ylabel('Voltage [p.u.]')

% Energy -----------------------------------------------------------------
subplot(3,3,3)
E0 = max(results.e_workers(:))/0.9;
E = results.e_workers/E0;
plot(E); xlim([0,47]);
xlabel('time [30 min]') 
ylim([0,1]);
ylabel('SOC [-]')

subplot(3,3,6)
E0 = max(results.q1.e_workers(:))/0.9;
E = results.q1.e_workers/E0;
plot(E); xlim([0,47]);
xlabel('time [30 min]') 
ylim([0,1]);
ylabel('SOC [-]')

subplot(3,3,9)
E0 = max(results.q1.q1.e_workers(:))/0.9;
E = results.q1.q1.e_workers/E0;
plot(E); xlim([0,47]);
xlabel('time [30 min]') 
ylim([0,1]);
ylabel('SOC [-]')

