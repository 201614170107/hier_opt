function P=generate_load_and_pv(P_Load,P_PV, shuffle_interval, num_loads, num_pv, PV_power_intervals,day)
%shuffle load
ix=floor((0:length(P_Load)-1)/96)+1;
P_shuffled=P_Load;
for p=1:shuffle_interval:floor(ix(end)/(shuffle_interval))*shuffle_interval-shuffle_interval
    ix2 = p:min(p+shuffle_interval,ix(end));
    ix3 = ix2(randperm(length(ix2)));
    for q=1:length(ix2)
        try
        P_shuffled(find(ix==ix2(q)),:)=P_Load(find(ix==ix3(q)),:);
        catch
            disp('err')
        end
    end
end

P=sum(P_shuffled(:,randi(size(P_Load,2),num_loads,1)),2);
P=P-sum(P_PV(:,randi(size(P_PV,2),num_pv,1)),2)/num_pv*mean(P)*random('unif',PV_power_intervals(1),PV_power_intervals(2),1);
% plot(P)
% hold on
% plot((day-1)*96+1:day*96,P((day-1)*96+1:day*96),'linewidth',2);
P=P((day-1)*96+1:day*96);