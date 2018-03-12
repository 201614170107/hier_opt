function violPerc=get_violation_perc(P,viol)
[f,x] = ecdf(abs(sum(P,2)));

violPerc = 1 - interp1(x(2:end),f(2:end),viol,'pchip');
plot(x,f)
hold on
% hline(0.9)
% vline(viol)
hold off
drawnow;