function save_pic(h,name,painters)
if nargin<3
    painters = true;
end
set(h,'Units','inches');
h.set('position',[0,0,10,8]);
screenposition = get(gcf,'Position');
set(h,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
set(gca,'FontSize',8)
if painters
print(name, '-dpdf', '-painters')
else
print(name, '-dpdf')
end