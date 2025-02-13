function save_plot(plotname,txtsz)
set(findobj(gcf,'type','axes'),'FontSize',txtsz)
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 25 25]); 
exportgraphics(gcf,join([plotname,".png"],""),'Resolution',300) 
saveas(gcf,join([plotname,".png"],""));
ax = gcf;
ax.RendererMode = 'manual';  % use a set with older versions of Matlab
hgexport(ax, join([plotname,".eps"],""),...
    hgexport('factorystyle'), 'Format', 'eps');
savefig(gcf,plotname)
end
