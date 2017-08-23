function f = boundedline_plot_one(x,y1,ci1,col1,alph,p_vals,xlab,ylab,xlim_in,ylim_in,fsize,fsize2,fig_path,fig_name)

% Function to plot one line with confidence intervals and p-values.
%
% Input:
% x         x-axis coordinates of data
% y1        y-axis coordinates of first dataset
% ci1       confidence intervals of first dataset (as absolute deviations from y1, [lower upper])
% col1      color for first dataset
% alph      alpha value for confidence interval transparency
% p_vals    p-values indicating significance of group differences (eg: from ranksum, or ttest2)
% xlab      x-axis label
% ylab      y-axis label
% fsize     font size for axis text labels
% fsize2    font size for axis numbers
% fig_path  path where figure is exported
% fig_name  figure name (for exporting as .png)
%
% Frantisek Vasa, % fv247@cam.ac.uk

f = figure; hold on;

ind1 = ~isnan(y1);

% plot line
[hls,hp] = boundedline(x(ind1),y1(ind1),ci1(ind1));
set(hls,'Color',col1,'LineWidth',2); set(hp,'FaceColor',col1,'FaceAlpha',alph); 

% axis limits
xlim(xlim_in); ylim(ylim_in);

% y-coordinate of p-value markers
yl = ylim_in; sig_mrk = yl(2)-0.05*(yl(2)-yl(1)); 

% p-values
s = scatter(find(p_vals<0.05),sig_mrk+zeros(sum(p_vals<0.05),1),20,'MarkerFaceColor',rgb('black'),'MarkerEdgeColor',rgb('black')); 

plot(x',zeros(length(x),1),'k--') % line at y=0

hold off;

ylim(yl)

% axis size and labels
set(gca,'FontSize',fsize2);
xlabel(xlab,'FontSize',fsize,'FontName','Arial');
ylabel(ylab,'FontSize',fsize,'FontName','Arial');

% adjust position
pos = get(gcf,'Position');
set(gcf,'Position', [pos(1) pos(2) 1.25*pos(3) 0.70*pos(4)]); set(gcf,'color','w');

% adjust limits
set(gca,'YTick',[ylim_in(1),ylim_in(1)+(ylim_in(2)-ylim_in(1))/2,ylim_in(2)])
if strcmp(ylab,'efficiency_{norm}')
    set(gca,'YTick',[ylim_in(1),ylim_in(1)+(1-ylim_in(1))/2,1])
end

movegui(f,'southwest')

% export figure
print([fig_path fig_name '.png'],'-dpng')

end