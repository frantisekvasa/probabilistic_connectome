function f = boundedline_plot_2p_drop(x,y1,ci1,y2,ci2,drop1,drop2,col1,col2,alph,r3to1,p_vals,p_vals2,xlab,ylab,ylab2,leg1,leg2,xlim_in,ylim_in,ylim_in2,fsize,fsize2,fig_path,fig_name)

% Function to plot two lines with confidence intervals, two sets of p-values and an additional two lines (indicating drop-off in group size).
%
% Input:
% x         x-axis coordinates of data
% y1        y-axis coordinates of first dataset
% ci1       confidence intervals of first dataset (as absolute deviations from y1, [lower upper])
% y2        y-axis coordinates of second dataset
% ci2       confidence intervals of second dataset (as absolute deviations from y2, [lower upper])
% drop1     "drop-off" curve one (size of y1 as a function of x)
% drop2     "drop-off" curve two (size of y2 as a function of x)
% col1      color for first dataset
% col2      color for second dataset
% alph      alpha value for confidence interval transparency
% r3to1     horizontal lines of group size ratios (for specific applications only) - if unused set as NaN
% p_vals    p-values indicating significance of group differences (eg: from ranksum, or ttest2)
% p_vals2   a second set of p-values indicating significance of group differences (eg: from permutation testing)
% xlab      x-axis label
% ylab      y-axis label
% ylab2     second y-axis label (for right-hand y-axis)
% leg1      legend label for first dataset
% leg2      legend label for second dataset
% fsize     font size for axis text labels
% fsize2    font size for axis numbers
% fig_path  path where figure is exported
% fig_name  figure name (for exporting as .png)
%
% Frantisek Vasa, % fv247@cam.ac.uk

f = figure; hold on;

ind1 = ~isnan(y1);

% first line
[hls,hp] = boundedline(x(ind1),y1(ind1),ci1(ind1));
set(hls,'Color',col1,'LineWidth',2); set(hp,'FaceColor',col1,'FaceAlpha',alph); 

ind2 = ~isnan(y2);

% second line
[hln,hp] = boundedline(x(ind2),y2(ind2),ci2(ind2));
set(hln,'Color',col2,'LineWidth',2); set(hp,'FaceColor',col2,'FaceAlpha',alph); 

% horizontal line for group size ratios (1:3, 3:1)
if ~isnan(r3to1)
    yl = ylim_in;
    p = plot(r3to1(1)*ones(10,1),linspace(yl(1),yl(2),10),'--','Color',rgb('black'));
    p = plot(r3to1(2)*ones(10,1),linspace(yl(1),yl(2),10),'--','Color',rgb('black'));
end

% axis limits
xlim(xlim_in); ylim(ylim_in);

% y-coordinate of p-value markers
yl = ylim_in; sig_mrk = yl(2)-0.05*(yl(2)-yl(1)); 

% p-values
s = scatter(find(p_vals2<0.05),sig_mrk+zeros(sum(p_vals2<0.05),1),70,'MarkerFaceColor',rgb('white'),'MarkerEdgeColor',rgb('black')); 
s = scatter(find(p_vals<0.05),sig_mrk+zeros(sum(p_vals<0.05),1),10,'MarkerFaceColor',rgb('black'),'MarkerEdgeColor',rgb('black')); 

% legend
if ~isnan(r3to1)
    legend([hls,hln,s,p],leg1,leg2,'P < 0.05','3:1 / 1:3 ratios','Location','SouthEast','FontSize',fsize)
else
    legend([hls,hln,s],leg1,leg2,'P < 0.05','Location','SouthEast','FontSize',fsize)
end

hold off;

ylim(yl)

% axis size and labels
set(gca,'FontSize',fsize2);
xlabel(xlab,'FontSize',fsize,'FontName','Arial');
ylabel(ylab,'FontSize',fsize,'FontName','Arial');

pos = get(gcf,'Position');
set(gcf,'Position', [pos(1) pos(2) 1.25*pos(3) 0.70*pos(4)]); set(gcf,'color','w');

set(gca,'YTick',[ylim_in(1),ylim_in(1)+(ylim_in(2)-ylim_in(1))/2,ylim_in(2)])
if strcmp(ylab,'efficiency_{norm}')
    set(gca,'YTick',[ylim_in(1),ylim_in(1)+(1-ylim_in(1))/2,1])
end

% drop off curves
yyaxis right
hold on
plot(x,drop1,'--','Color',col1,'LineWidth',1);
plot(x,drop2,'--','Color',col2,'LineWidth',1);
ylb2 = ylabel(ylab2,'rot',-90,'FontSize',fsize,'FontName','Arial')
ylab2pos = get(ylb2,'Position');
ylab2pos(1) = 1.065*ylab2pos(1);
set(ylb2,'Position',ylab2pos)
hold off

% adjust position
gca; ylim(ylim_in2)
set(gca,'YColor',rgb('black'),'YTick',[ylim_in2(1),ylim_in2(1)+(ylim_in2(2)-ylim_in2(1))/2,ylim_in2(2)])

movegui(f,'southwest')

% export figure
print([fig_path fig_name '.png'],'-dpng')

end