function hist_box_plot(data1,data2,ylims,yincr,col1,col2,xlab,ylab,fsize,fsize2)

% Function to plot two box-whisker plots, along with histograms on either
% side.

% Input:
% data1     vector of data for first boxplot and histogram (if data are of unequal size, this is the longer vector)
% data2     vector of data for second boxplot and histogram
% ylims     y-axis limits (as [lower upper])
% yincr     increment for y-axis labels
% col1      color corresponding to data1
% col2      color corresponding to data2
% xlab      x-axis label
% ylab      y-axis label
% fsize     font size for axis text labels
% fsize2    font size for axis numbers

% length(data1) > length(data2)
nd = length(data1)-length(data2); % difference in data length (need to pad data with nans for boxplot)

if nd<0; error('length(data1) < length(data2)'); end
if size(data1,1) ~= 1; data1 = data1'; end
if size(data2,1) ~= 1; data2 = data2'; end

figure;

bar_width = 0.9;                % factor to scale bar width

[p,h] = ranksum(data1,data2);   % p-value for group difference

% first histogram
s1 = subplot(1,3,1);
n1 = hist(data1,ylims(1):yincr:ylims(2)); b1 = barh(ylims(1):yincr:ylims(2),n1/sum(n1),bar_width,'FaceColor',col1,'EdgeColor','k','FaceAlpha',0.7); 
ylim([ylims(1)-yincr ylims(2)+yincr]);

set(gca,'FontSize',fsize2);

% yaxis settings
ylabel(ylab,'FontSize',fsize,'FontName','Arial')
ylabh = get(gca,'YLabel');
ylabpos = get(ylabh,'Position');
set(gca,'yaxislocation','right','xdir','reverse')
x1 = xlim; y1 = ylim;
box off

% xaxis settings
xlabel(xlab,'FontSize',fsize,'FontName','Arial');

% second histogram
s2 = subplot(1,3,3);
n2 = hist(data2,ylims(1):yincr:ylims(2)); b2 = barh(ylims(1):yincr:ylims(2),n2/sum(n2),bar_width,'FaceColor',col2,'EdgeColor','k','FaceAlpha',0.7); 
ylim([ylims(1)-yincr ylims(2)+yincr]);
x2 = xlim; y2 = ylim;
box off
set(gca,'FontSize',fsize2);

% axis limits
set([s1,s2],'XLim',[min(x1(1),x2(1)) max(x1(2),x2(2))],'YLim',[min(y1(1),y2(1)) max(y1(2),y2(2))]);
set([s1,s2],'YTick',[min(y1(1),y2(1)):2*yincr:max(y1(2),y2(2))]); 

% switch ylabel to LHS
ylabpos(1) = abs(ylabpos(1))+1.2*max(x1(2),x2(2));
set(ylabh,'Position',ylabpos)

s3 = subplot(1,3,2);
b = boxplot([data1; [data2 nan(1,nd)]]','labels',{'',''},'symbol','k.','widths',0.5); % blank labels easiest way to "remove" them
title(['p = ' num2str(roundn(p,7))],'FontSize',fsize,'FontName','Arial','FontWeight','normal')

% change colors of boxplot patches
box_col = [col2; col1];
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    set(h(j),'Color','w')
    patch(get(h(j),'XData'),get(h(j),'YData'),box_col(j,:),'FaceAlpha',0.7);
end

% set linewidth on boxplot
set(b,{'linew'},{1.3});

% change color of median lines
m = findobj(gca,'Tag','Median');

hold on;
plot(get(m(1),'XData')',get(m(1),'YData')','LineWidth',3,'Color',rgb('black'));
plot(get(m(2),'XData')',get(m(2),'YData')','LineWidth',3,'Color',rgb('black'));
hold off;

set(s3,'YLim',[min(y1(1),y2(1)) max(y1(2),y2(2))]);

% adjust position of boxplot
p1 = get(s1,'Position'); p3 = get(s3,'Position');
set(s3,'Position', [p3(1) p1(2)-0.01 p3(3) p1(4)+0.01],'XColor','w','YColor','w','XTick',[],'YTick',[]) % C/YColor is a hack to keep title; otherwise: 'Visible','off'
% 0.01 increment is additional adjustment, as boxplot is not aligned in the same way as histograms

set(gcf,'color','w'); 

end