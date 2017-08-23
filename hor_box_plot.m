function hor_box_plot(data1,data2,xlims,xt,col1,col2,xlab,lab1,lab2,fsize,fsize2)

% Function to plot a horizontal box-whisker plot.
%
% Input:
% data1     vector of data for first boxplot (if data are of unequal size, this is the longer vector)
% data2     vector of data for second boxplot 
% xlims     x-axis limits (as [lower upper])
% xt        x-tick positions (on x-axis)
% col1      color corresponding to data1
% col2      color corresponding to data2
% xlab      x-axis label
% lab1      label for first group
% labe2     label for second group
% fsize     font size for axis text labels
% fsize2    font size for axis numbers
%
% Frantisek Vasa, % fv247@cam.ac.uk

% in case length(data1) > length(data2)
nd = length(data1)-length(data2); % difference in data length (need to pad data with nans for boxplot)

figure;
b = boxplot([[data2 nan(1,nd)]; data1]','labels',{lab1,lab2},'symbol','k.','widths',0.5,'Orientation','horizontal'); % blank labels easiest way to "remove" them

% change colors of boxplot patches
box_col = [col1; col2];
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    set(h(j),'Color','w')
    patch(get(h(j),'XData'),get(h(j),'YData'),box_col(j,:),'FaceAlpha',0.7);
end

% set linewidth on boxplot
set(b,{'linew'},{1.3});

% change color of median lines
m = findobj(gca,'Tag','Median');

set(gca,'FontSize',fsize2); 

hold on;
plot(get(m(1),'XData')',get(m(1),'YData')','LineWidth',3,'Color',rgb('black'));
plot(get(m(2),'XData')',get(m(2),'YData')','LineWidth',3,'Color',rgb('black'));
hold off;

xlabel(xlab,'FontSize',fsize,'FontName','Arial');

set(gca,'XLim',xlims,'XTick',xt);

box off
set(gcf,'color','w')%,'YTick',[]); 

% adjust size
pos = get(gcf,'Position');
set(gcf,'Position', [pos(1) pos(2) pos(3) 0.40*pos(4)]); set(gcf,'color','w');

end