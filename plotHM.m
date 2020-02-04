function ax = plotHM(ax, X, args)
%Plots an M xN matrix as a heatmap
% -------------------------------------------------------------------------
% ax = plotHM(ax, X, GraphLimits, xLabel, yLabel, cLabel)
% -------------------------------------------------------------------------
%     ax matlab.graphics.axis.Axes;
%     X double;
%     GraphLimits double = [floor(prctile(X(:),5)),ceil( prctile(X(:),95))];
%     XLabel string = "Time (h)"
%     YLabel string = "Single Cells"
%     CLabel string = "Nuclear NF\kappaB"
arguments 
    ax 
    X double;
    args.GraphLimits double = [floor(prctile(X(:),5)),ceil( prctile(X(:),95))];
    args.XLabel string = "Time (h)"
    args.YLabel string = "Single Cells"
    args.CLabel string = "Nuclear NF\kappaB"
end
%% Initialize variables
t=0:1/12:((size(X,2)-1)/12);
opt= maketicks(t,args.GraphLimits,0); 
%% Plot heatmap
mod_colormap = divergingmap(0:1/1023:1,[14 28 77]/255,[158 27 2]/255);
mod_colormap(1,:) = [0.1 0.1 0.1];

% Make plot, setting axes and other opt
XDOWN = [zeros(1,size(X,2));X(1:end-1,:)];
XUP = [X(2:end,:);zeros(1,size(X,2))];
imagesc(ax,opt.TimeBounds, [1 size(X,1)],nan(size(X)),opt.MeasurementBounds);
hold on
imagesc(ax, opt.TimeBounds, [1 size(X,1)],XUP,opt.MeasurementBounds);
imagesc(ax,opt.TimeBounds, [1 size(X,1)],XDOWN,opt.MeasurementBounds);
imagesc(ax,opt.TimeBounds, [1 size(X,1)],X,opt.MeasurementBounds); 
hold off
% set(h, 'Parent', ax);
colormap(mod_colormap)
ylabel(ax, args.YLabel ,'interpreter','tex','Rotation',90, 'FontWeight', "Bold");
xlabel(ax, args.XLabel, 'FontWeight', "Bold");
set(ax, 'YTickLabel',[], 'YColor','k');
if ~isempty(args.CLabel)
  yyaxis(ax,'right')
  ylabel(ax, args.CLabel ,'interpreter','tex','Rotation',90, 'Color','k', 'FontWeight', "Bold");
  cTick =linspace(0,max(args.GraphLimits),max([2,min([args.GraphLimits(2)+1,5])]));
 c= colorbar(ax,'YTick',cTick,'YTickLabel',string(cTick));
 c.AxisLocation = "out";
end
set(ax,'TickDir','in','LineWidth',1, 'box','on', 'YTick',[]);

end

