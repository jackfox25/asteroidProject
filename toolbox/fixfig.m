function fixfig(f, varargin)

    fontname = 'Times New Roman'; % leave alone.
    fontsize = 16;                % scale up or down as necessary.
    legloc = 'northeast';         % 'northwest', 'southeast', etc. legend can also be dragged after plot is created.
    linewidth = 2.5;              % scale up or down as necessary.
    marker = [];                  % replace brackets with '.', 'o', 'x' etc. to change the marker for each point. brackets leaves it alone.
    markersize = [];              % scale up or down as necessary.

    if ~isempty(varargin)
        linewidth = [];
    end    

    for i=1:length(f.Children)
        if isa(f.Children(i),'matlab.graphics.illustration.Legend')
            leg = f.Children(i);
            leg.FontName = fontname;
            leg.Location = legloc;
        elseif isa(f.Children(i),'matlab.graphics.axis.Axes')
            ax = f.Children(i);
            ax.XGrid = 'on'; ax.YGrid = 'on';
            ax.FontSize = fontsize; ax.FontName = fontname;
            elems = ax.Children;
            for j=1:length(elems)
                if isa(elems(j),'matlab.graphics.chart.primitive.Line')
                    l = elems(j);
                    if ~isempty(linewidth)
                        l.LineWidth = linewidth;
                    end
                    if ~isempty(marker)
                        l.Marker = marker;
                    end
                    if ~isempty(markersize)
                        l.MarkerSize = markersize;
                    end
                end
            end
        end
    end
    
end