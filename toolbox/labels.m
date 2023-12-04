% ============================================================
% Function: labels
% Author: Jack Fox
% Date: 11/01/23
%     Puts axis labels on plots
% ============================================================
% Inputs:
%   ax         axis handle (i.e. gca)
%   axStrArr   array of string handles: {"XAxisLabel", "YAxisLabel",
%                                                      "ZAxisLabel"}
%   titleStr   title string
% ============================================================
function labels(ax, axStrArr, titleStr)
    l(1) = xlabel(ax, axStrArr{1});
    l(2) = ylabel(ax, axStrArr{2});
    l(3) = title(ax, titleStr);
    if length(axStrArr) == 3
        l(4) = zlabel(ax, axStrArr{3});
    end
    set(l,'interpreter','latex');
end