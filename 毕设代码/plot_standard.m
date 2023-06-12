function plot_standard(h,MarkerSize,YLabelFontSize,FontSize,LineWidth,LegendFontSize,TitleFontSize)
%% Default value for figures.
% MarkerSize=9;  LineWidth =2; LineMethod =1; PlotMethod =1; FontSize=22; FontSize= 24 ; LineWidth = 3;
% TitleFontSize =20; LegendFontSize =24; axis_ratio=1.5; YLabelFontSize =24;
% myboldify: make lines and text bold;  boldifies the current figure; applies to the figure with the handle h
if nargin < 1
    h = gcf;
    MarkerSize=9; YLabelFontSize =24; FontSize= 24 ;
    LineWidth = 2;TitleFontSize =24; LegendFontSize =36; axis_ratio=1.5; %myboldify(h,FontSize,LineWidth,LegendFontSize,TitleFontSize);
end
ha = get(h, 'Children'); % the handle of each axis
for i = 1:length(ha)
    if strcmp(get(ha(i),'Type'), 'axes') % axis format
        set(ha(i), 'FontSize', LegendFontSize);      % tick mark and frame format
        set(ha(i), 'LineWidth', LineWidth);
        set(get(ha(i),'XLabel'), 'FontSize', YLabelFontSize);
        %set(get(ha(i),'XLabel'), 'VerticalAlignment', 'top');
        set(get(ha(i),'YLabel'), 'FontSize', YLabelFontSize);
        %set(get(ha(i),'YLabel'), 'VerticalAlignment', 'baseline');
        set(get(ha(i),'ZLabel'), 'FontSize', FontSize);
        %set(get(ha(i),'ZLabel'), 'VerticalAlignment', 'baseline');
        set(get(ha(i),'Title'), 'FontSize', TitleFontSize);
        %set(get(ha(i),'Title'), 'FontWeight', 'Bold');
    end
    hc = get(ha(i), 'Children'); % the objects within an axis
    for j = 1:length(hc)
        chtype = get(hc(j), 'Type');
        if strcmp(chtype(1:length(chtype)), 'text')
            set(hc(j), 'FontSize', LegendFontSize); % 14 pt descriptive labels
        elseif strcmp(chtype(1:length(chtype)), 'line')
            set(hc(j), 'LineWidth', LineWidth);
            set(hc(j), 'MarkerSize', MarkerSize);
        elseif strcmp(chtype, 'hggroup')
            hcc = get(hc(j), 'Children');
            if strcmp(get(hcc, 'Type'), 'hggroup')
                hcc = get(hcc, 'Children');
            end
            for k = 1:length(hcc) % all elements are 'line'
                set(hcc(k), 'LineWidth', LineWidth);
                set(hcc(k), 'MarkerSize', LegendFontSize);
            end
        end
    end
end
end

