function [enrichment, I] = GSEA(feature, max_i, draw_figure)
if nargin < 2
    draw_figure = false;
end

N = length(feature); % total number of qualifying genes
M = sum(feature);    % number of qualifying in-group genes
if M == 0
    enrichment = 1;
    I = [];
    return;
end

%A = cumsum(feature);
PValues = ones(max_i, 1);
for i = 1:max_i
    A = sum(feature(1:i));
    B = sum(~feature(1:i));
    C = sum(feature(i+1:end));
    D = sum(~feature(i+1:end));
    [~, PValues(i)] = fishertest([A B; C D], 'Tail', 'right');
    %PValues(i) = sum(hygepdf(A(i):i, N, M, i));
end

[enrichment, I] = min(PValues);

if draw_figure
    figure; 
    hold on;
    bar(1:N, -0.5*feature, 0.4);
    ylim([-0.5, -log10(enrichment)]);
    xlim([1, max_i]);
    plot(1:max_i, -log10(PValues), '-');
    xlabel('Set #');
    ylabel('-log10(p-value)');
    plot(I, -log10(enrichment), 'or');
    text(I + 0.02*max_i, -log10(enrichment), ['p = ' sprintf('%.1e', enrichment)]);
    hold off;
end
