% This code couples GO terms to proteins and plot them in a graph.

% INPUT
load('core_names')
load('GO_MF')
load('interact_matrix_small_phys')

% OUTPUT
% GEF
% GAP
% GTPase

% GAPs and GEFs
GEF = unique(GO_MF(contains(GO_MF(:,2),'guanyl'),1));
GAP = unique(GO_MF(contains(GO_MF(:,2),'GTPase activator'),1));
GTPase = unique(GO_MF(contains(GO_MF(:,2),'GTPase activity'),1));

figure
hold on
plot([NaN NaN], [NaN NaN], 'r')
plot([NaN NaN], [NaN NaN], 'g')
plot([NaN NaN], [NaN NaN], 'm')
plot([NaN NaN], [NaN NaN], 'b')
set(gca, 'ColorOrderIndex', 1)
set(gca,'XColor', 'none','YColor','none')
h = plot(graph(interact_matrix_small,core_names));
highlight(h,GEF,'NodeColor','r')
highlight(h,GAP,'NodeColor','g')
highlight(h,GTPase,'NodeColor','m')
legend('GEF','GAP','GTPase','other')