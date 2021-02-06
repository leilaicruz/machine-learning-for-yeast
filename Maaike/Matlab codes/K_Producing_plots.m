% This code produces various plots.

load('all_variables_e_s')
figure
hold on
histogram(all_variables_e_s(6,:))
xline(23.5,'--r')
xlabel('#proteins')
ylabel('#networks')
legend('network size distribution','cut-off for small networks','Location','northwest')
hold off

figure
hold on
histogram(sum(existing,2))
histogram(sum(non_existing,2))
title('All proteins')
legend('Existing','Non-existing')
xlabel('#proteins')
ylabel('#observations')
hold off