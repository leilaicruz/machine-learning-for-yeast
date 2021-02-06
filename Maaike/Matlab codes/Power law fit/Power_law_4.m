% distribution
x = degree_per_node{1};
%x_t = cell2mat(table2cell(readtable('Testdata1.txt')));

[alpha2, xmin2, L2] = plfit(x);            % plfit
h = plplot(x, xmin, alpha);             % plplot
[alpha_var, xmin_var, nn] = plvar(x);	% plvar
[pvalue2,gof2] = plpva(x, xmin);          % plpva