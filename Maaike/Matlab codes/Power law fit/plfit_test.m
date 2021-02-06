x = degree_per_node;
[alpha, xmin, L] = plfit(x);            % plfit
%%
h = plplot(x, xmin, alpha);             % plplot
[alpha_var, xmin_var, nn] = plvar(x);	% plvar
%%
[pvalue,gof] = plpva(x, xmin);          % plpva          