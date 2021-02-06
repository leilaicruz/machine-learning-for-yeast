% This code makes matrix scatterplots of network characteristics.

% INPUT
load('all_variables_e')
load('all_variables_ne')

% OUTPUT
% dataplots

all_var = all_variables_e;       % all_variables_(n)e(_s)

all_var_names = {'average cc','average degree','average pathlength','density','edges','nodes','alpha','p value'};
var_name = {'density'};

index = [repmat("existing",length(all_var),1);repmat("non-existing",length(all_var),1)];

figure
gplotmatrix([transpose(all_variables_e_s);transpose(all_variables_ne_s)],[],index,['b','r'],[],[],'off',[],all_var_names)