% INPUT
load('non_existing')
load('core_names')
load('GO')

% OUTPUT
% network_terms
% missing_terms
% unique_terms
% set
% GO_MF

network_terms = cell(1,length(existing));
missing_terms = cell(1,length(existing));

% which GO terms are present in which networks
for ii=1:length(network_terms)
    GO_names = core_names(existing(ii,:));
    for jj=1:length(GO_names)
        index = strcmpi(GO(:,1),GO_names(jj));
        network_terms{ii} = [network_terms{ii};GO(index,:)];
    end
end

% which GO terms are missing in which networks
for ii=1:length(network_terms)
    if isempty(network_terms{ii})
        missing_terms{ii} = unique(GO(:,3));
    else
        missing_terms{ii} = setdiff(GO(:,3),network_terms{ii}(:,3));
    end
end

unique_terms = cell(size(network_terms));
for ii=1:length(network_terms)
    if ~isempty(network_terms{ii})
        unique_terms{ii} = unique(network_terms{ii}(:,3));
    end
end

set = unique_terms{2};
for ii=3:length(unique_terms)
    if ~isempty(unique_terms)
        set = intersect(set,unique_terms{ii});
    end
end

% GAPs and GEFs
GO_MF = GO(strcmp(GO(:,2),'Molecular function'),:);