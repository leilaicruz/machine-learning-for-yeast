% This code calculates the matching index.

% INPUT
load('network')
load('network_graph')
load('names')

% remove empty networks
network(26) = [];
network(1) = [];
network_graph(26) = [];
network_graph(1) = [];
names(26) = [];
names(1) = [];

% calculate matching indices
% Matching index = distinct common neighbours / total number of neighbours
matching_index = cell(1,length(network_graph));
index = zeros(1,length(names));
max_ind = zeros(1,length(network_graph));
match_ind = cell(1,length(network_graph));
match_CDC = cell(1,length(network_graph));

ind_rac = zeros(1,length(network_graph));
ind_boi = zeros(1,length(network_graph));
match_rac = zeros(1,length(network_graph));
match_boi = zeros(1,length(network_graph));

for ii=1:length(network_graph)

    matching_index{ii} = zeros(numnodes(network_graph{ii}));
    
    % calculate matching index
    for jj=1:numnodes(network_graph{ii})
        added = network{ii} + network{ii}(jj,:);
        nonzero = sum((added~=0),2);
        overlap = sum((added==2),2);
        matching_index{ii}(:,jj) = overlap./nonzero;
        matching_index{ii}(jj,jj) = NaN;
    end
    
    % find the most resembling protein
    if sum(strcmpi(names{ii},'CDC42')) ~= 0
        index(ii) = find(strcmpi(names{ii},'CDC42'));
        max_ind(ii) = max(matching_index{ii}(:,index(ii)));
        match_ind{ii} = find(matching_index{ii}(:,index(ii)) == max_ind(ii));
        match_CDC{ii} = names{ii}(match_ind{ii});
        
        % find the matching indices of the most resembling proteins
        if sum(strcmpi(names{ii},'RAC1')) ~= 0
            ind_rac(ii) = find(strcmpi(names{ii},'RAC1'));
            match_rac(ii) = matching_index{ii}(index(ii),ind_rac(ii));
        end
        
        if sum(strcmpi(names{ii},'BOI2')) ~= 0
            ind_boi(ii) = find(strcmpi(names{ii},'BOI2'));
            match_boi(ii) = matching_index{ii}(index(ii),ind_boi(ii));
        end
    end
   
end

avg_rac = average(nonzeros(match_rac));
avg_boi = average(nonzeros(match_boi));

clear added
clear nonzero
clear overlap
clear ii
clear jj