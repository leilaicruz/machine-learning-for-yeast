% This code generates small non-existing networks, either with a random
% probability distribution or one equal to that of core or small existing networks.

% INPUT
load('existing')
load('core_names')
choice = 2;                 % way of constructing the non-existing networks: 0 for 'equal', 1 for 'representable'   

% OUTPUT
% existing_s
% non_existing_s: networks of equal size and optionally probability distribution

% select smaller networks
lower_cutoff = 2;
upper_cutoff = 23;
existing_s = existing(lower_cutoff <= sum(existing,2) & sum(existing,2) <= upper_cutoff,:);

% define probabilities
if choice == 0          % uniform probability, small size
    prob = sum(existing_s,'all')/(size(existing_s,1)*size(existing_s,2));   
elseif choice == 1      % probability of small networks
    prob = sum(existing_s,1)/size(existing_s,1);
elseif choice == 2      % probability of core networks, corrected to create small networks
    prob = sum(existing,1)/(size(existing,1)) .* sum(sum(existing_s,1)/size(existing_s,1)) ./ sum(sum(existing,1)/(size(existing,1)));
elseif choice == 3      % uniform probability, large size
    prob = sum(existing,'all')/(size(existing,1)*size(existing,2)); 
end

% generate non-existing networks
non_existing_s = rand(size(existing_s));
non_existing_s = non_existing_s < prob;

% make sure that non_existing doesn't overlap with existing
for ii=1:size(existing_s,1)
    for jj=1:size(non_existing_s,1)
        if isequal(existing_s(ii,:),non_existing_s(jj,:))
            non_existing_s(jj,:) = rand(1,size(existing_s,2));
            non_existing_s = non_existing_s < prob;
            jj=jj-1;
        end
    end
end

fraction_s = sum(existing_s/size(existing_s,1));
fraction = sum(existing/size(existing,1));
fraction_eq = repmat(sum(existing,'all')/(size(existing,1)*size(existing,2)),1,length(core_names));
fraction_e = sum(existing,'all')/(size(existing,1)*size(existing,2));

figure
hold on
bar(categorical(core_names),[fraction;fraction_s])
yline(fraction_e,'--g')
legend('Distribution of all networks','Distribution of small networks','Uniform distribution (P = 0.57)')
ylabel('Probability of presence')
hold off

clear choice
clear ii
clear jj
clear prob
clear index
clear fraction
clear fraction_s
clear lower_cutoff
clear upper_cutoff