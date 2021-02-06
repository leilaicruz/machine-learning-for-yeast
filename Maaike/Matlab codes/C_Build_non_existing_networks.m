% This code generates non-existing networks, either with a random
% probability distribution or one equal to that of existing networks.

% INPUT
%load('existing')           % from Repertoire_Edited
choice = 0;                 % way of constructing the non-existing networks: 0 for 'uniform', 1 for 'representable'

% OUTPUT
% non_existing: networks of equal size and optionally probability distribution

% define probabilities
if choice == 0
    prob = sum(existing,'all')/(size(existing,1)*size(existing,2));   
elseif choice == 1 
    prob = sum(existing,1)/size(existing,1);
end

% generate non-existing networks
non_existing = rand(size(existing));
non_existing = non_existing < prob;

% make sure that non_existing doesn't overlap with existing
for ii=1:size(existing,1)
    for jj=1:size(non_existing,2)
        if isequal(existing(ii,:),non_existing(jj,:))
            non_existing(jj,:) = rand(size(existing,1));
            non_existing = non_existing < prob;
            jj=jj-1;
        end
    end
end

% plot size distribution to check manually for stochastic deviations
figure
hold on
histogram(sum(existing,2))
histogram(sum(non_existing,2))
legend('existing','non-existing')
hold off

clear choice
clear ii
clear jj
clear prob