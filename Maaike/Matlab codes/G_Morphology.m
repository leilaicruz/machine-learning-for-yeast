% This code sorts morphological data and normalises them with respect to
% the wild-types of essential and non-essential genes.

% INPUT
load('ess')             % retrieved from SCMD2 datasets (select)
load('non_ess')
load('wt_e')
load('wt_ne')
load('name_ess')
load('name_non_ess')
load('core_names')
load('core_codes')

% OUTPUT
% corrected_sorted
% sorted

sorted = zeros(length(core_names),7);
sorted_ess = zeros(length(core_names),7);
sorted_non_ess = zeros(length(core_names),7);

% select information from morphology data
for ii=1:length(core_names)
    index = find(strcmpi(name_ess,core_codes(ii)));
    if ~isempty(index)
        sorted(ii,:) = ess(index,:);
        sorted_ess(ii,:) = ess(index,:);
    end
    index = find(strcmpi(name_non_ess,core_codes(ii)));
    if ~isempty(index)
        sorted(ii,:) = non_ess(index,:);
        sorted_non_ess(ii,:) = non_ess(index,:);
    end
end

% correct with respect to wt for both non-essential and essential
corrected_sorted = sorted_ess./mean(wt_e) + sorted_non_ess./mean(wt_ne);

% genes not represented in database are assumed to have no influence
corrected_sorted(corrected_sorted==0) = 1;

clear index
clear ii
clear sorted_ess
clear sorted_non_ess

% plot morphological influences
for ii=1:size(corrected_sorted_alt,2)
    
    figure
    hold on
    bar(categorical(core_names),corrected_sorted_alt(:,ii));
    
    ind = zeros(length(GAP),1);
    for jj=1:length(GAP)
        ind(jj) = find(strcmpi(core_names,GAP(jj)));
    end
    bar(categorical(GAP),corrected_sorted_alt(ind,ii),'g');
    
    ind = zeros(length(GEF),1);
    for jj=1:length(GEF)
        ind(jj) = find(strcmpi(core_names,GEF(jj)));
    end
    bar(categorical(GEF),corrected_sorted_alt(ind,ii),'r');
    
    ind = zeros(length(GTPase),1);
    for jj=1:length(GTPase)
        ind(jj) = find(strcmpi(core_names,GTPase(jj)));
    end
    bar(categorical(GTPase),corrected_sorted_alt(ind,ii),'m');

    title(morph_title(ii))
    if ii==1 || ii==2
        legend('other','GAP','GEF','GTPase','Location','northwest')
    else
        legend('other','GAP','GEF','GTPase','Location','northeast')
    end
    hold off
end

% delete proteins without data for plotting
corrected_sorted_alt(corrected_sorted_alt == 0) = NaN;
corrected_sorted_alt_abs = abs(corrected_sorted_alt);

% plot morphological influence against centralities
figure
gplotmatrix(corrected_sorted_alt,avg_cent_e,index_e,['b' 'g' 'r' 'm'],'o',[],[],[],morph_title,cent_names)