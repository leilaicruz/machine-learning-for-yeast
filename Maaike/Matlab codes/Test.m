% waarom is de matrix niet symmetrisch?
% er gaat iets mis met de 'Y' proteins en IMP2'

yy=cell(1,size_XL);
for ii=1:size_XL
    xx=find(transpose(interact_matrix_XL(ii,:))~=interact_matrix_XL(:,ii));
    yy{ii}=gene_names_XL(xx);
end

size_yy = zeros(1,size_XL);
for ii=1:size_XL
    [~,size_yy(ii)] = size(yy{ii});
end
missing = cell(size_XL,1);
index = 0;
for ii=1:size_XL
    for jj=1:size_yy(ii)
        missing(index+jj) = yy{ii}(1,jj);
    end
    index=index+size_yy(ii);
end

zz = unique(missing,'stable');

%%
AA = readtable('results.csv');
%%
BB = table2cell(AA);
%%
pp=cellfun(@isempty,BB);
[qq,rr]=find(pp(2,:));
%%
CC=BB;
CC(qq,rr)=CC(qq,rr-1);