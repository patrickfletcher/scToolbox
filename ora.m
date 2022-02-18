function results = ora(queries, terms, domainsize, options)
arguments
    queries
    terms
    domainsize
    options.min_K=0
    options.str_intersections=false
    options.sortby='none'
end

%query, references: array of structs with fields "name" and "genes"

nQ=length(queries);
nRef=length(terms);

qnames=cat(1,queries(:).name);
term_names=cat(1,terms(:).name);
term_size=arrayfun(@(x)length(x.genes),terms);

results=table;
results.term_name=term_names(:);
results.term_size=term_size(:);

keep=true(size(term_size));
keep=keep & term_size>options.min_K;

terms=terms(keep);
results=results(keep,:);

resnames=["n_int","expected","fc","p","fdr","int"];

for q = 1:nQ
    thisName=queries(q).name;
    thisG=queries(q).genes;
    querysize=length(thisG);

    intersections=arrayfun(@(x)intersect(thisG,x.genes(:)),terms,'UniformOutput',false)';
    n_int=cellfun(@length,intersections(:));

%     overlaps=cell(nRef,1);
%     K=zeros(nRef,1);
%     X=zeros(nRef,1);
%     for i=1:nRef
%         overlaps{i}=intersect(thisG,references(i).genes);
%         K(i)=length(references(i).genes);
%         X(i)=length(overlaps{i});
%     end

    expectedfrac=results.term_size./domainsize;
    expected=expectedfrac.*querysize;
    fc=n_int./expected;
    
    p = hygecdf(n_int,domainsize,results.term_size,querysize,'upper'); %P(at least n_int)
%     p = hygepdf(n_int,domainsize,results.term_size,querysize); %P(exactly n_int)
    [~,~,~,fdr]=fdr_bh(p);

    thisResnames=resnames;
    if nQ>1
        thisResnames=thisResnames+"_"+thisName;
    end
    A=[n_int(:),expected(:),fc(:),p(:),fdr(:)];
    results=[results,array2table(A,VariableNames=thisResnames(1:end-1))];

    if options.str_intersections
        intersections=cellfun(@(x) strjoin(x,', '), intersections);
    end
    results.(thisResnames(end))=intersections(:);

end

resvars=results.Properties.VariableNames;
if nQ>1
    results.max_int=max(results{:,contains(resvars,"n_int")},[],2);
    results.min_fdr=min(results{:,contains(resvars,"fdr")},[],2);
    results.max_fc=max(results{:,contains(resvars,"fc")},[],2);
%     results=sortrows(results,'min_fdr','ascend');
else
%     results=sortrows(results,'fdr','ascend');
end