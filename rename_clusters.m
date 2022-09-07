function clustid = rename_clusters(old_clustid, newids, options)
arguments
    old_clustid
    newids=[]
    options.sortby='size'
    options.datavals=[]
    options.sortdir='descend'
    options.do_numbers=false
end

oldids=unique(old_clustid);

clustid=old_clustid;
if isempty(newids)
    newids=oldids;
    if options.do_numbers
        newids = 1:length(oldids);
    end
end

datavals=options.datavals;
switch options.sortby
    case "size"
        vals=arrayfun(@(x)nnz(old_clustid==x),oldids);
    case "meanval"
        vals=arrayfun(@(x)mean(datavals(old_clustid==x)),oldids);
    case "medval"
        vals=arrayfun(@(x)median(datavals(old_clustid==x)),oldids);
end
[~,ixs]=sort(vals,options.sortdir);
for i=1:length(oldids)
    clustid(old_clustid==oldids(ixs(i)))=newids(i);
end