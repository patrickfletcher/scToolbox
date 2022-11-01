function [clustid, newids] = rename_clusters(old_clustid, newids, options)
arguments
    old_clustid
    newids=[]
    options.sortby='none'
    options.datavals=[]
    options.sortdir='descend'
    options.do_numbers=false
    options.do_categorical=false
end

if ~iscategorical(old_clustid)
    old_clustid=categorical(old_clustid);
end
oldids=categories(old_clustid);
% oldids=unique(old_clustid,'stable');
ixs=1:length(old_clustid);

datavals=options.datavals;
if options.sortby~="none"
    switch options.sortby
        case "size"
            vals=arrayfun(@(x)nnz(old_clustid==x),oldids);
        case "meanval"
            vals=arrayfun(@(x)mean(datavals(old_clustid==x)),oldids);
        case "medval"
            vals=arrayfun(@(x)median(datavals(old_clustid==x)),oldids);
    end
    [~,ixs]=sort(vals,options.sortdir);
end

clustid=old_clustid;
if isempty(newids)
    newids=oldids(ixs);
    if options.do_numbers
        newids = 1:length(oldids);
        clustid=zeros(size(clustid));
    end
end

for i=1:length(oldids)
    clustid(old_clustid==oldids(ixs(i)))=newids(i);
end

if options.do_categorical
    clustid=categorical(clustid,newids,newids);
end