function vars = make_vars_unique(vars)
vars=string(vars);
[~, ind] = unique(vars);
duplicate_ind = setdiff(1:length(vars), ind);
duplicate_vars = vars(duplicate_ind);

unique_dupes = unique(duplicate_vars);
for i=1:length(unique_dupes)
    this_dupe=duplicate_vars==unique_dupes(i);
    n=sum(this_dupe);
    new_names=strcat(duplicate_vars(this_dupe),"-",string((1:n)'));
    vars(duplicate_ind(this_dupe))=new_names;
end