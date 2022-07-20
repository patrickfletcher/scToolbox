function vars = make_vars_unique(vars, options)
arguments
    vars
    options.sep="_"
end
vars=string(vars);
[~, ind] = unique(vars); %indices of first occurrence of each unique
duplicate_ind = setdiff(1:length(vars), ind); %missing indices are dupes
duplicate_vars = vars(duplicate_ind);

unique_dupes = unique(duplicate_vars);
for i=1:length(unique_dupes)
    this_dupe=duplicate_vars==unique_dupes(i);
    n=sum(this_dupe);
    new_names=strcat(duplicate_vars(this_dupe),options.sep,string((1:n)'));
    vars(duplicate_ind(this_dupe))=new_names;
end