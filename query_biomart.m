function genetab = query_biomart(dataset, attributes, filter_name, filter_values)
arguments
    dataset
    attributes
    filter_name = []
    filter_values = []
end
    
attributes_file = 'tmp_attributes.csv';
vals_file = 'tmp_filter_values.csv';
result_file = 'tmp_result.csv';

%if nargin==1 print list of attributes?

if ~isstring(attributes)
    attributes=string(attributes(:));
end
if ~isstring(filter_values)
    filter_values=string(filter_values(:));
end

writematrix(attributes,attributes_file)
writematrix(filter_values,vals_file)


Rpath = [FindRpath, filesep, 'Rscript.exe', '" "', '--vanilla '];
scriptfile = "C:\Users\fletcherpa\Documents\GitHub\scToolbox\query_biomaRt.R";
command = char(strjoin([scriptfile,dataset,attributes_file,filter_name,vals_file,result_file]," "));

commandline=['"', Rpath, command '"'];

[status, cmdout]=system(commandline);
if status~=0
    error(cmdout)
end
disp(cmdout)

genetab=readtable(result_file);
[~,Locb]=ismember(filter_values,genetab.ensembl_gene_id);
genetab=genetab(Locb,:);

% genetab=sortrows(genetab);