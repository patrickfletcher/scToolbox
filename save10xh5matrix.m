function save10xh5matrix(filename, X, gene_ids, gene_names, barcodes, options)
arguments
    filename, 
    X,  
    gene_ids, 
    gene_names, 
    barcodes, 
    options.cellchunks=16
    options.datachunks=128
    options.genechunks=16
    options.deflate=6
    options.shuffle=true
    options.datatype="double"
    options.indextype="int32"
    options.overwrite=true
    options.featuretype=""
    options.genome=""
end
% Write count data to an h5 file in 10X v3 format

%TODO: attributes etc..

shape = size(X);

if ~issparse(X)
    tic
    disp("convert to sparse")
    X=sparse(X);
    toc
end

[indices,c,data]=find(X);
indices = indices - 1; %zero indexed in 10X h5's

%construct cell-wise indptr
[~,ic]=unique(c,'last');
indptr = [0;ic];

% data = cast(data,options.datatype);
% indices = cast(indices,options.indextype);
% indptr = cast(indptr,options.indextype);
% shape = int32(shape);

featuretype=strings(size(gene_names));
genome=strings(size(gene_names));

cellchunksize = ceil(length(barcodes)/options.cellchunks);
datachunksize = ceil(length(data)/options.datachunks);
genechunksize = ceil(length(gene_ids)/options.genechunks);

if exist(filename,"file") && options.overwrite
    delete(filename)
end

tic
h5create(filename,"/matrix/barcodes", length(barcodes), Datatype="string", ...
    Chunksize=cellchunksize, Deflate=options.deflate, Shuffle=options.shuffle)
h5write(filename,"/matrix/barcodes", barcodes)

h5create(filename,"/matrix/data", length(data), Datatype=options.datatype, ...
    Chunksize=datachunksize, Deflate=options.deflate, Shuffle=options.shuffle)
h5write(filename,"/matrix/data", data)

h5create(filename,"/matrix/indices", length(indices), Datatype=options.indextype, ...
    Chunksize=datachunksize, Deflate=options.deflate, Shuffle=options.shuffle)
h5write(filename,"/matrix/indices", indices)

h5create(filename,"/matrix/indptr", length(indptr), Datatype=options.indextype, ...
    Chunksize=cellchunksize, Deflate=options.deflate, Shuffle=options.shuffle)
h5write(filename,"/matrix/indptr", indptr)

h5create(filename,"/matrix/shape", 2, Datatype="int32", ...
    Chunksize=2, Deflate=options.deflate, Shuffle=options.shuffle)
h5write(filename,"/matrix/shape", shape)

h5create(filename,"/matrix/features/feature_type", length(featuretype), Datatype="string", ...
    Chunksize=genechunksize, Deflate=options.deflate, Shuffle=options.shuffle)
h5write(filename,"/matrix/features/feature_type", featuretype)

h5create(filename,"/matrix/features/genome", length(genome), Datatype="string", ...
    Chunksize=genechunksize, Deflate=options.deflate, Shuffle=options.shuffle)
h5write(filename,"/matrix/features/genome", genome)

h5create(filename,"/matrix/features/id", length(gene_ids), Datatype="string", ...
    Chunksize=genechunksize, Deflate=options.deflate, Shuffle=options.shuffle)
h5write(filename,"/matrix/features/id", gene_ids)

h5create(filename,"/matrix/features/name", length(gene_names), Datatype="string", ...
    Chunksize=genechunksize, Deflate=options.deflate, Shuffle=options.shuffle)
h5write(filename,"/matrix/features/name", gene_names)

toc

% matlab string default:
% Datatype:   H5T_STRING
%     String Length: variable      <------------
%     Padding: H5T_STR_NULLTERM    <------------       
%     Character Set: H5T_CSET_UTF8
%     Character Type: H5T_C_S1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% example V3 chem:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h5disp(f3)
% HDF5 filtered_feature_bc_matrix.h5 
% Group '/' 
%     Attributes:
%         'filetype':  'matrix'
%         'version':  2
%         'software_version':  'cellranger-6.0.0'
%         'library_ids':  'E19A'
%         'original_gem_groups':  1 
%         'chemistry_description':  'Single Cell 3' v3'
%     Group '/matrix' 
%         Dataset 'barcodes' 
%             Size:  7143
%             MaxSize:  7143
%             Datatype:   H5T_STRING
%                 String Length: 18
%                 Padding: H5T_STR_NULLPAD
%                 Character Set: H5T_CSET_ASCII
%                 Character Type: H5T_C_S1
%             ChunkSize:  447
%             Filters:  deflate(4)
%         Dataset 'data' 
%             Size:  15867120
%             MaxSize:  Inf
%             Datatype:   H5T_STD_I32LE (int32)
%             ChunkSize:  80000
%             Filters:  shuffle, deflate(4)
%         Dataset 'indices' 
%             Size:  15867120
%             MaxSize:  Inf
%             Datatype:   H5T_STD_I64LE (int64)
%             ChunkSize:  80000
%             Filters:  shuffle, deflate(4)
%         Dataset 'indptr' 
%             Size:  7144
%             MaxSize:  Inf
%             Datatype:   H5T_STD_I64LE (int64)
%             ChunkSize:  80000
%             Filters:  shuffle, deflate(4)
%         Dataset 'shape' 
%             Size:  2
%             MaxSize:  Inf
%             Datatype:   H5T_STD_I32LE (int32)
%             ChunkSize:  80000
%             Filters:  shuffle, deflate(4)
%         Group '/matrix/features' 
%             Dataset '_all_tag_keys' 
%                 Size:  1
%                 MaxSize:  1
%                 Datatype:   H5T_STRING
%                     String Length: 6
%                     Padding: H5T_STR_NULLPAD
%                     Character Set: H5T_CSET_ASCII
%                     Character Type: H5T_C_S1
%                 ChunkSize:  []
%                 Filters:  none
%                 FillValue:  '      '
%             Dataset 'feature_type' 
%                 Size:  32883
%                 MaxSize:  32883
%                 Datatype:   H5T_STRING
%                     String Length: 15
%                     Padding: H5T_STR_NULLPAD
%                     Character Set: H5T_CSET_ASCII
%                     Character Type: H5T_C_S1
%                 ChunkSize:  1028
%                 Filters:  deflate(4)
%             Dataset 'genome' 
%                 Size:  32883
%                 MaxSize:  32883
%                 Datatype:   H5T_STRING
%                     String Length: 19
%                     Padding: H5T_STR_NULLPAD
%                     Character Set: H5T_CSET_ASCII
%                     Character Type: H5T_C_S1
%                 ChunkSize:  1028
%                 Filters:  deflate(4)
%             Dataset 'id' 
%                 Size:  32883
%                 MaxSize:  32883
%                 Datatype:   H5T_STRING
%                     String Length: 18
%                     Padding: H5T_STR_NULLPAD
%                     Character Set: H5T_CSET_ASCII
%                     Character Type: H5T_C_S1
%                 ChunkSize:  1028
%                 Filters:  deflate(4)
%             Dataset 'name' 
%                 Size:  32883
%                 MaxSize:  32883
%                 Datatype:   H5T_STRING
%                     String Length: 18
%                     Padding: H5T_STR_NULLPAD
%                     Character Set: H5T_CSET_ASCII
%                     Character Type: H5T_C_S1
%                 ChunkSize:  1028
%                 Filters:  deflate(4)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% example V2 chem:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f2="D:/scRNAseq/rat/pituitary/nov18/E6/F1/outs/filtered_gene_bc_matrices_h5.h5";
% h5disp(f2)
% HDF5 filtered_gene_bc_matrices_h5.h5 
% Group '/' 
%     Attributes:
%         'TITLE':  ''
%         'CLASS':  'GROUP'
%         'VERSION':  '1.0'
%         'FILTERS':  65793
%         'PYTABLES_FORMAT_VERSION':  '2.1'
%         'filetype':  'matrix'
%         'chemistry_description':  'Single Cell 3' v2'
%         'library_ids':  'F1'
%         'original_gem_groups':  1 
%     Group '/Ens_6' 
%         Attributes:
%             'TITLE':  ''
%             'CLASS':  'GROUP'
%             'VERSION':  '1.0'
%             'FILTERS':  65793
%         Dataset 'barcodes' 
%             Size:  1902
%             MaxSize:  1902
%             Datatype:   H5T_STRING
%                 String Length: 18
%                 Padding: H5T_STR_NULLTERM
%                 Character Set: H5T_CSET_ASCII
%                 Character Type: H5T_C_S1
%             ChunkSize:  3640
%             Filters:  shuffle, deflate(1)
%             Attributes:
%                 'CLASS':  'CARRAY'
%                 'VERSION':  '1.1'
%                 'TITLE':  ''
%         Dataset 'data' 
%             Size:  4813280
%             MaxSize:  4813280
%             Datatype:   H5T_STD_I32LE (int32)
%             ChunkSize:  32768
%             Filters:  shuffle, deflate(1)
%             Attributes:
%                 'CLASS':  'CARRAY'
%                 'VERSION':  '1.1'
%                 'TITLE':  ''
%         Dataset 'gene_names' 
%             Size:  32883
%             MaxSize:  32883
%             Datatype:   H5T_STRING
%                 String Length: 15
%                 Padding: H5T_STR_NULLTERM
%                 Character Set: H5T_CSET_ASCII
%                 Character Type: H5T_C_S1
%             ChunkSize:  4369
%             Filters:  shuffle, deflate(1)
%             Attributes:
%                 'CLASS':  'CARRAY'
%                 'VERSION':  '1.1'
%                 'TITLE':  ''
%         Dataset 'genes' 
%             Size:  32883
%             MaxSize:  32883
%             Datatype:   H5T_STRING
%                 String Length: 18
%                 Padding: H5T_STR_NULLTERM
%                 Character Set: H5T_CSET_ASCII
%                 Character Type: H5T_C_S1
%             ChunkSize:  3640
%             Filters:  shuffle, deflate(1)
%             Attributes:
%                 'CLASS':  'CARRAY'
%                 'VERSION':  '1.1'
%                 'TITLE':  ''
%         Dataset 'indices' 
%             Size:  4813280
%             MaxSize:  4813280
%             Datatype:   H5T_STD_I64LE (int64)
%             ChunkSize:  16384
%             Filters:  shuffle, deflate(1)
%             Attributes:
%                 'CLASS':  'CARRAY'
%                 'VERSION':  '1.1'
%                 'TITLE':  ''
%         Dataset 'indptr' 
%             Size:  1903
%             MaxSize:  1903
%             Datatype:   H5T_STD_I64LE (int64)
%             ChunkSize:  8192
%             Filters:  shuffle, deflate(1)
%             Attributes:
%                 'CLASS':  'CARRAY'
%                 'VERSION':  '1.1'
%                 'TITLE':  ''
%         Dataset 'shape' 
%             Size:  2
%             MaxSize:  2
%             Datatype:   H5T_STD_I32LE (int32)
%             ChunkSize:  16384
%             Filters:  shuffle, deflate(1)
%             Attributes:
%                 'CLASS':  'CARRAY'
%                 'VERSION':  '1.1'
%                 'TITLE':  ''

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% example dropletutils
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HDF5 dev_all5.h5 
% Group '/' 
%     Attributes:
%         'chemistry_description':  'Single Cell 3' v3'
%         'filetype':  'matrix'
%         'library_ids':  'custom'
%         'original_gem_groups':  1 
%         'version':  2 
%     Group '/matrix' 
%         Dataset 'barcodes' 
%             Size:  60773
%             MaxSize:  60773
%             Datatype:   H5T_STRING
%                 String Length: 21
%                 Padding: H5T_STR_NULLPAD
%                 Character Set: H5T_CSET_ASCII
%                 Character Type: H5T_C_S1
%             ChunkSize:  60773
%             Filters:  shuffle, deflate(6)
%             Attributes:
%                 'rhdf5-NA.OK':  1 
%         Dataset 'data' 
%             Size:  156114576
%             MaxSize:  156114576
%             Datatype:   H5T_IEEE_F64LE (double)
%             ChunkSize:  156114576
%             Filters:  shuffle, deflate(6)
%             Attributes:
%                 'rhdf5-NA.OK':  1 
%         Dataset 'indices' 
%             Size:  156114576
%             MaxSize:  156114576
%             Datatype:   H5T_STD_I32LE (int32)
%             ChunkSize:  156114576
%             Filters:  shuffle, deflate(6)
%             Attributes:
%                 'rhdf5-NA.OK':  1 
%         Dataset 'indptr' 
%             Size:  60774
%             MaxSize:  60774
%             Datatype:   H5T_STD_I32LE (int32)
%             ChunkSize:  60774
%             Filters:  shuffle, deflate(6)
%             Attributes:
%                 'rhdf5-NA.OK':  1 
%         Dataset 'shape' 
%             Size:  2
%             MaxSize:  2
%             Datatype:   H5T_STD_I32LE (int32)
%             ChunkSize:  2
%             Filters:  shuffle, deflate(6)
%             Attributes:
%                 'rhdf5-NA.OK':  1 
%         Group '/matrix/features' 
%             Dataset '_all_tag_keys' 
%                 Size:  1
%                 MaxSize:  1
%                 Datatype:   H5T_STRING
%                     String Length: 6
%                     Padding: H5T_STR_NULLPAD
%                     Character Set: H5T_CSET_ASCII
%                     Character Type: H5T_C_S1
%                 ChunkSize:  1
%                 Filters:  shuffle, deflate(6)
%                 Attributes:
%                     'rhdf5-NA.OK':  1 
%             Dataset 'feature_type' 
%                 Size:  32883
%                 MaxSize:  32883
%                 Datatype:   H5T_STRING
%                     String Length: 15
%                     Padding: H5T_STR_NULLPAD
%                     Character Set: H5T_CSET_ASCII
%                     Character Type: H5T_C_S1
%                 ChunkSize:  32883
%                 Filters:  shuffle, deflate(6)
%                 Attributes:
%                     'rhdf5-NA.OK':  1 
%             Dataset 'genome' 
%                 Size:  32883
%                 MaxSize:  32883
%                 Datatype:   H5T_STRING
%                     String Length: 7
%                     Padding: H5T_STR_NULLPAD
%                     Character Set: H5T_CSET_ASCII
%                     Character Type: H5T_C_S1
%                 ChunkSize:  32883
%                 Filters:  shuffle, deflate(6)
%                 Attributes:
%                     'rhdf5-NA.OK':  1 
%             Dataset 'id' 
%                 Size:  32883
%                 MaxSize:  32883
%                 Datatype:   H5T_STRING
%                     String Length: 18
%                     Padding: H5T_STR_NULLPAD
%                     Character Set: H5T_CSET_ASCII
%                     Character Type: H5T_C_S1
%                 ChunkSize:  32883
%                 Filters:  shuffle, deflate(6)
%                 Attributes:
%                     'rhdf5-NA.OK':  1 
%             Dataset 'name' 
%                 Size:  32883
%                 MaxSize:  32883
%                 Datatype:   H5T_STRING
%                     String Length: 18
%                     Padding: H5T_STR_NULLPAD
%                     Character Set: H5T_CSET_ASCII
%                     Character Type: H5T_C_S1
%                 ChunkSize:  32883
%                 Filters:  shuffle, deflate(6)
%                 Attributes:
%                     'rhdf5-NA.OK':  1 