function result = build_sparse_matrix(vals, indptr, indices, nrow, ncol, format)

% CSR(C):  
% - data, row_ind and col_ind satisfy the relationship a[row_ind[k], col_ind[k]] = data[k]
% - column(row) indices for row(col) i are stored in    indices[indptr[i]:indptr[i+1]] 
% and their corresponding values are stored in   vals[indptr[i]:indptr[i+1]] 

switch lower(format)
    case 'csr' 
        row_ind = indices;
        cdiff=diff(indptr);
        col_ind=zeros(size(row_ind), 'int64');
        for i=1:length(cdiff)
            col_ind(indptr(i)+1:indptr(i+1))=i*ones(cdiff(i),1,'like',row_ind);
        end

    case 'csc' 
        col_ind = indices;
        rdiff=diff(indptr);
        row_ind=zeros(size(col_ind), 'int64');
        for i=1:length(rdiff)
            row_ind(indptr(i)+1:indptr(i+1))=i*ones(rdiff(i),1,'like',col_ind);
        end
end

result=sparse(row_ind, col_ind, vals, nrow, ncol);