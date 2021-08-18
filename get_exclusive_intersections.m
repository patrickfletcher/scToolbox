function [exclusive_intersections, intersection_ix, int_combo_names]=get_exclusive_intersections(sets, setnames, options)
arguments
    sets
    setnames
    options.sortOrder = 'stable'
    options.max_k=Inf %really need to add "other combos" bucket
end
% computes the "exclusive intserctions" of the input sets. That is, the
% partitions of all elements into all possible combinations of
% intersections.
% input: one or more sets (cell array, each element is a set, elements of sets are set-operation compatible)
% outputs:
% - exclusive_intersections: cell array with intersection data
% - intersection_ix: which sets 

%TODO: output a nice struct instead. 

%TODO: intersection_ix should contain the indices of gene in each set

% sets=varargin;
n_sets=length(sets);

max_k=min(n_sets,options.max_k); %max k for nchoosek

n_exclusive_intersections=0;
for k=1:max_k
    n_exclusive_intersections=n_exclusive_intersections+nchoosek(n_sets,k);
end

set_ix_list=1:n_sets;
exclusive_intersections={};
intersection_ix={};
for k=1:max_k
    in_ix=nchoosek(set_ix_list,k);
%     if k>1, in_ix=flipud(in_ix); end %put combos with first ix's first
    for i=1:size(in_ix,1)
        %compute the intersection of in_sets
        this_in_ix=in_ix(i,:);
        this_excl_int=sets{in_ix(i,1)};
        for j=2:length(this_in_ix)
            [this_excl_int,ia,ib]=intersect(this_excl_int, sets{this_in_ix(j)},options.sortOrder);
        end
        %subtract away the other (out) sets
        out_ix=setdiff(set_ix_list, this_in_ix);
        for j=1:length(out_ix)
            [this_excl_int,ia]=setdiff(this_excl_int, sets{out_ix(j)},options.sortOrder);
        end
        exclusive_intersections(end+1)={this_excl_int};
        intersection_ix(end+1)={this_in_ix};
    end
end

for i=1:length(intersection_ix)
    int_combo_names{i}=strjoin(setnames(intersection_ix{i}),'_');
end




% % sets=varargin;
% n_sets=length(sets);
% 
% % n_exclusive_intersections=0;
% % for k=1:n_sets
% %     n_exclusive_intersections=n_exclusive_intersections+nchoosek(n_sets,k);
% % end
% 
% set_ix_list=1:n_sets;
% % exclusive_intersections={};
% intersection_ix={};
% for k=1:n_sets
%     in_ix=combnk(set_ix_list,k);
%     for i=1:size(in_ix,1)
%         %compute the intersection of in_sets
%         this_in_ix=in_ix(i,:);
%         this_excl_int=sets{in_ix(i,1)};
%         for j=2:length(this_in_ix)
%             this_excl_int=intersect(this_excl_int, sets{this_in_ix(j)});
%         end
%         %subtract away the other (out) sets
%         out_ix=setdiff(set_ix_list, this_in_ix);
%         for j=1:length(out_ix)
%             this_excl_int=setdiff(this_excl_int, sets{out_ix(j)});
%         end
% %         exclusive_intersections(end+1)={this_excl_int};
%         intersection_ix(end+1)={this_in_ix};
%         this_fieldname=strjoin(setnames(this_in_ix),"_");
%         exclusive_intersections.(this_fieldname)=this_excl_int;
%     end
% end