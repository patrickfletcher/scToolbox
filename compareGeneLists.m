function [unique1,unique2,common,iu1,iu2,ic1,ic2]=compareGeneLists(genes1,genes2)
genes1(strcmp(genes1,''))=[];
genes2(strcmp(genes2,''))=[];
[unique1,iu1]=setdiff(genes1,genes2);
[unique2,iu2]=setdiff(genes2,genes1);
[common,ic1,ic2]=intersect(genes1,genes2);



% [~,iu1,iu2]=setxor(genes1,genes2);
% unique1=genes1(iu1);
% unique2=genes2(iu2);