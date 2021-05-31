function [wc, ax]= plot_wordcloud(docs_raw)

commonstop=importdata('stopwords.txt');
biostop=importdata('bio-stopwords.txt');
stopwords=[commonstop;biostop;{'positive';'negative';'regulation'}];

% docs_raw - cell array of text
nDocs=length(docs_raw);
doc_terms={};
raw_counts={};
for i=1:nDocs
    %process the doc
    doc=docs_raw{i};
    doc = strjoin(doc,' '); 
    %here could replace some sequences 
    
    doc = split(doc); %would be nice to keep certain terms as a unit
    doc(ismember(doc,stopwords))=[];
    
    %porter stemmer?
%     docstem=cellfun(@(x)porterStemmer(x),doc,'UniformOutput',false);
%     [u_stem, iu, ia] = unique(docstem);
%     stem_counts= accumarray(ia, 1);
%     multi_stem=find(stem_counts>1);
%     for k=1:length(multi_stem)
%         thissstem=ia==ia(multi_stem(k));
%         stemset=doc(thissstem);
%         [~,smallest]=min(cellfun(@length,stemset));
%         doc(thissstem)=stemset(smallest);
%     end
    
    [u_term, ~, idxU] = unique(doc);
    doc_terms{i}=u_term;
    raw_counts{i} = accumarray(idxU, 1); 
end

all_terms=unique(cat(1,doc_terms{:}));
nTerms=length(all_terms);

%tf(t,d)=term_counts/length(d)
tf=zeros(nTerms,nDocs);
for i=1:nDocs
    [termsInDoc,locInDoc]=ismember(all_terms,doc_terms{i});
    locInDoc(locInDoc==0)=[];
    thisTf=raw_counts{i};
%     thisTf=raw_counts{i}/length(docs{i});
%     thisTf=raw_counts{i}/max(raw_counts{i});
    tf(termsInDoc,i)=thisTf(locInDoc);
end


docsWithTerm=cellfun(@(x)ismember(all_terms,x),doc_terms,'UniformOutput',false);
docsWithTerm=sum(cat(2,docsWithTerm{:}),2);
% idf=log10(nDocs./docsWithTerm);
idf=log10(nDocs./(docsWithTerm+1))+1;

%tf-idf 
tfidf=tf.*idf;

% words(ismember(words,stopwords))=[];
%
fh=figure();clf
% fh.Units='inches';
% fh.Position=[3,4,8,4];
ht=tiledlayout(1,nDocs,'TileSpacing','compact','Padding','compact');
for i=1:nDocs
    ax(i)=nexttile(ht);

%     fh=figure(i);clf
    term_ix=ismember(all_terms,doc_terms{i});
    terms=all_terms(term_ix);
    
    size_data=tfidf(term_ix,i);

    [size_data,ixs]=sort(size_data,'descend');
    terms=terms(ixs);

%     size_data=rescale(size_data,6,10);
    
    wc{i}=wordcloud(terms, size_data,'MaxDisplayWords',50,'LayoutNum',2,'SizePower',.5);
    
end
