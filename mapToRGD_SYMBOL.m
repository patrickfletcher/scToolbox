%get the RGD SYMBOL for each gene in different references (for commmon mapping)
%-should match more than name: eg. Bex1 name and Refseq ID in UCSC is incorrect?
clear

doSave=0;

%% get RGD data
load RGDtable.mat

%%
N=height(ref);

id_ix=nan(N,1);
ens_id_ix=nan(N,1);
name_ix=nan(N,1);
ens_ix=nan(N,1);
start_ix=nan(N,1);
stop_ix=nan(N,1);

old_name_ix=nan(N,1);
refseq_ix=nan(N,1);

startEns=find(contains(ref.gene_id,'ENS'),1,'first');

D_start=nan(N,1);
D_stop=nan(N,1);
missing=[];
many=[];

%% Find genomic position matches
for i=1:N
    
%     for i=1:length(ref.chr
    
    %assume there is match in chr, strand
    chr_match=RGD.CHROMOSOME_6_0==ref.chr(i);
    strand_match=RGD.STRAND_6_0==ref.strand(i);
    thisRegion=chr_match&strand_match;
    regix=find(thisRegion);
    
    if ~isempty(regix)
        
    %find best start/stop position matches
    % check all positions listed for that gene?? or just most extreme?
    start_dist=RGD.START_POS_6_0(thisRegion)-ref.start(i);
    [min_start_dist,this_start]=min(abs(start_dist));
    start_ix(i)=regix(this_start);
    D_start(i)=start_dist(this_start);
    
    stop_dist=RGD.STOP_POS_6_0(thisRegion)-ref.stop(i);
    [min_stop_dist,this_stop]=min(abs(stop_dist));
    D_stop(i)=stop_dist(this_stop);
    stop_ix(i)=regix(this_stop);
    
    end
    
    if mod(i,round(N/10))==0
        fprintf('.')
    end
end

%%

%exact location match of both start and stop
match_start=D_start==0;
match_stop=D_stop==0;
match_startstop=match_start&match_stop;
match_startstopix=start_ix==stop_ix;
exactPosMatch=match_startstop&match_startstopix;
rgdix=start_ix(exactPosMatch);

exactPos=table();
exactPos.Ref_Name=ref.gene_name(exactPosMatch);
exactPos.Ref_ID=ref.gene_id(exactPosMatch);
% exactPos.Replace_Reason=repmat("Position",height(exactPos),1);
exactPos=[exactPos,rgd(rgdix,:)];

replaceableExactPos=exactPos.Ref_Name~=exactPos.RGD_Symbol;

posRename=exactPos(replaceableExactPos,:);


%LOC names = NCBI_GENE_ID
replaceableLOCNames=strcat('LOC',rgd.NCBI_ID)~=rgd.RGD_Symbol;
lmap=rgd(replaceableLOCNames,:);
locnames=strcat('LOC',lmap.NCBI_ID);
[~,ia,ib]=intersect(ref.gene_name,locnames);
lRename=table;
lRename.Ref_Name=ref.gene_name(ia);
lRename.Ref_ID=ref.gene_id(ia);
% lRename.Replace_Reason=repmat("NCBI_ID",height(lRename),1);
lRename=[lRename,lmap(ib,:)];

%RGD names = RGD_ID
replaceableRGDNames=strcat('RGD',rgd.RGD_ID)~=rgd.RGD_Symbol;
rmap=rgd(replaceableRGDNames,:);
rgdnames=strcat('RGD',rmap.RGD_ID);
[~,ia,ib]=intersect(ref.gene_name,rgdnames);
rRename=table;
rRename.Ref_Name=ref.gene_name(ia);
rRename.Ref_ID=ref.gene_id(ia);
% rRename.Replace_Reason=repmat("RGD_ID",height(rRename),1);
rRename=[rRename,rmap(ib,:)];

%ENS_IDs = ENSEMBL_ID
[~,ia,ib]=intersect(ref.gene_id,RGD.ENSEMBL_ID);
replaceableEnsNames=ref.gene_name(ia)~=RGD.SYMBOL(ib);
emap=rgd(replaceableEnsNames,:);

eRename=table;
eRename.Ref_Name=ref.gene_name(ia(replaceableEnsNames));
eRename.Ref_ID=ref.gene_id(ia(replaceableEnsNames));
% eRename.Replace_Reason=repmat("ENS_ID",height(eRename),1);
eRename=[eRename,rgd(ib(replaceableEnsNames),:)];

rgd_map=[posRename;lRename;rRename;eRename];
rgd_map=unique(rgd_map);

if doSave
save(reffile,'rgd_map','-append');
end





% D_stopOfStartMatch=RGD.STOP_POS_6_0(start_ix)-ncbi.stop;
% D_startOfStopMatch=RGD.START_POS_6_0(stop_ix)-ncbi.start;
% match_stopOfStart=D_stopOfStartMatch==0;
% match_startOfStop=D_startOfStopMatch==0;
% use_start_ix=(abs(D_start)+abs(D_stopOfStartMatch))<=(abs(D_stop)+abs(D_startOfStopMatch));


% % resolve matching
% match_id=~isnan(id_ix);
% match_ENS=~isnan(ens_id_ix);
% match_name=~isnan(name_ix);





    %     start_match=abs(start_dist)<tol;
    %     stop_match=abs(stop_dist)<tol2;
    
    
%     %first try exact symbol match
%     symbol_match_id=strcmp(RGD.SYMBOL(thisRegion),ncbi.gene_id(i));
%     if any(symbol_match_id)
%         id_ix(i)=regix(symbol_match_id);
%     end
%     %do more work to find a match
    
%     if i>=startEns
%         %exact match of ENSRNOG ID (for MT genes that were added) - doesn't work. original ncbi has the names as old_symbol
%         ensID_match=contains(RGD.ENSEMBL_ID(thisRegion),ncbi.gene_id(i));
%         if any(ensID_match)
%             ens_id_ix(i)=regix(ensID_match);
%         end
%     end
    
%     %this will return several matches when the NCBI "ID" has -2, -3 etc.
%     symbol_match_name=strcmp(RGD.SYMBOL(thisRegion),ncbi.gene_name(i));
%     if any(symbol_match_name)
%         name_ix(i)=regix(symbol_match_name);
%     end

    %always do the above, but if there's no symbol match do more work:
%     if ~any(symbol_match_id)
        
%         %old symbol match
%         old_symbol_match_name=false(size(regix));
%         for j=1:length(regix)
%             r=split(RGD.OLD_SYMBOL(regix(j)),';');
%             if any(r==ncbi.gene_name(i))
%                 old_symbol_match_name(j)=true;
%             end
%         end
%         if any(old_symbol_match_name)
%             matchix=find(old_symbol_match_name);
%             old_name_ix(i,1:length(matchix))=regix(matchix);
%         end
        
%         %genbank nucleotide match (transcript IDs)
%         thisTID=ncbi.transcript_id{i}';
%         refseq_match=false(size(regix));
%         for j=1:length(regix)
%             r=split(RGD.GENBANK_NUCLEOTIDE(regix(j)),';');
%             for k=1:length(thisTID)
%                 if any(r==thisTID(k))
%                     refseq_match(j)=true;
%                 end
%             end
%         end
%         if any(refseq_match)
%             matchix=find(refseq_match);
%             refseq_ix(i,1:length(matchix))=regix(matchix);
%         end
%     end




% match_namestart=name_ix==start_ix;
% match_namestop=name_ix==stop_ix;
% 
% match_refseq_start=refseq_ix(:,1)==start_ix;
% match_refseq_stop=refseq_ix(:,1)==stop_ix;
% % match_refseq_start=any(refseq_ix==start_ix,2);
% % match_refseq_stop=any(refseq_ix==stop_ix,2);
% 
% 
% tol1=1;
% tol2=1;
% start_lt_thr=abs(D_start)<tol1;
% stop_lt_thr=abs(D_stop)<tol2;
% one_lt_thr=start_lt_thr|stop_lt_thr;
% 
% % 1) ID matches symbol - this could fail in the cases of "name","name-2"
% match_ix=id_ix;
% 
% %add in the ENS_ID matches
% match_ix(match_ENS)=ens_id_ix(match_ENS);
% 
% % 2) if no ID match, check for position match & refseq match.
% match_ix(match_refseq_start&match_startstop&~match_id)=start_ix(match_refseq_start&match_startstop&~match_id);
% 
% %remaining: try position match alone - Dstart/stop 
% 
% rgdmap=table();
% rgdmap.nbci_id=ncbi.gene_id;
% rgdmap.nbci_name=ncbi.gene_name;
% rgdmap.posmatch=ncbi.gene_name;
% rgdmap.posmatch(~(start_lt_thr|stop_lt_thr))="";
% rgdmap.posmatch(start_lt_thr|stop_lt_thr)=RGD.SYMBOL(start_ix(start_lt_thr|stop_lt_thr));
% % rgdmap.startmatch=ncbi.gene_name;
% % rgdmap.startmatch(~start_lt_thr)="";
% % rgdmap.startmatch(start_lt_thr)=RGD.SYMBOL(start_ix(start_lt_thr));
% % rgdmap.stoptmatch=ncbi.gene_name;
% % rgdmap.stoptmatch(~stop_lt_thr)="";
% % rgdmap.stoptmatch(stop_lt_thr)=RGD.SYMBOL(start_ix(stop_lt_thr));
% 
% res=table();
% res.id_ix=id_ix;
% res.name_ix=name_ix;
% % sym.id_ix=rgd_ix_sym_id;
% res.start_ix=start_ix;
% res.stop_ix=stop_ix;
% res.match_startstop=match_startstop;
% res.D_start=D_start;
% res.D_stopOfStartMatch=D_stopOfStartMatch;
% res.D_stop=D_stop;
% res.D_startOfStopMatch=D_startOfStopMatch;
% res.use_start_ix=use_start_ix;
% res.

% ncbi=ncbi(:,1:13);
%
% ncbi.rgd_name=rgd_name;
% ncbi.rgd_id=rgd_id;
% ncbi.ens_id=ens_id;
% ncbi.rgd_chr=rgd_chr;
% ncbi.rgd_strand=rgd_strand;
% ncbi.rgd_start=rgd_start;
% ncbi.rgd_stop=rgd_stop;
% ncbi.rgd_D_start=D_start;
% ncbi.rgd_D_stop=D_stop;
%
% save ncbi_genetypes.mat ncbi







%     thismatch=false(height(RGD),1);
%     if any(symbol_match)
%         match_type(i,1)=true;
%         thismatch=symbol_match;
%     else
%         %next try exact match of RefSeq ID in GENBANK_NUCLEOTIDE
%         refseq_match=contains(RGD.GENBANK_NUCLEOTIDE,ncbi.transcript_id{i});
%         if any(refseq_match)
%             match_type(i,2)=true;
%             thismatch=refseq_match;
%         else
%             symbol_match_loose=contains(RGD.SYMBOL,ncbi.gene_name{i});
%             old_symbol_match=contains(RGD.OLD_SYMBOL,ncbi.gene_name{i}); %more than one
%             if any(symbol_match_loose)
%                 match_type(i,3)=true;
%                 thismatch=symbol_match_loose;
%             else
% %                 if any(old_symbol_match)
% %                     match_type(i,4)=true;
% %                     thismatch=old_symbol_match;
% %                 end
%             end
%         end
%     end

%     if any(thismatch)
%
%         ix=find(thismatch);
%
%
% %         thisGenbank=false(height(RGD),1);
% %         thisGenbank(ix(chr_match&strand_match&start_match))=true;
%
%         location_match=false(size(ix));
%         if (mix_start==mix_stop)
%             location_match(mix_start)=true;
%         end
%         matchix=ix(chr_match&strand_match&location_match); %match if name/chr/str all match
%
%         if ~isempty(matchix)
%             rgd_id{i}=RGD.GENE_RGD_ID(matchix);
%             ens_id{i}=RGD.ENSEMBL_ID(matchix);
%             rgd_name{i}=RGD.SYMBOL(matchix);
%             rgd_chr{i}=RGD.CHROMOSOME_60(matchix);
%             rgd_strand{i}=RGD.STRAND_60(matchix);
%             rgd_start(i)=str2double(RGD.START_POS_60(matchix));
%             rgd_stop(i)=str2double(RGD.STOP_POS_60(matchix));
%         end
%     else
%         matchix=[];
%
%     end
%
%     %sometimes still can't find it
%     if isempty(matchix)
%         missing(end+1)=i;
%     end
%
%     %sometimes still returns more than one
%     if sum(matchix)>1
%         many(end+1)=i;
%     end
