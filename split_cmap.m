function cmap=split_cmap(brewerseq1, brewerseq2, NameValueArgs)
arguments
    brewerseq1='Blues'
    brewerseq2='Reds'
    NameValueArgs.N=20
    NameValueArgs.Skip=5
    NameValueArgs.MidCol=[0.75,0.75,0.75]
end
nc=NameValueArgs.N;
nskip=NameValueArgs.Skip;
midcol=NameValueArgs.MidCol;

c_low=cbrewer('seq',brewerseq1,nc); c_low=flipud(c_low);
c_high=cbrewer('seq',brewerseq2,nc);
cmap=[c_low(1:end-nskip,:);midcol;c_high(1+nskip:end,:)];