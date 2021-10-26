function cmap=split_cmap(mapname_low, mapname_high, NameValueArgs)
arguments
    mapname_low='Blues'
    mapname_high='Reds'
    NameValueArgs.maptype_low='seq'
    NameValueArgs.flip_low=true
    NameValueArgs.maptype_high='seq'
    NameValueArgs.flip_high=false
    NameValueArgs.N=64
    NameValueArgs.Skip=16
    NameValueArgs.MidCol=[0.75,0.75,0.75]
end
nc=NameValueArgs.N;
nskip=NameValueArgs.Skip;
midcol=NameValueArgs.MidCol;

%final map to have N levels in each direction, with grey in the middle to indicate zero
maptype_low=NameValueArgs.maptype_low;
switch maptype_low
    case {'div', 'seq', 'qual'}
        c_low=cbrewer(maptype_low,mapname_low,nc+nskip); 
%     case 'matlab'
%         c_low=
end
if NameValueArgs.flip_low
    c_low=flipud(c_low);
end

maptype_high=NameValueArgs.maptype_high;
switch maptype_high
    case {'div', 'seq', 'qual'}
        c_high=cbrewer(maptype_high,mapname_high,nc+nskip);
%     case 'matlab'
%         c_low=
end
if NameValueArgs.flip_high
    c_high=flipud(c_high);
end

midcol=repmat(midcol,NameValueArgs.Skip,1);

cmap=[c_low(1:end-nskip,:);midcol;c_high(1+nskip:end,:)];