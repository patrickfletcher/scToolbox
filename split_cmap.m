function cmap=split_cmap(mapname_low, mapname_high, args)
arguments
    mapname_low='Blues'
    mapname_high='Reds'
    args.maptype_low='seq'
    args.flip_low=true
    args.maptype_high='seq'
    args.flip_high=false
    args.N=64
    args.Skip=16
    args.nMid=1
    args.MidCol=[0.75,0.75,0.75]
    args.RepMidSkip=false
end
nc=args.N;
nskip=args.Skip;
midcol=repmat(args.MidCol,args.nMid,1);

%final map to have N levels in each direction, with grey in the middle to indicate zero
maptype_low=args.maptype_low;
switch maptype_low
    case {'div', 'seq', 'qual'}
        c_low=cbrewer(maptype_low,mapname_low,nc+nskip); 
%     case 'matlab'
%         c_low=
end
if args.flip_low
    c_low=flipud(c_low);
end

maptype_high=args.maptype_high;
switch maptype_high
    case {'div', 'seq', 'qual'}
        c_high=cbrewer(maptype_high,mapname_high,nc+nskip);
%     case 'matlab'
%         c_low=
end
if args.flip_high
    c_high=flipud(c_high);
end

if args.RepMidSkip
    midcol=repmat(midcol,nskip,1);
end

cmap=[c_low(1:end-nskip,:);midcol;c_high(1+nskip:end,:)];