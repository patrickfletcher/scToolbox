% wrapper for Paul Tol's color schemes
% https://personal.sron.nl/~pault/ 
function cmap = tol_col(name, args)
arguments
    name = "rainbow_PuRd"
    args.N = 256
    args.list=false
end
N=args.N;


%add path to python script if needed...
Pypath = py.sys.path;
MLpath=string(path).split(';');
sctoolpath=MLpath(contains(MLpath,'scToolbox'));
if count(Pypath,sctoolpath) == 0 && ~isempty(sctoolpath)
    insert(Pypath,int32(0),sctoolpath);
end

cmap_names=cellfun(@(x) string(x), cell(py.tol_colors.tol_cmap()) );
cset_names=cellfun(@(x) string(x), cell(py.tol_colors.tol_cset()) );
Nset=get_N(cset_names);
color_counts=table(); color_counts.cset_name=cset_names(:); color_counts.N=Nset;

names=[cmap_names, cset_names];
if nargin==0 || args.list
    cmap_names
    cset_names
    disp(color_counts)
    return
end

mustBeMember(name, names)

cmap = [];

if ismember(name,cmap_names)
    cmap = py.tol_colors.tol_cmap(name, int32(N));
    cmap = double(cmap(linspace(0,1,N)));
    cmap=cmap(:,1:3);
else
    cset = py.tol_colors.tol_cset(name);
    csethex=cellfun(@(x) string(x), cell(cset) );
    cmap = hex2rgb(cellstr(csethex));
    i=cset_names==name;
    if N<=Nset(i)
        cmap=cmap(1:N,:);
    else
        cmap=repmat(cmap,ceil(N/Nset(i)),1);
        cmap=cmap(1:N,:);
        warning(name+" contains only "+string(Nset(i))+" colors, but "+string(N)+" were requested. Repeating colors...");
    end
end

end

function N = get_N(cset_names)
N = zeros(length(cset_names),1);
for i=1:length(cset_names)
    name=cset_names(i);
    cset = py.tol_colors.tol_cset(name);
    csethex=cellfun(@(x) string(x), cell(cset) );
    N(i) = length(csethex);
end
end