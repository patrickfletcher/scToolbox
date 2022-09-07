function ggalluvial(data,strat_var,alluv_var,fill_var,cols,figurefile,plotopts,figopts,options)
arguments
    data
    strat_var
    alluv_var
    fill_var
    cols
    figurefile
    plotopts.xvals="NULL"
    plotopts.yvals="proportions"
    plotopts.strat_levels="NULL"
    plotopts.alluv_levels="NULL"
    plotopts.fill_levels="NULL"
    figopts.width=6
    figopts.height=4
    figopts.resolution=150
    options.parfile="D:\tmp\tmp_scran_score_markers.mat"
    options.verbose=false
end 

if options.verbose
    disp("Running " + mfilename + "...")
end

%do the summary table here, write a smaller CSV. it is the ggplot df

typeage=stratify_factors(cells.age(:),type(:));
summary(typeage)

agecnt=countcats(cells.age);
typeagecnt=countcats(typeage); 
typeagecnt=reshape(typeagecnt,K,[])';

datafile='D:\tmp\tmp_ggalluvial.csv';
writetable(data,datafile)

figurefile=cellstr(figurefile);
strat_var=cellstr(strat_var);
alluv_var=cellstr(alluv_var);
fill_var=cellstr(fill_var);

% if ~isempty(cols)
    cols=cellstr(rgb2hex(cols));
% end

save(options.parfile,"datafile","figurefile","strat_var","alluv_var","fill_var","cols","figopts")

Rpath = [FindRpath, filesep, 'Rscript.exe', '" "', '--vanilla '];
scriptfile = mfilename('fullpath')+".R";
command = char(strjoin([scriptfile,options.parfile]," "));
commandline=['"', Rpath, command '"'];

tic
[status, cmdout]=system(commandline);

if status~=0
    error(cmdout)
end
if options.verbose
    disp(cmdout)
end
disp(mfilename+": ")
toc

%load and plot the image if it's a png
% figure
% I=imread()
