function hhm = heatmap_plot(cvals, opts)
arguments
    cvals
    opts.xorder=[]
    opts.yorder=[]
    
    opts.cmap=[]
    opts.isdiverging=false

    opts.title=[] %label to use as title
    opts.fig=[] %specify figure to plot in (if empty, new figure)
    opts.margins=[0.1,0.1,0.1,0.1] %[left, bottom, right, top]

    opts.cbgap=0.01 %space between ax and cb
    opts.cbdims=[0.3,0.01] %[length, width]
    opts.cblabel=[]
    opts.cbLoc {mustBeMember(opts.cbLoc,["east","west","north","south"])}='east'
    opts.cbJust {mustBeMember(opts.cbJust,["low","midlo","mid","midhi","high"])}='mid'
    opts.cbDigits=1
end