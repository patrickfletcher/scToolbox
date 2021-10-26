classdef upsetplot < handle %graphics object??
    
    %for set intersection visuals of >3 sets (or alternative to venn for 3)
    
    %make something like: https://github.com/ImSoErgodic/py-upset, https://github.com/hms-dbmi/UpSetR
    
    %exposed properties
    properties
        
        %data selection and display order options
        mincount=1
        order_method='size'
%         order_method='comb'
%         n_member_thr=Inf %max num members in a combo
        combs
        
        %data
        excl_inter
        inter_counts
        N
        
        %deviation from expected counts (%)
        expected_frac
        observed_frac
        dev
        
        %bar chart appearance shortcuts
        
        %combination plot appearance shortcuts
    end
    
    %internal handles and data (Access=private) 
    properties 
        tiles
        ax_intersect %bar showing size of intersections
        h_bar %the bar chart object
        ht_counts %text labels for counts above bars
        
        ax_combs %line/dot plot showing set combination for each bar
        hl_bg_dot %gray background line plots with wider marker/lines
        hl_combs %the set of line plots indicating combinations
        
        ax_setcount %bar showing size of each set - alt just text on right y-axis
    end
    
    methods 
        function hup=upsetplot(sets, setnames, excl_inter, setix, options)
            arguments
                sets
                setnames
                excl_inter=[]
                setix=[]
                options.ordermethod='size'
                options.mincount=1
                options.pooled_color=0.4*[1,1,1]
                options.pooled_marker='d'
                options.set_sortorder='none'
                options.Padding='compact'
                options.TileSpacing='tight';
                options.mincomb=0 %restrict to at least k-combs
                options.maxcomb=Inf %restrict to up to k-combs
                options.showremaining=true
                options.showdev=false
            end

            hup.order_method = options.ordermethod;
            hup.mincount = options.mincount;
            
            %prepare data
            for i=1:length(sets), sets{i}=sets{i}(:); end
            N=length(unique(cat(1,sets{:})));
            n_sets=length(sets);
            set_sizes=cellfun(@length,sets);
            switch options.set_sortorder
                case 'ascend'
                    [~,ixs]=sort(set_sizes,'ascend');
                case 'descend'
                    [~,ixs]=sort(set_sizes,'descend');
                case 'alpha'
                    [~,ixs]=natsort(cellstr(setnames));
                otherwise
                    ixs=1:n_sets;
            end
            sets=sets(ixs);
            set_sizes=set_sizes(ixs);
            setnames=setnames(ixs);
            if isempty(excl_inter)
                [excl_inter,setix]=get_exclusive_intersections(sets, setnames,'max_k',options.maxcomb);
%             else
%                 excl_inter=excl_inter(ixs);
%                 setix=setix(ixs);
            end
            inter_set_counts=cellfun(@length,setix);
            
            discard=inter_set_counts<options.mincomb;
            excl_inter(discard)=[];
            setix(discard)=[];
            
            inter_counts=cellfun(@length,excl_inter);
            
            Ninter=sum(inter_counts);
            %use N-Ninter to set discard_counts initially. N above is total
            
            discard=inter_counts<hup.mincount;
            discard_counts=sum(inter_counts(discard))+(N-Ninter);
            
            setix(discard)=[];
            inter_counts(discard)=[];
            
            %here need to remove any rows (set names) that have no
            %remaining items.
            remainsets=unique(cat(1,setix{:}));
            discardsets=setdiff((1:n_sets)',remainsets);
            if ~options.showremaining
                n_sets=length(remainsets);
                set_sizes(discardsets)=[];
                setnames(discardsets)=[];
                setix=cellfun(@(x) find(ismember(remainsets,x)),setix,'UniformOutput',false);
            end
            
            %prepare axis layout
            if options.showdev
                hup.tiles=tiledlayout(3,1);
            else
                hup.tiles=tiledlayout(2,1);
            end
            hup.tiles.Padding=options.Padding;
            hup.tiles.TileSpacing=options.TileSpacing;
            ax_intersect=nexttile(hup.tiles);
            ax_combs=nexttile(hup.tiles); 
            
            expected_frac=zeros(size(inter_counts));
            for i=1:length(inter_counts)
                expected_frac(i) = prod(set_sizes(setix{i})/N) * prod(1-set_sizes(setdiff(1:n_sets,setix{i}))/N);
            end
            observed_frac = inter_counts/N;
            dev= observed_frac - expected_frac;
            dev=dev*100;
            
            switch hup.order_method
                case 'size'
                    [~, ixs]=sort(inter_counts,'descend');
                case 'comb'
                    ixs=1:length(inter_counts);
                case 'combrev'
                    ixs=length(inter_counts):-1:1;
                case 'dev'
                    [~, ixs]=sort(dev,'descend');
                case 'absdev'
                    [~, ixs]=sort(abs(dev),'descend');
                case 'expected'
                    [~, ixs]=sort(expected_frac,'descend');
            end
            s_inter_counts=inter_counts(ixs);
            s_observed_frac=observed_frac(ixs);
            s_expected_frac=expected_frac(ixs);
            s_dev=dev(ixs);
            s_setix=setix(ixs);
            
            if discard_counts>0 && options.showremaining
                s_inter_counts(end+1)=discard_counts;
                s_setix{end+1}=1:n_sets;
            end
            
            n_combs=length(s_inter_counts);
            XLIM=[0.25,n_combs+0.75];
            [~,max_ix]=max(s_inter_counts);
            
            %%%plotting
            %bar chart for intersection counts
            comb_ix=1:n_combs;
            set_ix=1:n_sets;
            h_bar=bar(ax_intersect,comb_ix,s_inter_counts);
            h_bar.FaceColor=[0.8,0.8,0.8];
            
            ax_intersect.XLim=XLIM;
            ax_intersect.XTick=[];
            ax_intersect.XTickLabel=[];
            ax_intersect.Box='off';
            
            %text labels for counts
            t_gap=0.1*max(s_inter_counts);
            ht_counts=text(ax_intersect,comb_ix,h_bar.YEndPoints+t_gap,string(s_inter_counts));
            for i=1:length(ht_counts)
                ht_counts(i).Rotation=90;
            end
            YLIM=ax_intersect.YLim;
            YLIM(2)=ht_counts(max_ix).Extent(2)+ht_counts(max_ix).Extent(4)+t_gap;
            ax_intersect.YLim=YLIM;
            
            %%%
            %line plot for combinations
            yyaxis(ax_combs,'left')
            yback=repmat(set_ix(:), 1, n_combs);
            xback=repmat(comb_ix, n_sets, 1);
            hl_bg_dot=line(ax_combs,xback,yback,'color',[0.8,0.8,0.8],...
                'linestyle','none','marker','.','markersize',30);
            
            %preallocate with gobjects?
            for i=1:n_combs
                yfore=s_setix{i};
                xfore=comb_ix(i)*ones(size(yfore));
                hl_combs(i)=line(ax_combs,xfore,yfore,'color','k',...
                    'linewidth',1,'marker','.','markersize',20);
            end
            if discard_counts>0 && options.showremaining
                hl_combs(end).Marker=options.pooled_marker;
                hl_combs(end).MarkerSize=4;
                hl_combs(end).MarkerFaceColor=options.pooled_color;
                hl_combs(end).Color=options.pooled_color;
            end
            
            YLIMcombs=[0.5,n_sets+0.5];
            
            ax_combs.XLim=XLIM;
            ax_combs.YLim=YLIMcombs;
            ax_combs.XTick=[];
            ax_combs.XAxis.Visible='off';
%             ax_combs.XTickLabel=[];
            ax_combs.YTick=set_ix;
            ax_combs.YTickLabels=setnames;
            ax_combs.YDir='reverse';
            
            %a dummy axis to put set numbers on the right
            yyaxis(ax_combs,'right')
            ax_combs.YLim=YLIMcombs;
            ax_combs.YTick=set_ix;
            ax_combs.YTickLabels=string(set_sizes);
            ax_combs.YDir='reverse';
            ax_combs.YAxis(1).Color='k';
            ax_combs.YAxis(2).Color='k';
            
            all_ax=[ax_intersect,ax_combs];
            
            %%%
            %dev values
            if options.showdev
            end
            
            linkaxes(all_ax,'x')
            
            hup.N=N;
            hup.observed_frac=s_observed_frac;
            hup.expected_frac=s_expected_frac;
            hup.dev=s_dev;
            hup.combs=s_setix;
            hup.inter_counts=s_inter_counts;
            hup.excl_inter=excl_inter;
            
            hup.ax_combs=ax_combs;
            hup.ax_intersect=ax_intersect;
            hup.hl_combs=hl_combs;
            hup.h_bar=h_bar;
            hup.ht_counts=ht_counts;
        end
        
        %getters/setters
    end
end