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
        
%         ax_setcount %bar showing size of each set
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
                options.maxcomb=Inf %really need to add "other combos" bucket
            end

            hup.order_method = options.ordermethod;
            hup.mincount = options.mincount;
            
            %prepare data
            n_sets=length(sets);
            set_sizes=cellfun(@length,sets);
            if isempty(excl_inter)
                [excl_inter,setix]=get_exclusive_intersections(sets, setnames, 'max_k', options.maxcomb);
            end
            inter_counts=cellfun(@length,excl_inter);
                        
            %prepare axis layout
            hup.tiles=tiledlayout(2,1);
%             hup.tiles=tiledlayout(3,1);
            hup.tiles.TileSpacing='tight';
            hup.tiles.Padding='compact';
            ax_intersect=nexttile(hup.tiles);
            ax_combs=nexttile(hup.tiles);   
            
%             if ismember(hup.order_method,{'dev','absdev','expected'})
%                 subtiles=tiledlayout(hup.tiles,2,1);
%                 ax_combs=nexttile(subtiles);
%                 ax_dev=nexttile(subtiles);
%             else
%                 ax_combs=nexttile(hup.tiles);   
%             end
            
            N=sum(inter_counts);
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
                
            discard=s_inter_counts<hup.mincount;
            ixs(discard)=[];
            s_setix(discard)=[];
            s_inter_counts(discard)=[];
            s_observed_frac(discard)=[];
            s_expected_frac(discard)=[];
            s_dev(discard)=[];
            
            n_combs=length(ixs);
            XLIM=[0.25,n_combs+0.75];
            [~,max_ix]=max(s_inter_counts);
            
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
            
            %line plot for combinations
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
            ax_combs.XLim=XLIM;
            ax_combs.YLim=[0.5,n_sets+0.5];
            ax_combs.XAxis.Visible='off';
%             ax_combs.XTick=comb_ix;
%             ax_combs.XTickLabel=[];
            ax_combs.YTick=set_ix;
            ax_combs.YTickLabels=setnames;
            ax_combs.YDir='reverse';
            
            linkaxes([ax_intersect,ax_combs],'x')
            
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