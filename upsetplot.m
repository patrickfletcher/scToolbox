classdef upsetplot < handle %graphics object??
    
    %for set intersection visuals of >3 sets (or alternative to venn for 3)
    
    %make something like: https://github.com/ImSoErgodic/py-upset, https://github.com/hms-dbmi/UpSetR
    
    %exposed properties
    properties
        
        %data selection and display order options
        count_thr=1
%         order_method='size'
        order_method='comb'
        n_member_thr=Inf %max num members in a combo
        
        %bar chart look shortcuts
        
        %combination plot look shortcuts
    end
    
    %internal handles and data
    properties (Access=private) 
        ax_intersect %bar showing size of intersections
        h_bar %the bar chart object
        
        ax_combs %line/dot plot showing set combination for each bar
        hl_bg_dot %gray background line plots with wider marker/lines
        hl_combs %the set of line plots indicating combinations
        
        ax_setcount %bar showing size of each set
    end
    
    methods 
        function hup=upsetplot(sets, setnames)
            
            %prepare axis layout
            t=tiledlayout(2,1);
            t.TileSpacing='compact';
            t.Padding='compact';
            ax_intersect=nexttile();
            ax_combs=nexttile();
            
            if nargin==0 
                return
            end
            
            %prepare data
            n_sets=length(sets);
            [excl_inter,setix]=get_exclusive_intersections(sets);
            set_counts=cellfun(@length,sets);
            inter_counts=cellfun(@length,excl_inter);
            
            switch hup.order_method
                case 'size'
                    [s_inter_counts, ixs]=sort(inter_counts,'descend');
                case 'comb'
                    s_inter_counts=inter_counts;
                    ixs=1:length(inter_counts);
            end
            s_setix=setix(ixs);
            
            discard=s_inter_counts<hup.count_thr;
            ixs(discard)=[];
            s_setix(discard)=[];
            s_inter_counts(discard)=[];
            
            n_combs=length(ixs);
            XLIM=[0.5,n_combs+0.5];
            [~,max_ix]=max(s_inter_counts);
            
            %bar chart for intersection counts
            comb_ix=1:n_combs;
            set_ix=1:n_sets;
            h_bar=bar(ax_intersect,comb_ix,s_inter_counts);
            h_bar.FaceColor=[0.8,0.8,0.8];
            
            ax_intersect.XLim=XLIM;
%             ax_intersect.XTick=comb_ix;
            ax_intersect.XTick=[];
            ax_intersect.XTickLabel=[];
            ax_intersect.Box='off';
%             ax_intersect.YScale='log';

            %text labels for counts
            t_gap=0.02*max(s_inter_counts);
            ht_counts=text(ax_intersect,comb_ix,s_inter_counts+t_gap,string(s_inter_counts));
            for i=1:length(ht_counts)
                ht_counts(i).Rotation=90;
            end
            YLIM=ax_intersect.YLim;
            YLIM(2)=ht_counts(max_ix).Extent(2)+ht_counts(max_ix).Extent(4)+t_gap;
            ax_intersect.YLim=YLIM;
            
            %line plot for combinations
%             yback=1:n_sets;
%             xback=x*ones(size(yback));
%             hl_bg_dot(i)=line(ax_combs,xx,yback,'color',[0.8,0.8,0.8],...
%                 'linestyle','none','marker','.','markersize',30);
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
        end
        
        %getters/setters
    end
end