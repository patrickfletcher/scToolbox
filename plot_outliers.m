function plot_outliers(QCdata,xname,yname,tf,ctnames,panels,gap,margh,margw)
%TODO: tf+ctnames vs cellID
%TODO: autopanels -> new tiledlayout?

msz=3;

amb=sum(tf,1)>1;
% plotAmb=false;

highlightOut=false;
if ~isempty(QCdata)
    highlightOut=true;
    Xstruct=QCdata.(xname);
    X=Xstruct.vals;
    xlo=Xstruct.lowthr; %xlo
    xhi=Xstruct.hithr; %xhi
%     xout=Xstruct.outsub; %x-grayed
    xout=Xstruct.outliers; %x-grayed

    Ystruct=QCdata.(yname);
    Y=Ystruct.vals;
    ylo=Ystruct.lowthr; %ylo
    yhi=Ystruct.hithr; %yhi
%     yout=Ystruct.outsub; %y-grayed
    yout=Ystruct.outliers; %y-grayed
end

isout=xout|yout;

Xrange=[min(X),max(X)]; 
XLIM=Xrange+.05*[-1,1].*Xrange; 
Yrange=[min(Y),max(Y)]; 
YLIM=Yrange+.05*[-1,1].*Yrange; 
XLIM(XLIM<0)=0;
YLIM(YLIM<0)=0;
if contains(xname,'frac')
    XLIM(XLIM>1)=1;
end
if contains(yname,'frac')
    YLIM(YLIM>1)=1;
end

if Xstruct.params.logval
    xlo=10.^(xlo);
    xhi=10.^(xhi);
end
if Ystruct.params.logval
    ylo=10.^(ylo);
    yhi=10.^(yhi);
end

nr=panels(1);
nc=panels(2);
for i=1:size(tf,1)
    ax=tight_subplot(nr,nc,i,gap,margh,margw);
    thissub=tf(i,:);
    x=X(thissub);
    y=Y(thissub);
    
    thisamb=amb(thissub);
    thisout=isout(thissub);
    
    xamb=x(thisamb);
    yamb=y(thisamb);
%     if ~plotAmb
%         xna=x(~thisamb);
%         yna=y(~thisamb);
%         thisxout=thisxout(~thisamb);
%         thisyout=thisyout(~thisamb);
%     end
    
    if ~isempty(x) && ~isempty(y)
        
%         if plotAmb
%             line(xamb, yamb,'color','r','linestyle','none','marker','o','markersize',msz)
            line(xamb, yamb,'color','r','linestyle','none','marker','x','markersize',msz)
%         end
        line(x,y,'color','k','linestyle','none','marker','.','markersize',msz);
        
        if highlightOut
%             line(x(thisxout), y(thisxout),'color','g','linestyle','none','marker','x','markersize',msz)
%             line(x(thisyout), y(thisyout),'color','b','linestyle','none','marker','+','markersize',msz)
            line(x(thisout), y(thisout),'color',0.5*[1,1,1],'linestyle','none','marker','.','markersize',msz)
        end
        axis tight
        
        if ismember(xname,{'n_genes','n_umi'})
            ax.XScale='log';
        end
        if ismember(yname,{'n_genes','n_umi'})
            ax.YScale='log';
        end

%         XLIM=xlim();
%         YLIM=ylim();
        if highlightOut
            if xlo(i)>min(Xrange)
                line(xlo(i)*[1,1], YLIM, 'color',[0,0,0.7]);
            end
            if xhi(i)<max(Xrange)
                line(xhi(i)*[1,1], YLIM, 'color',[0,0,0.7]);
            end
            if ylo(i)>min(Yrange)
                line(XLIM, ylo(i)*[1,1],'color',[0.7,0,0]); 
            end
            if yhi(i)<max(Yrange)
                line(XLIM, yhi(i)*[1,1],'color',[0.7,0,0]); 
            end
        end

    end
    
    
    title(ctnames(i), 'FontSize',8,'Interpreter','none');
    
    ax.FontSize=8;
    box on
    if i==size(tf,1)
        xlabel(xname,'Interpreter','none')
    end
    if i==1
        ylabel(yname,'Interpreter','none')
    end
end

% old stuff (?)

% density_resolution=64;
        
%         xx=x;
%         yy=y;
%         if ismember(xname,{'gpc','cpc'})
%             xx=log10(xx);
%         end
%         if ismember(yname,{'gpc','cpc'})
%             yy=log10(yy);
%         end
%         plotDensity([xx,yy],density_resolution);
%         line(xx(thisamb), yy(thisamb),'color','r','linestyle','none','marker','x','markersize',msz)
%         line(xx(thisyout), yy(thisyout),'color',[0.5,0.5,0.5],'linestyle','none','marker','.','markersize',msz)
%         line(xx(thisxout), yy(thisxout),'color',[0.5,0.5,0.5],'linestyle','none','marker','.','markersize',msz)
