% -----------------------------------------------------------------
%  graph_I_vs_C_raw.m
% ----------------------------------------------------------------- 
%  programmer: Americo Cunha Jr (UERJ)
%               
%  last update: Jan 17, 2021
% -----------------------------------------------------------------
function fig = graph_I_vs_C_raw(C,I,CMA,IMA,graphobj)
    
    % check number of arguments
    if nargin < 5
        error('Too few inputs.')
    elseif nargin > 5
        error('Too many inputs.')
    end

    % check arguments
    if length(C) ~= length(I)
        error('C and I vectors must be same length')
    end
    
    if length(CMA) ~= length(IMA)
        error('CMA and IMA vectors must be same length')
    end
    
    % convert to row vectors so fliplr can work
    if find( size(C) == max(size(C)) ) < 2
        C=C';
    end
    
    if find( size(I) == max(size(I)) ) < 2
        I=I';
    end
    
    if find( size(CMA) == max(size(CMA)) ) < 2
        CMA = CMA';
    end
    
    if find( size(IMA) == max(size(IMA)) ) < 2
        IMA = IMA';
    end
    
    fig = figure('Name',graphobj.gname,'NumberTitle','off');
    
    fh1 = loglog(C,I,'om','LineWidth',1);
    hold all
    fh2 = loglog(CMA,IMA  ,'.-g','LineWidth',3);
    set(gcf,'color','white');
    set(gca,'position',[0.2 0.2 0.7 0.7]);
    set(gca,'Box','on');
    set(gca,'TickDir','out','TickLength',[.02 .02]);
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'XGrid','off','YGrid','off');
    set(gca,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);
    set(gca,'FontName','Helvetica');
    set(gca,'FontSize',18);
    %set(gca,'XTick',xmin:xmax);
    %set(gca,'YTick',ymin:ymax);
    %axis([xmin xmax ymin ymax]);
    
    set(fh1,'DisplayName',graphobj.leg1);
    set(fh2,'DisplayName',graphobj.leg2);
    leg = [fh1; fh2];
    leg = legend(leg,'Location','Best');
    set(leg,'FontSize',12);

    if ( strcmp(graphobj.xmin,'auto') || strcmp(graphobj.xmax,'auto') )
        xlim('auto');
    else
        xlim([graphobj.xmin graphobj.xmax]);
    end
    
    if ( strcmp(graphobj.ymin,'auto') || strcmp(graphobj.ymax,'auto') )
        ylim('auto');
    else
        ylim([graphobj.ymin graphobj.ymax]);
    end
    
    labX = xlabel(graphobj.xlab,'FontSize',16,'FontName','Helvetica');
    labY = ylabel(graphobj.ylab,'FontSize',16,'FontName','Helvetica');
    
    hold off
    
	title(graphobj.gtitle,'FontSize',20,'FontName','Helvetica');
    
    if ( strcmp(graphobj.flag,'eps') )
        saveas(gcf,graphobj.gname,'epsc2');
        graphobj.gname = [graphobj.gname, '.eps'];
    end

end
% -----------------------------------------------------------------
