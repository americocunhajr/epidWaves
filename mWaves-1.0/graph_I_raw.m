% -----------------------------------------------------------------
%  graph_I_raw.m
% ----------------------------------------------------------------- 
%  programmer: Americo Cunha Jr (UERJ)
%               
%  last update: Jan 17, 2021
% -----------------------------------------------------------------
function fig = graph_I_raw(time,yraw,yMA,graphobj)
    
    % check number of arguments
    if nargin < 4
        error('Too few inputs.')
    elseif nargin > 4
        error('Too many inputs.')
    end

    % check arguments
    if length(time) ~= length(yraw)
        error('time and yraw vectors must be same length')
    end
    
    if length(time) ~= length(yMA)
        error('time and yMA vectors must be same length')
    end
    
    % convert to row vectors so fliplr can work
    if find( size(time) == max(size(time)) ) < 2
        time=time';
    end
    
    if find( size(yraw) == max(size(yraw)) ) < 2
        yraw=yraw';
    end
    
    if find( size(yMA) == max(size(yMA)) ) < 2
        yMA = yMA';
    end
    
    fig = figure('Name',graphobj.gname,'NumberTitle','off');
    
    fh1 = plot(time,yraw,'om','LineWidth',1);
    hold all
    fh2 = plot(time,yMA  ,'.-g','LineWidth',3);
    set(gcf,'color','white');
    set(gca,'position',[0.2 0.2 0.7 0.7]);
    set(gca,'Box','on');
    set(gca,'TickDir','out','TickLength',[.02 .02]);
    set(gca,'XMinorTick','off','YMinorTick','off');
    set(gca,'XGrid','off','YGrid','off');
    set(gca,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);
    set(gca,'FontName','Helvetica');
    set(gca,'FontSize',18);
    datetick('x',28,'keeplimits');
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
