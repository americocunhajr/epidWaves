% -----------------------------------------------------------------
%  graph_I_1w.m
% ----------------------------------------------------------------- 
%  programmer: Americo Cunha Jr (UERJ)
%               
%  last update: Jan 17, 2021
% -----------------------------------------------------------------
function fig = graph_I_1w(time,yraw,yMA,yupp,ylow,ypred,tau_ast,graphobj)
    
    % check number of arguments
    if nargin < 8
        error('Too few inputs.')
    elseif nargin > 8
        error('Too many inputs.')
    end

    % check arguments
    if length(time) ~= length(yraw)
        error('time and yraw vectors must be same length')
    end

    if length(time) ~= length(yMA)
        error('time and yMA vectors must be same length')
    end
    
    if length(ylow) ~= length(yupp)
        error('ylow and yupp vectors must be same length')
    end
    
    if length(time) ~= length(ylow)
        error('time and ylow vectors must be same length')
    end
    
    if length(time) ~= length(ypred)
        error('time and ypred vectors must be same length')
    end
    
    % convert to row vectors so fliplr can work
    if find( size(time) == max(size(time)) ) < 2
        time=time';
    end
    
    if find( size(yraw) == max(size(yraw)) ) < 2
        yraw=yraw';
    end
    
    if find( size(yMA) == max(size(yMA)) ) < 2
        yMA=yMA';
    end
    
    if find( size(ylow) == max(size(ylow)) ) < 2
        ylow = ylow';
    end
    
    if find( size(yupp) == max(size(yupp)) ) < 2
        yupp = yupp';
    end
    
    if find( size(ypred) == max(size(ypred)) ) < 2
        ypred = ypred';
    end
    
    % custom color
    MyGray = [0.8 0.8 0.8];
    
    fig = figure('Name',graphobj.gname,'NumberTitle','off');
    
    fh4 = fill([time fliplr(time)],[yupp fliplr(ylow)],MyGray);
    hold all
    fh1 = plot(time,yraw,'o-m','LineWidth',1);
    fh2 = plot(time,yMA  ,'.-g','LineWidth',3);
    fh3 = plot(time,ypred,' -b','LineWidth',3);
    fh5 = plot([tau_ast tau_ast],...
               [0 0.85*graphobj.tauy*graphobj.ymax],'--k','LineWidth',2);
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
    set(fh3,'DisplayName',graphobj.leg3);
    set(fh4,'DisplayName',graphobj.leg4);
    leg = [fh1; fh2; fh3; fh4];
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
    
    text(graphobj.taux,graphobj.tauy, datestr(tau_ast), ...
         'Units', 'normalized', ...
         'HorizontalAlignment',graphobj.tauh, ...
         'VerticalAlignment'  ,'top',...
         'FontSize', 16);
    
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
