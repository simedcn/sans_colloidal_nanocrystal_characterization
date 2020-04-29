function cornerplot(x,x_median,x_mle,left_CR,right_CR,names,plotvec)
    % Generate a 'corner plot'

    % Make the corner plot figure using the subplot command in MATLAB looping
    % over each plot element
    figure;
    ssize = get(groot,'Screensize');
    screen_h = ssize(4);
    screen_w = ssize(3);
    set(gcf,'position',[500,65,900,900])
    Nmap = 16;
    map = [linspace(1,0,Nmap)',linspace(1,0,Nmap)',ones(Nmap,1)]; % Colormap for shading the sample data

    % Number of parameters
    p = size(x,1);
    Nplot = sum(plotvec==1);
    plotinds = find(plotvec==1);
    haxes = zeros(Nplot^2,1);
    
    limmat = zeros(Nplot,2);
    tickcell = cell(Nplot,1);
    
    sp = 0.005;
    width = 0.9/Nplot - sp;
    height = 0.9/Nplot - sp;
    % Loop through each individual parameter to get plotting limits for 2D
    % marginal distributions.
    for i = 1:Nplot
        
        plotind = plotinds(i);

        % Make the (i,i)th subplot on the diagonal which depicts the PDF of x_i
        k = Nplot*(i-1) + i;
        haxes(k) = subplot(Nplot,Nplot,k);
        [ f, xi ] = ksdensity(x(plotind,:));
        ytop = max(f)/0.8;
        plot( xi, f, '-b' );
        hold on
        errorbar(x_median(plotind),ytop/2,x_median(plotind)-left_CR(plotind),right_CR(plotind)-x_median(plotind),'horizontal','-k','LineWidth',1.5)
        plot([x_median(plotind) x_median(plotind)],[0 ytop],'-k')
        plot([x_mle(plotind) x_mle(plotind)],[0 ytop],'--r')
        xlimits = [min(xi), max(xi)];
        ylimits = [0 ytop];
        limmat(i,:) = xlimits;
        set(haxes(k),'XLim',xlimits,'YLim',ylimits)
        tickcell{i} = get(gca,'XTick');
        
        pbaspect([1 1 1])
        set(gca,'TickLength',[0.06 0.06])
        set(gca,'YTickLabel',[])
        set(gca,'YTick',[])
        if i < Nplot
            set(gca,'XTickLabel',[])
        else
            xlabel(names{plotind},'FontWeight','normal','FontSize',10)
            set(gca,'FontWeight','normal','FontSize',10)
            xtickangle(45)
        end
        
    end
    
    for i = 1:Nplot
        
        plotind_i = plotinds(i);
        
        for j = i+1:Nplot
            
            plotind_j = plotinds(j);

            % Make the (i,j)th subplot which contains samples of x_i vs x_j and
            % a contour plot of the joint PDF of x_i and x_j
            k = Nplot * ( j - 1 ) + i;
            haxes(k) = subplot( Nplot, Nplot, k );

            % Make a histogram of the samples in the x_i, x_j projection as a
            % model of the marginalized distribution
            xdens = [x(plotind_i,:)',x(plotind_j,:)'];
            [f,xi] = ksdensity(xdens);
            Xplot = reshape(xi(:,1),30,30);
            Yplot = reshape(xi(:,2),30,30);
            Zplot = reshape(f,30,30);
            
            contourf(Xplot,Yplot,Zplot,'-b','LineWidth',0.25)
            colormap(map)
            hold on
            set(haxes(k),'Xlim',limmat(i,:),'YLim',limmat(j,:))
            
            pbaspect([1 1 1])
            set(gca,'TickLength',[0.06 0.06])
            if j < Nplot
                set(gca,'XTickLabel',[])
            else
                xlabel(names{plotind_i},'FontWeight','normal','FontSize',10)
                set(gca,'FontWeight','normal','FontSize',10)
                xtickangle(45)
            end
            if i > 1
                set(gca,'YTickLabel',[])
            else
                ylabel(names{plotind_j},'FontWeight','normal','FontSize',10)
                set(gca,'FontWeight','normal','FontSize',10)
                ytickangle(45)
            end

        end;

    end;
    
    for i = 1:Nplot
        k = Nplot*(i-1) + i;
        set(haxes(k),'position',[0.1-sp+(i-1)*(sp+width),0.1-sp+(Nplot-i)*(sp+height),width,height])
        set(haxes(k),'XLim',limmat(i,:),'XTick',tickcell{i})
    end
    for i = 1:Nplot
        for j = i+1:Nplot
            k = Nplot * ( j - 1 ) + i;
            set(haxes(k),'position',[0.1-sp + (i-1)*(sp+width),0.1-sp + (Nplot-j)*(sp+height),width,height])
            set(haxes(k),'XLim',limmat(i,:),'YLim',limmat(j,:),'XTick',tickcell{i},'YTick',tickcell{j})
        end

    end

end