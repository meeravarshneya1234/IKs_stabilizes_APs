function plotting(X,settings) 

% Colors Selection
if length(settings.Ca_scale)>1
    colors = hsv(length(settings.Ca_scale));
else
    colors = 'k';
end

for ii = 1:length(settings.Ks_scale) 

    str = settings.model_name;
    figure
    set(gcf,'Position',[25 25 500 900])
    for i = 1:length(settings.Ca_scale)
        subplot(5,1,1)
        p1(i) = plot(X.times{ii,i},X.V{ii,i},'color',colors(i,:),'linewidth',2);
        hold on
        xlabel('time (ms)')
        ylabel('V (mv)')
        
        subplot(5,1,2)
        p2 = plot(X.times{ii,i},X.IKs{ii,i},'color',colors(i,:),'linewidth',2);
        hold on
        xlabel('time (ms)')
        ylabel('IKs (A/F)')
        
        subplot(5,1,3)
        p3 = plot(X.times{ii,i},X.IKr{ii,i},'color',colors(i,:),'linewidth',2);
        hold on
        xlabel('time (ms)')
        ylabel('IKr (A/F)')
        
        subplot(5,1,4)
        p4 = plot(X.times{ii,i},X.ICaL{ii,i},'color',colors(i,:),'linewidth',2);
        hold on
        xlabel('time (ms)')
        ylabel('ICaL (A/F)')
        
        figurestrings{i} = ['P_C_a* ' num2str(settings.Ca_scale(i))]; % for legend
    end
    sh = subplot(5,1,5);
    p=get(sh,'position');
    lg_handles=p1;
    lh=legend(lg_handles,figurestrings);
    set(lh,'position',p);
    axis(sh,'off');
    
    %titlestring = [settings.model_name ' : G_K_s* ' num2str(settings.Ks_scale(ii))];
    titlestring = [str ' : G_K_s* ' num2str(settings.Ks_scale(ii))];
    mtit(titlestring,...
	     'fontsize',14);
    
end