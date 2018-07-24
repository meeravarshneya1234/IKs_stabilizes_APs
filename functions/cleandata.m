function [APfails,EADs] = cleandata(APDs,times,volts,t_cutoff,flag)

if ~exist('t_cutoff','var')
    t_cutoff = 5;
end

if ~exist('flag','var')
    flag =0;
end 

APfails = zeros(1,length(times));
EADs = zeros(1,length(times));

for i = 1:length(times)
    
    t = times{i};
    V = volts{i};
    
    if isnan(APDs(i)) %APD is NaN - cell fails
       % plot(t,V)
        APfails(i) = 1;
        
    elseif V(1) > -70 %AP does not begin at rest - cell fails
       % plot(t,V)        
        APfails(i) = 1;
        
    elseif max(V) < 0 % AP does not reach full potential - cell fails
       % plot(t,V)        
        APfails(i) = 1;
        
    else
        EADs(i) = findEAD(t,V,t_cutoff,flag);
%         slope_threshold = 0.1; %mV/ms, rate of change considered to be a significant depolarization or repolarization
%         t_cutoff = 10 ; %ms, the time between depolarization periods above which it is considered two separate depolarizations
%         V_slope = diff(V)./diff(t);
%         dices_above = find(V_slope > slope_threshold);
%         time_between_depols = diff(t(dices_above));
%         if flag==1
%             if sum(time_between_depols > t_cutoff) > 1
%                 %plot(t,V)
%                 APfails(i) = 1;
%             end
%         elseif any(time_between_depols > t_cutoff)
%            % plot(t,V)
%             APfails(i) = 1;
%         end
        
    end 

end

% cleanAPDs = APDs(~APfails);
% cleantimes = times(~APfails);
% cleanAPs = volts(~APfails);

%% Old Code 
% function [cleanAPDs,cleantimes,cleanAPs] = cleandata(APDs,times,actionpotentials)
% nans = isnan(APDs);
% cleanAPDs = APDs;
% cleanAPs = actionpotentials;
% cleantimes = times;
% 
% % remove APs whose APD = NaN
% cleanAPDs(nans) = [];
% cleanAPs(nans) = [];
% cleantimes(nans) = [];
% 
% % remove APs that do not repolarize 
% for i = 1:length(cleanAPs)
%     AP = cleanAPs{i};
%     if AP(1) > -20 
%         ind(i) = 1;
%     else 
%         ind(i) = 0;       
%     end 
% end
% cleantimes(ind==1) = [];
% cleanAPs(ind==1) = [];
% cleanAPDs(ind==1) = [];
% 
% % remove APs with EADs
% I =[];
% for i = 1:length(cleanAPs)
%     V = (cleanAPs{i});
%     t = (cleantimes{i});
%     [~,iV]=max(V);
%     indV = find(t>t(iV) & V<10);
%     indV = indV(1);
%     
%     deriv = diff(V(indV:end))./diff(t(indV:end));
%     negs = sign(deriv);
%     findzero = diff(negs);
%     if numel(find(~findzero==0))>1
%         I =[I;i];
%     end
% end
% cleantimes(I) = [];
% cleanAPs(I) = [];
% cleanAPDs(I) = [];
% 
% figure
% colors = repmat('rbcmkg',1,1000);
% for i = 1:length(cleanAPs)
%     plot(cleantimes{i},cleanAPs{i},'color',colors(i))
%     hold on   
% end
% 
% for i = 1:length(cleanAPs)
%     slope_threshold = 0.5; %mV/ms, rate of change considered to be a significant depolarization or repolarization
%     t_cutoff = 10 ; %ms, the time between depolarization periods above which it is considered two separate depolarizations
%     V_slope = diff(V)./diff(t);
%     dices_above = find(V_slope > slope_threshold);
%     dices_below = find(V_slope < slope_threshold*-1);
%     time_between_depols = diff(t(dices_above));
%     if any(time_between_depols > t_cutoff)
%         doesthiscellfail(i) = 1;
%     end
% end
% % less_100 = find(cleanAPDs<100);
% % cleantimes(less_100) = [];
% % cleanAPs(less_100) = [];
% % cleanAPDs(less_100) = [];
% % 
% % A =[];
% % for i = 1:length(cleanAPs)
% %     V = (cleanAPs{i});
% %     t = (cleantimes{i});
% %     [pks,locs] = findpeaks(V);
% % %     [~,iV]=max(V);
% % %     indV = find(t>t(iV) & V<1);
% % %     indV = indV(1);
% % %     
% % %     deriv = diff(V(indV:end))./diff(t(indV:end));
% % %     negs = sign(deriv);
% % %     findzero = diff(negs);
% %     if numel(pks)>2
% %         %find(~findzero==0)
% %         A =[A;i];
% % %         figure(i)
% % %         plot(t,V)
% %     end
% % end
% % cleantimes(A) = [];
% % cleanAPs(A) = [];
% % cleanAPDs(A) = [];
% % 
% % 
% % %%
% % % find_nan = find(isnan(APDs_new));
% % % actionpotentials_new(find_nan)= [];
% % % times_new(find_nan) = [];
% % % APDs_new(find_nan) = [];
% % % 
% % % APDs = [APDs APDs_new];
% % % actionpotentials = [actionpotentials;actionpotentials_new];
% % % times = [times;times_new];
% % % 
% % % 
% % % figure
% % % for i = 1:length(actionpotentials)
% % %     plot(times{i},actionpotentials{i},'linewidth',2)
% % %     hold on 
% % % end
% % % xlabel('Time (ms)')
% % % ylabel('Voltage (mV)')
% % % set(gca,'FontSize',12,'FontWeight','bold')
% % % 
% % % %%
% % % % AP end and start different
% % % for i = 1:length(actionpotentials)
% % %     V = (actionpotentials{i});
% % %     V1(i) = (V(1));
% % % end
% % % V1_find = find(V1>-70);
% % % 
% % % actionpotentials(V1_find)= [];
% % % times(V1_find) = [];
% % % APDs(V1_find) = [];
% % % 
% % % for i = 1:length(actionpotentials)
% % %     figure(2)
% % %     plot(times{i},actionpotentials{i})
% % %     hold on
% % %     
% % % end
% % % 
% % % I =[];
% % % for i = 1:length(actionpotentials)
% % %     V = (actionpotentials{i});
% % %     t = (times{i});
% % %     [~,iV]=max(V);
% % %     indV = find(t>t(iV) & V<1);
% % %     indV = indV(1);
% % %     
% % %     deriv = diff(V(indV:end))./diff(t(indV:end));
% % %     negs = sign(deriv);
% % %     findzero = diff(negs);
% % %     if numel(find(~findzero==0))>1
% % %         %find(~findzero==0)
% % %         I =[I;i];
% % % %         figure(i)
% % % %         plot(t,V)
% % %     end
% % % end
% % % 
% % % actionpotentials(I)= [];
% % % times(I) = [];
% % % APDs(I) = [];
% % % 
% % % for i = 1:length(actionpotentials)
% % %     figure(3)
% % %     plot(times{i},actionpotentials{i},'linewidth',2)
% % %     hold on 
% % % end
% % % 
% % % less_100 = find(APDs<100);
% % % actionpotentials(less_100)= [];
% % % times(less_100) = [];
% % % APDs(less_100) = [];
% % % 
% % % for i = 1:length(actionpotentials)
% % %     figure(4)
% % %     plot(times{i},actionpotentials{i},'linewidth',2)
% % %     hold on 
% % % end
% % % 
% % % find_nan = find(isnan(APDs));
% % % actionpotentials(find_nan)= [];
% % % times(find_nan) = [];
% % % APDs(find_nan) = [];
% % % 
% % % for i = 1:length(actionpotentials)
% % %     figure(6)
% % %     plot(times{i},actionpotentials{i},'linewidth',2)
% % %     hold on 
% % % end
% % % xlabel('Time (ms)')
% % % ylabel('Voltage (mV)')
