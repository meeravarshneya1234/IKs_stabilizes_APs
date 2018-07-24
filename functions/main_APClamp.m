function datatable= main_APClamp(data)

%%
[p,c] = parameters(data.model_name,data.celltype);
start_point = find(data.volts>-75,1);
[~,ind_V] = max(data.volts);
end_point = find( data.times>data.times(ind_V) & data.volts<-75,1);

[APD_base,repol_index,tinit] = find_APD(data.times,data.volts);

odefcn = str2func(strcat('dydt_',char(data.model_name),'APclamp'));

% code to include a parfor loop 
tic 

parfor ind = 1:length(data.repol_change)
    
    statevar_i_original = data.statevar_i_original;

    State = zeros(0,length(data.statevar_i_original));
    Ti = zeros(0,1);
    Voltages = zeros(0,1);
    IKr = zeros(0,1);
    IKs = zeros(0,1);
    Area_Kr = zeros(0,1);
    Area_Ks = zeros(0,1);
    Ratio = zeros(0,1);
    
    new_APD = (APD_base*data.repol_change(ind))+tinit;
    new_end = find(new_APD>data.times,1,'last');
    
    % new V & t values between -75 mV back to -75 mV
    %l1 = (repol_index-start_point)+1;
    l1 = 350;
    t_interp1 = linspace(data.times(start_point),data.times(end_point),l1); % original data
    t_repol = linspace(data.times(start_point),data.times(new_end),l1);
    v_repol = interp1(data.times(start_point:end_point),data.volts(start_point:end_point),t_interp1);
    
    % new V & t during relaxation
    %l2 = (length(times)- (length(v_repol)+start_point))+1;
    l2 = 350;
    t_interp2 = linspace(data.times(end_point+1),data.times(end),l2);
    t_relax = linspace((data.times(new_end)+1),data.times(end),l2);
    v_relax = interp1(data.times(end_point+1:end),data.volts(end_point+1:end),t_interp2);
    
    t = [data.times(1:start_point-1); t_repol';t_relax'];
    V = [data.volts(1:start_point-1); v_repol';v_relax'];
    
    %set up segements for ode solver
    simints = length(t)-1 ;
    intervals = cell(simints,1) ;
    voltage = zeros(simints,1) ;
    for i=1:length(intervals)
        intervals{i} = [t(i),t(i+1)] ;
        voltage(i) = V(i) ;
    end
    
    index = 0;
    for ii = 1:data.numberofAPs
        time = 0 ;
        Voltage =voltage(1) ;
        statevar_i = statevar_i_original;
        statevars = statevar_i ;
        
        for i=1:simints
            [post,posstatevars] = ode15s(odefcn,intervals{i},statevar_i,[],voltage(i),p,c) ;
            time = [time;post(2:end)] ;
            statevars = [statevars;posstatevars(2:end,:)] ;
            statevar_i = posstatevars(end,:) ;
            Voltage = [Voltage;voltage(i)*ones(length(post)-1,1)] ;
        end
        
        true_vals=ismember(time,t);
        true_statevars = statevars(true_vals,:);
        outputcell = num2cell(true_statevars,1) ;
        statevar_i_original = cellfun(@(v) v(end), outputcell); % inital conditions for next AP
        
        disp(['AP # ' num2str(ii)])
        
        if ii > (data.numberofAPs - data.numberkeep) 
            [~,Vind] = ICs(data.model_name,1,1000);
            index = index + 1;
            [nRows,~] = size(true_statevars);
            State(end+1:end+nRows,:)  = true_statevars;
            Ti(end+1:end+nRows,1) = time(true_vals)+data.PCL*(index-1);
            Voltages(end+1:end+nRows,1) = Voltage(true_vals);
            sv_matrix = [true_statevars(:,1:(Vind-1)) Voltage(true_vals) true_statevars(:,(Vind:end))];
            [~,iKr,iKs]=currents(data.model_name,num2cell(sv_matrix,1),p,c);
            IKr(end+1:end+nRows,1)= iKr;
            IKs(end+1:end+nRows,1)= iKs;
            Area_Kr(end+1:end+index,1) = trapz(time(true_vals),iKr);
            Area_Ks(end+1:end+index,1) = trapz(time(true_vals),iKs);
            Ratio(end+1:end+index,1) =  trapz(time(true_vals),iKs)/(trapz(time(true_vals),iKr) + trapz(time(true_vals),iKs));
        end        
    end       
    datatable(ind).statevars = State;
    datatable(ind).Voltage = Voltages;
    datatable(ind).time = Ti;
    datatable(ind).IKs = IKs;
    datatable(ind).IKr = IKr;
    datatable(ind).Area_Kr = Area_Kr;
    datatable(ind).Area_Ks = Area_Ks;
    datatable(ind).Ratio =  Ratio;
    
    disp(['APD ' num2str(data.repol_change(ind))])
end
toc

