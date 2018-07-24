%% Determine the APD Spread for all the models 
load('population_data.mat')

modelnames = {'Fox','Hund','Shannon','Devenyi','Livshitz','Grandi','TT04','TT06','Ohara','Heijman_ISO_0'};

for i = 1:length(modelnames)
    APDs = new_datatable.(modelnames{i}).APDs;   
    pert1 = prctile(APDs,90);
    pert2 = prctile(APDs,10);
    spread =(pert1 - pert2 )/ median(APDs);
    new_datatable.(modelnames{i}).APDSpread = spread;   
end

