% Analyze protein usage
% Assumed that the models are already constructed by generateProtModels

organisms = {'Scer','Kmarx'};

% Load models
ecModel_Scer  = load('../../models/WTstudy_ecYeastGEM_protConstrained.mat');
ecModels{1}  = ecModel_Scer.ecModelP; 
ecModel_Kmarx = load('../../models/WTstudy_ecKmarx_protConstrained.mat');
ecModels{2} = ecModel_Kmarx.ecModelP;

% Get enzyme usages and reaction fluxes
clear sol absUsage capUsage UB
for i = 1:length(organisms)
   sol{i} = solveLP(ecModels{i},1);
   [absUsage{i}, capUsage{i}, UB{i}] = getEnzymeUsage(ecModels{i},sol{i}.x,true);
   printFluxes(ecModels{i},sol{i}.x,false,0, ...
       fullfile('../../','results','modelSimulation',['allFluxes_' organisms{i} '.txt']), ...
       '%rxnID\t%rxnName\t%eqn\t%flux\n');
end

% Prepare output
for i = 1:length(organisms)
    clear out
    out(:,1) = ecModels{i}.enzymes;
    out(:,2) = ecModels{i}.enzGenes;
    out(:,3) = ecModels{i}.enzNames;
    out(:,4) = strtrim(cellstr(num2str(capUsage{i},3)));
    out(:,5) = strtrim(cellstr(num2str(absUsage{i},3)));
    out(:,6) = strtrim(cellstr(num2str(UB{i},3)));
    % Write all usages to file
    header = {'protID','geneID','protName','capUsage','absUsage','UB'};
    out    = cell2table(out,'VariableNames',header);
    writetable(out,fullfile('..','Results','enzymeUsage',['enzymeUsage_' organisms{i} '.txt']),'Delimiter','\t')
end
    