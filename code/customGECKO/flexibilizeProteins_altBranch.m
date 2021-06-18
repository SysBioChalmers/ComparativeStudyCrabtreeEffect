function [model,enzUsages,modifications] = flexibilizeProteins_altBranch(model,gRate,c_UptakeExp,c_source)
% flexibilizeProteins
%   Function that takes an ecModel with proteomic constraints and, if it is
%   overconstrained with respect to the provided experimental growth rate,
%   iterates finding the top growth-limiting enzyme or enzyme complex. The
%   exchange rate upper bound for each of the identified enzyme or subunits 
%   is then set to infinity and after all iterations they are set equal to 
%   the usage value provided by an overall enzyme usage minimization 
%   simulation (subject to the provided growth rate and nutrient uptake 
%   constraints).
%
%   model           ecModel with proteomic constraints (individual enzyme
%                   levels)
%   gRate           Minimum growth rate the model should grow at [1/h]. For
%                   finding the growth reaction, GECKO will choose the
%                   non-zero coeff in the objective function.
%   c_UptakeExp     (Opt) Experimentally measured glucose uptake rate 
%                   [mmol/gDw h]
%	c_source        (Opt) The name of the exchange reaction that supplies
%                   the model with carbon.
%
%   model           ecModel with calibrated enzyme usage upper bounds
%   enzUsages       Calculated enzyme usages after final calibration 
%                   (enzyme_i demand/enzyme_i upper bound)
%   modifications   Table with all the modified values 
%                   (Protein ID/old value/Flexibilized value)
%
%   Usage: [model,enzUsages,modifications] = flexibilizeProteins(model,gRate,c_UptakeExp,c_source)
%
%   Benjamin J. Sanchez     2018-12-11
%   Ivan Domenzain          2020-02-24
%
current = pwd;
flexFactor    = 10;
flexProts     = {};
enzUsages     = [];
modifications = {};
%constrain glucose uptake if an experimental measurement is provided
if nargin > 2
    glucUptkIndx           = strcmp(model.rxnNames,c_source);
    model.ub(glucUptkIndx) = c_UptakeExp;
end
% get measured protein exchange rxns indexes
cd ../../../GECKO/geckomat/limit_proteins
measuredIndxs = getMeasuredProtsIndexes(model);
if ~isempty(measuredIndxs)
    abundances = model.ub(measuredIndxs);
    objIndex   = find(model.c==1);
    sol        = solveLP(model,1);
    growth     = sol.x(objIndex);
    difference = (gRate-growth)/gRate;
	tempModel  = model;
    % iterate while growth is underpredicted
    while growth<gRate && difference>(0.005)
        % Get a lower bound for the glucose uptake that matches the
        % experimentally measured value
        testModel = setParam(tempModel,'lb',tempModel.rxns(glucUptkIndx),0.95*5.5576);
        testSol   = solveLP(testModel);
        if ~isempty(testSol.x)
            tempModel = testModel;
        end
        [limIndex,flag] = findLimitingUBs(tempModel,measuredIndxs,flexFactor,1);
        if ~flag
            [limIndex,~] = findLimitingUBs(tempModel,measuredIndxs,flexFactor,2);
        end
        if ~isempty(limIndex)
            flexProts = [flexProts; tempModel.rxns(limIndex)];
            if length(tempModel.rxns(limIndex)) > 1
                prot = cell(size(tempModel.rxns(limIndex)))';
                proteins = tempModel.rxns(limIndex);
                for i = 1:length(proteins)
                    p = proteins(i);
                    p = split(p,'_');
                    p = p{2};
                    prot{i} = p;
                end
            else
                prot = split(tempModel.rxns{limIndex},'_');
                prot = prot{2};
            end
            % If subunit is part of complex, flexibilize all subunits
            complexIII_subunits = {'W0T4N4','W0TGH3','W0T7G1','W0TFW4','W0TG79',...
                'W0T4A3','W0T342','W0T8X6','W0TBV3'};
            complexIV_subunits = {'W0TI46','W0TIU0','W0TH97','W0TGI5','W0T9H3',...
                'W0TA28','W0T557','W0T9K0','W0TCG0','W0THT4','W0TA38'};
            ATPsynthase_subunits = {'W0TG03','W0TFY6','W0TDL5','W0T794','W0TH05','W0TA81',...
                'W0TKQ2','W0T470','A0A1L7LMC3','W0TCT8','W0TC70','W0T9V6',...
                'W0T6Z4','W0TFB3','W0TAU7','W0TEE5','W0TFI9','W0TI50'};
            PDH_subunits = {'W0T5D1','W0T7E5','W0THU7','W0TFW6','W0TER0'};
            SDH_subunits = {'W0TA03','W0TBF6','W0TGF8','W0TKL8'};
            KGD_subunits = {'W0TDY2','W0TBX1'};
            if sum(ismember(ATPsynthase_subunits,prot)) > 0
                ATPsynthase_subunits = cellfun(@(c)['prot_' c '_exchange'],ATPsynthase_subunits,'uni',false)';
                limIndex = find(ismember(model.rxns,ATPsynthase_subunits));
                flexProts = [flexProts;ATPsynthase_subunits];
                flexProts = unique(flexProts);
            elseif sum(ismember(complexIII_subunits,prot)) > 0
                complexIII_subunits = cellfun(@(c)['prot_' c '_exchange'],complexIII_subunits,'uni',false)';
                limIndex = find(ismember(model.rxns,complexIII_subunits));
                flexProts = [flexProts;complexIII_subunits];
                flexProts = unique(flexProts);
            elseif sum(ismember(complexIV_subunits,prot)) > 0
                complexIV_subunits = cellfun(@(c)['prot_' c '_exchange'],complexIV_subunits,'uni',false)';
                limIndex = find(ismember(model.rxns,complexIV_subunits));
                flexProts = [flexProts;complexIV_subunits];
                flexProts = unique(flexProts);
            elseif sum(ismember(PDH_subunits,prot)) > 0
                PDH_subunits = cellfun(@(c)['prot_' c '_exchange'],PDH_subunits,'uni',false)';
                limIndex = find(ismember(model.rxns,PDH_subunits));
                flexProts = [flexProts;PDH_subunits];
                flexProts = unique(flexProts);
            elseif sum(ismember(PDH_subunits,prot)) > 0
                SDH_subunits = cellfun(@(c)['prot_' c '_exchange'],SDH_subunits,'uni',false)';
                limIndex = find(ismember(model.rxns,SDH_subunits));
                flexProts = [flexProts;SDH_subunits];
                flexProts = unique(flexProts);
            elseif sum(ismember(KGD_subunits,prot)) > 0
                KGD_subunits = cellfun(@(c)['prot_' c '_exchange'],KGD_subunits,'uni',false)';
                limIndex = find(ismember(model.rxns,KGD_subunits));
                flexProts = [flexProts;KGD_subunits];
                flexProts = unique(flexProts);
            end
            %Flexibilize the top growth limiting protein on the original ecModel
            tempModel.ub(limIndex) = 1000;
            sol = solveLP(tempModel,1);
            if ~isempty(sol.x)
            	growth = sol.x(objIndex);
                for j=1:length(limIndex)
                    indx = limIndex(j);
                    disp(['Modified ub for: ' tempModel.rxnNames{indx} ' gRate: ' num2str(growth)])
                end
            else
               %In case that the resulting model is a non-functional one
               %then proceed with a suboptimal growth rate (this makes the 
               %while loop to break)
               warning(['Unfeasible flexibilization of ' model.rxnNames{indx} ' UB'])
               gRate = growth;
            end
        else
            %In case that no limiting enzymes have been found then proceed 
            %with a suboptimal growth rate (this makes the while loop to break)
            warning('No limiting individual enzyme was found. All fully saturated enzymes are flexibilized by 1%')
            ratios = sol.x(measuredIndxs)./model.ub(measuredIndxs);
            sat_enz = find(ratios>0.9999);
            idxs = measuredIndxs(sat_enz);
            tempModel.ub(idxs) = tempModel.ub(idxs)*1.01;
        end
    end
    [model,enzUsages]  = getNewBounds(tempModel,gRate,measuredIndxs,flexProts,objIndex,abundances);
    modifiedAbundances = model.ub(measuredIndxs);
    exchangedProteins  = model.rxnNames(measuredIndxs);
    modifications      = getDifferences(abundances,modifiedAbundances,exchangedProteins,model);
    %Update model.concs field, taking flexibilized protein mass into account
    model.concs = nan(size(model.enzymes));
    for i = 1:length(model.enzymes)
        rxn_name = ['prot_' model.enzymes{i} '_exchange'];
        index    = find(strcmpi(rxn_name,model.rxns));
        if ~isempty(index)
            model.concs(i) = model.ub(index)*model.MWs(i); %g/gDW
        end
    end
end
cd(current)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function diffTable = getDifferences(originalBounds,newBounds,exchangedProteins,model)
protein_IDs     = {};
previous_values = [];
modified_values = [];
Mweights        = [];
for i=1:length(originalBounds)
    if newBounds(i)>originalBounds(i)
        proteinID   = exchangedProteins{i};
        proteinID   = proteinID(1:((strfind(proteinID,'_exchange'))-1));
        protein_IDs = [protein_IDs; proteinID];
        %Get protein MW
        index = find(strcmpi(model.enzymes,strrep(proteinID,'prot_','')));
        if ~isempty(index)
            Mweights = [Mweights; model.MWs(index)];
        end
        %Get previous and modified usage values [mmol/gDw h]
        previous_values = [previous_values; originalBounds(i)];
        if newBounds(i)~=Inf
            modified_values = [modified_values; newBounds(i)];
        else
            modified_values = [modified_values; originalBounds(i)];
        end
    end
end
%Calculate flexibilized mass for each protein [mmol/gDw h]*[g/mmol]>g/gDwh
flex_Mass = (modified_values-previous_values).*Mweights;
diffTable = table(protein_IDs,previous_values,modified_values,flex_Mass);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function measuredIndxs = getMeasuredProtsIndexes(model)
measuredIndxs  = intersect(find(contains(model.rxnNames,'prot_')),find(contains(model.rxnNames,'_exchange')));
measuredIndxs  = measuredIndxs(1:end-1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [model,usagesTable] = getNewBounds(model,gRate,protIndxs,flexProts,gPos,abundances)
%Now that the model is growing at least at the specified dilution rate
%lets fix the growth rate and minimize enzymes usage
objectiveVector = model.c;
model.lb(gPos)  = 0.9999*gRate;
model.c(:)      = 0;
protNames       = model.rxnNames(protIndxs);
pool_Index      = contains(model.rxnNames,'prot_pool_');
%Forces flux to split over measured enzymes
model.c(pool_Index) = -1;
optSolution         = solveLP(model,1);
optSolution         = optSolution.x;
enzUsages           = zeros(length(protIndxs),1);
for i=1:length(protIndxs)
    index = protIndxs(i);
    name  = model.rxns(index);
    %If protein was flexibilized set its upper bound to the simulated
    %concentration
    if ismember(name,flexProts)
        if optSolution(index)>abundances(i)
            model.ub(index) = optSolution(index);
        else
            model.ub(index) = abundances(i);
        end
    end
    enzUsages(i) = optSolution(index)/model.ub(index);
end
model.c     = objectiveVector;
usagesTable = table(protNames,enzUsages,'VariableNames',{'prot_IDs' 'usage'});
end