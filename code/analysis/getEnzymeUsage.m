function [absUsage,capUsage,UB,protID] = getEnzymeUsage(ecModel,fluxes,zero)
% getEnzymeUsage
%
% collects enzyme usages based on a provided flux distribution
% 1) Absolute usage: the specific enzyme usage in mmol/gDW/h, which can be
%    given for enzymes both with- and without measured abundances
% 2) Capacity usage: the ratio of the available enzyme that is used,
%    calculated by (absUsage/UB). If the enzyme was not constrained by
%    experimental data, the capacity usage will be 0
% 3) UB: the upper bound of each enzyme exchange reaction
% 4) protID: the protein identifiers for each enzyme
%
% Input:
% ecModel       enzyme-constrained model (1x1 struct)
% fluxes        vector of fluxes
% zero          logical whether enzymes with zero absolute usage should be
%               included (opt, default true)
%
% Output:
% capUsage      vector of enzyme capacity usages
% absUsage      vector of absolute enzyme usages
% UB            vector of enzyme exchange reaction upper bounds
% protID        string array of protein IDs
%
% Usage: [absUsage,capUsage,UB,protID] = getEnzymeUsage(ecModel,fluxes,zero)

if nargin < 3
    zero = true;
end

if isfield(ecModel,'enzymes')
    protInd    = find(contains(ecModel.rxnNames,ecModel.enzymes));
    matchProt  = regexprep(ecModel.rxnNames(protInd),'(draw_)?prot_','');
    matchProt  = regexprep(matchProt,'_exchange.*','');
    [~,b]      = ismember(ecModel.enzymes,matchProt);
    
    absUsage   = fluxes(protInd);
    UB         = ecModel.ub(protInd);
    capUsage   = absUsage./UB;
    
    absUsage   = absUsage(b);
    capUsage   = capUsage(b);
    UB         = UB(b);
    protID      = ecModel.enzymes;
else
    expression = '^(draw_prot_.*)|(prot_.*_exchange$)';
    protRxn    = regexp(ecModel.rxnNames,expression);
    protRxn    = find(~cellfun(@isempty,protRxn));
    expression = '(^draw_prot_)|(^prot_)|(_exchange$)';
    protID     = regexprrep(ecModel.rxnNames(protRxn),expression,'');
    
    absUsage   = fluxes(protRxn);
    UB         = ecModel.ub(protRxn);
    capUsage   = absUsage./UB;
end

if true(~zero)
    nonZero    = absUsage > 0;
    absUsage   = absUsage(nonZero);
    capUsage   = capUsage(nonZero);
    UB         = UB(nonZero);
    protID     = protID(nonZero);
end
    
end
    
    
    
    