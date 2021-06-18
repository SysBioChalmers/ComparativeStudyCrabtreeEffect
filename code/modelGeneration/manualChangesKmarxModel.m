function model = manualChangesKmarxModel(model)
%manualKcatChanges
%
%   Make changes to kcat values that are limiting at maximum growth rate
%
% Carl Malina. Last edited: 2021-03-16
%

current = pwd;
% The external and internal NADH dehydrogenase proteins are not included in
% model. 
% Correct the grRules in model
model.grRules{strcmp(model.rxns,'r_1192')} = 'KLMA_10695';
model.grRules{strcmp(model.rxns,'r_1193')} = 'KLMA_10279';

% load data from UniProt and KEGG
data = load('../../data/databases/ProtDatabase_Kmarx.mat');
swissprot = data.swissprot;
kegg      = data.kegg;

% Add proteins to reactions
proteins = {'W0T415';'W0T527'}; % NDE;NDI
rxnIDs = {'r_1192';'r_1193'};
rxnNames = model.rxnNames(ismember(model.rxns,rxnIDs));
cd ../../../GECKO/geckomat/change_model/
for i = 1:length(rxnIDs)
    rxnID   = rxnIDs{i};
    newID   = [rxnID 'No1'];
    newName = [rxnNames{i} ' (No1)'];
    kvalues = 500*3600; % s-1 -> h-1
    newMet  = ['prot_' proteins{i}];
    grRule  = model.grRules{strcmp(model.rxns,rxnID)};
    model   = addEnzymesToRxn(model,kvalues,rxnID,newMet,{newID,newName},grRule);
    % remove original reaction
    model = removeReactions(model,{rxnID});
end
% add additional information for proteins
for i = 1:length(proteins)
    model = addProtein(model,proteins{i},kegg,swissprot);
end

% r_0006 - Phosphopyruvate hydratase (ENO)
proteins = {'W0T7K9'}; % ENO
rxnID    = 'r_0006';
rxnName  = model.rxnNames{strcmp(model.rxns,rxnID)};
newID    = [rxnID 'No1'];
newName  = [rxnName ' (No1)'];
kvalues  = 71*3600; % s-1 -> h-1, use value from S. cerevisiae
newMet  = 'prot_W0T7K9';
grRule   = model.grRules{strcmp(model.rxns,rxnID)};
model    = addEnzymesToRxn(model,kvalues,rxnID,newMet,{newID,newName},grRule);
% remove original reaction
model = removeReactions(model,{rxnID});
% add additional information for proteins
for i = 1:length(proteins)
    model = addProtein(model,proteins{i},kegg,swissprot);
end

% r_1202 - Pyruvate decarboxylase (PDC1)
proteins = {'W0TEZ4'}; % ENO
rxnID    = 'r_1202';
rxnName  = model.rxnNames{strcmp(model.rxns,rxnID)};
newID    = [rxnID 'No1'];
newName  = [rxnName ' (No1)'];
kvalues  = 145*3600; % s-1 -> h-1, use value from S. cerevisiae PDC1 (PMID: 23423327)
newMet  = 'prot_W0TEZ4';
grRule   = model.grRules{strcmp(model.rxns,rxnID)};
model    = addEnzymesToRxn(model,kvalues,rxnID,newMet,{newID,newName},grRule);
% remove original reaction
model = removeReactions(model,{rxnID});

% r_0007 - Aldehyde dehydrogenase (NAD-dependent)
proteins = {'W0T8I5';'W0TCX8';'W0TBN1';'W0T8F8'}; % ALD2, ALD4, ALD5, ALD6
rxnID    = 'r_0007';
model.grRules{strcmp(model.rxns,rxnID)} = 'KLMA_20673 or KLMA_50012 or KLMA_40404 or KLMA_10742';
grR = model.grRules{strcmp(model.rxns,rxnID)};
rxnName  = model.rxnNames{strcmp(model.rxns,rxnID)};
% Add arm reaction
model  = addArmReaction(model,rxnID);
for i = 1:length(proteins)
    protein = proteins{i};
    newID   = [rxnID 'No' num2str(i)];
    newName = [rxnName ' (No' num2str(i) ')'];
    kvalues = 31.7*3600; % s-1 -> h-1, use value from S. cerevisiae
    newMet  = ['prot_' protein];
    grRule  = split(grR,' or ');
    grRule  = grRule{i};
    model   = addEnzymesToRxn(model,kvalues,rxnID,newMet,{newID,newName},grRule);
end
% remove original reaction
model = removeReactions(model,{rxnID});
% add additional information for proteins
for i = 1:length(proteins)
    model = addProtein(model,proteins{i},kegg,swissprot);
end

% r_0008 - Aldehyde dehydrogenase (NAD-dependent), mitochondrial
% Since cellular localization not well established, assume the same enzymes
% as cytosolic equivalent
proteins = {'W0T8I5';'W0TCX8';'W0TBN1';'W0T8F8'}; % ALD2, ALD4, ALD5, ALD6
rxnID    = 'r_0008';
model.grRules{strcmp(model.rxns,rxnID)} = 'KLMA_20673 or KLMA_50012 or KLMA_40404 or KLMA_10742';
grR = model.grRules{strcmp(model.rxns,rxnID)};
rxnName  = model.rxnNames{strcmp(model.rxns,rxnID)};
% Add arm reaction
model  = addArmReaction(model,rxnID);
for i = 1:length(proteins)
    protein = proteins{i};
    newID   = [rxnID 'No' num2str(i)];
    newName = [rxnName ' (No' num2str(i) ')'];
    kvalues = 31.7*3600; % s-1 -> h-1, use value from S. cerevisiae
    newMet  = ['prot_' protein];
    grRule  = split(grR,' or ');
    grRule  = grRule{i};
    model   = addEnzymesToRxn(model,kvalues,rxnID,newMet,{newID,newName},grRule);
end
% remove original reaction
model = removeReactions(model,{rxnID});

% r_0009 - Aldehyde dehydrogenase (NADP-dependent)
proteins = {'W0T8I5';'W0T8F8'}; % ALD2, ALD6
rxnID = 'r_0009';
model.grRules{strcmp(model.rxns,rxnID)} = 'KLMA_20673 or KLMA_10742';
grR = model.grRules{strcmp(model.rxns,rxnID)};
rxnName = model.rxnNames{strcmp(model.rxns,rxnID)};
% Add arm reaction
model  = addArmReaction(model,rxnID);
for i = 1:length(proteins)
    protein = proteins{i};
    newID   = [rxnID 'No' num2str(i)];
    newName = [rxnName ' (No' num2str(i) ')'];
    kvalues = 2.75*3600; % s-1 -> h-1, use value from S. cerevisiae calculated from specific activity PMID:23454351 
    newMet  = ['prot_' protein];
    grRule  = split(grR,' or ');
    grRule  = grRule{i};
    model   = addEnzymesToRxn(model,kvalues,rxnID,newMet,{newID,newName},grRule);
end
% remove original reaction
model = removeReactions(model,{rxnID});

% r_0010 - Aldehyde dehydrogenase (NADP-dependent), mitochondrial
% Since cellular localization not well established, assume the same enzymes
% as cytosolic equivalent
proteins = {'W0T8I5';'W0T8F8'}; % ALD2, ALD6
rxnID = 'r_0010';
model.grRules{strcmp(model.rxns,rxnID)} = 'KLMA_20673 or KLMA_10742';
grR = model.grRules{strcmp(model.rxns,rxnID)};
rxnName = model.rxnNames{strcmp(model.rxns,rxnID)};
% Add arm reaction
model  = addArmReaction(model,rxnID);
for i = 1:length(proteins)
    protein = proteins{i};
    newID   = [rxnID 'No' num2str(i)];
    newName = [rxnName ' (No' num2str(i) ')'];
    kvalues = 2.75*3600; % s-1 -> h-1, use value from S. cerevisiae calculated from specific activity PMID:23454351 
    newMet  = ['prot_' protein];
    grRule  = split(grR,' or ');
    grRule  = grRule{i};
    model   = addEnzymesToRxn(model,kvalues,rxnID,newMet,{newID,newName},grRule);
end
% remove original reaction
model = removeReactions(model,{rxnID});

% r_0067 - 6-phosphogluconolactonase
proteins = {'W0TF97';'W0T4I6'}; 
rxnID    = 'r_0067';
rxnName  = model.rxnNames{strcmp(model.rxns,rxnID)};
grR      = model.grRules{strcmp(model.rxns,rxnID)};
% Add arm reaction
model  = addArmReaction(model,rxnID);
for i = 1:length(proteins)
    protein = proteins{i};
    newID   = [rxnID 'No' num2str(i)];
    newName = [rxnName ' (No' num2str(i) ')'];
    kvalues = 27.6*3600; % s-1 -> h-1, specific activity from Bos taurus used (PMID: 6852020)
    newMet  = ['prot_' protein];
    grRule  = split(grR,' or ');
    grRule  = grRule{i};
    model   = addEnzymesToRxn(model,kvalues,rxnID,newMet,{newID,newName},grRule);
end
% remove original reaction
model = removeReactions(model,{rxnID});
% add additional information for proteins
for i = 1:length(proteins)
    model = addProtein(model,proteins{i},kegg,swissprot);
end

% r_0587 - NADP-specific glutamate dehydrogenase GDH3
% Remove reaction where GDH3 catalyses NADH-specific reaction
model = removeReactions(model,{'r_0586No2'});
rxnID    = 'r_0587';
% Change direction of reaction, since it's most likely synthesizing
% L-glutamate
rxnInd = strcmp(model.rxns,rxnID);
subsPos = find(model.S(:,rxnInd)<0);
prodPos = find(model.S(:,rxnInd)>0);
model.S(subsPos,rxnInd) = model.S(subsPos,rxnInd)*-1;
model.S(prodPos,rxnInd) = model.S(prodPos,rxnInd)*-1;
rxnName  = model.rxnNames{strcmp(model.rxns,rxnID)};
newID    = [rxnID 'No1'];
newName  = [rxnName ' (No1)'];
kvalues  = 14*3600; % s-1 -> h-1, use value from S. cerevisiae GDH3
newMet  = 'prot_W0TCB9';
grRule   = model.grRules{strcmp(model.rxns,rxnID)};
model    = addEnzymesToRxn(model,kvalues,rxnID,newMet,{newID,newName},grRule);
% remove original reaction
model = removeReactions(model,{rxnID});

% The directionality of succinyl-CoA ligase reaction (r_0042No1), favors
% ATP hydrolys and production of succinyl-CoA from succinate. Change
% direction to correct this. Move substrates to products side and products
% to substrate side
rxnInd = strcmp(model.rxns,'r_0042No1');
% Switch substrates and products
subsPos = find(model.S(:,rxnInd) == -1);
prodPos = find(model.S(:,rxnInd) == 1);
model.S(subsPos,rxnInd) = 1;
model.S(prodPos,rxnInd) = -1;

% W0TF60 (KLMA_60417, EC1.14.19.41) in r_0326No1 was found as the most 
% growth-limiting enzyme by modifyKcats. Change kcat according to the 
% identified value (0.17 s-1)
rxnInd                     = strcmp(model.rxns,'r_0326No1');
proteinInd                 = strcmp(model.mets,'prot_W0TF60');
newCoeff                   = -1/(0.17*3600);
model.S(proteinInd,rxnInd) = newCoeff;

% W0TBC7 (KLMA_50283, EC6.3.5.3) in r_0510 was also found by modifyKcats 
% to be growth-limiting. Change kcat according to the modified value
% (5.06748 s-1)
rxnInd                     = strcmp(model.rxns,'r_0510No1');
proteinInd                 = strcmp(model.mets,'prot_W0TBC7');
newCoeff                   = -1/(5.06748*3600);
model.S(proteinInd,rxnInd) = newCoeff;

% W0TFW7 (KLMA_80221, EC5.3.1.6), previous value (2 s-1) was found to be
% growth-limiting. Replace with kcat used in ecYeastGEM (9.6 s-1)
rxnInd                     = strcmp(model.rxns,'r_0059_REVNo1');
proteinInd                 = strcmp(model.mets,'prot_W0TFW7');
newCoeff                   = -1/(9.6*3600);
model.S(proteinInd,rxnInd) = newCoeff;

% W0TBM3 (KLMA_30204, EC1.14.14.17) identified as growth-limiting by
% modifyKcats. Use modifeid value (0.181333 s-1)
rxnInd                     = strcmp(model.rxns,'r_0316No1');
proteinInd                 = strcmp(model.mets,'prot_W0TBM3');
newCoeff                   = -1/(0.181333*3600);
model.S(proteinInd,rxnInd) = newCoeff;

% W0T2P5 (KLMA_10220, EC2.6.1.52) identified as growth-limiting by
% modifyKcats. Use modifeid value (22.675 s-1)
rxnInd                     = strcmp(model.rxns,'r_0638No1');
proteinInd                 = strcmp(model.mets,'prot_W0T2P5');
newCoeff                   = -1/(22.675*3600);
model.S(proteinInd,rxnInd) = newCoeff;

% W0T9E2 (KLMA_20256, EC2.4.1.1) in r_1151_REV has very low kcat 
% (~0.001 s-1), resulting in a very high requirement of the protein at max 
% growth rate. Use same kcat at forward direction (50.5 s-1).
rxnInd                     = strcmp(model.rxns,'r_1151_REVNo1');
proteinInd                 = strcmp(model.mets,'prot_W0T9E2');
newCoeff                   = -1/(50.5*3600);
model.S(proteinInd,rxnInd) = newCoeff;

% W0TCS1 (KLMA_40587, EC2.4.2.7) is growth-limiting (previous kcat = 0.0095
% s-1. Use value from giardia intestinalis (2.8 s-1)
rxnInd                     = ismember(model.rxns,{'r_0450No1';'r_0475No1';'r_0509No1'});
proteinInd                 = strcmp(model.mets,'prot_W0TCS1');
newCoeff                   = -1/(2.8*3600);
model.S(proteinInd,rxnInd) = newCoeff;

% W0T2Z0/W0T7K6 (KLMA_10153/KLMA_30312, EC4.2.1.3) growth-limiting. Use
% cerevisiae kcat (143.3 s-1) and account for two-step reaction
% ACO2a
rxnInd                     = ismember(model.rxns,{'r_0048No1';'r_0051_REVNo1'});
proteinInd                 = strcmp(model.mets,'prot_W0T2Z0');
newCoeff                   = -1/(2*143.3*3600);
model.S(proteinInd,rxnInd) = newCoeff;
% ACO2b
rxnInd                     = ismember(model.rxns,{'r_0048No2';'r_0051_REVNo2'});
proteinInd                 = strcmp(model.mets,'prot_W0T7K6');
newCoeff                   = -1/(2*143.3*3600);
model.S(proteinInd,rxnInd) = newCoeff;

% W0TBX9 (KLMA_50503, EC1.2.1.11) growth-limiting. Enzyme shows 82 %
% sequence similarity to S. cerevisiae ortholog -> use kcat from cerevisiae
% calculated from specific activity and Mw (189.8 s-1)
rxnInd                     = strcmp(model.rxns,'r_0634No1');
proteinInd                 = strcmp(model.mets,'prot_W0TBX9');
newCoeff                   = -1/(189.8*3600);
model.S(proteinInd,rxnInd) = newCoeff;

% W0TDY2 (KLMA_60541 (KGD1), EC1.2.4.2) is limiting. Update kcat to value 
% used for S. cerevisiae (27.8 s-1). Involved in a 2-step reaction 
% -> kcat value should be multiplied by 2
rxnInd                     = ismember(model.rxns,{'r_0044No1';'r_0055No1'});
proteinInd                 = strcmp(model.mets,'prot_W0TDY2');
newCoeff                   = -1/(2*27.8*3600);
model.S(proteinInd,rxnInd) = newCoeff;

% W0T9G3/W0TF38 (KLMA_30163/KLMA_60125 (gabD), EC1.2.1.16) 
rxnInd                     = strcmp(model.rxns,'r_0601No1');
proteinInd                 = strcmp(model.mets,'prot_W0T9G3');
newCoeff                   = -1/(0.066*3600);
model.S(proteinInd,rxnInd) = newCoeff;
rxnInd                     = strcmp(model.rxns,'r_0601No2');
proteinInd                 = strcmp(model.mets,'prot_W0TF38');
newCoeff                   = -1/(0.066*3600);
model.S(proteinInd,rxnInd) = newCoeff;

% W0TEG9 (KLMA_60167 (MDH1), EC1.1.1.37) growth-limiting. Use kcat value
% for cerevisiae (190.4 s-1)
rxnInd                     = strcmp(model.rxns,'r_0038No1');
proteinInd                 = strcmp(model.mets,'prot_W0TEG9');
newCoeff                   = -1/(190.4*3600);
model.S(proteinInd,rxnInd) = newCoeff;

% Change the kcat value for ATP synthase to match the maximum value
% measured experimentally in BRENDA (390 s-1)
rxnIndxs = find(contains(model.rxns,'r_0177No'));
newCoeff = -1/(200*3600);
% Update kcat value for each reaction
for i = 1:length(rxnIndxs)
    rxnInd      = rxnIndxs(i);
    subsPos     = find(model.S(:,rxnInd)<0);
    subsMets    = model.mets(subsPos);
    protSubsInd = subsPos(contains(subsMets,'prot_'));
    % Update coefficients for proteins
    model.S(protSubsInd,rxnInd) = newCoeff;
end

% Change the kcat value for complex III to maximum value measured in S.
% cerevisiae
rxnIndxs = find(contains(model.rxns,'r_0176No'));
newCoeff = -1/((240/2)*3600); % Assume functional complex is dimeric, as in S. cerevisiae
% Update kcat value for each reaction
for i = 1:length(rxnIndxs)
    rxnInd      = rxnIndxs(i);
    subsPos     = find(model.S(:,rxnInd)<0);
    subsMets    = model.mets(subsPos);
    protSubsInd = subsPos(contains(subsMets,'prot_'));
    % Update coefficients for proteins
    model.S(protSubsInd,rxnInd) = newCoeff;
end

% Change the kcat value for complex IV to maximum value measured in S.
% cerevisiae
rxnIndxs = find(contains(model.rxns,'r_0175No'));
newCoeff = -1/((693.5/2)*3600); % Assume functional complex is dimeric, as in S. cerevisiae
% Update kcat value for each reaction
for i = 1:length(rxnIndxs)
    rxnInd      = rxnIndxs(i);
    subsPos     = find(model.S(:,rxnInd)<0);
    subsMets    = model.mets(subsPos);
    protSubsInd = subsPos(contains(subsMets,'prot_'));
    % Update coefficients for proteins
    model.S(protSubsInd,rxnInd) = newCoeff;
end

% kcat of W0T8F1 (KLMA_30602 (ACS1)) and W0TC09 (KLMA_50533 (ACS2)) is
% orders of magnitude higher than the value for other yeasts (9400 s-1). 
% Replace with kcat from S. cerevisiae (62.4 s-1)
rxnInd                     = strcmp(model.rxns,'r_0004No1');
proteinInd                 = strcmp(model.mets,'prot_W0T8F1');
newCoeff                   = -1/(62.4*3600);
model.S(proteinInd,rxnInd) = newCoeff;

rxnInd                     = strcmp(model.rxns,'r_0004No2');
proteinInd                 = strcmp(model.mets,'prot_W0TC09');
newCoeff                   = -1/(62.4*3600);
model.S(proteinInd,rxnInd) = newCoeff;

% kcat of W0T4R5 (KLMA_10763, RAG2) is an order of magnitude higher 
% (3300 s-1) than in S. cerevisiae (685 s-1). Replace with value from 
% S. cerevisiae 
rxnInd                     = strcmp(model.rxns,'r_0030No1');
proteinInd                 = strcmp(model.mets,'prot_W0T4R5');
newCoeff                   = -1/(685*3600);
model.S(proteinInd,rxnInd) = newCoeff;

% kcat of W0TDE9 (KLMA_40125, TPI1) is twice as high (16700 s-1) as in 
% S. cerevisiae (8662.5 s-1). Replace with value from S. cerevisiae 
rxnInd                     = strcmp(model.rxns,'r_0020No1');
proteinInd                 = strcmp(model.mets,'prot_W0TDE9');
newCoeff                   = -1/(8662.5*3600);
model.S(proteinInd,rxnInd) = newCoeff;

% kcat of W0T3V9, replace with value from S. cerevisiae 
rxnInd                     = strcmp(model.rxns,'r_0061No1');
proteinInd                 = strcmp(model.mets,'prot_W0T3V9');
newCoeff                   = -1/(36.4*3600);
model.S(proteinInd,rxnInd) = newCoeff;

% kcat of W0TD99, replace with value from S. cerevisiae 
rxnInd                     = strcmp(model.rxns,'r_0062No1');
proteinInd                 = strcmp(model.mets,'prot_W0T3V9');
newCoeff                   = -1/(36.4*3600);
model.S(proteinInd,rxnInd) = newCoeff;
rxnInd                     = strcmp(model.rxns,'r_0062_REVNo1');
model.S(proteinInd,rxnInd) = newCoeff;

% kcat of W0THR5, asssume same kcat in both directions
rxnInd                     = strcmp(model.rxns,'r_0063No1');
proteinInd                 = strcmp(model.mets,'prot_W0THR5');
newCoeff                   = -1/(25.31*3600);
model.S(proteinInd,rxnInd) = newCoeff;
rxnInd                     = strcmp(model.rxns,'r_0063_REVNo1');
newCoeff                   = -1/(25.31*3600); % value from S. stipitis with correct substrates
model.S(proteinInd,rxnInd) = newCoeff;

% kcat of W0THR5 
rxnInd                     = strcmp(model.rxns,'r_0066No1');
proteinInd                 = strcmp(model.mets,'prot_W0THR5');
newCoeff                   = -1/(26.87*3600); % value from S. stipitis with fructose-6P as substrate
model.S(proteinInd,rxnInd) = newCoeff;
rxnInd                     = strcmp(model.rxns,'r_0066_REVNo1');
newCoeff                   = -1/(18.02*3600); % value from S. stipitis with correct substrates
model.S(proteinInd,rxnInd) = newCoeff;

% kcat of W0THR5 
rxnInd                     = strcmp(model.rxns,'r_0159No1');
proteinInd                 = strcmp(model.mets,'prot_W0THR5');
newCoeff                   = -1/(26.87*3600); % value from S. stipitis with fructose-6P as substrate
model.S(proteinInd,rxnInd) = newCoeff;
rxnInd                     = strcmp(model.rxns,'r_0159_REVNo1');
newCoeff                   = -1/(18.02*3600); % value from S. stipitis with correct substrates
model.S(proteinInd,rxnInd) = newCoeff;

% set kcat of W0TCY6 (KLMA_50631, ZWF) to value from S. cerevisiae
rxnInd                     = strcmp(model.rxns,'r_0068No1');
proteinInd                 = strcmp(model.mets,'prot_W0TCY6');
newCoeff                   = -1/(632.8*3600);
model.S(proteinInd,rxnInd) = newCoeff;

% Update rxn r_1204 - Pyruvate dehydrogenase (add all subunits required +
% update the kcat value)
% rename original reaction to enable removing it afterwards
rxnID    = 'r_1204No1';
proteins = {'W0T5D1';'W0T7E5';'W0THU7';'W0TFW6';'W0TER0'}; 
rxnName  = model.rxnNames{strcmp(model.rxns,rxnID)};
grRule   = 'KLMA_20172 and KLMA_20278 and KLMA_80065 and KLMA_60420 and KLMA_60529';
newMets  = cell(size(proteins'));
for i = 1:length(proteins)
    protein    = proteins{i};
    newMet     = ['prot_' protein];
    newMets{i} = newMet;
end
kvalues = 83.3*3600*ones(size(proteins')); % s-1 -> h-1, use value from S. cerevisiae
newID   = rxnID;
newName = rxnName;
model   = addEnzymesToRxn(model,kvalues,rxnID,newMets,{newID,newName},grRule);
% add additional information for proteins
for i = 1:length(proteins)
    if sum(strcmp(model.enzymes,proteins{i})) == 0
        model = addProtein(model,proteins{i},kegg,swissprot);
    end
end

% r_1010 - 5,10-methylenetetrahydrofolate:NADP+ oxidoreductase
proteins = {'W0TE99';'W0TEY9'}; % ADE3, MIS1
rxnID    = 'r_1011';
grR = model.grRules{strcmp(model.rxns,rxnID)};
rxnName  = model.rxnNames{strcmp(model.rxns,rxnID)};
% Add arm reaction
model  = addArmReaction(model,rxnID);
for i = 1:length(proteins)
    protein = proteins{i};
    newID   = [rxnID 'No' num2str(i)];
    newName = [rxnName ' (No' num2str(i) ')'];
    kvalues = 9.4*3600; % s-1 -> h-1, use value from Human
    newMet  = ['prot_' protein];
    grRule  = split(grR,' or ');
    grRule  = grRule{i};
    model   = addEnzymesToRxn(model,kvalues,rxnID,newMet,{newID,newName},grRule);
end
% remove original reaction
model = removeReactions(model,{rxnID});
% add additional information for proteins
for i = 1:length(proteins)
    model = addProtein(model,proteins{i},kegg,swissprot);
end
cd(current)
end