%% Flexibilize kcats in ecYeastGEM to allow growth at experimentally measured rate
ecYeastGEM = load('../../models/ecYeastGEM.mat');
ecYeastGEM = ecYeastGEM.ecModel;
cd ../customGECKO/
[~,modifications] = manualModifications(ecYeastGEM);
% Don't allow further modification of FBA1 since already modified manually
modifications{end+1} = 'P14540_3378';
% Remove incorrect usage of biotin protein ligase BPL1 (P48445) from model
ecYeastGEM.S(strcmp(ecYeastGEM.mets,'prot_P48445'),strcmp(ecYeastGEM.rxns,'r_0109No1')) = 0;
% Block use of r_4264No1
ecYeastGEM = setParam(ecYeastGEM,'ub','r_4264No1',0);
% Block use of r_0714_REVNo1
ecYeastGEM = setParam(ecYeastGEM,'ub','r_0714_REVNo1',0);

% Flexibilize kcats to allow experimentally observed growth rate
cd ../customGECKO
[ecYeastGEM_batch_modified,OptSigma] = getConstrainedModel(ecYeastGEM,modifications,'ecYeastGEM');
disp(['Sigma factor (fitted for growth on glucose): ' num2str(OptSigma)])
protPos = strcmp(ecYeastGEM_batch_modified.rxns,'prot_pool_exchange');
ecYeastGEM_batch_modified.ub(protPos) = ecYeastGEM_batch_modified.ub(protPos)*(0.53/0.52);
%% Modify kcats identified as limiting by getConstrained model
% P38604 (ERG7) in lanosterol synthase (No1)
ecYeastGEM_modified = ecYeastGEM;
ecYeastGEM_modified.S(strcmp(ecYeastGEM_modified.mets,'prot_P38604'),3657) = -1/(4.0763*3600);
% P38972 (ADE6) in 5'-phosphoribosylformyl glycinamidine synthetase (No1)
ecYeastGEM_modified.S(strcmp(ecYeastGEM_modified.mets,'prot_P38972'),2898) = -1/(5.0675*3600);
% P00931 (TRP5) in tryptophan synthase (indoleglycerol phosphate) (No1)
ecYeastGEM_modified.S(strcmp(ecYeastGEM_modified.mets,'prot_P00931'),4061) = -1/(775.75*3600);
% P05694 (MET6) in methionine synthase (No1)
ecYeastGEM_modified.S(strcmp(ecYeastGEM_modified.mets,'prot_P05694'),3687) = -1/(3.5*3600);
% P39006 (PSD1) in PS decarboxylase (1-16:1, 2-16:1), mitochondrial membrane (No1)
ecYeastGEM_modified.S(strcmp(ecYeastGEM_modified.mets,'prot_P39006'),[4565;4566;4567;4568;4569;4570;4571;4572]) = -1/(366.667*3600);
% P36148 (GPT2) in glycerol-3-phosphate acyltransferase (16:1), ER membrane (No2)
ecYeastGEM_modified.S(strcmp(ecYeastGEM_modified.mets,'prot_P36148'),4290) = -1/(9.9733*3600);
% P54839 (ERG13) in hydroxymethylglutaryl CoA synthase (No1)
ecYeastGEM_modified.S(strcmp(ecYeastGEM_modified.mets,'prot_P54839'),3527) = -1/(4.55*3600);
% Q12122 (LYS21) in lysine biosynthesis
ecYeastGEM_modified.S(strcmp(ecYeastGEM_modified.mets,'prot_Q12122'),4106) = -1/(5.133*3600);
% P32476 (ERG1) in ergosterol biosynthesis
ecYeastGEM_modified.S(strcmp(ecYeastGEM_modified.mets,'prot_P32476'),4009) = -1/(0.181333*3600);
% P39533 (ACO2) in lysine biosynthesis
ecYeastGEM_modified.S(strcmp(ecYeastGEM_modified.mets,'prot_P39533'),3501) = -1/(2.5*3600);
% P28777 (ARO2) in chorismate biosynthesis
ecYeastGEM_modified.S(strcmp(ecYeastGEM_modified.mets,'prot_P28777'),3203) = -1/(105.93*3600);
% Correct kcat for PDC5
% P16467 (PDC5) in pyruvate decarboxylase (No3)
ecYeastGEM_modified.S(strcmp(ecYeastGEM_modified.mets,'prot_P16467'),strcmp(ecYeastGEM_modified.rxns,'r_0959No3')) = -1/(207*3600);

% Limit secretion of unmeasured metabolites
% limit secretion of 2,3-butanediol
ecYeastGEM_modified = setParam(ecYeastGEM_modified,'ub','r_1549',1e-5);
% limit isobutanol secretion
ecYeastGEM_modified = setParam(ecYeastGEM_modified,'ub','r_1866',1e-5);
% limit secretion of isobutyraldehyde
ecYeastGEM_modified = setParam(ecYeastGEM_modified,'ub','r_1870',1e-5);
% limit secretion of L-valine
ecYeastGEM_modified = setParam(ecYeastGEM_modified,'ub','r_1914',1e-5);
% limit secretion of 6-phospho-D-gluconate
ecYeastGEM_modified = setParam(ecYeastGEM_modified,'ub','r_4538',1e-5);
% limit secretion of glycerone
ecYeastGEM_modified = setParam(ecYeastGEM_modified,'ub','r_4507',1e-5);
% limit secretion of palmitoleate
ecYeastGEM_modified = setParam(ecYeastGEM_modified,'ub','r_1994',1e-5);
% limit secretion of L-phenylalanine
ecYeastGEM_modified = setParam(ecYeastGEM_modified,'ub','r_1903',1e-5);
% limit secretion of hypxanthine
ecYeastGEM_modified = setParam(ecYeastGEM_modified,'ub','r_1841',1e-5);
% limit secretion of AMP
ecYeastGEM_modified = setParam(ecYeastGEM_modified,'ub','r_4548',1e-5);
% limit secretion of acetaldehyde
ecYeastGEM_modified = setParam(ecYeastGEM_modified,'ub','r_4548',1e-5);
% Limit secretion of L-alanine
ecYeastGEM_modified = setParam(ecYeastGEM_modified,'ub','r_1873',0);

%% Generate a protein constrained model from experimental data
generate_protModels(ecYeastGEM_modified,3,'WTstudy_ecYeastGEM_protConstrained',ecYeastGEM_batch_modified);
%% Test performance of protein constrained model

%% Flexibilize kcats in eciSM966 to allow growth at experimentally measured rate
% Load the model
ecKmarx_mod = load('../../models/eciSM996.mat');
ecKmarx = ecKmarx_mod.ecModel;
% Update the amino acid composition in the protein pseudoreaction
% to sum up to 1 g amino acids/g protein (currently 0.45 g/g)
% get index of protein pseudoreaction
protRxnPos = strcmpi(ecKmarx.rxnNames,'protein pseudoreaction');
aaPos = find(full(ecKmarx.S(:,protRxnPos))<0);
aas   = ecKmarx.metNames(aaPos);
aaCoeffs = full(ecKmarx.S(aaPos,protRxnPos));
% Read file with corrected abundances for amino acids (mmol/g protein)
% calculated from proteome composition
fileName = '../../data/databases/Kmarx_AA_composition.txt';
fid = fopen(fileName);
aaData = textscan(fid,'%s%f','delimiter','\t','HeaderLines',1);
aaName = aaData{1};
aaAbundance = aaData{2};
fclose(fid);
% Update coefficients in protein pseudoreaction
for i = 1:length(aaName)
    aa = aaName{i};
    abundance = aaAbundance(i);
    % Get metabolite index of amino acid
    ind = aaPos(strcmp(aas,aa));
    % Update coefficient in S matrix
    disp([aa ' coeffcient before: ' num2str(ecKmarx.S(ind,protRxnPos))])
    ecKmarx.S(ind,protRxnPos) = -abundance;
    disp([aa ' coefficient after: ' num2str(ecKmarx.S(ind,protRxnPos))])
end

%% Generate a model constrained by experimental measurements
% Limit secretion of unmeasured metabolites
ecKmarx_modified = setParam(ecKmarx,'ub',ecKmarx.rxns(strcmpi(ecKmarx.rxnNames,'acetate exchange')),0.4783);
ecKmarx_modified = setParam(ecKmarx_modified,'ub',ecKmarx_modified.rxns(strcmpi(ecKmarx_modified.rxnNames,'pyruvate exchange')),1e-5);
ecKmarx_modified = setParam(ecKmarx_modified,'ub',ecKmarx_modified.rxns(strcmpi(ecKmarx_modified.rxnNames,'acetaldehyde exchange')),1e-5);
ecKmarx_modified = setParam(ecKmarx_modified,'ub',ecKmarx_modified.rxns(strcmpi(ecKmarx_modified.rxnNames,'ethanol exchange')),1e-5);
ecKmarx_modified = setParam(ecKmarx_modified,'ub',ecKmarx_modified.rxns(strcmpi(ecKmarx_modified.rxnNames,'L-serine exchange')),1e-5);
ecKmarx_modified = setParam(ecKmarx_modified,'ub',ecKmarx_modified.rxns(strcmpi(ecKmarx_modified.rxnNames,'citrate exchange')),1e-5);
ecKmarx_modified = setParam(ecKmarx_modified,'ub',ecKmarx_modified.rxns(strcmpi(ecKmarx_modified.rxnNames,'ethyl acetate exchange')),1e-5);
ecKmarx_modified = setParam(ecKmarx_modified,'ub',ecKmarx_modified.rxns(strcmpi(ecKmarx_modified.rxnNames,'2-isopropylmalate exchange')),1e-5);
ecKmarx_modified = setParam(ecKmarx_modified,'ub',ecKmarx_modified.rxns(strcmpi(ecKmarx_modified.rxnNames,'L-valine exchange')),1e-5);
ecKmarx_modified = setParam(ecKmarx_modified,'ub',ecKmarx_modified.rxns(strcmpi(ecKmarx_modified.rxnNames,'L-glutamate exchange')),1e-5);
ecKmarx_modified = setParam(ecKmarx_modified,'ub',ecKmarx_modified.rxns(strcmpi(ecKmarx_modified.rxnNames,'glyoxylate exchange')),1e-5);
ecKmarx_modified = setParam(ecKmarx_modified,'ub',ecKmarx_modified.rxns(strcmpi(ecKmarx_modified.rxnNames,'formate exchange')),1e-5);
% Block incorrect NADPH-producing reactions
% cytosolic NADP-dependent isocitrate dehydrogeanase, repressed by glucose
ecKmarx_modified = setParam(ecKmarx_modified,'ub','arm_r_0861',0);
ecKmarx_modified = setParam(ecKmarx_modified,'ub','arm_r_0049',0);
ecKmarx_modified = setParam(ecKmarx_modified,'ub','arm_r_0035',0);
% Avoid cycles where NADH is converted to NADPH by ALD genes
ecKmarx_modified = setParam(ecKmarx_modified,'ub','r_0778',0);
ecKmarx_modified = setParam(ecKmarx_modified,'ub','r_0782',0);
ecKmarx_modified = setParam(ecKmarx_modified,'ub','r_0788',0);
ecKmarx_modified = setParam(ecKmarx_modified,'ub','r_0792',0);
ecKmarx_modified = setParam(ecKmarx_modified,'ub','r_0794',0);
ecKmarx_modified = setParam(ecKmarx_modified,'ub','r_0796',0);
ecKmarx_modified = setParam(ecKmarx_modified,'ub','r_0803',0);
ecKmarx_modified = setParam(ecKmarx_modified,'ub','r_1093',0);
ecKmarx_modified = setParam(ecKmarx_modified,'ub','r_0624',0);
% Avoid TCA bypass
% Block r_0599 Succinate-semialdehyde:NAD+ oxidoreductase to avoid TCA
% cycle bypass
ecKmarx_modified = setParam(ecKmarx_modified,'ub','arm_r_0599',0);
ecKmarx_modified = setParam(ecKmarx_modified,'ub','arm_r_0601',0);

% Run pipeline to constrain model
[ecKmarx_batch,OptSigma] = getConstrainedModel_Kmarx(ecKmarx_modified,{},'ecKmarx_batch12');
disp(['Sigma factor (fitted for growth on glucose): ' num2str(OptSigma)])

% Perform kcat modifications on batch model and regular model according to getConstrainedModel
cd ../modelGeneration
ecKmarx_batch_modified = manualChangesKmarxModel(ecKmarx_batch);
ecKmarx_modified = manualChangesKmarxModel(ecKmarx_modified);
% Curate bounds on reactions involving NADPH/NADH to prevent cycling between the two
%ecKmarx_modified = setParam(ecKmarx_modified,'ub','arm_r_0009',0.1);
%ecKmarx_batch_modified = setParam(ecKmarx_batch_modified,'ub','arm_r_0009',0.1);
%ecKmarx_modified = setParam(ecKmarx_modified,'ub','r_0009No1',0);
%ecKmarx_batch_modified = setParam(ecKmarx_batch_modified,'ub','r_0009No1',0);
ecKmarx_modified = setParam(ecKmarx_modified,'ub','arm_r_0609_REV',0);
ecKmarx_batch_modified = setParam(ecKmarx_batch_modified,'ub','arm_r_0609_REV',0);
ecKmarx_modified = setParam(ecKmarx_modified,'ub','r_0624',0);
ecKmarx_batch_modified = setParam(ecKmarx_batch_modified,'ub','r_0624',0);
ecKmarx_modified = setParam(ecKmarx_modified,'ub','r_0624_REV',0);
ecKmarx_batch_modified = setParam(ecKmarx_batch_modified,'ub','r_0624_REV',0);
ecKmarx_modified = setParam(ecKmarx_modified,'ub','r_0778_REV',0);
ecKmarx_batch_modified = setParam(ecKmarx_batch_modified,'ub','r_0778_REV',0);
ecKmarx_modified = setParam(ecKmarx_modified,'ub','r_0782_REV',0);
ecKmarx_batch_modified = setParam(ecKmarx_batch_modified,'ub','r_0782_REV',0);
ecKmarx_modified = setParam(ecKmarx_modified,'ub','r_0788_REV',0);
ecKmarx_batch_modified = setParam(ecKmarx_batch_modified,'ub','r_0788_REV',0);
ecKmarx_modified = setParam(ecKmarx_modified,'ub','r_0792_REV',0);
ecKmarx_batch_modified = setParam(ecKmarx_batch_modified,'ub','r_0792_REV',0);
ecKmarx_modified = setParam(ecKmarx_modified,'ub','r_0794_REV',0);
ecKmarx_batch_modified = setParam(ecKmarx_batch_modified,'ub','r_0794_REV',0);
ecKmarx_modified = setParam(ecKmarx_modified,'ub','r_0796_REV',0);
ecKmarx_batch_modified = setParam(ecKmarx_batch_modified,'ub','r_0796_REV',0);
ecKmarx_modified = setParam(ecKmarx_modified,'ub','r_0803_REV',0);
ecKmarx_batch_modified = setParam(ecKmarx_batch_modified,'ub','r_0803_REV',0);
ecKmarx_modified = setParam(ecKmarx_modified,'ub','r_1093_REV',0);
ecKmarx_batch_modified = setParam(ecKmarx_batch_modified,'ub','r_1093_REV',0);
ecKmarx_modified = setParam(ecKmarx_modified,'ub','r_0072No1',0);
ecKmarx_batch_modified = setParam(ecKmarx_batch_modified,'ub','r_0072No1',0);
ecKmarx_modified = setParam(ecKmarx_modified,'ub','r_0630No1',0);
ecKmarx_batch_modified = setParam(ecKmarx_batch_modified,'ub','r_0630No1',0);
ecKmarx_modified = setParam(ecKmarx_modified,'ub','r_0727No1',0);
ecKmarx_batch_modified = setParam(ecKmarx_batch_modified,'ub','r_0727No1',0);
ecKmarx_modified = setParam(ecKmarx_modified,'ub','r_1168',0);
ecKmarx_batch_modified = setParam(ecKmarx_batch_modified,'ub','r_1168',0);
ecKmarx_modified = setParam(ecKmarx_modified,'ub','r_1016No1',0);
ecKmarx_batch_modified = setParam(ecKmarx_batch_modified,'ub','r_1016No1',0);
ecKmarx_modified = setParam(ecKmarx_modified,'ub','arm_r_0450',0.5);
ecKmarx_batch_modified = setParam(ecKmarx_batch_modified,'ub','arm_r_0450',0.5);
%ecKmarx_modified = setParam(ecKmarx_modified,'ub',ecKmarx_modified.rxns(contains(ecKmarx_modified.rxns,'prot_W0T8I5')),0);
%ecKmarx_batch_modified = setParam(ecKmarx_batch_modified,'ub',ecKmarx_batch_modified.rxns(contains(ecKmarx_batch_modified.rxns,'prot_W0T8I5')),0);

% Generate a model constrained with experimental measurements
generate_protModels_Kmarx(ecKmarx_modified,3,'WTstudy_ecKmarx_protConstrained11',ecKmarx_batch_modified);
