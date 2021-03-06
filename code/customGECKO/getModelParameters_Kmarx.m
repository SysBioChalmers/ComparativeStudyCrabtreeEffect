function parameters = getModelParameters_Kmarx
% getModelParameters
%
%   Set model and organism specific parameters that are used by the
%   ecModel generation pipeline.
%
%   Ivan Domenzain. Last edited: 2019-07-29
%

%Average enzyme saturation factor
parameters.sigma          = 0.5;
%Total protein content in the cell [g protein/gDw]
parameters.Ptot           = 0.39163;      %Assumed constant
%Minimum growth rate the model should grow at [1/h]
parameters.gR_exp         = 0.4411;     %[g/gDw h] 
%Provide your organism scientific name
parameters.org_name       = 'kluyveromyces marxianus';
%Provide your organism KEGG ID
parameters.keggID         = 'kmx';
%The name of the exchange reaction that supplies the model with carbon (rxnNames)
parameters.c_source       = 'D-glucose exchange (reversible)'; 
%Rxn Id for biomass pseudoreaction
parameters.bioRxn         = 'r_1912';
%Compartment name in which the added enzymes should be located
parameters.enzyme_comp    = 'Cytoplasm';
%Rxn names for the most common experimentally measured "exchange" fluxes
parameters.exch_names{1}  = 'biomass exchange';
parameters.exch_names{2}  = 'D-glucose exchange (reversible)';
parameters.exch_names{3}  = 'oxygen exchange (reversible)';
parameters.exch_names{4}  = 'carbon dioxide exchange';
%biomass components pseudoreactions (proteins, carbs and lipids lumped pools)
parameters.bio_comp{1}  = 'protein';
parameters.bio_comp{2}  = 'carbohydrate';
parameters.bio_comp{3}  = 'lipid';
%Rxn IDs for reactions in the oxidative phosphorylation pathway (optional)
parameters.oxPhos{1} = 'r_0052';
parameters.oxPhos{2} = 'r_0175';
parameters.oxPhos{3} = 'r_0176';
parameters.oxPhos{4} = 'r_0177';
%
parameters.NGAM = 'r_1905';
end