function [ModelAdjusted, RxnBalancingList] = AdjustOriginalAraCore
%==========================================================================
% This function introduces several adjustmens to the AraCore model
% published by Arnold and Nikoloski (2014) Plant Physiol. 165(3):1380-1391.
% doi: 10.1104/pp.114.235358, as described below, in order to guarantee a
% complete integrattion with the LipidModule.
%
% USAGE:
%
%    [ModelAdjusted, RxnBalancingList] = AdjustOriginalAraCore
%
%                          
% OUTPUTS:
%    ModelAdjusted:     COBRA model structure that includes the same set 
%                       of reactions as the source model, but with the 
%                       modifications described below:
%                       
%                     * Harmonization of KEGGID for metabolites which are
%                       equivalent in the LipidModule but appear with
%                       different identifier.
%
%                     * Addition of missig KEGG identifiers for
%                       metabolites.
%
%                     * Harmonization of metabolites charges and molecular
%                       formulas.
%
%                     * Adjustment of name for selected metabolites to
%                       avoid confusion among metabolites which are not
%                       equivalent in LipidModule.
%
%   RxnBalancingList:   List containing all the reaction which are mass- 
%                       and/or charge-imbalanced reaction, excluding
%                       reactions naturally imbalanced such as exchange,
%                       import, siphon, demand, biomass and SLIME reactions.                     
%==========================================================================
%% 1). Load 'AraCore Model'

 modelPath = fullfile('InputFiles', 'AraCorewKEGG.xls'); 
 model = readCbModel(modelPath);

%% 2). Unification of KEGG identifiers: 
% The KEGGID for several metabolites differ between AraCore and the 
% LipidModule, hence, they are harmonized for the reasons explained next:
%
% * FBP[h] = the isomer that is metabolized in the Glycolysis /
%            Gluconeogenesis pathway correspond to beta-D-Fructose 
%            1,6-bisphosphate, which has the KEGGID C05378
%
% * F6P[h] = the isomer that is metabolized in the Glycolysis /
%            Gluconeogenesis pathway correspond to beta-D-Fructose 
%            6-phosphate which has the KEGGID C05345
%
% * G6P[c] = the isomer that act as intermediate in the PPP pathway is
%            beta-D-Glucose 6-phosphate with KEGGID C01172

MetaboliteList = {'G6P[h]' 'G6P[c]' 'FBP[h]' 'FBP[c]' 'F6P[h]' 'F6P[c]'};

KEGGIDs = {'C01172' 'C01172' 'C05378' 'C05378' 'C05345' 'C05345'};

for i = 1:size(MetaboliteList,2)
    MetID = MetaboliteList{1,i};
    KEGGID_i = KEGGIDs{1,i};
    MetIdx = strcmp(model.mets, MetID);

    if ~isempty(MetIdx)
        model.metKEGGID{MetIdx} = KEGGID_i;
    end
end

%% 3). Addition of missig KEGG identifiers: 
% There are a set of metabolites whose corresponding KEGGIDs are available
% but are missing in AraCore.

MetaboliteList = {'ADPG[h]' 'starch1[h]' 'starch2[h]' 'starch3[h]' 'starch5[h]'...
    'Mas[h]' 'Mas[c]' 'Glc[h]' 'Glc[c]' 'starch1[c]' 'starch2[c]'...
    'cellulose1[c]' 'cellulose2[c]' 'cellulose3[c]' 'A-DHL[m]' 'DHL[m]'...
    'S-DHL[m]' 'hnu[h]' 'DHP[m]' 'LPA[m]' 'ORO[m]' 'ORO[c]' 'CAIR[h]' 'DHO[c]'...
    'Ala[c]' 'Arg[c]' 'Asn[c]' 'dCTP[c]' 'dGTP[c]' 'dTTP[c]' 'Frc[c]' 'GABA[c]'...
    'His[c]' 'Ile[c]' 'Leu[c]' 'Orn[h]' 'Pro[c]' 'Suc[c]' 'Thr[c]' 'Tre[c]'...
    'Trp[c]' 'Tyr[c]' 'urea[m]' 'Val[c]'};

KEGGIDs = {'C00498' 'G10495' 'G10495' 'G10495' 'G10495'...
    'C00208' 'C00208' 'C00031' 'C00031' 'G10495' 'G10495'...
    'C00760' 'C00760' 'C00760' 'C16255' 'C15973'...
    'C16254' 'C00205' 'C15973' 'C15972' 'C00295' 'C00295' 'C04751' 'C00337'...
    'C00041' 'C00062' 'C00152' 'C00458' 'C00286' 'C00459' 'C02336' 'C00334'...
    'C00135' 'C00407' 'C00123' 'C00077' 'C00148' 'C00089' 'C00188' 'C01083'...
    'C00078' 'C00082' 'C00086' 'C00183'};

for i = 1:size(MetaboliteList,2)
    MetID = MetaboliteList{1,i};
    KEGGID_i = KEGGIDs{1,i};
    MetIdx = strcmp(model.mets, MetID);

    if ~isempty(MetIdx)
        model.metKEGGID{MetIdx} = KEGGID_i;
    end
end

%% 4). Harmonization of metabolites charges and molecular formulas: 
%     There are a set of metabolites whose charge and/or formulas differ 
%     from the equivalent metabolite in LipidModule, which posses a problem 
%     when trying to connect the respective metabolites between the models. 
%     To prevent later incorporation of duplicated metabolites, the charges 
%     and formulas are harmonized. There are other metabolites whose charges
%     and/or formulas were misassigned. An explanation is included case by 
%     case next:
%
% * DHL[m] = The reaction EC:2.3.1.12 catalyzed by the pyruvate dehydrogenase
%            complex, where the component E2 (dihydrolipoamide acetyltransferase) is 
%            represented in https://www.uniprot.org/uniprot/Q0WQF7, where 
%            (R)-N6-dihydrolipoyl-L-lysyl-[protein] has the molecular formula
%            C14H26N2O2S2 
%
% * A-DHL[m] = This metabolite is generated as a product in the same reaction
%            2.3.1.12, where (R)-N6-(S8-acetyldihydrolipoyl)-L-lysine residue has the
%            molecular formula C16H28N2O3S2
%
% * S-DHL[m] = The reaction EC:2.3.1.61 is catalyzed by 2-oxoglutarate
%            dehydrogenase complex, where the 2-oxoglutarate dehydrogenase complex 
%            component E2-1 (dihydrolipoamide succinyltransferase) is represented in
%            https://www.uniprot.org/uniprot/Q9FLQ4, where 
%            (R)-N6-(S8-succinyldihydrolipoyl)-L-lysine residue has the molecular 
%            formula C18H30N2O5S2
%
% * cplx4_m = The reaction EC 1.9.3.1 Transferred to EC:7.1.1.9 is catalyzed
%             by cytochrome c oxidase, the last enzyme in the mitochondrial electron
%             transport chain is represented in https://www.uniprot.org/uniprot/P93285,
%             where the respective charges for [Fe(II)cytochrome c] and
%             [Fe(III)cytochrome c] are +2 and +3, respectively. These charges will be
%             corrected in AraCore model since a -2 charge was assigned to both 
%             metabolites. A generic molecular formula (FeX) is used for these 
%             metabolites to facilitate the mass balancing.
%
% * amDHP[m] = The reaction EC:1.4.4.2 is catalyzed by The glycine
%            decarboxylase (GDC) or glycine cleavage system where
%            (R)-N6-(S8-aminomethyldihydrolipoyl)-L-lysine residue has a molecular
%            formula C15H30N3O2S2
%
% * ACP[h] = The order of the elements in the molecular formula is adjusted
%            alphabetically to HRS. The charge of the compound is changed
%            from o to -1.       
%
% * M-ACP[h] = The charge ia adjusted to take into account the extra -1
%            contributed by ACP
%
% * NAD+ = The charge of the oxidized form of this metabolite is adjusted to
%          be -1 since when this cofactor is oxidized the negative charge of one
%          phosphate group is neutralized by the cation derived from the 
%          tetravalent state of the nitrogen atom of the nicotinamide group.
%
% * NADH[c] = The charge originally was set to -1, and was corrected to -2,
%           since the reduced form of this cofactor has two phosphate groups and the
%           nitrogen atom of the nicotinamide group has a trivalent state.
%
% * NADP+/NADPH = The charge of this cofactor was also adjusted for the
%               reason explained before. In the original model the charge for the
%               oxidixed fors was -4, but it is -3, and the reduced form in the cytosol
%               has a charge -3, which should be -4. The required adjustments are made to
%               have ozidized form -3 and reduced form -4.
%
% * GDP[h]/GDP[c] = The charge and molecular form of this cofactor is
%               adjusted, since is has two phosphate groups the charge correspinds to -3
%               and the formula according to the charge is C10H12N5O11P2
%               https://www.rhea-db.org/rhea/15753.
%
% * GMP[c] = The charge and formula in the original model is -1 and
%           C10H13N5O8P, and are adjusted to -2 and C10H12N5O8P, resepctively, 
%           according to https://www.uniprot.org/uniprot/P93757,
%           https://www.uniprot.org/uniprot/Q9CAD1.
%
% * GTP[h]/GTP[c] = The charge and formula are adjusted according to
%               https://www.rhea-db.org/rhea/15753 to -4 and C10H12N5O14P3
%
% * UDP[c] = The charge and formula are adjusted according to
%           https://www.rhea-db.org/rhea/19929 to -3 C9H11N2O12P2.
%
% * AICAR[h] = Adjust charge and formula to -2 and C9H13N4O8P, according to
%           https://www.uniprot.org/uniprot/Q8RY94,
%           https://www.rhea-db.org/rhea/23920 
%
% * AIR[h] = Adjust charge and formula to -1 and C8H13N3O7P, according to
%           https://www.rhea-db.org/rhea/10792,
%           https://www.uniprot.org/uniprot/Q84TI2.
%
% * A-CoA[c] = Adjust charge and formula to -4 and C23H34N7O17P3S, according
%           to https://www.rhea-db.org/rhea/16845,
%           https://www.uniprot.org/uniprot/Q9LXS6.
%
% * H[h]/H[l]/H[c]/H[m]/H[i]/H[p] = Hydrogen atoms have a charge of +1. 
%           The adjustments will be made for all species in the model.

% First adjust the metabolites charges:
MetaboliteList = {'Fdox[h]' 'Fdrd[h]' 'Cytcox[m]' 'Cytcrd[m]' 'ACP[h]'...
    'M-ACP[h]' 'NAD[c]' 'NAD[h]' 'NAD[m]' 'NAD[p]' 'NADH[c]'...
    'NADP[h]' 'NADP[c]' 'NADP[m]' 'NADPH[c]' 'GDP[h]' 'GDP[c]'...
    'GMP[c]' 'GTP[h]' 'GTP[c]' 'UDP[c]' 'AICAR[h]' 'AIR[h]' 'A-CoA[c]'...
    'H[h]' 'H[l]' 'H[c]' 'H[m]' 'H[i]' 'H[p]' 'PCox[h]' 'LPA[m]' 'CoA[c]'...
    'PPi[h]' 'H2S[h]' 'H2S[c]' 'H2S[m]' 'SO3[h]' 'CAIR[h]' 'PRPP[c]'...
    'DHO[c]' 'DHO[h]' 'DHO[m]' 'aMet[c]' 'FGAM[h]' 'CDP[c]'};

Charges = [2 1 3 2 -1 -2 -1 -1 -1 -1 -2 -3 -3 -3 -4 -3 -3 -2 -4 -4 -3 -2 -1 -4,...
    1 1 1 1 1 1 1 0 -4 -3 -1 -1 -1 -2 -2 -5 -1 -1 -1 1 -1 -3];

for i = 1:size(MetaboliteList,2)
    MetID = MetaboliteList{1,i};
    Charge_i = Charges(1,i);
    MetIdx = strcmp(model.mets, MetID);

    if ~isempty(MetIdx)
        model.metCharges(MetIdx) = Charge_i;
    end
end

% Next, adjust the molecular formulas:
MetaboliteList =  {'Fdox[h]' 'Fdrd[h]' 'DHL[m]' 'A-DHL[m]' 'S-DHL[m]'...
    'Cytcox[m]' 'Cytcrd[m]' 'amDHP[m]' 'ACP[h]' 'GDP[h]' 'GDP[c]' 'GMP[c]'...
     'GTP[h]' 'GTP[c]' 'UDP[c]' 'AICAR[h]' 'AIR[h]' 'A-CoA[c]' 'NADH[c]'...
     'M-ACP[h]' 'NADPH[c]' 'DHP[m]' 'LPA[m]' 'CoA[c]' 'PPi[h]' 'H2S[h]' 'H2S[c]' 'H2S[m]'...
     'SO3[h]' 'CAIR[h]' 'PRPP[c]' 'DHO[c]' 'DHO[h]' 'DHO[m]' 'aMet[c]' 'FGAM[h]'...
     'CDP[c]'};

MolecularFormulas = {'FeS8X' 'FeS8X' 'C14H26N2O2S2' 'C16H28N2O3S2' 'C18H29N2O5S2'...
    'FeX' 'FeX' 'C15H30N3O2S2' 'HRS' 'C10H12N5O11P2' 'C10H12N5O11P2' 'C10H12N5O8P'...
    'C10H12N5O14P3' 'C10H12N5O14P3' 'C9H11N2O12P2' 'C9H13N4O8P' 'C8H13N3O7P'...
    'C23H34N7O17P3S' 'C21H27N7O14P2' 'C3H2O3RS' 'C21H26N7O17P3' 'C14H26N2O2S2'...
    'C14H24N2O2S2' 'C21H32N7O16P3S' 'HO7P2' 'HS' 'HS' 'HS' 'O3S' 'C9H12N3O9P'...
    'C5H8O14P3' 'C5H5N2O4' 'C5H5N2O4' 'C5H5N2O4' 'C15H23N6O5S' 'C8H15N3O8P'...
    'C9H12N3O11P2'};

for i = 1:size(MetaboliteList,2)
    MetID = MetaboliteList{1,i};
    Form_i = MolecularFormulas{1,i};
    MetIdx = strcmp(model.mets, MetID);

    if ~isempty(MetIdx)
        model.metFormulas{MetIdx} = Form_i;
    end
end

%% 5). Adjust name of selected metabolites: 
% In reactions EC1.2.4.1/1.8.1.4/1.2.4.2 the intermediate that participates
% correspond to a protein N6-lipoyl-L-lysine. In the AraCore model there
% are two metabolites which represent this intermediate: LPA[m] (Lipoamide)
% and LPL[m] (Lipoylprotein). Since, there´s little informaion about the
% origin of these metabolites and to facilitate the integration
% of the LipidModule with the AraCore model, it is assumed that the
% N6-lipoyl-L-lysine corresponds to the same species derived from the
% Lipoic acid metabolism named DLipoyl generated in the mitochondria.
% Hence, the metabolites LPA[m] and LPL[m] in the AraCore model will be
% integrated as one single metabolite corresponding to DLipoyl[m].

% Replace LPA[m] ID: 
OldMetID = 'LPA[m]';
MetIdx = strcmp(model.mets, OldMetID);
model.mets{MetIdx} = 'DLipoyl[m]'; % New Met ID
model.metNames{MetIdx} = '(R)-N6-lipoyl-L-lysyl-[protein]'; %New Met Name

% The metabolite LPL[m] (Lipoylprotein) will be also replaced by DLipoyl[m]:
% Find metabolite index:
MetIndx = strcmp(model.mets, 'LPL[m]');
GetSindx = model.S(MetIndx, :);

% Find non-zero entries in S:
GetCoefficients = model.S(MetIndx, GetSindx ~= 0);
model.S(MetIdx, GetSindx ~= 0) = GetCoefficients;

% Set all coefficients to zero and remove from model:
model.S(MetIndx, :) = 0;
model = removeMetabolites(model, 'LPL[m]');

% The metabolite DHP[m] (Dihydrolipolprotein) will be also replaced by DHL[m]
% ((R)-N6-dihydrolipoyl-L-lysyl-(protein)):
% Find metabolite index:
MetIndx = strcmp(model.mets, 'DHP[m]');
GetSindx = model.S(MetIndx, :);

% Find non-zero entries in S and add to DHL[m]:
GetCoefficients = model.S(MetIndx, GetSindx ~= 0);
MetIdx = strcmp(model.mets, 'DHL[m]'); % Find DHL[m] index
model.S(MetIdx, GetSindx ~= 0) = GetCoefficients;

% Set all coefficients to zero and remove from model:
model.S(MetIndx, :) = 0;
model = removeMetabolites(model, 'DHP[m]');

% Replace DHL[m] Name: 
model.metNames{MetIdx} = '(R)-N6-dihydrolipoyl-L-lysyl-(protein)';

%% 6). Adjust coefficients of ferredoxin in reaction SO3 reductase EC1.8.7.1 R00859
% OldRxn : 7 H[h] + 3 Fdrd[h] + SO3[h]  -> 3 H2O[h] + 3 Fdox[h] + H2S[h]
% NewRxn : 7 H[h] + 6 Fdrd[h] + SO3[h]  -> 3 H2O[h] + 6 Fdox[h] + H2S[h]

MetsList = {'Fdrd[h]', 'Fdox[h]'};
GetRxnID = strcmp(model.rxns, 'SO3R_h');

for j = 1:size(MetsList,2)
    MetID = MetsList{1,j};
    GetMetIdx = strcmp(model.mets,  MetID);
    
    switch MetID
        case {'Fdrd[h]'}
            % Add 1 pho_loss[h] to the products side of the reaction:
            model.S(GetMetIdx, GetRxnID) = -6;
            
        case {'Fdox[h]'}
            % Add 4 pho_loss[h] in the products side of the reaction:
            model.S(GetMetIdx, GetRxnID) = 6;
    end
end

%% 6). Adjust selected reactions of model:
%
% a). Adjust photon balance in photosystems reactions:

% Create new metabolite:
model = addMetabolite(model, 'pho_loss[h]', 'metName','photon_loss', 'metFormula','X', 'Charge',0);

% Add exchange reaction for 'Ex_pho_loss' metabolite: 
model = addExchangeRxn(model, 'pho_loss[h]', 0, 1000);

% Add 'pho_loss[h]' metabolite to photosystems reactions:
% PSI_h, PSI_h  hnu[h] + PCrd[h] + Fdox[h]  -> PCox[h] + Fdrd[h] + pho_loss[h]
% PSII_h,   PSII_h 4 hnu[h] + 2 PQ[h] + 2 H2O[h] + 4 H[h]  -> 2 PQH2[h] + O2[h] + 4 H[l] + 4 pho_loss[h]

RxnsList = {'PSI_h', 'PSII_h'};
GetMetID = strcmp(model.mets, 'pho_loss[h]');

for j = 1:size(RxnsList,2)
    RxnID = RxnsList{1,j};
    GetRxnIdx = strcmp(model.rxns,  RxnID);
    
    switch RxnID
        case {'PSI_h'}
            % Add 1 pho_loss[h] to the products side of the reaction:
            model.S(GetMetID, GetRxnIdx) = 1;
            
        case {'PSII_h'}
            % Add 4 pho_loss[h] in the products side of the reaction:
            model.S(GetMetID, GetRxnIdx) = 4;
    end
end

% b). Adjust stoichiometry of H atoms in ATPAse reaction in chloroplast:
MetsList = {'H[l]', 'H[h]'};
GetRxnID = strcmp(model.rxns, 'ATPase_h');

for j = 1:size(MetsList,2)
    MetID = MetsList{1,j};
    GetMetIdx = strcmp(model.mets,  MetID);
    
    switch MetID
        case {'H[l]'}
            % Add 1 pho_loss[h] to the products side of the reaction:
            model.S(GetMetIdx, GetRxnID) = -13;
            
        case {'H[h]'}
            % Add 4 pho_loss[h] in the products side of the reaction:
            model.S(GetMetIdx, GetRxnID) = 10;
    end
end

%% 7). Change Orotate name to harmonize it with ID used in LipidModule:

ListOldIDs = {'ORO[m]' 'ORO[c]' 'CAIR[h]'};
ListNewIDs = {'OA[m]' 'OA[c]' 'APRC[h]'};

for j = 1:size(ListOldIDs,2)
    OldID_i = ListOldIDs{1,j};
    NewID_i = ListNewIDs{1,j};
    GetMetIdx = strcmp(model.mets,  OldID_i);
    model.mets{GetMetIdx} = NewID_i;
end 

%% 8). Run the 'checkMassChargeBalance' function to obtain the list of imbalancedMass reactions:

RxnBalancingList = VerifyMassChargeBalance(model);

% Remove reactions which are only charge imbalanced:
RxnBalancingList = RxnBalancingList(~cellfun(@isempty,RxnBalancingList(:,2)),:); 

% Remove reactions imbalanced in atoms different to H and O:
IdxAtoms = contains(RxnBalancingList(:,2), {'C','S','N','Mg','Fe', 'X'});
RxnBalancingList = RxnBalancingList(IdxAtoms == 0,:);

RxnBalancingList = RxnBalancingList(:,1:2);

%% 9). Run 'AdjustUnbalancedAtoms' function to balance H- and O-atoms in the model:

ModelBalanced = AdjustUnbalancedAtoms(model,RxnBalancingList);

%%
ModelAdjusted = ModelBalanced;

end