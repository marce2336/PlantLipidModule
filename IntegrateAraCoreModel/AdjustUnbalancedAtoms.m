function ModelBalanced = AdjustUnbalancedAtoms(model,RxnBalancing)
%==========================================================================
% 'AdjustUnbalancedAtoms' function adjustes the reactions which present 
% mass imbalance due to H- or O-atoms, by correcting the stoichiometric
% coefficients for these elements.
%
% USAGE:
%       ModelBalanced = AdjustUnbalancedAtoms(model,RxnBalancing)
%
% INPUTS:
%    
%       model:          Standard COBRA model structure obtained after
%                       running the "LipidModuleIntegration" function
%
%       RxnBalancing:   List with the H- and O-atoms-imbalanced reactions
%                       obtained after running the "checkMassChargeBalance"
%                       function of the CobraToolbox. The list does not
%                       include the reactions naturally imbalanced:
%                       exchange, import, siphon, demand, biomass, SLIME
%                       reactions and reactions where metabolites without
%                       known molecular formula are participating.
%
% OUTPUTS:                  
%       ModelBalanced:  COBRA model structure were the mass- and charge- 
%                       imbalanced reactions included in the RxnBalancing 
%                       list are corrected accordingly.
%==========================================================================
%% 1). Balance mass of reactions

for j = 1:size(RxnBalancing,1)
    ImRxnID = RxnBalancing(j,1);
    ImRxnFormula = printRxnFormula(model, 'rxnAbbrList', cellstr(ImRxnID));
    SplitFormula = split(ImRxnFormula, '>');
    ImAtoms = RxnBalancing(j,2); 
    NbAtoms = strlength(ImAtoms);
    Formula = "";
    
    % Check if imbalanced atoms correspond to H2O:
    if NbAtoms > 6
        AtomH = extractBetween(ImAtoms, strfind(ImAtoms,"H"), strfind(ImAtoms,"H"));
        AtomO2 = extractBetween(ImAtoms, strfind(ImAtoms,"O"), strfind(ImAtoms,"O"));
        ListAtoms = [AtomH, AtomO2; "H[", "O2["]';
        for viii = 1:length(ListAtoms)
            Coeff = double(extractBetween(ImAtoms, strfind(ImAtoms,ListAtoms(viii,1))-3, ListAtoms(viii,1)));
            AtomH = extractBetween(ImAtoms, strfind(ImAtoms,ListAtoms(viii,1)), strfind(ImAtoms,ListAtoms(viii,1)));
            Formula = [char(Formula),char(AtomH),char(num2str(abs(Coeff)))];
        end
        if Formula == "H2O1"
            ListAtoms = [Formula,"H2O["]; %ID of water
        end
        
    % Check if imbalanced atoms correspond to H:    
    else
        % Extract atom abbreviation:
        AtomH = extractBetween(ImAtoms, strfind(ImAtoms,"H"), strfind(ImAtoms,"H"));
        ListAtoms = [AtomH, "H["];
        
        %Indentify what side of the reaction is imbalanced
        ImAtoms = RxnBalancing(j,2);
        LenImAtoms = strlength(ImAtoms);

        % Imbalance is on the left hand side of the reaction
        if LenImAtoms <= 3
            Coeff = double(extractBetween(ImAtoms, strfind(ImAtoms,ListAtoms(1,1))-2, ListAtoms(1,1)));

        % Imbalance is on the right hand side of the reaction    
        else
            Coeff = double(extractBetween(ImAtoms, strfind(ImAtoms,ListAtoms(1,1))-3, ListAtoms(1,1)));
        end
    
        % First verify if individual H atom is present in Reactants side:
        CompsProducts = '';
        CompsReactants = '';
        CheckIsAtom = SearchExistAtoms(SplitFormula, 'R');

        % If individual H atom is not present in Reactants side, search in Products side:
        if sum(~cellfun(@isempty, CheckIsAtom)) == 0
            CheckIsAtom = SearchExistAtoms(SplitFormula, 'P');
            
            % If individual H atom is neither present in Products side, act
            % accordingly:
            if sum(~cellfun(@isempty, CheckIsAtom)) == 0
                CheckIsAtom = '';
            else
                CompsProducts = extractBetween(CheckIsAtom, '[', ']');
                CompsProducts = [ListAtoms{1,2},CompsProducts{1},']'];
            end
            
        else
            CompsReactants = extractBetween(CheckIsAtom, '[', ']');
            CompsReactants = [ListAtoms{1,2},CompsReactants{1},']'];
        end        
    end
    
    % Get index of imbalanced reaction
    ImRxnIdx = find(strcmp(model.rxns,ImRxnID));
    
%% CASE 1: The rxn is imbalanced in H molecules, but H+ as individual species is not present in the formula      
    if isempty(CheckIsAtom) && ListAtoms(1,2) == "H[" 
        
        %Find compartment location of reactants
        CompsReactants = unique(extractBetween(SplitFormula(1,1),'[',']'));
        CompsReactants = CompsReactants(~cellfun(@isempty, CompsReactants(:,1)),:);

        %Find compartment location of products
        CompsProducts = unique(extractBetween(SplitFormula(2,1),'[',']'));
        CompsProducts = CompsProducts(~cellfun(@isempty, CompsProducts(:,1)),:);  
            
        % If metabolites belonging to different compartments are
        % encountered on the same side of reaction act accordingly:
        if size(CompsReactants,1) > 1 || size(CompsProducts,1) > 1

            % Create list with species in reactants and products:
            ReactantsList = split(SplitFormula(1,1), '+');
            ReactantsList = extractBetween(ReactantsList, 1, ']');
            ReactantsList = strrep(ReactantsList, ' ', '');
            ProductsList = split(SplitFormula(2,1), '+');
            ProductsList = extractBetween(ProductsList, 1, ']');
            ProductsList = strrep(ProductsList, ' ', '');

            % Obtain KEGGIDs and molecular formulas for species in reactants:
            [KEGG_reactants, Formula_reactants] = GetKEGGIDs(ReactantsList, model);

            % Obtain KEGGIDs and molecular formulas for species in products:
            [KEGG_products, Formula_products] = GetKEGGIDs(ProductsList, model);

            % Obtain stoichiometric coefficients for H atoms:
            H_Reactants = extractBetween(Formula_reactants, 'H', 'O');
            H_Products = extractBetween(Formula_products, 'H', 'O');
            
            % Verify if empty cells correspond to a coefficient for H equal
            % to 1. If empty cell is encountered assign coefficient:
            CheckEmptyReactants = cellfun(@isempty, H_Reactants(:)).*contains(Formula_reactants, 'H') == 1;
            CheckEmptyProducts = cellfun(@isempty, H_Products(:)).*contains(Formula_products, 'H') == 1;
            
            if sum(CheckEmptyReactants) > 0
                H_Reactants{CheckEmptyReactants} = '1';
                
            elseif sum(CheckEmptyProducts) > 0
                H_Products{CheckEmptyProducts} = '1';
            end

            % Check if N atoms may be included:
            IsN_Reactants = (sum(strlength(H_Reactants)) > 3).*contains(H_Reactants(:), 'N') == 1;
            IsN_Products = (sum(strlength(H_Products)) > 3).*contains(H_Products(:), 'N') == 1;
            
            if sum(IsN_Reactants) >= 1 || sum(IsN_Products) >= 1
                H_Reactants(IsN_Reactants) = extractBefore(H_Reactants(IsN_Reactants), 'N');
                H_Products(IsN_Products) = extractBefore(H_Products(IsN_Products), 'N');
            end
            
            % Identify reaction type, whether it is a transport or a
            % different type:
            RxnType = sum(ismember(KEGG_reactants, KEGG_products));
            
            switch RxnType
                case 0 % All metabolites at both sides of the reaction have different KEGGIDs, hence is different from transport type
                    % Get the compartments for all the H species included in the model:
                    ListH_atoms = model.mets(strcmp(model.metKEGGID, 'C00080'));
                    ListH_atoms = extractBetween(ListH_atoms, '[', ']');
                    
                    if Coeff > 0 %There is surplus of H atoms at the right hand side of the rxn
                        % Obtain list of compartments for all the H atoms in model:
                        CommonComps = ismember(CompsReactants, ListH_atoms);
                        getCompartment = CompsReactants{CommonComps};% Get compartmen ID from reactant
                        
                        %Find index of H atom:
                        H_idx = strcmp(model.mets,['H','[',getCompartment,']']);

                        %Add H atoms in the reactants side of the reaction
                        model.S(H_idx,ImRxnIdx) = -Coeff;
                        
                    elseif Coeff < 0 %There is surpluss of H atoms at the left hand side of the reaction
                        CommonComps = ismember(CompsProducts, ListH_atoms);
                        getCompartment = CompsProducts{CommonComps}; % Get compartmen ID from product

                        %Find index of H atom:
                        H_idx = strcmp(model.mets,['H','[',getCompartment,']']);

                        %Add H atoms in the products side of the reaction
                        model.S(H_idx,ImRxnIdx) = abs(Coeff);
                    end
                    
                otherwise % The same KEGGIDs are encountered at both sides of the reaction, hence it is a transport type
                    % Find the metabolite pair with different number of H atoms:
                    for k = 1:size(ReactantsList,1)
                        KEGG_i = KEGG_reactants{k,1};

                        % Find metabolite on the products side of the reaction:
                        MatchingKEGG = strcmp(KEGG_products, KEGG_i);

                        % If matching metabolite is found, compare Number of H
                        % atoms
                        if sum(MatchingKEGG) > 0
                            Numb_H_reactants = H_Reactants{k,1};
                            Numb_H_products  = H_Products{MatchingKEGG};

                            if sum(strcmp(Numb_H_reactants, Numb_H_products)) == 0

                                % If number of H atoms is different, obtain compartment information 
                                % according to where is the H surplus:

                                if Coeff > 0 %There is surpluss of H atoms at the right hand side of the rxn
                                    getCompartment = CompsReactants{k,1};% Get compartmen ID from reactant

                                    %Find index of H atom:
                                    H_idx = strcmp(model.mets,['H','[',getCompartment,']']);

                                    %Add H atoms in the reactants side of the reaction
                                    model.S(H_idx,ImRxnIdx) = -Coeff;

                                elseif Coeff < 0 %There is surpluss of H atoms at the left hand side of the reaction
                                    getCompartment = CompsProducts{MatchingKEGG}; % Get compartmen ID from product

                                    %Find index of H atom:
                                    H_idx = strcmp(model.mets,['H','[',getCompartment,']']);

                                    %Add H atoms in the products side of the reaction
                                    model.S(H_idx,ImRxnIdx) = abs(Coeff);
                                break
                                end
                            end
                        end
                    end
            end
        
        % All metabolites in the same side of the reactions belongs to the
        % same compartment!
        else
            if Coeff > 0 %There is surpluss of H atoms at the right hand side of the rxn
                getCompartment = CompsReactants{:};% Get compartmen ID from reactant
                            
                %Find index of H atom:
                H_idx = strcmp(model.mets,['H','[',getCompartment,']']);

                %Add H atoms in the reactants side of the reaction
                model.S(H_idx,ImRxnIdx) = -Coeff;

            elseif Coeff < 0 %There is surpluss of H atoms at the left hand side of the reaction
                getCompartment = CompsProducts{:}; % Get compartmen ID from product

                %Find index of H atom:
                H_idx = strcmp(model.mets,['H','[',getCompartment,']']);

                %Add H atoms in the products side of the reaction
                model.S(H_idx,ImRxnIdx) = abs(Coeff);
            end
        end
                
%% CASE 2: The rxn is imbalanced in H atoms, and H+ is present as individual species    
    elseif ~isempty(CheckIsAtom) && ListAtoms(1,2) == "H[" 
        
        % Adjust coefficients to balance H atoms:
        if Coeff > 0 && ~isempty(CompsProducts) %There is surpluss of H atoms at the right hand side of the rxn
            Hidx = find(strcmp(model.mets,CompsProducts));
            CurrCoeff = model.S(Hidx,ImRxnIdx);
            model.S(Hidx,ImRxnIdx) = CurrCoeff-Coeff; %Add H atoms in the reactants side of the Rxn

        elseif Coeff > 0 && ~isempty(CompsReactants) %There is surpluss of H atoms at the left hand side of the rxn
            Hidx = find(strcmp(model.mets,CompsReactants));
            CurrCoeff = model.S(Hidx,ImRxnIdx);
            model.S(Hidx,ImRxnIdx) = CurrCoeff-Coeff; %Add H atoms in the reactants side of the Rxn

        elseif Coeff < 0 && ~isempty(CompsProducts) %There is surpluss of H atoms at the left hand side of the rxn
            Hidx = find(strcmp(model.mets,CompsProducts));
            CurrCoeff = model.S(Hidx,ImRxnIdx);
            model.S(Hidx,ImRxnIdx) = abs(Coeff) + CurrCoeff; %Add H atoms in the products side of the Rxn

        elseif Coeff < 0 && ~isempty(CompsReactants) %There is surpluss of H atoms at the left hand side of the rxn
            Hidx = find(strcmp(model.mets,CompsReactants));
            CurrCoeff = model.S(Hidx,ImRxnIdx);
            model.S(Hidx,ImRxnIdx) = abs(Coeff) + CurrCoeff; %Add H atoms in the products side of the Rxn
        end
        
%% CASE 3: The rxn is imbalanced in H2O molecules          
    elseif isempty(CheckIsAtom) && ListAtoms(1,2) == "H2O" 
        
        % Indentify which side of the reaction is imbalanced:
        ImAtoms = RxnBalancing(j,2);
        Coeff = double(extractBetween(ImAtoms, strfind(ImAtoms,"H")-3, strfind(ImAtoms,"H")-1));
        
        %Find metabolite location:
        CompsReactants = unique(extractBetween(SplitFormula(1,1),'[',']'));
        CompsProducts = unique(extractBetween(SplitFormula(2,1),'[',']'));
        
        if length(CompsReactants) > 1 || length(CompsProducts) > 1
            Xr = ~cellfun(@isempty,(strfind(CompsReactants,'0')));
            Xp = ~cellfun(@isempty,(strfind(CompsProducts,'0')));
            
            % Extract compartment abbreviation:
            txtR = cellstr(['[',char(CompsReactants(Xr,1)),']']);
            txtP = cellstr(['[',char(CompsProducts(Xp,1)),']']);
        
        else
            % Extract compartment abbreviation:
            txtR = cellstr(['[',char(CompsReactants),']']);
            txtP = cellstr(['[',char(CompsProducts),']']);
        end

        if ~isempty(txtR)
            ImMetIDR = cellstr(['H2O',char(txtR)]);
        end
        
        if ~isempty(txtP)
            ImMetIDP = cellstr(['H2O',char(txtP)]);
        end

        if Coeff > 0 % There is surplus of H2O molecules at the right hand side of the reaction
            %Find index of H2O molecules in model
            Hidx = find(strcmp(model.mets,ImMetIDR));
            CurrCoeff = model.S(Hidx,ImRxnIdx);
            
            % Add H2O molecules in the reactants side of the Rxn
            model.S(Hidx,ImRxnIdx) = CurrCoeff-(Coeff/2);

        elseif Coeff < 0 %There is surplus of H2O molecules at the left hand side of the reaction
            %Find Idx of H2O molecules
            Hidx = find(strcmp(model.mets,ImMetIDP));
            CurrCoeff = model.S(Hidx,ImRxnIdx);
            
            %Add H2O molecules in the products side of the reaction
            model.S(Hidx,ImRxnIdx) = (abs(Coeff)/2) + CurrCoeff;

        end
    end
end

%% Generate output variable:

ModelBalanced = model;

end