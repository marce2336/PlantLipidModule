function RxnStoichiometry = VerifyUnbalancedRxns(UnbalancedRxnsList,model)
%==========================================================================
% 'VerifyUnbalancedRxns' function retrieves all the information for the
% metabolites participating in reactions which are mass-
% and/or-charge-imbalanced.
%
% USAGE:
%
%       RxnStoichiometry = VerifyUnbalancedRxns(UnbalancedRxnsList,model)
%
% INPUTS:       
%       model:              COBRA model structure containing the AraCore 
%                           model published by Arnold and Nikoloski (2014)
%
%       UnbalancedRxnsList: List containing the reactions which are mass-
%                           and/or charge-imbalanced, and the
%                           corresponding details. 
%
% OUTPUTS:      
%       RxnStoichiometry:   List containing the information for each 
%                           metabolite participating in a mass- and/or
%                           charge-imbalanced reaction, including:
%                           metabolites ID, formula, charge, KEGGID, and
%                           stoichiometric coefficient.
%==========================================================================
%% Retrieve metabolites information

ChargesList = cellstr(string(model.metCharges));
Metslist = [model.mets, model.metNames, model.metFormulas, ChargesList, model.metKEGGID];

%% Extract stoichiometries for each Rxn

RxnFormula = cell(length(UnbalancedRxnsList),2);

for i = 1:size(UnbalancedRxnsList,1)
    IDrxnUnb = cellstr(UnbalancedRxnsList{i,1});
    Rxn = printRxnFormula(model, 'rxnAbbrList', IDrxnUnb);
    RxnFormula(i,1) = IDrxnUnb;
    RxnFormula(i,2) = Rxn;
end

% Eliminate empty cells:
RxnFormula = RxnFormula(~cellfun(@isempty, RxnFormula(:,1)),:);

 % Split Rxns
RxnStoichiometry = cell(length(UnbalancedRxnsList)*13+length(UnbalancedRxnsList),12);
Count1 = 0;
for ii = 1:size(RxnFormula,1)
    Formula = RxnFormula(ii,2);
    Count1 = Count1 + 1;
    RxnStoichiometry(Count1,1:2) = RxnFormula(ii,:);
    
    %Split formula into reactants and products:
    SplitForm = split(Formula, '>');
    
    %Split reactants into individual metabolites:
    SplitRctants = split(SplitForm(1,1),'+');
    
    %Get reaction ID
    RxnID = find(strcmp(model.rxns,RxnFormula(ii,1)));
    
    for j = 1:size(SplitRctants,1)
        Reactant_i = extractBetween(SplitRctants{j,1}, 1, ']');
        
        if ~isempty(Reactant_i{:})
            %Eliminate empty spaces:
            Reactant_i       = strrep(Reactant_i, ' ', '');
            
            % Find position first letter and extract metabolite ID:
            TF = isstrprop(Reactant_i,'alpha');
            IsLetter = find(TF{:});
            IdxName = IsLetter(1,1);
            GetReactantID    = extractAfter(Reactant_i, IdxName-1);
       
            %Find metabolite:
            IdxMets = strcmp(Metslist(:,1),[GetReactantID{:},']']);
            
            % Obtain metabolite information:
            MetsData = Metslist(IdxMets,:);
            RxnStoichiometry(Count1+j,2:6) = cellstr(MetsData);
            
            %Get metabolite stoichiometric coefficient:
            SmetIdx = strcmp(model.mets,[GetReactantID{:},']']);
            Scoeff = full(model.S(SmetIdx,RxnID));
            RxnStoichiometry(Count1+j,1) = cellstr(string(Scoeff));
                        
        end
    end
    
    %Split products into individual metabolites:
    SplitPducts = split(SplitForm(2,1),'+');
    
    for j = 1:size(SplitPducts,1)
        Product_i = extractBetween(SplitPducts{j,1}, 1, ']');
        
        if ~isempty(Product_i{:})
            %Eliminate empty spaces:
            Product_i = strrep(Product_i, ' ', '');
            
            % Find position first letter and extract metabolite ID:
            TF = isstrprop(Product_i,'alpha');
            IsLetter = find(TF{:});
            IdxName = IsLetter(1,1);
            GetProductID    = extractAfter(Product_i, IdxName-1);
            
            %Find metabolite:
            IdxMets = strcmp(Metslist(:,1),[GetProductID{:},']']);
            
            % Obtain metabolite information:
            MetsData = cellstr(Metslist(IdxMets,:));
            RxnStoichiometry(Count1+j,8:12) = cellstr(MetsData);
            
            %Get metabolite stoichiometric coefficient
            SmetIdx = strcmp(model.mets,[GetProductID{:},']']);
            Scoeff = full(model.S(SmetIdx,RxnID));
            RxnStoichiometry(Count1+j,7) = cellstr(string(Scoeff));
        end
    end
    Count1 = Count1 + 5;
end

RxnStoichiometry = splitvars(table(RxnStoichiometry));
RxnStoichiometry.Properties.VariableNames = {'ReactStoichCoeff','ReactantsID',...
    'ReactantsName','ReactantsFormula','ReactantsCharge','ReactantsKEGG','ProdStoichCoeff',...
    'ProductsID','ProductsName','ProductsFormula','ProductsCharge','ProductsKEGG'};

writetable(RxnStoichiometry, 'RxnStoichiometry.txt');
pathDestination = fullfile('OutputFiles');
movefile('RxnStoichiometry.txt', pathDestination)

end
