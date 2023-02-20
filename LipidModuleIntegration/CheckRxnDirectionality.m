function RxnsToAdjust = CheckRxnDirectionality(RxnList1, RxnList2, model)
%% ========================================================================
% CheckRxnDirectionality is run after the usage of the function
% "checkDuplicateRxn" of the CobraToolbox. When non-unique reactions are
% identified and removed from the model, the function compares the location
% of reactants and products for each kept and removed reaction pair. When 
% differences are identified in the directionality of the reaction pair the 
% corresponding reaction ID is included in a list which later is used to
% decide whether or not it is necessary to adjust reaction direction.
%
% USAGE:
%
%    RxnsToAdjust = CheckRxnDirectionality(RxnList1, RnxList2, model)
%
% INPUTS:
%    RxnList1:         List with the indexes of the duplicated
%                      reactions kept in the model.
%                           
%    RnxList2:         List with the indexes of the duplicated
%                      reactions removed from the model.
%
%    model:            Standard COBRA model structure that contains
%                      non-unique reactions.
%
% OUTPUTS:             RxnsToAdjust: list of reactions to which
%                      directionality must be adjusted.
%
% =========================================================================
%%
RxnsToAdjust = cell(size(RxnList1,1),1);

for i = 1:size(RxnList1,1)
    IdKeptRxn_i = model.rxns{RxnList1(i)};
    IdRmvdRxn_i = model.rxns{RxnList2(i)};
    
    % Obtain reaction formulas:
    FormulaKeptRxn_i = printRxnFormula(model, IdKeptRxn_i);
    FormulaKeptRxn_i = strsplit(string(FormulaKeptRxn_i), '>');
    
    % Check that reaction is not exchange rxn:
    if strlength(FormulaKeptRxn_i(1,1)) <= 2
        FormulaKeptRxn_i(1,1) = extractBetween(FormulaKeptRxn_i(1,1), 1, strlength(FormulaKeptRxn_i(1,1)));
    else
        FormulaKeptRxn_i(1,1) = extractBefore(FormulaKeptRxn_i(1,1), strlength(FormulaKeptRxn_i(1,1))-2);
    end
    
    FormulaIdRmvdRxn_i = printRxnFormula(model, IdRmvdRxn_i);
    FormulaIdRmvdRxn_i = strsplit(string(FormulaIdRmvdRxn_i), '>');
    
    % Check that reaction is not exchange rxn:
    if strlength(FormulaIdRmvdRxn_i(1,1)) <= 2
        FormulaIdRmvdRxn_i(1,1) = extractBetween(FormulaIdRmvdRxn_i(1,1), 1, strlength(FormulaIdRmvdRxn_i(1,1)));
    else
        FormulaIdRmvdRxn_i(1,1) = extractBefore(FormulaIdRmvdRxn_i(1,1), strlength(FormulaIdRmvdRxn_i(1,1))-2);
    end
    
    % Check metabolites orientation:
    LocateReactants = contains(FormulaKeptRxn_i(1,1), FormulaIdRmvdRxn_i(1,1));
    LocateProducts = contains(FormulaKeptRxn_i(1,2), FormulaIdRmvdRxn_i(1,2));
    
    if LocateReactants == 1 && LocateProducts == 0
        LocateProducts = contains(FormulaIdRmvdRxn_i(1,2), FormulaKeptRxn_i(1,2));
    elseif LocateReactants == 0 && LocateProducts == 1
        LocateReactants = contains(FormulaIdRmvdRxn_i(1,1), FormulaKeptRxn_i(1,1));
    end
    
    % Get lower bounds of reactions:
    LBKeptRxn_i = model.lb(RxnList1(i));
    LBRmvdRxn_i = model.lb(RxnList2(i));
    
   % If metabolites are located in the same side of the reactions, but the
   % reaction directionality is different add to the list:
    if LocateReactants == 1 && LocateProducts == 1 && LBKeptRxn_i ~= LBRmvdRxn_i
        RxnsToAdjust{i,1}= IdKeptRxn_i;
    
    % If metabolites are located in opposite sides of the reaction and the
    % reactions are irreversible, add to the list:    
    elseif LocateReactants ~= 1 && LocateProducts ~= 1 && LBKeptRxn_i == 0 && LBRmvdRxn_i == 0
        RxnsToAdjust{i,1}= IdKeptRxn_i;
    
    % If metabolites are located in opposite sides of the reaction and the
    % reaction removed is irreversible, add to the list:
    elseif LocateReactants ~= 1 && LocateProducts ~= 1 && LBKeptRxn_i == -1000 && LBRmvdRxn_i == 0
        RxnsToAdjust{i,1}= IdKeptRxn_i;
        
    % If metabolites are located in opposite sides of the reaction and the
    % reaction kept is irreversible, add to the list:
    elseif LocateReactants ~= 1 && LocateProducts ~= 1 && LBKeptRxn_i == 0 && LBRmvdRxn_i == -1000
        RxnsToAdjust{i,1}= IdKeptRxn_i;
    end
    
end

RxnsToAdjust = RxnsToAdjust(~cellfun(@isempty, RxnsToAdjust(:,1)));

end
