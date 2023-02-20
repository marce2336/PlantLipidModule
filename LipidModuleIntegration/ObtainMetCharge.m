function MatchingMets = ObtainMetCharge(MatchingMets, LmNeutralForm)
%% ========================================================================
% The 'ObtainMetCharge' function searches for the charges of metabolites
% in the ModelTemplate, followed by obtaining the respective neutral 
% molecular formulas.
%
% USAGE: MatchingMets = ObtainMetCharge(MatchingMets, LmNeutralForm)
%
% INPUTS:   
%   MatchingMets:      List with equivalent metabolites encountered in ModelTemplate.
%
%   LmNeutralForm:     Neutral moelcular formula of LipidModule metabolite.
%
% OUTPUTS:
%   MatchingMets:      List containing the Neutral Formulas of the
%                      equivalent metabolites encountered in ModelTemplate
%                      after comparing the available model identifiers 
%                      against a defined LipidModule metabolite.
%
% =========================================================================                       
%% 
for i = 1:size(MatchingMets,1)
    
    GetCharge = MatchingMets(i,2);
    
    if ismissing(GetCharge)
        GetCharge = "NaN";
        TmCharge = str2double(GetCharge);
    else
        TmCharge = str2double(GetCharge);
    end
    
    % The metabolite is not charged, then retrieve neutral formula:
    if TmCharge == 0
        MatchingMets(i,4) = MatchingMets(i,3);

    % The metabolite has not charge assigned:    
    elseif isnan(TmCharge)
        MatchingMets(i,4) = MatchingMets(i,3); % Assume this is the neutral formula

    % The metabolite is charged, then retrieve charged formula:    
    else
        TmChargedForm = MatchingMets(i,3);

        % If charged formula and neutral formula for ModelTemplate and
        % LipidModule mets exist,respectively, then calculate neutral formula:
        if ~isempty(TmChargedForm{1,1}) && ~isempty(LmNeutralForm)
            TmNeutralForm = ObtainNeutralFormula(TmChargedForm, LmNeutralForm, TmCharge);
            MatchingMets(i,4) = cellstr(TmNeutralForm);
        end
    end
end

end