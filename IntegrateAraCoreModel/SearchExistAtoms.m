function CheckIsAtom = SearchExistAtoms(SplitFormula, RxnSide)
%==========================================================================
switch RxnSide
    case {'P'}
        MetForms   = split(SplitFormula{2}, '+');
        IdxMets     = contains(MetForms, 'H[');
        
    case {'R'}
        MetForms  = split(SplitFormula{1}, '+');
        IdxMets = contains(MetForms, 'H[');      
end
        
% Get metabolites with H atoms:
FindMets_w_H = MetForms(IdxMets > 0);

% Verify that H is not part of NADH/NADPH/GSH/CTH:
IsInCofactor = contains(FindMets_w_H, {'NAD','NADP','GS','CT'});

% Exclude metabolites were H is part of the molecule:
FindMets_w_H = FindMets_w_H(IsInCofactor == 0);

CheckIsAtom = FindMets_w_H;

end