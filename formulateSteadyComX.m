function modelCom = formulateSteadyComX(models, comGrowthRate, fluxBudget, metLossPercent, infValue)
%% formulateSteadyComX
% The code sources to run the SteadyComX simulations with RACs
% Kim et al. (2021) "Resource-allocation constraint governs structure and function of microbial communities in metabolic modeling" Submitted
%
% Generate a community model structure by combining individual metabolic models and impose resource-allocation constraints (i.e., flux budget constraints)
%
% Inputs
%   models          Cell array contains individual species models
%                   (e.g., models{1} = model_species_1, model{2} = model_species_2, and so on)
%   comGrowthRate   Community growth rate (default = 0.1 h-1)
%   fluxBudget      Intracellular flux budget (default = 100 mmol/g DCW/h)
%   metLossPercent  Minimum percentage of metabolites that are lost immediately after they are secreted by a microbe (default = 0 %)
%   infValue        Upper and lower limits of individual metabolic fluxes (default = 1000 mmol/g DCW/h)
%
% Outputs
%   modelCom        SteadyComX model in COBRA model structure
%                   Compatible w/ basic COBRA Toolbox functions (e.g., changeRxnBounds, optimizeCbModel, and so on)
%                   Note: There are some non-standard fields in the modelCom structure for downstream analyses

%% Process inputs
if nargin < 2
    comGrowthRate = 0.1;
end

if nargin < 3
    fluxBudget = 100;
end

if nargin < 4
    metLossPercent = 0;
end

if nargin < 5
    infValue = 1000;
end

%% Initialize modelCom structure
% Initialize modelCom structure (LP)
modelCom.S = [];
modelCom.b = [];
modelCom.c = [];
modelCom.lb = [];
modelCom.ub = [];
modelCom.csense = '';
modelCom.osense = -1;

% Initialize modelCom structure (row/col info)
modelCom.gr_pos = [];
modelCom.bg_pos = [];
modelCom.mets = {};
modelCom.metOrigNames = {};
modelCom.rxns = {};
modelCom.rxnOrigNames = {};

% Initialize modelCom structure (exchange met/rxn info)
modelCom.SEx = [];
modelCom.metExs = {};

%% Read individual models and add to the modelCom structure
for i = 1:length(models)
    model = models{i};
    model = splitExcRxns(model);
    [nMets, nRxns] = size(model.S);
    
    % Polish the model
    model.lb(model.lb == -1000) = -infValue;
    model.lb(model.lb ==  1000) =  infValue;
    model.ub(model.ub == -1000) = -infValue;
    model.ub(model.ub ==  1000) =  infValue;
    
    % Identify internal metabolites and reactions
    metEx = strcmp(getCompartment(model.mets),'e');
    metPp = strcmp(getCompartment(model.mets),'p');
    metInt = ~logical(metEx + metPp);
    SConnectivity = model.S;
    SConnectivity(SConnectivity ~= 0) = 1;
    rxnConnectivityMetAll = full(sum(SConnectivity))';
    SConnectivity(~metInt,:) = 0;
    rxnConnectivityMetInt = full(sum(SConnectivity))';
    rxnInt = find((rxnConnectivityMetAll - rxnConnectivityMetInt) == 0);
    
    % Identify exchange metabolites and reactions
    metEx = find(metEx);
    rxnExAll = find(sum(model.S ~= 0, 1) == 1)';
    [rxnEx, metExFound] = find(model.S(metEx, rxnExAll)');
    metEx = metEx(unique(metExFound));
    rxnEx = rxnExAll(rxnEx);
    model.mets(metEx) = strrep(model.mets(metEx), '[e]', '[u]');
    model.rxns(rxnEx) = strrep(model.rxns(rxnEx), '(e)', '(u)');
    
    % Identify reversible and irreversible reactions
    rxnIrrF = find(model.lb >= 0);
    rxnIrrR = find(model.ub <= 0);
    rxnBlck = intersect(rxnIrrF, rxnIrrR);
    rxnIrrF = setdiff(rxnIrrF, rxnBlck);
    rxnIrrR = setdiff(rxnIrrR, rxnBlck);
    rxnRev  = 1:length(model.rxns);
    rxnRev  = setdiff(rxnRev, union(rxnIrrF, rxnIrrR))';
    
    % Define budget constraint
    rxnMatExp1 = intersect(rxnInt, rxnRev);
    budgetRow = sparse(1,nRxns);
    budgetRow(intersect(rxnInt, rxnIrrF)) =  1;
    budgetRow(intersect(rxnInt, rxnIrrR)) = -1;
    budgetRow = [budgetRow, sparse(ones(1,2*length(rxnMatExp1)))];
    budgetRow = [-budgetRow, fluxBudget];

    % Expand the model to incorporate the budget constraint
    matSparse = sparse(nMets,nRxns);
    matSpeye  = speye(nRxns,nRxns);
    modelExp1.S = [model.S, matSparse(:,rxnMatExp1), matSparse(:,rxnMatExp1), sparse(nMets,1);
        matSpeye(rxnMatExp1,:), -matSpeye(rxnMatExp1,rxnMatExp1), matSpeye(rxnMatExp1,rxnMatExp1), sparse(length(rxnMatExp1),1);
        budgetRow];
    modelExp1.mets = [model.mets; strcat('AbsFlux_', model.rxns(rxnMatExp1)); 'Budget'];
    modelExp1.mets = modelExp1.mets(1:size(modelExp1.S,1));
    modelExp1.rxns = [model.rxns; strcat(model.rxns(rxnMatExp1), '_F'); strcat(model.rxns(rxnMatExp1), '_R'); 'BiomassVar'];
    lb = - infValue * ones(nRxns,1);
    lb(model.lb >= 0) = 0;
    modelExp1.lb = [lb; zeros(2*length(rxnMatExp1),1); 0];
    ub = infValue * ones(nRxns,1);
    ub(model.ub <= 0) = 0;
    modelExp1.ub = [ub; infValue * ones(2*length(rxnMatExp1),1); infValue];
    modelExp1.c = zeros(length(modelExp1.rxns),1);
    modelExp1.c(end) = 1;
    modelExp1.b = [model.b; zeros(length(rxnMatExp1),1); 0];
    modelExp1.b = modelExp1.b(1:size(modelExp1.S,1));
    modelExp1.csense = [model.csense; repmat('E', length(rxnMatExp1), 1); 'G'];
    modelExp1.csense = modelExp1.csense(1:size(modelExp1.S,1));
    bg_pos = sparse(nMets+length(rxnMatExp1)+1, size(modelExp1.S,2));
    bg_pos(end,end) = 1;
    bg_pos = bg_pos(1:size(modelExp1.S,1),:);
    
    % Reduce problem size by allowing SBCs for non-growing species if requested
    % Note: this part is related to the SteadyCom constraints (LB * X <= V <= UB * X)
    rxnMatExp2_1 = find(((model.lb > -infValue) .* (model.lb < infValue)) .* (model.lb ~= 0));
    rxnMatExp2_2 = intersect(rxnEx, find((model.lb == -infValue) + (model.lb == infValue)));
    rxnMatExp2_LB = union(rxnMatExp2_1, rxnMatExp2_2);
    rxnMatExp2_1 = find(((model.ub > -infValue) .* (model.ub < infValue)) .* (model.ub ~= 0));
    rxnMatExp2_2 = intersect(rxnEx, find((model.ub == -infValue) + (model.ub == infValue)));
    rxnMatExp2_UB = union(rxnMatExp2_1, rxnMatExp2_2);
    
    % Expand the model to incorporate the SteadyCom mass balance constraints
    % Step 1: introduce new rows for the SteadyCom constraints (LB * X <= V <= UB * X)
    modelExp2.S_LB = [-speye(nRxns,nRxns), sparse(nRxns,2*length(rxnMatExp1)), sparse(model.lb)];
    modelExp2.S_LB = modelExp2.S_LB(rxnMatExp2_LB,:);
    modelExp2.S_UB = [speye(nRxns,nRxns), sparse(nRxns,2*length(rxnMatExp1)), -sparse(model.ub)];
    modelExp2.S_UB = modelExp2.S_UB(rxnMatExp2_UB,:);
    modelExp2.S = [modelExp2.S_LB; modelExp2.S_UB];
    modelExp2.mets_LB = strcat('LB_', model.rxns);
    modelExp2.mets_LB = modelExp2.mets_LB(rxnMatExp2_LB,:);
    modelExp2.mets_UB = strcat('UB_', model.rxns);
    modelExp2.mets_UB = modelExp2.mets_UB(rxnMatExp2_UB,:);
    modelExp2.mets = [modelExp2.mets_LB; modelExp2.mets_UB];
    modelExp2.b = zeros(length(modelExp2.mets),1);
    modelExp2.csense = repmat('L', length(modelExp2.mets), 1);
    % Step 2: introduce a new row for the SteadyCom constraint (V_biomass = D * X)
    modelExp2.S = [modelExp2.S; sparse(1,size(modelExp2.S,2))];
    modelExp2.S(end, model.c == 1) = -1;
    modelExp2.S(end,end) = comGrowthRate;
    modelExp2.mets = [modelExp2.mets; 'Biomass'];
    modelExp2.b = [modelExp2.b; 0];
    modelExp2.csense = [modelExp2.csense; 'E'];
    % Step 3: merge with the existing model structure
    modelExp1.S = [modelExp1.S; modelExp2.S];
    modelExp1.mets = [modelExp1.mets; modelExp2.mets];
    modelExp1.b = [modelExp1.b; modelExp2.b];
    modelExp1.csense = [modelExp1.csense; modelExp2.csense];
    bg_pos = [bg_pos; sparse(size(modelExp2.S,1), size(modelExp2.S,2))];
    gr_pos = sparse(size(modelExp1.S,1), size(modelExp1.S,2));
    gr_pos(end,end) = 1;
    
    % Add the model to the modelCom structure
    % Step 1: expand LP matrix
    [nMets, nRxns] = size(modelExp1.S);
    [nRows, nCols] = size(modelCom.S);
    modelCom.S = [modelCom.S, sparse(nRows, nRxns);
        sparse(nMets, nCols), modelExp1.S];
    % Step 2: generate vectors for LP (b, c, lb, ub, csense)
    modelCom.b = [modelCom.b; modelExp1.b];
    modelCom.c = [modelCom.c; modelExp1.c];
    modelCom.lb = [modelCom.lb; modelExp1.lb];
    modelCom.ub = [modelCom.ub; modelExp1.ub];
    modelCom.csense = [modelCom.csense; modelExp1.csense];
    % Step 3: annotate rows/columns
    modelCom.mets = [modelCom.mets; strcat('M', num2str(i), '_', modelExp1.mets)];
    modelCom.metOrigNames = [modelCom.metOrigNames; modelExp1.mets];
    modelCom.rxns = [modelCom.rxns; strcat('M', num2str(i), '_', modelExp1.rxns)];
    modelCom.rxnOrigNames = [modelCom.rxnOrigNames; modelExp1.rxns];
    % Step 4: track the positions of comGrowthRate and fluxBudget values in S matrix
    modelCom.gr_pos = [modelCom.gr_pos, sparse(nRows, nRxns);
        sparse(nMets, nCols), gr_pos];
    modelCom.bg_pos = [modelCom.bg_pos, sparse(nRows, nRxns);
        sparse(nMets, nCols), bg_pos];
    
    % Set mass balance equations in common compartment
    model.mets(metEx) = strrep(model.mets(metEx), '[u]', '[e]');
    metExAbbr = setdiff(model.mets(metEx), modelCom.metExs);
    [nRows, nCols] = size(modelCom.SEx);
    modelCom.SEx = [modelCom.SEx, sparse(nRows, size(modelExp1.S,2));
        sparse(length(metExAbbr), nCols + size(modelExp1.S,2))];
    modelCom.metExs = [modelCom.metExs; metExAbbr];
    for j = 1:length(metEx)
        SEx_row_id = find(strcmp(model.mets(metEx(j)), modelCom.metExs));
        SEx_col_ids = rxnEx(find(model.S(metEx(j), rxnEx)'));
        if length(SEx_col_ids) == 1
            SEx_col_id = SEx_col_ids(1) + nCols;
            modelCom.SEx(SEx_row_id, SEx_col_id) = i;
        elseif length(SEx_col_ids) == 2
            for k = 1:length(SEx_col_ids)
                SEx_col_id = SEx_col_ids(k) + nCols;
                rxnExName = model.rxns{SEx_col_ids(k)};
                if strcmp(rxnExName(end), 'F')
                    modelCom.SEx(SEx_row_id, SEx_col_id) =  i;
                elseif strcmp(rxnExName(end), 'R')
                    modelCom.SEx(SEx_row_id, SEx_col_id) = -i;
                end
            end
        end
    end
end

%% Finalize the modelCom structure
% Sort exchange metabolites in alphabetical order
[modelCom.metExs, metExs_order] = sort(modelCom.metExs);
modelCom.SEx = modelCom.SEx(metExs_order,:);

% Generate exchange metabolite-model map
metExsToModels_F = zeros(length(modelCom.metExs), length(models));
metExsToModels_R = zeros(length(modelCom.metExs), length(models));
for i = 1:length(modelCom.metExs)
    rxnEx_row_id = find(modelCom.SEx(i,:));
    for j = 1:length(rxnEx_row_id)
        if modelCom.SEx(i,rxnEx_row_id(j)) > 0
            metExsToModels_F(i,  modelCom.SEx(i,rxnEx_row_id(j))) = rxnEx_row_id(j);
        elseif modelCom.SEx(i,rxnEx_row_id(j)) < 0
            metExsToModels_R(i, -modelCom.SEx(i,rxnEx_row_id(j))) = rxnEx_row_id(j);
        end
    end
end
modelCom.metExsToModels_F = metExsToModels_F;
modelCom.metExsToModels_R = metExsToModels_R;
modelCom.SEx = (1 - metLossPercent/100) * double(modelCom.SEx > 0) - double(modelCom.SEx < 0);

% Finalize the modelCom structure
modelCom.S = [modelCom.S, sparse(size(modelCom.S,1), length(modelCom.metExs));
    modelCom.SEx, -speye(length(modelCom.metExs), length(modelCom.metExs))];
modelCom.b = [modelCom.b; zeros(length(modelCom.metExs), 1)];
modelCom.c = [modelCom.c; zeros(length(modelCom.metExs), 1)];
modelCom.lb = [modelCom.lb; zeros(length(modelCom.metExs), 1)];
modelCom.ub = [modelCom.ub; infValue * ones(length(modelCom.metExs), 1);];
modelCom.csense = [modelCom.csense; repmat('E',length(modelCom.metExs),1)];
modelCom.mets = [modelCom.mets; strcat('C_', strrep(modelCom.metExs, '[e]', '[u]'))];
modelCom.metOrigNames = [modelCom.metOrigNames; strrep(modelCom.metExs, '[e]', '[u]')];
modelCom.rxns = [modelCom.rxns; strcat('EX_', strrep(modelCom.metExs, '[e]', '(e)'))];
modelCom.rxnOrigNames = [modelCom.rxnOrigNames; strcat('EX_', strrep(modelCom.metExs, '[e]', '(e)'))];
modelCom.gr_pos = [modelCom.gr_pos, sparse(size(modelCom.gr_pos,1), length(modelCom.metExs));
    sparse(length(modelCom.metExs), size(modelCom.gr_pos,2) + length(modelCom.metExs))];
modelCom.bg_pos = [modelCom.bg_pos, sparse(size(modelCom.bg_pos,1), length(modelCom.metExs));
    sparse(length(modelCom.metExs), size(modelCom.bg_pos,2) + length(modelCom.metExs))];
modelCom = rmfield(modelCom, 'SEx');
modelCom.comGrowthRate  = comGrowthRate;
modelCom.fluxBudget     = fluxBudget;
modelCom.metLossPercent = metLossPercent;
modelCom.infValue       = infValue;
modelCom = orderfields(modelCom,{'rxns', 'rxnOrigNames', 'mets', 'metOrigNames', ...
    'S', 'b', 'c', 'lb', 'ub', 'csense', 'osense', 'gr_pos', 'bg_pos', ...
    'metExs', 'metExsToModels_F', 'metExsToModels_R', 'comGrowthRate', ...
    'fluxBudget', 'metLossPercent', 'infValue'});

end

function model = splitExcRxns(model)
    % Identify exchange metabolites and reactions
    metEx = find(strcmp(getCompartment(model.mets),'e'));
    rxnExAll = find(sum(model.S ~= 0, 1) == 1)';
    [rxnEx, ~] = find(model.S(metEx, rxnExAll)');
    rxnEx = sort(rxnExAll(rxnEx));
    
    % Split exchange reactions into forward and reverse irreversible reactions
    for i = 1:length(rxnEx)
        rxnId = rxnEx(i);
        model.S = [model.S(:,1:rxnId), -model.S(:,rxnId), model.S(:,(rxnId+1):end)];
        model.rxns = [model.rxns(1:(rxnId-1)); strcat(model.rxns(rxnId), '_F'); ...
            strcat(model.rxns(rxnId), '_R'); model.rxns((rxnId+1):end)];
        lb = model.lb(rxnId);
        ub = model.ub(rxnId);
        if lb >= 0
            lbF =  lb; ubF =  ub;
            lbR =   0; ubR =   0;
        elseif ub <= 0
            lbF =   0; ubF =   0;
            lbR = -ub; ubR = -lb;
        else
            lbF =   0; ubF =  ub;
            lbR =   0; ubR = -lb;
        end
        model.lb = [model.lb(1:(rxnId-1)); lbF; lbR; model.lb((rxnId+1):end)];
        model.ub = [model.ub(1:(rxnId-1)); ubF; ubR; model.ub((rxnId+1):end)];
        model.c = [model.c(1:rxnId); -model.c(rxnId); model.c((rxnId+1):end)];
        rxnEx = rxnEx + 1;
    end
end