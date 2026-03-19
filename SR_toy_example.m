%% SR_general_v2.m
% General-purpose symbolic regression via genetic programming.
% v2 fixes: operator frequency weighting, ramped half-and-half init,
%           tighter protected operators, full-reset on repeated stagnation.
%
% ---- QUICK START --------------------------------------------------------
%  1. Set your data + VAR_NAMES
%  2. Trim OPS_WEIGHTS if you know some operators don't belong
%  3. Run — works for polynomial, trig, exponential, rational, mixed targets
% -------------------------------------------------------------------------

rng(42);

%% =========================================================
%% CONFIGURATION
%% =========================================================

% --- Data ----------------------------------------------------------------
X_data    = linspace(-3, 3, 60)';
VAR_NAMES = {'x'};
Y_data    = sin(X_data) + 0.5*X_data.^2 + 0.05*randn(size(X_data));

% Multi-variable example (uncomment):
% [X1,X2]  = meshgrid(linspace(-2,2,10), linspace(-2,2,10));
% X_data    = [X1(:), X2(:)];
% VAR_NAMES = {'x1','x2'};
% Y_data    = sin(X_data(:,1)) .* exp(-X_data(:,2).^2) + 0.05*randn(100,1);

% --- Operator set + weights ----------------------------------------------
% Weight controls how often each op appears in random trees.
% Arithmetic ops should dominate; exotic ops should be rare.
% Set weight=0 to disable an operator entirely.
OPS      = {'+',  '-',  '*',  '/',   ...
            'sin','cos','tan',        ...
            'exp','log','sqrt',       ...
            'pow2','abs','tanh'};
ARITY    = [ 2,    2,    2,    2,    ...
             1,    1,    1,          ...
             1,    1,    1,          ...
             1,    1,    1   ];
% Relative sampling weight per operator (need not sum to 1)
OPS_W    = [ 8,    8,    8,    4,    ...   % arithmetic: high
             5,    5,    1,          ...   % trig: medium, tan: low
             3,    2,    3,          ...   % exp/log/sqrt: low-medium
             4,    2,    3   ];            % pow2: medium, abs/tanh: low

% --- GP parameters -------------------------------------------------------
MAX_D      = 5;     % max tree depth
POP_SZ     = 600;   % population size
N_GEN      = 400;   % max generations
TOURN_K    = 7;     % tournament size
P_CROSS    = 0.80;  % crossover probability
P_MUT      = 0.20;  % mutation probability
TOP_K      = 40;    % individuals to constant-tune each generation
ALPHA      = 0.006; % parsimony pressure weight
STAG_LIM   = 20;    % gens without improvement → partial restart
FULL_RESET = 3;     % partial restarts before full population reset
RMSE_TOL   = 1e-4;  % stop if RMSE drops below this

%% =========================================================
%% INITIALISE — ramped half-and-half
%% =========================================================
% Vary tree depth across population for better initial diversity
pop       = cell(POP_SZ, 1);
depths    = repmat(2:MAX_D, 1, ceil(POP_SZ/(MAX_D-1)));
depths    = depths(randperm(numel(depths), POP_SZ));
for i = 1:POP_SZ
    if rand < 0.5
        pop{i} = grow_tree(depths(i), OPS, ARITY, OPS_W, VAR_NAMES);
    else
        pop{i} = full_tree(depths(i), OPS, ARITY, OPS_W, VAR_NAMES);
    end
end

best_fitness  = inf;
best_rmse     = inf;
best_tree     = [];
stagnation    = 0;
n_restarts    = 0;
history       = nan(N_GEN, 2);

%% =========================================================
%% MAIN LOOP
%% =========================================================
for gen = 1:N_GEN

    %% Constant-tune top-K -----------------------------------------
    fit_raw = arrayfun(@(i) rmse_tree(pop{i}, X_data, Y_data, VAR_NAMES), ...
                       1:POP_SZ);
    [~, sidx] = sort(fit_raw);
    for k = 1:min(TOP_K, POP_SZ)
        pop{sidx(k)} = tune_constants(pop{sidx(k)}, X_data, Y_data, VAR_NAMES);
    end

    %% Fitness with parsimony --------------------------------------
    fit = arrayfun(@(i) fitness(pop{i}, X_data, Y_data, VAR_NAMES, ALPHA), ...
                   1:POP_SZ);

    [f_best, idx] = min(fit);
    r_best        = rmse_tree(pop{idx}, X_data, Y_data, VAR_NAMES);

    if f_best < best_fitness
        best_fitness = f_best;
        best_rmse    = r_best;
        best_tree    = pop{idx};
        stagnation   = 0;
    else
        stagnation = stagnation + 1;
    end
    history(gen,:) = [best_fitness, best_rmse];

    fprintf('Gen %3d | fitness: %.5f | RMSE: %.5f | size: %2d | %s\n', ...
            gen, best_fitness, best_rmse, ...
            tree_size(best_tree), tree2str(best_tree, VAR_NAMES));

    %% Early stop --------------------------------------------------
    if best_rmse < RMSE_TOL
        fprintf('\n  >> RMSE tolerance reached — stopping.\n');
        break;
    end

    %% Restart logic -----------------------------------------------
    if stagnation >= STAG_LIM
        n_restarts = n_restarts + 1;
        if n_restarts >= FULL_RESET
            % Full reset — keep only the single best
            fprintf('  >> FULL RESET (restart #%d)\n', n_restarts);
            for i = 2:POP_SZ
                d = randi([2, MAX_D]);
                if rand < 0.5
                    pop{i} = grow_tree(d, OPS, ARITY, OPS_W, VAR_NAMES);
                else
                    pop{i} = full_tree(d, OPS, ARITY, OPS_W, VAR_NAMES);
                end
            end
            n_restarts = 0;
        else
            % Partial restart — reseed bottom 60%
            fprintf('  >> Partial restart #%d — reseeding 60%%\n', n_restarts);
            rs = floor(POP_SZ*0.4)+1;
            for i = rs:POP_SZ
                d = randi([2, MAX_D]);
                pop{i} = grow_tree(d, OPS, ARITY, OPS_W, VAR_NAMES);
            end
        end
        stagnation = 0;
    end

    %% Next generation ---------------------------------------------
    new_pop    = cell(POP_SZ, 1);
    new_pop{1} = best_tree;   % elitism

    for i = 2:POP_SZ
        p1 = tournament(pop, fit, TOURN_K);
        if rand < P_CROSS
            p2    = tournament(pop, fit, TOURN_K);
            child = crossover(p1, p2, MAX_D);
        else
            child = p1;
        end
        if rand < P_MUT
            child = mutate(child, MAX_D, OPS, ARITY, OPS_W, VAR_NAMES);
        end
        new_pop{i} = child;
    end
    pop = new_pop;
end

%% =========================================================
%% RESULTS
%% =========================================================
fprintf('\n================================================\n');
fprintf('Best expression : %s\n',    tree2str(best_tree, VAR_NAMES));
fprintf('RMSE            : %.6f\n',  best_rmse);
fprintf('Tree size       : %d\n',    tree_size(best_tree));
fprintf('================================================\n');

%% Plots
figure('Name','SR Results','Position',[80 80 1300 460]);
n_cols = 1 + (size(X_data,2) == 1);

if size(X_data, 2) == 1
    subplot(1, n_cols, 1);
    Yhat = eval_tree(best_tree, X_data, VAR_NAMES);
    plot(X_data, Y_data, 'k.',  'MarkerSize', 10, 'DisplayName', 'Noisy data'); hold on;
    plot(X_data, Yhat,   'r-',  'LineWidth', 2.5, ...
         'DisplayName', ['SR: ' tree2str(best_tree, VAR_NAMES)]);
    legend('Interpreter','none','Location','best');
    grid on; xlabel('x'); ylabel('y');
    title(sprintf('Best fit  (RMSE=%.4f)', best_rmse), 'Interpreter','none');
end

subplot(1, n_cols, n_cols);
gens_run = find(~isnan(history(:,2)), 1, 'last');
semilogy(1:gens_run, history(1:gens_run,2), 'b-',  'LineWidth',1.5); hold on;
semilogy(1:gens_run, history(1:gens_run,1), 'r--', 'LineWidth',1.0);
legend('RMSE','Fitness','Interpreter','none','Location','best');
grid on; xlabel('Generation'); ylabel('Error (log scale)');
title('Convergence history');


%% =========================================================
%% LOCAL FUNCTIONS
%% =========================================================

%% --- Tree generation -----------------------------------------------------

% GROW: choose terminal with increasing probability near leaves
function node = grow_tree(max_d, ops, arity, ops_w, var_names)
    p_term = max(0.1, 1 - max_d/5);
    if max_d == 0 || rand < p_term
        node = terminal_node(var_names);
    else
        k             = weighted_choice(ops_w);
        node.type     = 'op';
        node.op       = ops{k};
        node.ar       = arity(k);
        node.children = cell(arity(k),1);
        for c = 1:arity(k)
            node.children{c} = grow_tree(max_d-1, ops, arity, ops_w, var_names);
        end
    end
end

% FULL: all leaves at exactly max_d
function node = full_tree(max_d, ops, arity, ops_w, var_names)
    if max_d == 0
        node = terminal_node(var_names);
    else
        k             = weighted_choice(ops_w);
        node.type     = 'op';
        node.op       = ops{k};
        node.ar       = arity(k);
        node.children = cell(arity(k),1);
        for c = 1:arity(k)
            node.children{c} = full_tree(max_d-1, ops, arity, ops_w, var_names);
        end
    end
end

function k = weighted_choice(w)
    cdf = cumsum(w) / sum(w);
    k   = find(rand <= cdf, 1, 'first');
    if isempty(k), k = numel(w); end
end

function node = terminal_node(var_names)
    if rand < 0.55
        node.type = 'var';
        node.op   = var_names{randi(numel(var_names))};
    else
        node.type = 'const';
        node.op   = round(randn * 2, 3);
    end
    node.ar = 0; node.children = {};
end

%% --- Tree evaluation -----------------------------------------------------
function y = eval_tree(node, X, var_names)
    N = size(X,1);
    switch node.type
        case 'var'
            col = find(strcmp(var_names, node.op),1);
            y   = X(:,col);
        case 'const'
            y = repmat(double(node.op), N, 1);
        case 'op'
            c = cellfun(@(ch) eval_tree(ch,X,var_names), ...
                        node.children,'UniformOutput',false);
            switch node.op
                case '+',    y = c{1} + c{2};
                case '-',    y = c{1} - c{2};
                case '*',    y = c{1} .* c{2};
                case '/',    y = pdiv(c{1}, c{2});
                case 'sin',  y = sin(c{1});
                case 'cos',  y = cos(c{1});
                case 'tan',  y = ptan(c{1});
                case 'exp',  y = exp(clamp(c{1},-20,20));
                case 'log',  y = plog(c{1});
                case 'sqrt', y = sqrt(abs(c{1}));
                case 'pow2', y = c{1}.^2;
                case 'abs',  y = abs(c{1});
                case 'tanh', y = tanh(c{1});
                otherwise,   y = zeros(N,1);
            end
            y(~isfinite(y)) = 1e6;
    end
end

% Protected operators
function y = pdiv(a,b)
    y = a ./ b;
    y(abs(b) < 1e-8) = 1e6;
end
function y = plog(a)
    y = log(abs(a) + 1e-10);
end
function y = ptan(a)
    % Clamp strictly away from ±π/2 (1.5708) with margin
    a = clamp(a, -1.4, 1.4);
    y = tan(a);
    y(abs(y) > 50) = sign(y(abs(y)>50)) * 50;   % hard cap
end
function v = clamp(v, lo, hi)
    v = max(lo, min(hi, v));
end

%% --- Fitness -------------------------------------------------------------
function f = fitness(node, X, Y, var_names, alpha)
    r = rmse_tree(node, X, Y, var_names);
    f = r + alpha * tree_size(node);
end

function r = rmse_tree(node, X, Y, var_names)
    try
        Yhat = eval_tree(node, X, var_names);
        r    = sqrt(mean((Yhat - Y).^2));
        if ~isfinite(r), r = 1e9; end
    catch
        r = 1e9;
    end
end

%% --- Constant optimisation -----------------------------------------------
function node = tune_constants(node, X, Y, var_names)
    [vals, locs] = extract_constants(node);
    if isempty(vals), return; end
    obj  = @(c) eval_tree(insert_constants(node,locs,c), X, var_names) - Y;
    opts = optimoptions('lsqnonlin','Display','off', ...
                        'MaxIter',100,'FunctionTolerance',1e-7);
    try
        vo   = lsqnonlin(obj, vals, [], [], opts);
        node = insert_constants(node, locs, vo);
    catch, end
end

function [vals,locs] = extract_constants(node)
    node      = tag_nodes(node);
    all_nodes = list_nodes(node);
    vals = []; locs = [];
    for k = 1:numel(all_nodes)
        if strcmp(all_nodes{k}.type,'const')
            vals(end+1) = all_nodes{k}.op; %#ok<AGROW>
            locs(end+1) = all_nodes{k}.id; %#ok<AGROW>
        end
    end
end

function node = insert_constants(node,locs,vals)
    node = tag_nodes(node);
    for k = 1:numel(locs)
        node = set_const(node,locs(k),vals(k));
    end
end

function node = set_const(node,tid,val)
    if node.id == tid, node.op = val; return; end
    for c = 1:numel(node.children)
        node.children{c} = set_const(node.children{c},tid,val);
    end
end

%% --- Selection -----------------------------------------------------------
function w = tournament(pop,fit,k)
    idx     = randperm(numel(pop),k);
    [~,b]   = min(fit(idx));
    w       = pop{idx(b)};
end

%% --- Crossover -----------------------------------------------------------
function child = crossover(p1,p2,max_d)
    n1    = list_nodes(p1);
    n2    = list_nodes(p2);
    pt    = n1{randi(numel(n1))};
    sub   = n2{randi(numel(n2))};
    child = replace_node(p1,pt.id,sub);
    if tree_depth(child) > max_d, child = p1; end
end

%% --- Mutation ------------------------------------------------------------
function node = mutate(node,max_d,ops,arity,ops_w,var_names)
    nodes  = list_nodes(node);
    target = nodes{randi(numel(nodes))};
    d_left = max(0,max_d - target.depth);
    new_sub = grow_tree(d_left,ops,arity,ops_w,var_names);
    node    = replace_node(node,target.id,new_sub);
end

%% --- Tree utilities ------------------------------------------------------
function [node,counter] = tag_nodes(node,depth,counter)
    if nargin < 2, depth=0;    end
    if nargin < 3, counter={0}; end
    counter{1} = counter{1}+1;
    node.id    = counter{1};
    node.depth = depth;
    for c = 1:numel(node.children)
        [node.children{c},counter] = tag_nodes(node.children{c},depth+1,counter);
    end
end

function nodes = list_nodes(root)
    root  = tag_nodes(root);
    nodes = {};
    queue = {root};
    while ~isempty(queue)
        n=queue{1}; queue(1)=[];
        nodes{end+1}=n; %#ok<AGROW>
        queue=[queue, n.children(:)'];
    end
end

function node = replace_node(node,tid,rep)
    node = tag_nodes(node);
    node = do_replace(node,tid,rep);
end

function node = do_replace(node,tid,rep)
    if node.id==tid, node=rep; return; end
    for c = 1:numel(node.children)
        node.children{c} = do_replace(node.children{c},tid,rep);
    end
end

function d = tree_depth(node)
    if isempty(node.children), d=0; return; end
    d = 1 + max(cellfun(@tree_depth,node.children));
end

function n = tree_size(node)
    n = 1 + sum(cellfun(@tree_size,node.children));
end

function s = tree2str(node,var_names)
    switch node.type
        case 'var',   s = node.op;
        case 'const', s = num2str(node.op,'%.4g');
        case 'op'
            switch node.op
                case {'+','-','*','/'}
                    s = sprintf('(%s %s %s)', ...
                        tree2str(node.children{1},var_names), node.op, ...
                        tree2str(node.children{2},var_names));
                otherwise
                    s = sprintf('%s(%s)',node.op, ...
                        tree2str(node.children{1},var_names));
            end
    end
end