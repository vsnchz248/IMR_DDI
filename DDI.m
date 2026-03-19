%% ============================================================
%  DDI Stress Reconstruction for N datasets (RAW + Padé derivatives)
%  - NO smoothing
%  - Derivatives via 4th-order compact Padé (tridiagonal solve)
%  - Subplot: R(t) + (R,V) colored by S
% =============================================================

clear; clc; close all;

%% =========================
% LOAD DATA
% =========================
load('synthetic_data.mat');

t_data = cell(numel(synthetic_data),1);
R_data = cell(numel(synthetic_data),1);
for i = 1:numel(synthetic_data)
    t_data{i} = synthetic_data{i}(:,1);
    R_data{i} = synthetic_data{i}(:,2);
end
N = numel(R_data);

%% =========================
% PARAMETERS
% =========================
rho   = 1000;        % kg/m^3
Pinf  = 101325;      % Pa
gamma = 0.072;       % N/m
k     = 1;           % pressure evolution parameter (tunable)

% Clustering (optional)
do_clustering = true;
n_clusters    = 30;

%% =========================
% STORAGE
% =========================
S_data  = cell(N,1);
V_data  = cell(N,1);
pb_data = cell(N,1);

Z_all = zeros(0,3);  % [R, V, S] stacked across experiments

%% =========================
% LOOP OVER EXPERIMENTS
% =========================
for alpha = 1:N

    t = t_data{alpha}(:);
    R = R_data{alpha}(:);

    if numel(t) ~= numel(R)
        error('Dataset %d: t and R must have same length.', alpha);
    end

    % Assume (approximately) uniform time step
    dt_vec = diff(t);
    dt = median(dt_vec);
    if any(abs(dt_vec - dt) > 1e-10*max(1,abs(dt)))
        warning('Dataset %d: t is not perfectly uniform; using median dt = %.6g.', alpha, dt);
    end

    % -------------------------
    % RAW derivatives via Padé
    % -------------------------
    V = pade_first_derivative(R, dt);
    A = pade_first_derivative(V, dt);

    % Inertia term (incompressible RP LHS)
    inertia = R .* A + 1.5 .* (V.^2);

    % Raw pressure from RP assuming S = 0
    pb_raw = rho .* inertia + Pinf + 2*gamma ./ R;

    % Pressure evolution (simple model): dpb/dt = (3k/R) * V * pb
    pb_consistent = zeros(size(R));
    pb_consistent(1) = pb_raw(1);

    for n = 2:numel(R)
        rhs = (3*k / R(n)) * V(n) * pb_consistent(n-1);
        pb_consistent(n) = pb_consistent(n-1) + dt * rhs;  % forward Euler
    end

    % Reconstructed stress
    S = pb_raw - pb_consistent;

    % Store
    S_data{alpha}  = S;
    V_data{alpha}  = V;
    pb_data{alpha} = pb_consistent;

    % Stack phase space
    Z_all = [Z_all; R(:), V(:), S(:)];
end

%% =========================
% GLOBAL CLUSTERING (optional)
% =========================
if do_clustering
    [idx, centroids] = kmeans(Z_all, n_clusters, 'Replicates', 5);
else
    idx = [];
    centroids = [];
end

%% =========================
% Visualization: R(t) + (R,V) colored by S
% =========================
figure('Color','w','Position',[100 100 1100 420])

% ---- LEFT: R(t)
subplot(1,2,1); hold on
for alpha = 1:N
    plot(t_data{alpha}, R_data{alpha}, 'LineWidth', 1.2)
end
xlabel('Time')
ylabel('R')
title('Radius vs Time')
grid on; box on

% ---- RIGHT: (R,V) colored by S
subplot(1,2,2)

Rvals = Z_all(:,1);
Vvals = Z_all(:,2);
Svals = Z_all(:,3);

% Robust color limits to avoid outliers dominating
cl = prctile(Svals(isfinite(Svals)), [2 98]);
mask = isfinite(Svals) & (Svals >= cl(1)) & (Svals <= cl(2));

scatter(Rvals(mask), Vvals(mask), 10, Svals(mask), 'filled')
xlabel('R')
ylabel('V')
title('Phase Space (R,V) colored by S')
grid on; box on; axis tight

colormap(turbo)
cb = colorbar;
cb.Label.String = 'S';
caxis(cl)
%%
%% =========================
% Visualization: R(t) + 3D phase space scatter (R,V,t)
% =========================
figure('Color','w','Position',[100 100 1200 450])

% ---- LEFT: R(t)
subplot(1,2,1); hold on
for alpha = 1:N
    plot(t_data{alpha}, R_data{alpha}, 'LineWidth', 1.2)
end
xlabel('Time')
ylabel('R')
title('Radius vs Time')
grid on; box on

% ---- RIGHT: 3D scatter (R,V,t) across all datasets
subplot(1,2,2); hold on

colors = lines(N);   % distinct color per dataset

for alpha = 1:N
    t = t_data{alpha}(:);
    R = R_data{alpha}(:);
    V = V_data{alpha}(:);

    % basic finite checks
    m = isfinite(t) & isfinite(R) & isfinite(V);

    scatter3(R(m), V(m), t(m), 10, ...
        'MarkerFaceColor', colors(alpha,:), ...
        'MarkerEdgeColor', 'none', ...
        'MarkerFaceAlpha', 0.6);
end

xlabel('R')
ylabel('V')
zlabel('t')
title('3D Phase Space (R, V, t) - colored by dataset')
grid on; box on; axis tight
view(45,30)


%% ============================================================
% Local function: 4th-order compact Padé first derivative
% ============================================================
function df = pade_first_derivative(f, h)
% pade_first_derivative: 4th-order compact derivative on a (nearly) uniform grid
%
%   (1/4) f'_{i-1} + f'_i + (1/4) f'_{i+1} = (3/(4h)) (f_{i+1} - f_{i-1})
%
% Boundaries: 2nd-order one-sided finite differences.

f = f(:);
N = numel(f);

% Tridiagonal coefficients
a = (1/4) * ones(N,1);   % sub-diagonal
b = 1     * ones(N,1);   % main diagonal
c = (1/4) * ones(N,1);   % super-diagonal

rhs = zeros(N,1);

% Interior RHS
rhs(2:N-1) = (3/(4*h)) * (f(3:N) - f(1:N-2));

% Boundary derivatives (2nd-order one-sided)
df1 = (-3*f(1) + 4*f(2) - f(3)) / (2*h);
dfN = ( 3*f(N) - 4*f(N-1) + f(N-2)) / (2*h);

% Build sparse system A * df = rhs
A = spdiags([a b c], -1:1, N, N);

% Impose boundary rows df(1)=df1, df(N)=dfN
A(1,:) = 0; A(1,1) = 1; rhs(1) = df1;
A(N,:) = 0; A(N,N) = 1; rhs(N) = dfN;

% Solve
df = A \ rhs;
end

