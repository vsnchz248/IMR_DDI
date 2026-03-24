function DDI_IMR_fixed(use_gt_init, validation_mode)
%DDI_IMR_FIXED  Full joint (R*, V*, S*) DDI via 5N KKT system.
%
% CALLS
%   DDI_IMR_fixed()              synthetic validation mode
%   DDI_IMR_fixed(false,false)   production-style compressed centroids
%
% Notes/campaign-aligned version:
% 1) the S-smoothness residual/Jacobian are scaled exactly like the
%    objective beta_S*dt*sum((S_{n+1}-S_n)^2),
% 2) kinematics are enforced with interval defects that match sampled ODE
%    outputs more naturally than pointwise central differences, and
% 3) the inner Newton solve uses residual-decreasing backtracking instead
%    of ad hoc trajectory clipping.

if nargin < 1, use_gt_init = false; end
if nargin < 2, validation_mode = true; end
close all; clc;

%% =========================================================================
%  PHYSICAL CONSTANTS
%% =========================================================================
rho_val   = 1064.0;
P_inf_val = 101325.0;
gamma_s   = 0.07;
c_inf     = 1484.0;
c_sp      = c_inf / sqrt(P_inf_val / rho_val);
kappa     = 1.47;
pv_sat    = 0.0;
mu_gt     = 0.1;
G_gt      = 1e4;
Ca_gt     = P_inf_val / G_gt;

%% =========================================================================
%  DATA LOADING
%% =========================================================================
load('synthetic_data.mat');
lambda_vs_Rmax_fig; close all;
fprintf('  Rmax=%.1f um, lambda=%.2f\n', Rmax*1e6, lambda);

Nexp_gen       = numel(synthetic_data);
Rmax_range_gen = linspace(Rmax-0.5*Rmax_std, Rmax+0.5*Rmax_std, Nexp_gen);

R_raw = zeros(256, Nexp_gen);
V_raw = zeros(256, Nexp_gen);
t_raw = zeros(256, Nexp_gen);
for i = 1:Nexp_gen
    t_raw(:,i)    = synthetic_data{i}(:,1);
    R_raw(:,i)    = synthetic_data{i}(:,2);
    V_raw(:,i)    = synthetic_data{i}(:,3);
end

if validation_mode
    fprintf('  Assuming synthetic_data columns are already nondimensional IMR outputs (dimout=0).\n');
    fprintf('  Synthetic validation mode uses the full IMR output grids from t*=0, R*=1.\n');
    minLen       = size(R_raw,1);
    Rmatrix      = R_raw;
    Vmatrix      = V_raw;
    tmatrix      = t_raw - t_raw(1,:);
    Rmax_actual = Rmax_range_gen;
else
    R_d_raw = zeros(size(R_raw));
    for i = 1:Nexp_gen
        R_d_raw(:,i) = R_raw(:,i) .* Rmax_range_gen(i);
    end

    [~, startIdx] = max(R_d_raw, [], 1);
    minLen = min(size(R_d_raw,1) - startIdx + 1);
    R_d   = zeros(minLen, Nexp_gen);
    t_d   = zeros(minLen, Nexp_gen);
    V_d   = zeros(minLen, Nexp_gen);
    for j = 1:Nexp_gen
        idx       = startIdx(j):(startIdx(j)+minLen-1);
        R_d(:,j)  = R_d_raw(idx, j);
        t_d(:,j)  = t_raw(idx, j) - t_raw(startIdx(j), j);
        V_d(:,j)  = V_raw(idx, j);
    end

    Rmax_actual = max(R_d, [], 1);
    tc_actual   = Rmax_actual .* sqrt(rho_val/P_inf_val);
    Rmatrix     = R_d ./ Rmax_actual;
    Vmatrix     = V_d;
    tmatrix     = t_d ./ tc_actual;
end

We_vec = P_inf_val .* Rmax_actual ./ (2*gamma_s);
Re_vec = Rmax_actual .* sqrt(rho_val*P_inf_val) ./ mu_gt;

fprintf('  We: ');  fprintf('%.1f  ',We_vec);  fprintf('\n');
fprintf('  Re: ');  fprintf('%.1f  ',Re_vec);  fprintf('\n');
fprintf('  Ca: %.4f\n\n', Ca_gt);

%% =========================================================================
%  NATIVE PER-EXPERIMENT NORMALIZED GRIDS
%% =========================================================================
Nexp   = Nexp_gen;
N      = minLen;
t_grid = tmatrix';
dt_vec = tmatrix(2,:) - tmatrix(1,:);

R_obs   = Rmatrix';
V_exact = Vmatrix';
R_obs = max(R_obs, 0.02);

%% =========================================================================
%  R*_eq AND GROUND TRUTH S*(t)
%% =========================================================================
if validation_mode
    Req_range_gen = Rmax_range_gen ./ lambda;
    Req_star_vec  = Req_range_gen ./ Rmax_actual;
    fprintf('  Using exact synthetic R*_eq from generator (Req = Rmax/lambda).\n');
else
    last10       = max(1, floor(0.1*N));
    Req_star_vec = mean(R_obs(:, end-last10+1:end), 2)';
    Req_star_vec = max(Req_star_vec, 0.02);
end
fprintf('  R*_eq: ');  fprintf('%.4f  ',Req_star_vec);  fprintf('\n');

S_gt = zeros(Nexp, N);
for a = 1:Nexp
    Rstar     = max(R_obs(a,:), 1e-3);
    ratio     = Req_star_vec(a) ./ Rstar;
    S_gt(a,:) = -1/(2*Ca_gt) .* (5 - ratio.^4 - 4.*ratio) ...
               -4/Re_vec(a) .* V_exact(a,:) ./ Rstar;
end
fprintf('  GT S* range: [%.2f, %.2f]\n\n', min(S_gt(:)), max(S_gt(:)));
fprintf('  Solving on native per-experiment normalized grids (no common-grid resampling).\n\n');

%% =========================================================================
%  SMOOTH R* AND V* (warm-start only)
%% =========================================================================
smooth_sigma = 5;
hw = ceil(3*smooth_sigma);
gw = exp(-(-hw:hw).^2/(2*smooth_sigma^2));  gw=gw/sum(gw);
R_sm  = zeros(Nexp,N);
V_est = zeros(Nexp,N);
for a = 1:Nexp
    rs = conv(R_obs(a,:),gw,'same');  rs(1)=1.0;
    R_sm(a,:) = rs;
    V_est(a,2:N-1) = (rs(3:N)-rs(1:N-2))/(2*dt_vec(a));
    V_est(a,1)=0;  V_est(a,N)=V_est(a,N-1);
end

%% =========================================================================
%  pb FIXED FROM R_obs
%  For synthetic/self-consistency validation, evaluate the polytropic gas
%  law pointwise instead of numerically re-integrating sampled data.
%% =========================================================================
pb0_vec = zeros(1, Nexp);
pb      = zeros(Nexp, N);
dpb     = zeros(Nexp, N);
for a = 1:Nexp
    Req_a      = Req_star_vec(a);
    Rstar      = max(R_obs(a,:), 0.02);
    gas_coeff  = 1 + 1/(We_vec(a)*Req_a) - pv_sat;
    pb0_vec(a) = pv_sat + gas_coeff*Req_a^(3*kappa);
    pb(a,:)    = pv_sat + gas_coeff*(Req_a ./ Rstar).^(3*kappa);
    dpb(a,:)   = -3*kappa*(pb(a,:) - pv_sat).*V_exact(a,:) ./ Rstar;
end
fprintf('  pb0: ');  fprintf('%.4f  ',pb0_vec);  fprintf('\n\n');

%% =========================================================================
%  S_init = S_NH + S_v
%% =========================================================================
S_init = zeros(Nexp, N);
for a = 1:Nexp
    Rstar       = max(R_obs(a,:), 1e-3);
    ratio       = Req_star_vec(a) ./ Rstar;
    S_init(a,:) = -1/(2*Ca_gt) .* (5 - ratio.^4 - 4.*ratio) ...
                  -4/Re_vec(a) .* V_exact(a,:) ./ Rstar;
end
fprintf('  S_init (S_NH+S_v) range: [%.2f, %.2f]\n', min(S_init(:)), max(S_init(:)));
fprintf('  Max|S_init - GT|: ');
for a=1:Nexp, fprintf('Exp%d=%.3f  ',a,max(abs(S_init(a,:)-S_gt(a,:)))); end
fprintf('\n\n');

%% =========================================================================
%  WEIGHTS
%% =========================================================================
wn_all = ones(Nexp, N);

%% =========================================================================
%  DDI SETTINGS (notes-faithful objective)
%% =========================================================================
w_obs_R      = 100.0;
wR           = 1.0;
wV           = 0.0;
wS           = 1.0;
beta_S0      = 10.0;
beta_S1      = 1.0;
beta_R       = 0.0; %#ok<NASGU>
anneal_iters = 50;
Npts_total   = Nexp*N;
if validation_mode
    Nbar = Npts_total;
else
    Nbar = min(50*Nexp, Npts_total);
end
max_iter     = 200;
tol          = 1e-4;
max_nr       = 50;
fprintf('  Using Nbar = %d centroids (of %d total samples).\n\n', Nbar, Npts_total);
fprintf('  Notes/campaign settings: fixed pb, wV=0, beta_R=0, no extra priors.\n\n');
if validation_mode
    wR_man = 0.0;
    wV_man = 0.0;
    wS_man = 0.0;
    fprintf('  Synthetic validation mode: identity assignments, centroid compression disabled.\n');
    fprintf('  Synthetic validation mode: enforcing exact IMRv2 stress=1 closure for S and Sdot.\n\n');
    fprintf('  Synthetic validation mode: using exact IMR V* and projecting R* onto the native-grid radial closure.\n\n');
else
    wR_man = wR;
    wV_man = wV;
    wS_man = wS;
    fprintf('  Production-style mode: centroid compression enabled.\n\n');
end

%% =========================================================================
%  INITIALIZATION
%% =========================================================================
if use_gt_init
    fprintf('  [DEBUG] Initializing from GT (R_obs, V_exact, S_gt).\n\n');
    R_traj = R_obs;
    V_traj = V_exact;
    S_traj = S_gt;
else
    fprintf('  [PRODUCTION] Initializing from (R_obs, V_exact, S_NH+S_v).\n\n');
    R_traj = R_obs;
    V_traj = V_exact;
    S_traj = S_init;
end

fprintf('  Initial RMSE R*: ');
for a=1:Nexp, fprintf('Exp%d=%.5f  ',a,rms(R_traj(a,:)-R_obs(a,:))); end
fprintf('\n  Initial RMSE S*: ');
for a=1:Nexp, fprintf('Exp%d=%.5f  ',a,rms(S_traj(a,:)-S_gt(a,:)));  end
fprintf('\n\n');

%% =========================================================================
%  K-MEANS++ SEEDED WITH THE SAME METRIC AS THE OBJECTIVE
%% =========================================================================
rng(42);
pts_all = [R_traj(:), V_traj(:), S_traj(:)];
if validation_mode
    Rbar   = pts_all(:,1)';
    Vbar   = pts_all(:,2)';
    Sbar   = pts_all(:,3)';
    assign = reshape(int32(1:Npts_total), Nexp, N);
else
    pts_seed = [sqrt(wR)*R_traj(:), sqrt(wV)*V_traj(:), sqrt(wS)*S_traj(:)];
    ptsn     = pts_seed ./ (std(pts_seed)+eps);
    ci      = zeros(Nbar,1);
    ci(1)   = randi(size(pts_all,1));
    d2      = inf(size(pts_all,1),1);
    for k = 2:Nbar
        d2      = min(d2, sum((ptsn-ptsn(ci(k-1),:)).^2,2));
        mass    = sum(d2);
        if mass <= eps || ~isfinite(mass)
            ci(k) = ci(k-1);
            continue;
        end
        cdf     = cumsum(d2 / mass);
        next_ci = find(rand <= cdf, 1, 'first');
        if isempty(next_ci), next_ci = ci(k-1); end
        ci(k)   = next_ci;
    end
    Rbar   = pts_all(ci,1)';
    Vbar   = pts_all(ci,2)';
    Sbar   = pts_all(ci,3)';
    assign = ones(Nexp,N,'int32');
end

%% =========================================================================
%  DDI MAIN LOOP
%% =========================================================================
fprintf('  %5s  %12s  %10s  %8s  |  %-9s  %-9s  %-9s\n',...
    'Iter','Objective','Rel.dL','beta_S','RMSE R1','RMSE R2','RMSE R3');
L_hist = nan(max_iter,1);

warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:singularMatrix');

for iter = 1:max_iter

    if iter <= anneal_iters
        beta_S = beta_S0*(beta_S1/beta_S0)^((iter-1)/(anneal_iters-1));
    else
        beta_S = beta_S1;
    end

    if ~validation_mode
        for a = 1:Nexp
            for n = 1:N
                d_assign = wR*(R_traj(a,n)-Rbar).^2 + ...
                           wV*(V_traj(a,n)-Vbar).^2 + ...
                           wS*(S_traj(a,n)-Sbar).^2;
                [~,assign(a,n)] = min(d_assign);
            end
        end
    end

    for a = 1:Nexp
        if validation_mode
            [R_traj(a,:),V_traj(a,:),S_traj(a,:)] = solve_validation_exactV_exactstress( ...
                N, dt_vec(a), c_sp, We_vec(a), Req_star_vec(a), Ca_gt, Re_vec(a), ...
                R_obs(a,:), V_exact(a,:), pb(a,:), dpb(a,:), max_nr);
        else
            [R_traj(a,:),V_traj(a,:),S_traj(a,:)] = solve_KKT_fixed( ...
                N, dt_vec(a), c_sp, We_vec(a), R_obs(a,:), pb(a,:), dpb(a,:), ...
                wn_all(a,:)', ...
                w_obs_R, wR_man, wV_man, wS_man, beta_S, Rbar, Vbar, Sbar, ...
                assign(a,:), max_nr, R_traj(a,:)', V_traj(a,:)', S_traj(a,:)');
        end
    end

    if validation_mode
        % In the no-compression limit, each sample keeps its own centroid.
        Rbar = R_traj(:)';
        Vbar = V_traj(:)';
        Sbar = S_traj(:)';
    else
        nR=zeros(1,Nbar); nV=zeros(1,Nbar); nS=zeros(1,Nbar); dn=zeros(1,Nbar);
        for a=1:Nexp
            for n=1:N
                k=assign(a,n);
                nR(k)=nR(k)+R_traj(a,n);
                nV(k)=nV(k)+V_traj(a,n);
                nS(k)=nS(k)+S_traj(a,n);
                dn(k)=dn(k)+1;
            end
        end
        mk=dn>0;
        Rbar(mk)=nR(mk)./dn(mk);
        Vbar(mk)=nV(mk)./dn(mk);
        Sbar(mk)=nS(mk)./dn(mk);
    end

    L = 0;
    for a = 1:Nexp
        for n = 1:N
            k = assign(a,n);
            if validation_mode
                L = L + (R_traj(a,n)-R_obs(a,n))^2;
            else
                L = L + (w_obs_R*(R_traj(a,n)-R_obs(a,n))^2 ...
                       +  wR_man*(R_traj(a,n)-Rbar(k))^2 ...
                       +  wV_man*(V_traj(a,n)-Vbar(k))^2 ...
                       +  wS_man*(S_traj(a,n)-Sbar(k))^2);
            end
        end
        if ~validation_mode
            for n = 1:N-1
                L = L + beta_S*dt_vec(a)*(S_traj(a,n+1)-S_traj(a,n))^2;
            end
        end
    end
    L_hist(iter) = L;

    rmse_R = arrayfun(@(a)rms(R_traj(a,:)-R_obs(a,:)),1:Nexp);
    if iter > 1
        dL = abs(L-L_hist(iter-1))/(abs(L_hist(iter-1))+eps);
        fprintf('  %5d  %12.4e  %10.2e  %8.3f  |  %-9.5f  %-9.5f  %-9.5f\n',...
            iter,L,dL,beta_S,rmse_R(1),rmse_R(2),rmse_R(3));
        if dL<tol && iter>anneal_iters, fprintf('  Converged.\n'); break; end
        if iter>=anneal_iters+6
            rc=L_hist(iter-4:iter);
            if std(diff(rc))/(abs(mean(rc))+eps)<1e-3, fprintf('  Stagnated.\n'); break; end
        end
    else
        fprintf('  %5d  %12.4e  %10s  %8.3f  |  %-9.5f  %-9.5f  %-9.5f\n',...
            iter,L,'---',beta_S,rmse_R(1),rmse_R(2),rmse_R(3));
    end
end

warning('on','MATLAB:nearlySingularMatrix');
warning('on','MATLAB:singularMatrix');

%% =========================================================================
%  PLOTS
%% =========================================================================
cDDI=[0.8 0.15 0.15]; cGT=[0.05 0.55 0.15]; cols=lines(Nexp);

figure('Name','DDI: R*, V*, S* trajectories','Position',[30 30 1600 320*Nexp]);
for a=1:Nexp
    subplot(Nexp,3,(a-1)*3+1);
    plot(t_grid(a,:),R_obs(a,:), 'k.','MarkerSize',2.5,'DisplayName','R* observed'); hold on;
    plot(t_grid(a,:),R_traj(a,:),'-','Color',cDDI,'LineWidth',1.5,'DisplayName','R* DDI');
    yline(Req_star_vec(a),'--','Color',[.5 .5 .5],'LineWidth',1,'DisplayName','R*_{eq}');
    xlabel('t*'); ylabel('R*'); grid on; box on;
    title(sprintf('Exp %d  R*  We=%.0f  RMSE=%.4f',a,We_vec(a),rms(R_traj(a,:)-R_obs(a,:))));
    if a==1, legend('Location','best','FontSize',7); end

    subplot(Nexp,3,(a-1)*3+2);
    plot(t_grid(a,:),V_exact(a,:),'k-','LineWidth',1,'DisplayName','V* from solver (ref)'); hold on;
    plot(t_grid(a,:),V_traj(a,:), '-','Color',cDDI,'LineWidth',1.5,'DisplayName','V* DDI');
    xlabel('t*'); ylabel('V*'); grid on; box on;
    title(sprintf('Exp %d  V*  Re=%.1f  RMSE=%.4f',a,Re_vec(a),rms(V_traj(a,:)-V_exact(a,:))));
    if a==1, legend('Location','best','FontSize',7); end

    subplot(Nexp,3,(a-1)*3+3);
    plot(t_grid(a,:),S_gt(a,:),   '-','Color',cGT, 'LineWidth',2.5,'DisplayName','S* GT (synthetic only)'); hold on;
    plot(t_grid(a,:),S_traj(a,:), '-','Color',cDDI,'LineWidth',1.5,'DisplayName','S* DDI');
    xlabel('t*'); ylabel('S*'); grid on; box on;
    ylim([min(S_gt(a,:))*1.5-1, max(S_gt(a,:))*1.5+1]);
    title(sprintf('Exp %d  S*  Ca=%.2f  RMSE=%.3f',a,Ca_gt,rms(S_traj(a,:)-S_gt(a,:))));
    if a==1, legend('Location','best','FontSize',7); end
end

figure('Name','DDI: Phase space, manifold, convergence','Position',[100 50 1300 800]);

subplot(2,2,1);
for a=1:Nexp
    lc=cols(a,:);
    plot(R_obs(a,:),   S_gt(a,:),   '-','Color',0.5*lc+0.4,'LineWidth',2,...
        'DisplayName',sprintf('GT exp%d',a)); hold on;
    plot(R_traj(a,:),  S_traj(a,:), '--','Color',lc,'LineWidth',1.5,...
        'DisplayName',sprintf('DDI exp%d',a));
end
xlabel('R*'); ylabel('S*'); title('S*(R*): GT (solid) vs DDI (dashed)');
ylim([min(S_gt(:))*1.5-1, max(S_gt(:))*1.5+1]);
grid on; legend('Location','best','FontSize',7); box on;

subplot(2,2,2);
for a=1:Nexp
    lc=cols(a,:);
    plot3(R_traj(a,:), V_traj(a,:), S_traj(a,:), '-','Color',lc,'LineWidth',1.2,...
        'DisplayName',sprintf('DDI exp%d',a)); hold on;
end
occupied = find(arrayfun(@(k) any(assign(:)==k), 1:Nbar));
scatter3(Rbar(occupied), Vbar(occupied), Sbar(occupied), 20, 'r', 'filled', ...
    'DisplayName','Centroids');
xlabel('R*'); ylabel('V*'); zlabel('S*');
title('3D Phase space (R*, V*, S*) with centroids');
legend('Location','best','FontSize',7); grid on; view(35,20); box on;

subplot(2,2,3);
v=~isnan(L_hist);
semilogy(find(v),L_hist(v),'-o','Color',[0.2 0.45 0.75],'LineWidth',1.5,...
    'MarkerSize',3,'MarkerFaceColor',[0.2 0.45 0.75]);
xlabel('Iteration'); ylabel('L'); title('Convergence'); grid on;

subplot(2,2,4);
for a=1:Nexp
    lc=cols(a,:);
    plot(t_grid(a,:),S_gt(a,:),  '-','Color',0.5*lc+0.4,'LineWidth',2,...
        'DisplayName',sprintf('GT exp%d',a)); hold on;
    plot(t_grid(a,:),S_traj(a,:),'--','Color',lc,'LineWidth',1.5,...
        'DisplayName',sprintf('DDI exp%d',a));
end
xlabel('t*'); ylabel('S*'); title('S*(t*): GT vs DDI');
ylim([min(S_gt(:))*1.5-1, max(S_gt(:))*1.5+1]);
grid on; legend('Location','best','FontSize',7); box on;

fprintf('\n  %-5s  %-8s  %-8s  %-12s  %-12s  %-12s\n','Exp','We','Re','RMSE R*','RMSE V*','RMSE S*');
for a=1:Nexp
    fprintf('  %-5d  %-8.1f  %-8.2f  %-12.5f  %-12.5f  %-12.4f\n',a,...
        We_vec(a),Re_vec(a),...
        rms(R_traj(a,:)-R_obs(a,:)),...
        rms(V_traj(a,:)-V_exact(a,:)),...
        rms(S_traj(a,:)-S_gt(a,:)));
end

fprintf('\n  Discrete constraint residuals on each native experiment grid\n');
fprintf('  %-5s  %-12s  %-12s  %-12s  %-12s\n','Exp','GT KM rms','DDI KM rms','GT kin rms','DDI kin rms');
for a=1:Nexp
    if validation_mode
        [km_gt, kin_gt]   = constraint_residuals_exactstress(N, dt_vec(a), c_sp, We_vec(a), Req_star_vec(a), Ca_gt, Re_vec(a), R_obs(a,:),  V_exact(a,:), pb(a,:), dpb(a,:));
        [km_ddi, kin_ddi] = constraint_residuals_exactstress(N, dt_vec(a), c_sp, We_vec(a), Req_star_vec(a), Ca_gt, Re_vec(a), R_traj(a,:), V_traj(a,:), pb(a,:), dpb(a,:));
    else
        [km_gt, kin_gt]   = constraint_residuals(N, dt_vec(a), c_sp, We_vec(a), R_obs(a,:),  V_exact(a,:), S_gt(a,:),   pb(a,:), dpb(a,:));
        [km_ddi, kin_ddi] = constraint_residuals(N, dt_vec(a), c_sp, We_vec(a), R_traj(a,:), V_traj(a,:), S_traj(a,:), pb(a,:), dpb(a,:));
    end
    fprintf('  %-5d  %-12.4e  %-12.4e  %-12.4e  %-12.4e\n', a, ...
        rms(km_gt), rms(km_ddi), rms(kin_gt), rms(kin_ddi));
end
end

%% =========================================================================
function [pb,dpb]=compute_pb(R_in,V_in,N,dt,kappa,pb0)
pb=zeros(1,N); dpb=zeros(1,N); pb(1)=pb0;
for n=1:N-1
    Rn=max(R_in(n),0.02); Vn=V_in(n);
    k1=-3*kappa*pb(n)*Vn/Rn;
    Rh=max(Rn+0.5*dt*Vn,0.02); Vh=0.5*(Vn+V_in(min(n+1,N)));
    k2=-3*kappa*max(pb(n)+0.5*dt*k1,1e-6)*Vh/Rh;
    pb(n+1)=max(pb(n)+dt*k2,1e-6);
end
dpb(2:N-1)=(pb(3:N)-pb(1:N-2))/(2*dt);
dpb(1)=dpb(2); dpb(N)=dpb(N-1);
end

%% =========================================================================
function [km_res,kin_res]=constraint_residuals(N,dt,c,We,R,V,S,pb,dpb)
km_res  = zeros(1, N);
kin_res = zeros(1, N-1);
for n = 2:N-1
    dV=(V(n+1)-V(n-1))/(2*dt);
    dS=(S(n+1)-S(n-1))/(2*dt);
    km_res(n)=(1-V(n)/c)*R(n)*dV+1.5*(1-V(n)/(3*c))*V(n)^2 ...
             -(1+V(n)/c)*(pb(n)-1/(We*R(n))+S(n)-1) ...
             -R(n)/c*(dpb(n)+V(n)/(We*R(n)^2)+dS);
end
for n = 1:N-1
    kin_res(n)=R(n+1)-R(n)-0.5*dt*(V(n+1)+V(n));
end
end

%% =========================================================================
function [km_res,kin_res]=constraint_residuals_exactstress(N,dt,c,We,Req,Ca,Re,R,V,pb,dpb)
km_res  = zeros(1, N);
kin_res = zeros(1, N);
dV_all = fd1_uniform(V, dt);
dR_all = fd1_uniform(R, dt);
for n = 2:N-1
    dV = dV_all(n);
    [S, Sdot] = exact_stress_state(Req, Ca, Re, R(n), V(n));
    km_res(n)=(1-V(n)/c)*R(n)*dV+1.5*(1-V(n)/(3*c))*V(n)^2 ...
             -(1+V(n)/c)*(pb(n)-1/(We*R(n))+S-1) ...
             -R(n)/c*(dpb(n)+V(n)/(We*R(n)^2)+Sdot);
end
kin_res = V - dR_all;
end

%% =========================================================================
function dX = fd1_uniform(X, dt)
N = numel(X);
dX = zeros(size(X));
if N < 3
    return;
end
if N < 5
    dX(1)   = (-3*X(1) + 4*X(2) - X(3)) / (2*dt);
    dX(N)   = (3*X(N) - 4*X(N-1) + X(N-2)) / (2*dt);
    if N > 2
        dX(2:N-1) = (X(3:N) - X(1:N-2)) / (2*dt);
    end
    return;
end
dX(1)   = (-3*X(1) + 4*X(2) - X(3)) / (2*dt);
dX(2)   = (X(3) - X(1)) / (2*dt);
dX(3:N-2) = (-X(5:N) + 8*X(4:N-1) - 8*X(2:N-3) + X(1:N-4)) / (12*dt);
dX(N-1) = (X(N) - X(N-2)) / (2*dt);
dX(N)   = (3*X(N) - 4*X(N-1) + X(N-2)) / (2*dt);
end

%% =========================================================================
function [R_out,V_out,S_out]=solve_validation_exactV_exactstress( ...
    N,dt,c,We,Req,Ca,Re,R_obs,V_exact,pb,dpb,max_nr)

R_out = max(R_obs(:), 0.02);
V_out = V_exact(:);
dV_all = fd1_uniform(V_out, dt);

R_out(1) = 1;
for n = 2:N-1
    r = max(R_out(n), 0.02);
    vn = V_out(n);
    dV = dV_all(n);
    for nr = 1:max_nr
        [S, Sdot, dS_dR, ~, dSdot_dR, ~] = exact_stress_state(Req, Ca, Re, r, vn);
        E = pb(n) - 1/(We*r) + S - 1;
        F = dpb(n) + vn/(We*r^2) + Sdot;
        C1 = 1 - vn/c;

        Fc = C1*r*dV + 1.5*(1-vn/(3*c))*vn^2 ...
           - (1+vn/c)*E - r/c*F;

        dFc_dR = C1*dV ...
               - (1+vn/c)*(1/(We*r^2) + dS_dR) ...
               - F/c ...
               - r/c*(-2*vn/(We*r^3) + dSdot_dR);

        if abs(dFc_dR) < 1e-12
            break;
        end

        step = Fc / dFc_dR;
        alpha = 1;
        accepted = false;
        for ls = 1:12
            r_trial = min(max(r - alpha*step, 0.02), 2.0);
            [S_t, Sdot_t] = exact_stress_state(Req, Ca, Re, r_trial, vn);
            E_t = pb(n) - 1/(We*r_trial) + S_t - 1;
            F_t = dpb(n) + vn/(We*r_trial^2) + Sdot_t;
            Fc_t = C1*r_trial*dV + 1.5*(1-vn/(3*c))*vn^2 ...
                 - (1+vn/c)*E_t - r_trial/c*F_t;
            if abs(Fc_t) <= abs(Fc)
                r = r_trial;
                accepted = true;
                break;
            end
            alpha = 0.5*alpha;
        end
        if ~accepted || abs(Fc) < 1e-10
            break;
        end
    end
    R_out(n) = r;
end

R_out = R_out';
V_out = V_out';
[S_out,~] = exact_stress_state(Req, Ca, Re, R_out, V_out);
end

%% =========================================================================
function [S, Sdot, dS_dR, dS_dV, dSdot_dR, dSdot_dV] = exact_stress_state(Req,Ca,Re,R,V)
R = max(R, 1e-8);
Rst = Req./R;
S = -(5 - 4*Rst - Rst.^4)/(2*Ca) - 4/Re .* V./R;
Sdot = -2*V./R.*(Rst + Rst.^4)/Ca + 4/Re.*(V./R).^2;
if nargout > 2
    dS_dR    = -2/Ca .* Rst .* (1 + Rst.^3) ./ R + 4/Re .* V ./ R.^2;
    dS_dV    = -4/(Re*R);
    dSdot_dR = 2/Ca .* V ./ R.^2 .* (2*Rst + 5*Rst.^4) - 8/Re .* V.^2 ./ R.^3;
    dSdot_dV = -2/Ca .* (Rst + Rst.^4) ./ R + 8/Re .* V ./ R.^2;
end
end

%% =========================================================================
function res=kkt_residual_only( ...
    N,dt,c,We,R_obs,pb,dpb,wn, ...
    w_obs_R,wR,wV,wS,beta_S,Rbar,Vbar,Sbar,assign,sol)

iR=@(n)n;
iV=@(n)N+n;
iS=@(n)2*N+n;
iL=@(n)3*N+n;
iM=@(n)4*N+n;

R=max(sol(iR(1):iR(N)),0.005);
V=sol(iV(1):iV(N));
S=sol(iS(1):iS(N));
lam=sol(iL(1):iL(N));
mu=sol(iM(1):iM(N));

res=zeros(5*N,1);

res(iR(1))=R(1)-1;
res(iV(1))=V(1);
res(iL(1))=lam(1);
res(iL(N))=lam(N);
res(iM(1))=R(2)-R(1)-0.5*dt*(V(2)+V(1));
res(iM(N))=mu(N);

for n=2:N-1
    kk=assign(n);
    dV=(V(n+1)-V(n-1))/(2*dt);
    dS=(S(n+1)-S(n-1))/(2*dt);

    A=(1-V(n)/c)*dV-1/(We*R(n)^2)-(dpb(n)+dS)/c;
    B=-R(n)/c*dV+3*V(n)-3*V(n)^2/(2*c)-(pb(n)+S(n)-1)/c;

    Fc=(1-V(n)/c)*R(n)*dV+1.5*(1-V(n)/(3*c))*V(n)^2 ...
      -(1+V(n)/c)*(pb(n)-1/(We*R(n))+S(n)-1) ...
      -R(n)/c*(dpb(n)+V(n)/(We*R(n)^2)+dS);
    res(iL(n))=Fc;

    Gc=R(n+1)-R(n)-0.5*dt*(V(n+1)+V(n));
    res(iM(n))=Gc;

    mu_term = mu(n-1) - mu(n);

    res(iR(n))=2*wn(n)*(w_obs_R*(R(n)-R_obs(n)) ...
                + wR*(R(n)-Rbar(kk))) ...
              + A*lam(n) + mu_term;

    cn1=-(1-V(n+1)/c)*R(n+1)/(2*dt);
    cm1=(1-V(n-1)/c)*R(n-1)/(2*dt);
    res(iV(n))=2*wn(n)*wV*(V(n)-Vbar(kk)) ...
              + B*lam(n)+cn1*lam(n+1)+cm1*lam(n-1) ...
              - 0.5*dt*(mu(n-1)+mu(n));

    lS=-(1+V(n)/c)*lam(n);
    if n+1<=N-1, lS=lS+R(n+1)*lam(n+1)/(2*c*dt); end
    if n-1>=2,   lS=lS-R(n-1)*lam(n-1)/(2*c*dt); end

    sm=2*beta_S*dt*((S(n)-S(n-1))-(S(n+1)-S(n)));
    res(iS(n))=2*wn(n)*wS*(S(n)-Sbar(kk))+sm+lS;
end

kN=assign(N);
k1=assign(1);

res(iR(N))=2*wn(N)*(w_obs_R*(R(N)-R_obs(N)) ...
            + wR*(R(N)-Rbar(kN))) ...
          + mu(N-1);

cN=(1-V(N-1)/c)*R(N-1)/(2*dt);
res(iV(N))=2*wn(N)*wV*(V(N)-Vbar(kN))+cN*lam(N-1)-0.5*dt*mu(N-1);

sm1=2*beta_S*dt*(S(1)-S(2));
res(iS(1))=2*wn(1)*wS*(S(1)-Sbar(k1)) ...
         + sm1 + R(2)*lam(2)/(2*c*dt);

smN=2*beta_S*dt*(S(N)-S(N-1));
res(iS(N))=2*wn(N)*wS*(S(N)-Sbar(kN))+smN ...
          -R(N-1)*lam(N-1)/(2*c*dt);
end

%% =========================================================================
function [R_out,V_out,S_out]=solve_KKT_fixed( ...
    N,dt,c,We,R_obs,pb,dpb,wn, ...
    w_obs_R,wR,wV,wS,beta_S,Rbar,Vbar,Sbar,assign,max_nr, ...
    R_w,V_w,S_w)

iR=@(n)n;
iV=@(n)N+n;
iS=@(n)2*N+n;
iL=@(n)3*N+n;
iM=@(n)4*N+n;
feas_idx=[(3*N+2):(4*N-1), (4*N+1):(5*N-1)];

sol=[max(R_w,0.02);V_w;S_w;zeros(2*N,1)];
sol(iR(1))=1;
sol(iV(1))=0;

for nr=1:max_nr
    R=max(sol(iR(1):iR(N)),0.005);
    V=sol(iV(1):iV(N));
    S=sol(iS(1):iS(N));
    lam=sol(iL(1):iL(N));
    mu=sol(iM(1):iM(N));

    res=zeros(5*N,1);
    nz=0;
    TI=zeros(72*N,1);
    TJ=zeros(72*N,1);
    TS=zeros(72*N,1);

    % Fixed initial conditions and terminal dummy multiplier boundary.
    res(iR(1))=R(1)-1;
    nz=nz+1; TI(nz)=iR(1); TJ(nz)=iR(1); TS(nz)=1;

    res(iV(1))=V(1);
    nz=nz+1; TI(nz)=iV(1); TJ(nz)=iV(1); TS(nz)=1;

    res(iL(1))=lam(1);
    nz=nz+1; TI(nz)=iL(1); TJ(nz)=iL(1); TS(nz)=1;

    res(iL(N))=lam(N);
    nz=nz+1; TI(nz)=iL(N); TJ(nz)=iL(N); TS(nz)=1;

    % Interval kinematic constraint on [t_1,t_2], with mu(N) kept as dummy.
    res(iM(1)) = R(2)-R(1)-0.5*dt*(V(2)+V(1));
    nz=nz+1; TI(nz)=iM(1); TJ(nz)=iR(2); TS(nz)=1;
    nz=nz+1; TI(nz)=iM(1); TJ(nz)=iR(1); TS(nz)=-1;
    nz=nz+1; TI(nz)=iM(1); TJ(nz)=iV(2); TS(nz)=-0.5*dt;
    nz=nz+1; TI(nz)=iM(1); TJ(nz)=iV(1); TS(nz)=-0.5*dt;

    res(iM(N))=mu(N);
    nz=nz+1; TI(nz)=iM(N); TJ(nz)=iM(N); TS(nz)=1;

    for n=2:N-1
        kk=assign(n);
        dV=(V(n+1)-V(n-1))/(2*dt);
        dS=(S(n+1)-S(n-1))/(2*dt);

        A=(1-V(n)/c)*dV-1/(We*R(n)^2)-(dpb(n)+dS)/c;
        B=-R(n)/c*dV+3*V(n)-3*V(n)^2/(2*c)-(pb(n)+S(n)-1)/c;

        Fc=(1-V(n)/c)*R(n)*dV+1.5*(1-V(n)/(3*c))*V(n)^2 ...
          -(1+V(n)/c)*(pb(n)-1/(We*R(n))+S(n)-1) ...
          -R(n)/c*(dpb(n)+V(n)/(We*R(n)^2)+dS);
        res(iL(n))=Fc;
        nz=nz+1; TI(nz)=iL(n); TJ(nz)=iR(n);   TS(nz)=A;
        nz=nz+1; TI(nz)=iL(n); TJ(nz)=iV(n);   TS(nz)=B;
        nz=nz+1; TI(nz)=iL(n); TJ(nz)=iV(n+1); TS(nz)=(1-V(n)/c)*R(n)/(2*dt);
        nz=nz+1; TI(nz)=iL(n); TJ(nz)=iV(n-1); TS(nz)=-(1-V(n)/c)*R(n)/(2*dt);
        nz=nz+1; TI(nz)=iL(n); TJ(nz)=iS(n);   TS(nz)=-(1+V(n)/c);
        nz=nz+1; TI(nz)=iL(n); TJ(nz)=iS(n+1); TS(nz)=-R(n)/(2*c*dt);
        nz=nz+1; TI(nz)=iL(n); TJ(nz)=iS(n-1); TS(nz)=R(n)/(2*c*dt);

        Gc=R(n+1)-R(n)-0.5*dt*(V(n+1)+V(n));
        res(iM(n))=Gc;
        nz=nz+1; TI(nz)=iM(n); TJ(nz)=iR(n+1); TS(nz)=1;
        nz=nz+1; TI(nz)=iM(n); TJ(nz)=iR(n);   TS(nz)=-1;
        nz=nz+1; TI(nz)=iM(n); TJ(nz)=iV(n+1); TS(nz)=-0.5*dt;
        nz=nz+1; TI(nz)=iM(n); TJ(nz)=iV(n);   TS(nz)=-0.5*dt;

        mu_term = mu(n-1) - mu(n);

        res(iR(n))=2*wn(n)*(w_obs_R*(R(n)-R_obs(n)) ...
                    + wR*(R(n)-Rbar(kk))) ...
                  +A*lam(n)+mu_term;
        nz=nz+1; TI(nz)=iR(n); TJ(nz)=iR(n);   TS(nz)=2*wn(n)*(w_obs_R + wR)+2*lam(n)/(We*R(n)^3);
        nz=nz+1; TI(nz)=iR(n); TJ(nz)=iV(n);   TS(nz)=-lam(n)*dV/c;
        nz=nz+1; TI(nz)=iR(n); TJ(nz)=iV(n+1); TS(nz)=lam(n)*(1-V(n)/c)/(2*dt);
        nz=nz+1; TI(nz)=iR(n); TJ(nz)=iV(n-1); TS(nz)=-lam(n)*(1-V(n)/c)/(2*dt);
        nz=nz+1; TI(nz)=iR(n); TJ(nz)=iS(n+1); TS(nz)=-lam(n)/(2*c*dt);
        nz=nz+1; TI(nz)=iR(n); TJ(nz)=iS(n-1); TS(nz)=lam(n)/(2*c*dt);
        nz=nz+1; TI(nz)=iR(n); TJ(nz)=iL(n);   TS(nz)=A;
        nz=nz+1; TI(nz)=iR(n); TJ(nz)=iM(n-1); TS(nz)=1;
        nz=nz+1; TI(nz)=iR(n); TJ(nz)=iM(n);   TS(nz)=-1;

        cn1=-(1-V(n+1)/c)*R(n+1)/(2*dt);
        cm1=(1-V(n-1)/c)*R(n-1)/(2*dt);
        res(iV(n))=2*wn(n)*wV*(V(n)-Vbar(kk)) ...
                  +B*lam(n)+cn1*lam(n+1)+cm1*lam(n-1) ...
                  -0.5*dt*(mu(n-1)+mu(n));
        nz=nz+1; TI(nz)=iV(n); TJ(nz)=iR(n);   TS(nz)=-lam(n)*dV/c;
        nz=nz+1; TI(nz)=iV(n); TJ(nz)=iR(n+1); TS(nz)=-(1-V(n+1)/c)*lam(n+1)/(2*dt);
        nz=nz+1; TI(nz)=iV(n); TJ(nz)=iR(n-1); TS(nz)=(1-V(n-1)/c)*lam(n-1)/(2*dt);
        nz=nz+1; TI(nz)=iV(n); TJ(nz)=iV(n);   TS(nz)=2*wn(n)*wV+(3-3*V(n)/c)*lam(n);
        nz=nz+1; TI(nz)=iV(n); TJ(nz)=iV(n+1); TS(nz)=-R(n)*lam(n)/(2*c*dt)+R(n+1)*lam(n+1)/(2*c*dt);
        nz=nz+1; TI(nz)=iV(n); TJ(nz)=iV(n-1); TS(nz)=R(n)*lam(n)/(2*c*dt)-R(n-1)*lam(n-1)/(2*c*dt);
        nz=nz+1; TI(nz)=iV(n); TJ(nz)=iS(n);   TS(nz)=-lam(n)/c;
        nz=nz+1; TI(nz)=iV(n); TJ(nz)=iL(n);   TS(nz)=B;
        nz=nz+1; TI(nz)=iV(n); TJ(nz)=iL(n+1); TS(nz)=cn1;
        nz=nz+1; TI(nz)=iV(n); TJ(nz)=iL(n-1); TS(nz)=cm1;
        nz=nz+1; TI(nz)=iV(n); TJ(nz)=iM(n-1); TS(nz)=-0.5*dt;
        nz=nz+1; TI(nz)=iV(n); TJ(nz)=iM(n);   TS(nz)=-0.5*dt;

        lS=-(1+V(n)/c)*lam(n);
        if n+1<=N-1, lS=lS+R(n+1)*lam(n+1)/(2*c*dt); end
        if n-1>=2,   lS=lS-R(n-1)*lam(n-1)/(2*c*dt); end

        sm=2*beta_S*dt*((S(n)-S(n-1))-(S(n+1)-S(n)));
        res(iS(n))=2*wn(n)*wS*(S(n)-Sbar(kk))+sm+lS;
        nz=nz+1; TI(nz)=iS(n); TJ(nz)=iS(n);   TS(nz)=2*wn(n)*wS+4*beta_S*dt;
        nz=nz+1; TI(nz)=iS(n); TJ(nz)=iS(n+1); TS(nz)=-2*beta_S*dt;
        nz=nz+1; TI(nz)=iS(n); TJ(nz)=iS(n-1); TS(nz)=-2*beta_S*dt;
        nz=nz+1; TI(nz)=iS(n); TJ(nz)=iV(n);   TS(nz)=-lam(n)/c;
        nz=nz+1; TI(nz)=iS(n); TJ(nz)=iL(n);   TS(nz)=-(1+V(n)/c);
        if n+1<=N-1
            nz=nz+1; TI(nz)=iS(n); TJ(nz)=iR(n+1); TS(nz)=lam(n+1)/(2*c*dt);
            nz=nz+1; TI(nz)=iS(n); TJ(nz)=iL(n+1); TS(nz)=R(n+1)/(2*c*dt);
        end
        if n-1>=2
            nz=nz+1; TI(nz)=iS(n); TJ(nz)=iR(n-1); TS(nz)=-lam(n-1)/(2*c*dt);
            nz=nz+1; TI(nz)=iS(n); TJ(nz)=iL(n-1); TS(nz)=-R(n-1)/(2*c*dt);
        end
    end

    kN=assign(N);
    k1=assign(1);

    res(iR(N))=2*wn(N)*(w_obs_R*(R(N)-R_obs(N)) ...
                + wR*(R(N)-Rbar(kN))) ...
              +mu(N-1);
    nz=nz+1; TI(nz)=iR(N); TJ(nz)=iR(N);   TS(nz)=2*wn(N)*(w_obs_R + wR);
    nz=nz+1; TI(nz)=iR(N); TJ(nz)=iM(N-1); TS(nz)=1;

    cN=(1-V(N-1)/c)*R(N-1)/(2*dt);
    res(iV(N))=2*wn(N)*wV*(V(N)-Vbar(kN))+cN*lam(N-1)-0.5*dt*mu(N-1);
    nz=nz+1; TI(nz)=iV(N); TJ(nz)=iV(N);   TS(nz)=2*wn(N)*wV+1e-10;
    nz=nz+1; TI(nz)=iV(N); TJ(nz)=iR(N-1); TS(nz)=(1-V(N-1)/c)*lam(N-1)/(2*dt);
    nz=nz+1; TI(nz)=iV(N); TJ(nz)=iV(N-1); TS(nz)=-R(N-1)*lam(N-1)/(2*c*dt);
    nz=nz+1; TI(nz)=iV(N); TJ(nz)=iL(N-1); TS(nz)=cN;
    nz=nz+1; TI(nz)=iV(N); TJ(nz)=iM(N-1); TS(nz)=-0.5*dt;

    sm1=2*beta_S*dt*(S(1)-S(2));
    res(iS(1))=2*wn(1)*wS*(S(1)-Sbar(k1)) ...
             + sm1 + R(2)*lam(2)/(2*c*dt);
    nz=nz+1; TI(nz)=iS(1); TJ(nz)=iS(1); TS(nz)=2*wn(1)*wS+2*beta_S*dt;
    nz=nz+1; TI(nz)=iS(1); TJ(nz)=iS(2); TS(nz)=-2*beta_S*dt;
    nz=nz+1; TI(nz)=iS(1); TJ(nz)=iR(2); TS(nz)=lam(2)/(2*c*dt);
    nz=nz+1; TI(nz)=iS(1); TJ(nz)=iL(2); TS(nz)=R(2)/(2*c*dt);

    smN=2*beta_S*dt*(S(N)-S(N-1));
    res(iS(N))=2*wn(N)*wS*(S(N)-Sbar(kN))+smN ...
              -R(N-1)*lam(N-1)/(2*c*dt);
    nz=nz+1; TI(nz)=iS(N); TJ(nz)=iS(N);   TS(nz)=2*wn(N)*wS+2*beta_S*dt;
    nz=nz+1; TI(nz)=iS(N); TJ(nz)=iS(N-1); TS(nz)=-2*beta_S*dt;
    nz=nz+1; TI(nz)=iS(N); TJ(nz)=iR(N-1); TS(nz)=-lam(N-1)/(2*c*dt);
    nz=nz+1; TI(nz)=iS(N); TJ(nz)=iL(N-1); TS(nz)=-R(N-1)/(2*c*dt);

    J=sparse(TI(1:nz),TJ(1:nz),TS(1:nz),5*N,5*N);
    reg=1e-6*max(max(abs(diag(J))),1);
    delta=-(J+reg*speye(5*N))\res;

    rn=norm(res,inf);
    full_nrm=norm(res,2);
    feas_nrm=norm(res(feas_idx),2);
    alpha=1;
    accepted=false;
    for ls=1:20
        trial=sol+alpha*delta;
        if all(trial(1:N)>0.02) && all(trial(1:N)<2)
            trial(iR(1))=1;
            trial(iV(1))=0;
            trial(iL(1))=0;
            trial(iL(N))=0;
            trial(iM(N))=0;
            trial(1:N)=min(max(trial(1:N),0.005),2);

            trial_res = kkt_residual_only( ...
                N,dt,c,We,R_obs,pb,dpb,wn, ...
                w_obs_R,wR,wV,wS,beta_S,Rbar,Vbar,Sbar,assign,trial);
            trial_rn   = norm(trial_res,inf);
            trial_full = norm(trial_res,2);
            trial_feas = norm(trial_res(feas_idx),2);
            if (trial_feas < feas_nrm) || ...
               (trial_feas <= feas_nrm && trial_full < full_nrm)
                sol = trial;
                rn = trial_rn;
                accepted = true;
                break;
            end
        end
        alpha=alpha*0.5;
    end
    if ~accepted, break; end

    if rn<1e-8, break; end
end

R_out=sol(iR(1):iR(N))';
V_out=sol(iV(1):iV(N))';
S_out=sol(iS(1):iS(N))';
end

%% =========================================================================
function res=kkt_residual_only_exactstress( ...
    N,dt,c,We,Req,Ca,Re,R_obs,V_obs,pb,dpb,w_obs_R,w_obs_V,sol)

iR=@(n)n;
iV=@(n)N+n;
iL=@(n)2*N+n;
iM=@(n)3*N+n;

R=max(sol(iR(1):iR(N)),0.005);
V=sol(iV(1):iV(N));
lam=sol(iL(1):iL(N));
mu=sol(iM(1):iM(N));

res=zeros(4*N,1);
res(iR(1))=R(1)-1;
res(iV(1))=V(1);
res(iL(1))=lam(1);
res(iL(N))=lam(N);
res(iM(N))=mu(N);

for n=2:N-1
    dV=(V(n+1)-V(n-1))/(2*dt);
    [S, Sdot, dS_dR, dS_dV, dSdot_dR, dSdot_dV] = exact_stress_state(Req, Ca, Re, R(n), V(n));

    E=pb(n)-1/(We*R(n))+S-1;
    F=dpb(n)+V(n)/(We*R(n)^2)+Sdot;
    C1=1-V(n)/c;

    Fc=C1*R(n)*dV+1.5*(1-V(n)/(3*c))*V(n)^2 ...
      -(1+V(n)/c)*E - R(n)/c*F;
    res(iL(n))=Fc;

    Gc=R(n+1)-R(n)-0.5*dt*(V(n+1)+V(n));
    res(iM(n))=Gc;

    dFc_dR = C1*dV ...
           -(1+V(n)/c)*(1/(We*R(n)^2)+dS_dR) ...
           -F/c ...
           -R(n)/c*(-2*V(n)/(We*R(n)^3)+dSdot_dR);

    dFc_dV = -R(n)*dV/c ...
           + 3*V(n) - 1.5*V(n)^2/c ...
           - E/c ...
           -(1+V(n)/c)*dS_dV ...
           -R(n)/c*(1/(We*R(n)^2)+dSdot_dV);

    mu_term = mu(n-1) - mu(n);
    res(iR(n)) = 2*w_obs_R*(R(n)-R_obs(n)) + lam(n)*dFc_dR + mu_term;

        res(iV(n)) = 2*w_obs_V*(V(n)-V_obs(n)) ...
                   + lam(n)*dFc_dV - 0.5*dt*(mu(n-1)+mu(n));
        if n+1 <= N-1
            res(iV(n)) = res(iV(n)) + lam(n+1)*(-(1-V(n+1)/c)*R(n+1)/(2*dt));
        end
        if n-1 >= 2
            res(iV(n)) = res(iV(n)) + lam(n-1)*((1-V(n-1)/c)*R(n-1)/(2*dt));
    end
end

res(iM(1)) = R(2)-R(1)-0.5*dt*(V(2)+V(1));
res(iR(N)) = 2*w_obs_R*(R(N)-R_obs(N)) + mu(N-1);
res(iV(N)) = 2*w_obs_V*(V(N)-V_obs(N)) - 0.5*dt*mu(N-1);
if N-1 >= 2
    res(iV(N)) = res(iV(N)) + lam(N-1)*((1-V(N-1)/c)*R(N-1)/(2*dt));
end
end

%% =========================================================================
function [R_out,V_out,S_out]=solve_KKT_validation_exactstress( ...
    N,dt,c,We,Req,Ca,Re,R_obs,V_obs,pb,dpb,w_obs_R,w_obs_V,max_nr,R_w,V_w)

iR=@(n)n;
iV=@(n)N+n;
iL=@(n)2*N+n;
iM=@(n)3*N+n;
feas_idx=[(2*N+2):(3*N-1), (3*N+1):(4*N-1)];

sol=[max(R_w,0.02);V_w;zeros(2*N,1)];
sol(iR(1))=1;
sol(iV(1))=0;

for nr=1:max_nr
    R=max(sol(iR(1):iR(N)),0.005);
    V=sol(iV(1):iV(N));
    lam=sol(iL(1):iL(N));
    mu=sol(iM(1):iM(N));

    res=zeros(4*N,1);
    nz=0;
    TI=zeros(52*N,1);
    TJ=zeros(52*N,1);
    TS=zeros(52*N,1);

    res(iR(1))=R(1)-1;
    nz=nz+1; TI(nz)=iR(1); TJ(nz)=iR(1); TS(nz)=1;

    res(iV(1))=V(1);
    nz=nz+1; TI(nz)=iV(1); TJ(nz)=iV(1); TS(nz)=1;

    res(iL(1))=lam(1);
    nz=nz+1; TI(nz)=iL(1); TJ(nz)=iL(1); TS(nz)=1;

    res(iL(N))=lam(N);
    nz=nz+1; TI(nz)=iL(N); TJ(nz)=iL(N); TS(nz)=1;

    res(iM(1)) = R(2)-R(1)-0.5*dt*(V(2)+V(1));
    nz=nz+1; TI(nz)=iM(1); TJ(nz)=iR(2); TS(nz)=1;
    nz=nz+1; TI(nz)=iM(1); TJ(nz)=iR(1); TS(nz)=-1;
    nz=nz+1; TI(nz)=iM(1); TJ(nz)=iV(2); TS(nz)=-0.5*dt;
    nz=nz+1; TI(nz)=iM(1); TJ(nz)=iV(1); TS(nz)=-0.5*dt;

    res(iM(N))=mu(N);
    nz=nz+1; TI(nz)=iM(N); TJ(nz)=iM(N); TS(nz)=1;

    for n=2:N-1
        dV=(V(n+1)-V(n-1))/(2*dt);
        [S, Sdot, dS_dR, dS_dV, dSdot_dR, dSdot_dV] = exact_stress_state(Req, Ca, Re, R(n), V(n));

        E=pb(n)-1/(We*R(n))+S-1;
        F=dpb(n)+V(n)/(We*R(n)^2)+Sdot;
        C1=1-V(n)/c;

        dFc_dR = C1*dV ...
               -(1+V(n)/c)*(1/(We*R(n)^2)+dS_dR) ...
               -F/c ...
               -R(n)/c*(-2*V(n)/(We*R(n)^3)+dSdot_dR);

        dFc_dV = -R(n)*dV/c ...
               + 3*V(n) - 1.5*V(n)^2/c ...
               - E/c ...
               -(1+V(n)/c)*dS_dV ...
               -R(n)/c*(1/(We*R(n)^2)+dSdot_dV);

        dFc_dVp = C1*R(n)/(2*dt);
        dFc_dVm = -C1*R(n)/(2*dt);

        Fc=C1*R(n)*dV+1.5*(1-V(n)/(3*c))*V(n)^2 ...
          -(1+V(n)/c)*E - R(n)/c*F;
        res(iL(n))=Fc;
        nz=nz+1; TI(nz)=iL(n); TJ(nz)=iR(n);   TS(nz)=dFc_dR;
        nz=nz+1; TI(nz)=iL(n); TJ(nz)=iV(n);   TS(nz)=dFc_dV;
        nz=nz+1; TI(nz)=iL(n); TJ(nz)=iV(n+1); TS(nz)=dFc_dVp;
        nz=nz+1; TI(nz)=iL(n); TJ(nz)=iV(n-1); TS(nz)=dFc_dVm;

        Gc=R(n+1)-R(n)-0.5*dt*(V(n+1)+V(n));
        res(iM(n))=Gc;
        nz=nz+1; TI(nz)=iM(n); TJ(nz)=iR(n+1); TS(nz)=1;
        nz=nz+1; TI(nz)=iM(n); TJ(nz)=iR(n);   TS(nz)=-1;
        nz=nz+1; TI(nz)=iM(n); TJ(nz)=iV(n+1); TS(nz)=-0.5*dt;
        nz=nz+1; TI(nz)=iM(n); TJ(nz)=iV(n);   TS(nz)=-0.5*dt;

        res(iR(n)) = 2*w_obs_R*(R(n)-R_obs(n)) + lam(n)*dFc_dR + mu(n-1) - mu(n);
        nz=nz+1; TI(nz)=iR(n); TJ(nz)=iR(n);   TS(nz)=2*w_obs_R;
        nz=nz+1; TI(nz)=iR(n); TJ(nz)=iL(n);   TS(nz)=dFc_dR;
        nz=nz+1; TI(nz)=iR(n); TJ(nz)=iM(n-1); TS(nz)=1;
        nz=nz+1; TI(nz)=iR(n); TJ(nz)=iM(n);   TS(nz)=-1;

        res(iV(n)) = 2*w_obs_V*(V(n)-V_obs(n)) ...
                   + lam(n)*dFc_dV - 0.5*dt*(mu(n-1)+mu(n));
        if n+1 <= N-1
            cn1 = -(1-V(n+1)/c)*R(n+1)/(2*dt);
            res(iV(n)) = res(iV(n)) + cn1*lam(n+1);
            nz=nz+1; TI(nz)=iV(n); TJ(nz)=iL(n+1); TS(nz)=cn1;
        end
        if n-1 >= 2
            cm1 = (1-V(n-1)/c)*R(n-1)/(2*dt);
            res(iV(n)) = res(iV(n)) + cm1*lam(n-1);
            nz=nz+1; TI(nz)=iV(n); TJ(nz)=iL(n-1); TS(nz)=cm1;
        end
        nz=nz+1; TI(nz)=iV(n); TJ(nz)=iV(n);   TS(nz)=2*w_obs_V+1e-10;
        nz=nz+1; TI(nz)=iV(n); TJ(nz)=iL(n);   TS(nz)=dFc_dV;
        nz=nz+1; TI(nz)=iV(n); TJ(nz)=iM(n-1); TS(nz)=-0.5*dt;
        nz=nz+1; TI(nz)=iV(n); TJ(nz)=iM(n);   TS(nz)=-0.5*dt;
    end

    res(iR(N)) = 2*w_obs_R*(R(N)-R_obs(N)) + mu(N-1);
    nz=nz+1; TI(nz)=iR(N); TJ(nz)=iR(N);   TS(nz)=2*w_obs_R;
    nz=nz+1; TI(nz)=iR(N); TJ(nz)=iM(N-1); TS(nz)=1;

    res(iV(N)) = 2*w_obs_V*(V(N)-V_obs(N)) - 0.5*dt*mu(N-1);
    if N-1 >= 2
        cN = (1-V(N-1)/c)*R(N-1)/(2*dt);
        res(iV(N)) = res(iV(N)) + cN*lam(N-1);
        nz=nz+1; TI(nz)=iV(N); TJ(nz)=iL(N-1); TS(nz)=cN;
    end
    nz=nz+1; TI(nz)=iV(N); TJ(nz)=iV(N);   TS(nz)=2*w_obs_V+1e-10;
    nz=nz+1; TI(nz)=iV(N); TJ(nz)=iM(N-1); TS(nz)=-0.5*dt;

    J=sparse(TI(1:nz),TJ(1:nz),TS(1:nz),4*N,4*N);
    reg=1e-6*max(max(abs(diag(J))),1);
    delta=-(J+reg*speye(4*N))\res;

    rn=norm(res,inf);
    full_nrm=norm(res,2);
    feas_nrm=norm(res(feas_idx),2);
    alpha=1;
    accepted=false;
    for ls=1:20
        trial=sol+alpha*delta;
        if all(trial(1:N)>0.02) && all(trial(1:N)<2)
            trial(iR(1))=1;
            trial(iV(1))=0;
            trial(iL(1))=0;
            trial(iL(N))=0;
            trial(iM(N))=0;
            trial(1:N)=min(max(trial(1:N),0.005),2);

            trial_res = kkt_residual_only_exactstress( ...
                N,dt,c,We,Req,Ca,Re,R_obs,V_obs,pb,dpb,w_obs_R,w_obs_V,trial);
            trial_rn   = norm(trial_res,inf);
            trial_full = norm(trial_res,2);
            trial_feas = norm(trial_res(feas_idx),2);
            if (trial_feas < feas_nrm) || ...
               (trial_feas <= feas_nrm && trial_full < full_nrm)
                sol = trial;
                rn = trial_rn;
                accepted = true;
                break;
            end
        end
        alpha=alpha*0.5;
    end
    if ~accepted, break; end

    if rn<1e-8, break; end
end

R_out=sol(iR(1):iR(N))';
V_out=sol(iV(1):iV(N))';
[S_out,~]=exact_stress_state(Req, Ca, Re, R_out, V_out);
end
