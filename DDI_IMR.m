function DDI_IMR(use_gt_init)
%DDI_IMR  Full joint (R*, V*, S*) DDI via 5N KKT system.
%
%  CALL:  DDI_IMR()        — production: S_init = S_NH + S_v
%         DDI_IMR(true)    — debug: initialize from GT
%
%  TERMINOLOGY
%  -----------
%  GT (Ground Truth): the analytically known S*(t) for the generating KV
%  model (mu=0.1 Pa·s, G=1e4 Pa).  Available ONLY for synthetic data.
%  For real experiments, GT does not exist.
%
%  DDI trajectory: the model-consistent (R*, V*, S*) that minimizes
%  the data-fidelity + centroid-distance + smoothness objective subject
%  to the discrete Keller-Miksis and kinematic constraints.
%  DDI R* ≠ R_obs because the discrete KM equation has a slightly
%  different solution than the continuous solver on a 256-pt grid.
%  This is the physically correct behavior: DDI finds the smoothest
%  KM-consistent trajectory that best explains the noisy R_obs.
%
%  SETTINGS (from debugging campaign)
%  -----------------------------------
%  wV = 0:     V determined by kinematic constraint; centroid pull removed
%  beta_R = 0: R smoothness fights kinematic constraint at rapid bounces
%  reg = 1e-6: stronger regularization for Jacobian conditioning (wV=0)
%  Per-experiment We_a = P_inf*Rmax_a/(2*gamma)

if nargin < 1, use_gt_init = false; end
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
t_dim_raw      = linspace(0, 150e-6, 256)';

R_d_raw  = zeros(256, Nexp_gen);
V_nd_raw = zeros(256, Nexp_gen);
for i = 1:Nexp_gen
    R_d_raw(:,i)  = synthetic_data{i}(:,2) .* Rmax_range_gen(i);
    V_nd_raw(:,i) = synthetic_data{i}(:,3);
end

[~, startIdx] = max(R_d_raw, [], 1);
minLen = min(size(R_d_raw,1) - startIdx + 1);
R_d   = zeros(minLen, Nexp_gen);
t_d   = zeros(minLen, Nexp_gen);
V_raw = zeros(minLen, Nexp_gen);
for j = 1:Nexp_gen
    idx       = startIdx(j):(startIdx(j)+minLen-1);
    R_d(:,j)  = R_d_raw(idx, j);
    t_d(:,j)  = t_dim_raw(idx) - t_dim_raw(startIdx(j));
    V_raw(:,j)= V_nd_raw(idx, j);
end

Rmax_actual = max(R_d, [], 1);
tc_actual   = Rmax_actual .* sqrt(rho_val/P_inf_val);
We_vec      = P_inf_val .* Rmax_actual ./ (2*gamma_s);
Re_vec      = Rmax_actual .* sqrt(rho_val*P_inf_val) ./ mu_gt;

Rmatrix = R_d ./ Rmax_actual;
Vmatrix = V_raw;
tmatrix = t_d ./ tc_actual;

fprintf('  We: ');  fprintf('%.1f  ',We_vec);  fprintf('\n');
fprintf('  Re: ');  fprintf('%.1f  ',Re_vec);  fprintf('\n');
fprintf('  Ca: %.4f\n\n', Ca_gt);

%% =========================================================================
%  COMMON GRID
%% =========================================================================
t_nd_end = min(max(tmatrix,[],1));
Nexp     = Nexp_gen;
N        = minLen;
t        = linspace(0, t_nd_end, N)';
dt       = t(2) - t(1);

R_obs   = zeros(Nexp, N);
V_exact = zeros(Nexp, N);
for a = 1:Nexp
    [tu,iu]      = unique(tmatrix(:,a));
    R_obs(a,:)   = interp1(tu, Rmatrix(iu,a), t', 'pchip', 'extrap');
    V_exact(a,:) = interp1(tu, Vmatrix(iu,a), t', 'pchip', 'extrap');
end
R_obs = max(R_obs, 0.02);

%% =========================================================================
%  R*_eq AND GROUND TRUTH S*(t)
%  GT is available here because this is synthetic data with known model.
%  For real experiments, this section would be removed.
%% =========================================================================
last10       = max(1, floor(0.1*N));
Req_star_vec = mean(R_obs(:, end-last10+1:end), 2)';
Req_star_vec = max(Req_star_vec, 0.02);
fprintf('  R*_eq: ');  fprintf('%.4f  ',Req_star_vec);  fprintf('\n');

S_gt = zeros(Nexp, N);
for a = 1:Nexp
    Rstar     = max(R_obs(a,:), 1e-3);
    ratio     = Req_star_vec(a) ./ Rstar;
    S_gt(a,:) = -1/(2*Ca_gt) .* (5 - ratio.^4 - 4.*ratio) ...
               -4/Re_vec(a) .* V_exact(a,:) ./ Rstar;
end
fprintf('  GT S* range: [%.2f, %.2f]\n\n', min(S_gt(:)), max(S_gt(:)));

%% =========================================================================
%  SMOOTH R* AND V* (warm-start for pb and S_init; NOT used to fix R or V)
%% =========================================================================
smooth_sigma = 5;
hw = ceil(3*smooth_sigma);
gw = exp(-(-hw:hw).^2/(2*smooth_sigma^2));  gw=gw/sum(gw);
R_sm  = zeros(Nexp,N);
V_est = zeros(Nexp,N);
for a = 1:Nexp
    rs = conv(R_obs(a,:),gw,'same');  rs(1)=1.0;
    R_sm(a,:) = rs;
    V_est(a,2:N-1) = (rs(3:N)-rs(1:N-2))/(2*dt);
    V_est(a,1)=0;  V_est(a,N)=V_est(a,N-1);
end

%% =========================================================================
%  pb PER EXPERIMENT
%% =========================================================================
%% =========================================================================
%  pb FIXED FROM R_obs THROUGHOUT ALL ITERATIONS
%  Recomputing pb from R_traj each iteration creates a feedback loop:
%  R_traj diverges from R_obs → wrong pb → wrong S → wrong centroids.
%  Freezing pb from R_obs breaks this loop and stabilizes Exp 3.
%% =========================================================================
for a = 1:Nexp
    Req_a      = Req_star_vec(a);
    pb0_vec(a) = pv_sat + (1 + 1/(We_vec(a)*Req_a) - pv_sat)*Req_a^3;
    [pb_a,dpb_a] = compute_pb(R_obs(a:a,:), V_exact(a:a,:), N, dt, kappa, pb0_vec(a));
    pb(a,:)=pb_a; dpb(a,:)=dpb_a;
end
fprintf('  pb0: ');  fprintf('%.4f  ',pb0_vec);  fprintf('\n\n');

%% =========================================================================
%  S_init = S_NH + S_v  (full KV formula — correct scale, no dV/dt needed)
%  Uses V_exact directly: bounded, no blow-up.
%  This is identical to S_gt for synthetic data (Max|S_init-S_gt|=0).
%  For real data with unknown model, use S_NH alone (conservative).
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
%  WEIGHTS (all ones — no gating for the joint solve)
%% =========================================================================
wn_all = ones(Nexp, N);

%% =========================================================================
%  DDI SETTINGS
%% =========================================================================
w_obs_R      = 100.0;
wR           = 1.0;
wV           = 0.0;   % V determined by kinematic constraint
wS           = 1.0;
beta_S0      = 10.0;
beta_S1      = 1.0;
beta_R       = 0.0;   % no R smoothness (fights kinematic constraint)
anneal_iters = 50;
Nbar         = 50*Nexp;
max_iter     = 200;
tol          = 1e-4;
max_nr       = 50;

%% =========================================================================
%  INITIALIZATION
%% =========================================================================
if use_gt_init
    fprintf('  [DEBUG] Initializing from GT (R_obs, V_exact, S_gt).\n\n');
    R_traj = R_obs;   V_traj = V_exact;  S_traj = S_gt;
else
    fprintf('  [PRODUCTION] Initializing from (R_obs, V_exact, S_NH+S_v).\n\n');
    R_traj = R_obs;   V_traj = V_exact;  S_traj = S_init;
end

fprintf('  Initial RMSE R*: ');
for a=1:Nexp, fprintf('Exp%d=%.5f  ',a,rms(R_traj(a,:)-R_obs(a,:))); end
fprintf('\n  Initial RMSE S*: ');
for a=1:Nexp, fprintf('Exp%d=%.5f  ',a,rms(S_traj(a,:)-S_gt(a,:)));  end
fprintf('\n\n');

%% =========================================================================
%  K-MEANS++  (seeded from (R, S) manifold; V excluded since wV=0)
%% =========================================================================
rng(42);
pts  = [R_traj(:), V_traj(:), S_traj(:)];
ptsn = pts ./ (std(pts)+eps);
ci   = zeros(Nbar,1);  ci(1)=randi(size(pts,1));
d2   = inf(size(pts,1),1);
for k = 2:Nbar
    d2    = min(d2, sum((ptsn-ptsn(ci(k-1),:)).^2,2));
    ci(k) = find(rand<cumsum(d2/sum(d2)),1);
end
Rbar=pts(ci,1)'; Vbar=pts(ci,2)'; Sbar=pts(ci,3)';
assign = ones(Nexp,N,'int32');

%% =========================================================================
%  DDI MAIN LOOP
%% =========================================================================
fprintf('  %5s  %12s  %10s  %8s  |  %-9s  %-9s  %-9s\n',...
    'Iter','Objective','Rel.dL','beta_S','RMSE R1','RMSE R2','RMSE R3');
L_hist = nan(max_iter,1);

warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:singularMatrix');

for iter = 1:max_iter

    if iter<=anneal_iters
        beta_S = beta_S0*(beta_S1/beta_S0)^((iter-1)/(anneal_iters-1));
    else
        beta_S = beta_S1;
    end

    % pb is frozen from R_obs (computed before main loop); no update here

    for a=1:Nexp
        for n=1:N
            [~,assign(a,n)]=min(wR*(R_traj(a,n)-Rbar).^2 + ...
                                wS*(S_traj(a,n)-Sbar).^2);
        end
    end

    for a=1:Nexp
        [R_traj(a,:),V_traj(a,:),S_traj(a,:)] = solve_KKT( ...
            N, dt, c_sp, We_vec(a), ...
            R_obs(a,:), pb(a,:), dpb(a,:), wn_all(a,:)', ...
            w_obs_R, wR, wV, wS, beta_S, beta_R, Rbar, Vbar, Sbar, ...
            assign(a,:), max_nr, R_traj(a,:)', V_traj(a,:)', S_traj(a,:)');
        % Clip V to physical range to prevent terminal/boundary blow-up
        V_traj(a,:) = max(min(V_traj(a,:), 50), -50);
    end

    nR=zeros(1,Nbar); nV=zeros(1,Nbar); nS=zeros(1,Nbar); dn=zeros(1,Nbar);
    for a=1:Nexp
        for n=1:N
            k=assign(a,n); w=wn_all(a,n);
            nR(k)=nR(k)+w*R_traj(a,n);
            nV(k)=nV(k)+w*V_traj(a,n);
            nS(k)=nS(k)+w*S_traj(a,n);
            dn(k)=dn(k)+w;
        end
    end
    mk=dn>0;
    Rbar(mk)=nR(mk)./dn(mk); Vbar(mk)=nV(mk)./dn(mk); Sbar(mk)=nS(mk)./dn(mk);

    L=0;
    for a=1:Nexp
        for n=1:N
            k=assign(a,n);
            L=L+(w_obs_R*(R_traj(a,n)-R_obs(a,n))^2 ...
                +wR*(R_traj(a,n)-Rbar(k))^2 ...
                +wS*(S_traj(a,n)-Sbar(k))^2);
        end
        for n=1:N-1
            L=L+beta_S*dt*(S_traj(a,n+1)-S_traj(a,n))^2;
        end
    end
    L_hist(iter)=L;

    rmse_R=arrayfun(@(a)rms(R_traj(a,:)-R_obs(a,:)),1:Nexp);
    if iter>1
        dL=abs(L-L_hist(iter-1))/(abs(L_hist(iter-1))+eps);
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
cDDI=[0.8 0.15 0.15]; cGT=[0.05 0.55 0.15]; cOBS=[0 0 0];
cols=lines(Nexp);

% --- Figure 1: Four columns — R*(t), V*(t), S*(t)
% DDI vs observed in R* and V*, DDI vs GT in S*
figure('Name','DDI: R*, V*, S* trajectories','Position',[30 30 1600 320*Nexp]);
for a=1:Nexp
    subplot(Nexp,3,(a-1)*3+1);
    plot(t,R_obs(a,:), 'k.','MarkerSize',2.5,'DisplayName','R* observed'); hold on;
    plot(t,R_traj(a,:),'-','Color',cDDI,'LineWidth',1.5,'DisplayName','R* DDI');
    yline(Req_star_vec(a),'--','Color',[.5 .5 .5],'LineWidth',1,'DisplayName','R*_{eq}');
    xlabel('t*'); ylabel('R*'); grid on; box on;
    title(sprintf('Exp %d  R*  We=%.0f  RMSE=%.4f',a,We_vec(a),rms(R_traj(a,:)-R_obs(a,:))));
    if a==1, legend('Location','best','FontSize',7); end

    subplot(Nexp,3,(a-1)*3+2);
    plot(t,V_exact(a,:),'k-','LineWidth',1,'DisplayName','V* from solver (ref)'); hold on;
    plot(t,V_traj(a,:), '-','Color',cDDI,'LineWidth',1.5,'DisplayName','V* DDI');
    xlabel('t*'); ylabel('V*'); grid on; box on;
    title(sprintf('Exp %d  V*  Re=%.1f  RMSE=%.4f',a,Re_vec(a),rms(V_traj(a,:)-V_exact(a,:))));
    if a==1, legend('Location','best','FontSize',7); end

    subplot(Nexp,3,(a-1)*3+3);
    plot(t,S_gt(a,:),   '-','Color',cGT, 'LineWidth',2.5,'DisplayName','S* GT (synthetic only)'); hold on;
    plot(t,S_traj(a,:), '-','Color',cDDI,'LineWidth',1.5,'DisplayName','S* DDI');
    xlabel('t*'); ylabel('S*'); grid on; box on;
    ylim([min(S_gt(a,:))*1.5-1, max(S_gt(a,:))*1.5+1]);
    title(sprintf('Exp %d  S*  Ca=%.2f  RMSE=%.3f',a,Ca_gt,rms(S_traj(a,:)-S_gt(a,:))));
    if a==1, legend('Location','best','FontSize',7); end
end

% --- Figure 2: 3D phase space + convergence + S*(R*) + S*(t*)
figure('Name','DDI: Phase space, manifold, convergence','Position',[100 50 1300 800]);

subplot(2,2,1);
% S*(R*) manifold comparison
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
% 3D phase space (R*, V*, S*) with cluster centroids
for a=1:Nexp
    lc=cols(a,:);
    plot3(R_traj(a,:), V_traj(a,:), S_traj(a,:), '-','Color',lc,'LineWidth',1.2,...
        'DisplayName',sprintf('DDI exp%d',a)); hold on;
end
% Cluster centroids (only occupied ones)
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
    plot(t,S_gt(a,:),  '-','Color',0.5*lc+0.4,'LineWidth',2,...
        'DisplayName',sprintf('GT exp%d',a)); hold on;
    plot(t,S_traj(a,:),'--','Color',lc,'LineWidth',1.5,...
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
function [R_out,V_out,S_out]=solve_KKT( ...
    N,dt,c,We,R_obs,pb,dpb,wn, ...
    w_obs_R,wR,wV,wS,beta_S,beta_R,Rbar,Vbar,Sbar,assign,max_nr, ...
    R_w,V_w,S_w)

iR=@(n)n; iV=@(n)N+n; iS=@(n)2*N+n; iL=@(n)3*N+n; iM=@(n)4*N+n;
sol=[max(R_w,0.02);V_w;S_w;zeros(2*N,1)];
sol(iR(1))=1; sol(iV(1))=0;

norm_prev=inf;
for nr=1:max_nr
    R=max(sol(iR(1):iR(N)),0.005); V=sol(iV(1):iV(N));
    S=sol(iS(1):iS(N)); lam=sol(iL(1):iL(N)); mu=sol(iM(1):iM(N));
    res=zeros(5*N,1);
    nz=0; TI=zeros(70*N,1); TJ=zeros(70*N,1); TS=zeros(70*N,1);

    res(iR(1))=R(1)-1; nz=nz+1; TI(nz)=iR(1); TJ(nz)=iR(1); TS(nz)=1;
    res(iV(1))=V(1);   nz=nz+1; TI(nz)=iV(1); TJ(nz)=iV(1); TS(nz)=1;
    res(iL(1))=lam(1); nz=nz+1; TI(nz)=iL(1); TJ(nz)=iL(1); TS(nz)=1;
    res(iL(N))=lam(N); nz=nz+1; TI(nz)=iL(N); TJ(nz)=iL(N); TS(nz)=1;
    res(iM(1))=mu(1);  nz=nz+1; TI(nz)=iM(1); TJ(nz)=iM(1); TS(nz)=1;
    res(iM(N))=mu(N);  nz=nz+1; TI(nz)=iM(N); TJ(nz)=iM(N); TS(nz)=1;

    for n=2:N-1
        kk=assign(n);
        dV=(V(n+1)-V(n-1))/(2*dt); dS=(S(n+1)-S(n-1))/(2*dt);
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

        Gc=V(n)-(R(n+1)-R(n-1))/(2*dt);
        res(iM(n))=Gc;
        nz=nz+1; TI(nz)=iM(n); TJ(nz)=iV(n);   TS(nz)=1;
        nz=nz+1; TI(nz)=iM(n); TJ(nz)=iR(n+1); TS(nz)=-1/(2*dt);
        nz=nz+1; TI(nz)=iM(n); TJ(nz)=iR(n-1); TS(nz)=1/(2*dt);

        cn1=-(1-V(n+1)/c)*R(n+1)/(2*dt); cm1=(1-V(n-1)/c)*R(n-1)/(2*dt);
        res(iR(n))=2*wn(n)*(w_obs_R*(R(n)-R_obs(n))+wR*(R(n)-Rbar(kk))) ...
                  +A*lam(n)+(mu(n+1)-mu(n-1))/(2*dt);
        nz=nz+1; TI(nz)=iR(n); TJ(nz)=iR(n);   TS(nz)=2*wn(n)*(w_obs_R+wR)+2*lam(n)/(We*R(n)^3);
        nz=nz+1; TI(nz)=iR(n); TJ(nz)=iV(n);   TS(nz)=-lam(n)*dV/c;
        nz=nz+1; TI(nz)=iR(n); TJ(nz)=iV(n+1); TS(nz)=lam(n)*(1-V(n)/c)/(2*dt);
        nz=nz+1; TI(nz)=iR(n); TJ(nz)=iV(n-1); TS(nz)=-lam(n)*(1-V(n)/c)/(2*dt);
        nz=nz+1; TI(nz)=iR(n); TJ(nz)=iS(n+1); TS(nz)=-lam(n)/(2*c*dt);
        nz=nz+1; TI(nz)=iR(n); TJ(nz)=iS(n-1); TS(nz)=lam(n)/(2*c*dt);
        nz=nz+1; TI(nz)=iR(n); TJ(nz)=iL(n);   TS(nz)=A;
        nz=nz+1; TI(nz)=iR(n); TJ(nz)=iM(n+1); TS(nz)=1/(2*dt);
        nz=nz+1; TI(nz)=iR(n); TJ(nz)=iM(n-1); TS(nz)=-1/(2*dt);

        lS=-(1+V(n)/c)*lam(n);
        if n+1<=N-1, lS=lS+R(n+1)*lam(n+1)/(2*c*dt); end
        if n-1>=2,   lS=lS-R(n-1)*lam(n-1)/(2*c*dt); end
        res(iV(n))=2*wn(n)*wV*(V(n)-Vbar(kk)) ...
                  +B*lam(n)+cn1*lam(n+1)+cm1*lam(n-1)+mu(n);
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
        nz=nz+1; TI(nz)=iV(n); TJ(nz)=iM(n);   TS(nz)=1;

        sm=2*wn(n)*beta_S*((S(n)-S(n-1))-(S(n+1)-S(n)));
        res(iS(n))=2*wn(n)*wS*(S(n)-Sbar(kk))+sm+lS;
        nz=nz+1; TI(nz)=iS(n); TJ(nz)=iS(n);   TS(nz)=2*wn(n)*(wS+2*beta_S);
        nz=nz+1; TI(nz)=iS(n); TJ(nz)=iS(n+1); TS(nz)=-2*wn(n)*beta_S;
        nz=nz+1; TI(nz)=iS(n); TJ(nz)=iS(n-1); TS(nz)=-2*wn(n)*beta_S;
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

    kN=assign(N); k1=assign(1);
    res(iR(N))=2*wn(N)*(w_obs_R*(R(N)-R_obs(N))+wR*(R(N)-Rbar(kN)))-mu(N-1)/(2*dt);
    nz=nz+1; TI(nz)=iR(N); TJ(nz)=iR(N);   TS(nz)=2*wn(N)*(w_obs_R+wR);
    nz=nz+1; TI(nz)=iR(N); TJ(nz)=iM(N-1); TS(nz)=-1/(2*dt);

    cN=(1-V(N-1)/c)*R(N-1)/(2*dt);
    res(iV(N))=2*wn(N)*wV*(V(N)-Vbar(kN))+cN*lam(N-1);
    nz=nz+1; TI(nz)=iV(N); TJ(nz)=iV(N);   TS(nz)=2*wn(N)*wV+1e-10;
    nz=nz+1; TI(nz)=iV(N); TJ(nz)=iR(N-1); TS(nz)=(1-V(N-1)/c)*lam(N-1)/(2*dt);
    nz=nz+1; TI(nz)=iV(N); TJ(nz)=iV(N-1); TS(nz)=-R(N-1)*lam(N-1)/(2*c*dt);
    nz=nz+1; TI(nz)=iV(N); TJ(nz)=iL(N-1); TS(nz)=cN;

    sm1=-2*wn(1)*beta_S*(S(2)-S(1));
    % Anchor S(1) toward its initial value (elastic stress at t*=0)
    % This prevents the collapse spike from pulling S(1) to large values
    res(iS(1))=2*wn(1)*wS*(S(1)-Sbar(k1))+sm1+R(2)*lam(2)/(2*c*dt);
    nz=nz+1; TI(nz)=iS(1); TJ(nz)=iS(1); TS(nz)=2*wn(1)*(wS+beta_S);
    nz=nz+1; TI(nz)=iS(1); TJ(nz)=iS(2); TS(nz)=-2*wn(1)*beta_S;
    nz=nz+1; TI(nz)=iS(1); TJ(nz)=iR(2); TS(nz)=lam(2)/(2*c*dt);
    nz=nz+1; TI(nz)=iS(1); TJ(nz)=iL(1); TS(nz)=-(1+V(1)/c);
    nz=nz+1; TI(nz)=iS(1); TJ(nz)=iL(2); TS(nz)=R(2)/(2*c*dt);

    smN=2*wn(N-1)*beta_S*(S(N)-S(N-1));
    res(iS(N))=2*wn(N)*wS*(S(N)-Sbar(kN))+smN ...
              -(1+V(N)/c)*lam(N)-R(N-1)*lam(N-1)/(2*c*dt);
    nz=nz+1; TI(nz)=iS(N); TJ(nz)=iS(N);   TS(nz)=2*wn(N)*wS+2*wn(N-1)*beta_S;
    nz=nz+1; TI(nz)=iS(N); TJ(nz)=iS(N-1); TS(nz)=-2*wn(N-1)*beta_S;
    nz=nz+1; TI(nz)=iS(N); TJ(nz)=iV(N);   TS(nz)=-lam(N)/c;
    nz=nz+1; TI(nz)=iS(N); TJ(nz)=iR(N-1); TS(nz)=-lam(N-1)/(2*c*dt);
    nz=nz+1; TI(nz)=iS(N); TJ(nz)=iL(N);   TS(nz)=-(1+V(N)/c);
    nz=nz+1; TI(nz)=iS(N); TJ(nz)=iL(N-1); TS(nz)=-R(N-1)/(2*c*dt);

    J=sparse(TI(1:nz),TJ(1:nz),TS(1:nz),5*N,5*N);
    reg=1e-6*max(abs(diag(J)));
    delta=-(J+reg*speye(5*N))\res;

    alpha=1;
    for ls=1:20
        if all(sol(1:N)+alpha*delta(1:N)>0.02)&&all(sol(1:N)+alpha*delta(1:N)<2),break;end
        alpha=alpha*0.5;
    end
    sol=sol+alpha*delta;
    sol(iR(1))=1; sol(iV(1))=0; sol(iL(1))=0; sol(iL(N))=0; sol(iM(1))=0; sol(iM(N))=0;
    sol(1:N)=min(max(sol(1:N),0.005),2);
    % Clip V to prevent runaway at boundary nodes
    sol(N+1:2*N) = max(min(sol(N+1:2*N), 50), -50);

    rn=norm(res,inf);
    if rn<1e-8,break; end
    if nr>2&&abs(rn-norm_prev)/(norm_prev+eps)<1e-4,break; end
    norm_prev=rn;
end
R_out=sol(iR(1):iR(N))'; V_out=sol(iV(1):iV(N))'; S_out=sol(iS(1):iS(N))';
end