clc; clearvars;

%% -------------------- Load SI workbook --------------------
xlsPath = "Material_Database_SI.xlsx";

% Use detectImportOptions and PRESERVE the original (numeric) column names
optsL = detectImportOptions(xlsPath, 'Sheet','Lamina_SI');
optsL.VariableNamingRule = 'preserve';
LaminaSI = readtable(xlsPath, optsL);

optsF = detectImportOptions(xlsPath, 'Sheet','Fiber_SI');
optsF.VariableNamingRule = 'preserve';
FiberSI = readtable(xlsPath, optsF);

optsM = detectImportOptions(xlsPath, 'Sheet','Matrix_SI');
optsM.VariableNamingRule = 'preserve';
MatrixSI = readtable(xlsPath, optsM);

getProp = @(T, matId, propName) T{strcmpi(string(T.Property), string(propName)), string(matId)};

%% -------------------- Layup inputs --------------------
Theta_list = input('Enter ply angles (deg), e.g. [0 90 0 90]: ');
t_input_mm = input('Enter single-ply thickness in mm (scalar or vector): ');
if isscalar(t_input_mm)
    Ply_Thk_m = repmat(t_input_mm/1000, 1, numel(Theta_list));  % mm -> m
else
    Ply_Thk_m = (t_input_mm(:).'/1000);
    if numel(Ply_Thk_m) ~= numel(Theta_list)
        error('Length of thickness vector must match Theta_list.');
    end
end
n = numel(Theta_list);

%% -------------------- Material properties --------------------
useLamina = questdlg('Use lamina properties?', 'Material properties', 'Yes','No','Yes');
if isempty(useLamina), error('Selection cancelled.'); end

if strcmp(useLamina,'Yes')
    Lamina_ID = input('Enter Lamina numeric ID (e.g., 301): ');
    if ~ismember(string(Lamina_ID), string(LaminaSI.Properties.VariableNames))
        error('Lamina ID %s not found in Lamina_SI sheet.');
    end

    E1_GPa  = getProp(LaminaSI, Lamina_ID, 'E1_GPa');
    E2_GPa  = getProp(LaminaSI, Lamina_ID, 'E2_GPa');
    G12_GPa = getProp(LaminaSI, Lamina_ID, 'G12_GPa');
    v12     = getProp(LaminaSI, Lamina_ID, 'v12_-');
    a1      = getProp(LaminaSI, Lamina_ID, 'alpha1_perC');
    a2      = getProp(LaminaSI, Lamina_ID, 'alpha2_perC');
    v21     = (v12 * E2_GPa) / E1_GPa;

    % Strengths
    Xt_MPa  = getProp(LaminaSI, Lamina_ID, 'Xt_MPa');
    Xc_MPa  = getProp(LaminaSI, Lamina_ID, 'Xc_MPa');
    Yt_MPa  = getProp(LaminaSI, Lamina_ID, 'Yt_MPa');
    Yc_MPa  = getProp(LaminaSI, Lamina_ID, 'Yc_MPa');
    S12_MPa = getProp(LaminaSI, Lamina_ID, 'S12_MPa');

    %% PATCH: if alpha stored as microstrain/°C, convert to 1/°C
    a1 = a1 * 1e-6;
    a2 = a2 * 1e-6;

else
    % Example: Fiber 101, Matrix 201
    Fiber_ID  = input('Enter Fiber numeric ID (e.g., 101): ');
    Matrix_ID = input('Enter Matrix numeric ID (e.g., 201): ');

    if ~ismember(string(Fiber_ID),  string(FiberSI.Properties.VariableNames))
        error('Fiber ID %s not found in Fiber_SI sheet.', string(Fiber_ID));
    end
    if ~ismember(string(Matrix_ID), string(MatrixSI.Properties.VariableNames))
        error('Matrix ID %s not found in Matrix_SI sheet.', string(Matrix_ID));
    end

    Vf = input('Enter fiber volume fraction Vf (0-1): ');
    if ~(isscalar(Vf) && Vf>0 && Vf<1)
        error('Vf must be a scalar in (0,1).');
    end
    Vm = 1 - Vf;

    % Fiber props
    E1f = getProp(FiberSI,  Fiber_ID, 'EfL_GPa');
    E2f = getProp(FiberSI,  Fiber_ID, 'EfT_GPa');
    Gf  = getProp(FiberSI,  Fiber_ID, 'GfLT_GPa');
    vf  = getProp(FiberSI,  Fiber_ID, 'vfLT_-');
    a1f = getProp(FiberSI,  Fiber_ID, 'afL_perC');
    a2f = getProp(FiberSI,  Fiber_ID, 'afT_perC');

    % Matrix props
    Em  = getProp(MatrixSI, Matrix_ID, 'Em_GPa');
    Gm  = getProp(MatrixSI, Matrix_ID, 'Gm_GPa');
    vm  = getProp(MatrixSI, Matrix_ID, 'vm_-');
    am  = getProp(MatrixSI, Matrix_ID, 'am_perC');

    %% PATCH: convert constituent alphas if stored as microstrain/°C
    a1f = a1f * 1e-6;
    a2f = a2f * 1e-6;
    am  = am  * 1e-6;

    % Mixture rules (E2 uses E2f via inverse ROM)
    E1_GPa  = Em*Vm + E1f*Vf;
    E2_GPa  = (E2f*Em) / (Em*Vf + E2f*Vm);       % inverse ROM
    G12_GPa = (Gf * Gm) / (Gm*Vf + Gf*Vm);       % inverse ROM
    v12     = vm*Vm + vf*Vf;
    v21     = (v12 * E2_GPa) / E1_GPa;

    % Longitudinal/transverse CTEs (stiffness-weighted)
    a1 = (am*Em*Vm + a1f*E1f*Vf) / E1_GPa;
    a2 = (am*Em*Vm + a2f*E2f*Vf) / (Em*Vm + E2f*Vf);

    % Strengths (mapped simply)
    Xt_MPa  = getProp(FiberSI,  Fiber_ID,  'FFt_MPa');
    Xc_MPa  = getProp(FiberSI,  Fiber_ID,  'FFC_MPa');
    Yt_MPa  = getProp(MatrixSI, Matrix_ID, 'Fmt_MPa');
    Yc_MPa  = getProp(MatrixSI, Matrix_ID, 'Fmc_MPa');
    S12_MPa = getProp(MatrixSI, Matrix_ID, 'Fms_MPa');
end

%% -------------------- Units to SI base (Pa, m) --------------------
E1  = E1_GPa  * 1e9;     % GPa -> Pa
E2  = E2_GPa  * 1e9;
G12 = G12_GPa * 1e9;

Xt  = Xt_MPa  * 1e6;     % MPa -> Pa
Xc  = Xc_MPa  * 1e6;
Yt  = Yt_MPa  * 1e6;
Yc  = Yc_MPa  * 1e6;
S12 = S12_MPa * 1e6;

%% -------------------- Loads & ΔT --------------------
dT_C = input('Enter Temperature Change (°C): ');
force_Npm  = input('Enter in-plane resultants [Nx Ny Nxy] in N/m: ');
nx = force_Npm(1); ny = force_Npm(2); nxy = force_Npm(3);
moment_N   = input('Enter bending resultants [Mx My Mxy] in N: ');
mx = moment_N(1); my = moment_N(2); mxy = moment_N(3);

%% -------------------- Qbar for each ply --------------------
Q11 = E1 / (1 - v12*v21);
Q12 = v12*E2 / (1 - v12*v21);
Q22 = E2 / (1 - v12*v21);
Q66 = G12;

Qbar = zeros(3,3,n);
for i = 1:n
    th = Theta_list(i);
    m = cosd(th); n_ = sind(th);
    m2 = m^2; n2 = n_^2; mn = m*n_;

    Q11b = Q11*m2^2 + 2*(Q12 + 2*Q66)*m2*n2 + Q22*n2^2;
    Q22b = Q11*n2^2 + 2*(Q12 + 2*Q66)*m2*n2 + Q22*m2^2;
    Q12b = (Q11 + Q22 - 4*Q66)*m2*n2 + Q12*(m2^2 + n2^2);
    Q16b = (Q11 - Q12 - 2*Q66)*m2*mn - (Q22 - Q12 - 2*Q66)*n2*mn;
    Q26b = (Q11 - Q12 - 2*Q66)*n2*mn - (Q22 - Q12 - 2*Q66)*m2*mn;
    Q66b = (Q11 + Q22 - 2*Q12 - 2*Q66)*m2*n2 + Q66*(m2^2 + n2^2);

    Qbar(:,:,i) = [Q11b Q12b Q16b; Q12b Q22b Q26b; Q16b Q26b Q66b];
end

%% -------------------- z coordinates (m) --------------------
tTot = sum(Ply_Thk_m);
z_bot = -tTot/2;
Z_k = zeros(n,2);
for i = 1:n
    z_top = z_bot + Ply_Thk_m(i);
    Z_k(i,:) = [z_bot, z_top];
    z_bot = z_top;
end

%% -------------------- Thermal loads (global)
alpha_bar = zeros(3,1,n);
Nt_unsum  = zeros(3,n);
Mt_unsum  = zeros(3,n);

for i = 1:n
    th = Theta_list(i);
    m = cosd(th); n_ = sind(th);
    m2 = m^2; n2 = n_^2; mn = m*n_;

    alphax  = a1*m2 + a2*n2;
    alphay  = a2*m2 + a1*n2;
    alphaxy = 2*(a1 - a2)*mn;
    alpha_bar(:,:,i) = [alphax; alphay; alphaxy];

    Nt_unsum(:,i) = Qbar(:,:,i) * alpha_bar(:,:,i) * (Z_k(i,2) - Z_k(i,1));
    Mt_unsum(:,i) = Qbar(:,:,i) * alpha_bar(:,:,i) * (Z_k(i,2)^2 - Z_k(i,1)^2);
end

Nt = dT_C * sum(Nt_unsum, 2);         % N/m
Mt = 0.5 * dT_C * sum(Mt_unsum, 2);   % N

%% -------------------- ABD --------------------
A = zeros(3); B = zeros(3); D = zeros(3);
for i = 1:n
    dz  = (Z_k(i,2) - Z_k(i,1));
    dz2 = (Z_k(i,2)^2 - Z_k(i,1)^2);
    dz3 = (Z_k(i,2)^3 - Z_k(i,1)^3);
    A = A + Qbar(:,:,i) * dz;
    B = B + Qbar(:,:,i) * (0.5*dz2);
    D = D + Qbar(:,:,i) * (dz3/3);
end
ABD = [A B; B D];

%% -------------------- Strains/curvatures --------------------
% PATCH: free thermal expansion (N=-Nt, M=-Mt) or constrained (N=+Nt, M=+Mt)
isFreeThermal = False;   % set true for free thermal response; false for fully constrained thermal loads

if isFreeThermal
    N = [nx; ny; nxy] - Nt;    % Free thermal: total resultants = 0  => A*eps + ... + Nt = 0
    M = [mx; my; mxy] - Mt;
else
    N = [nx; ny; nxy] + Nt;    % Constrained thermal: apply Nt and Mt as external loads
    M = [mx; my; mxy] + Mt;
end

ek = ABD \ [N; M];
epsilon_0 = ek(1:3);
kappa     = ek(4:6);

%% -------------------- LOCAL stresses (σ1,σ2,τ12) --------------------
Q_local = [Q11, Q12, 0;  Q12, Q22, 0;  0, 0, Q66];
alpha_local = [a1; a2; 0];
T_eps = @(th) [ cosd(th)^2,            sind(th)^2,            cosd(th)*sind(th);
                sind(th)^2,            cosd(th)^2,           -cosd(th)*sind(th);
               -2*cosd(th)*sind(th),   2*cosd(th)*sind(th),   cosd(th)^2 - sind(th)^2 ];

sigma12_top    = zeros(3,1,n);
sigma12_bottom = zeros(3,1,n);
for i = 1:n
    th = Theta_list(i);
    Te = T_eps(th);
    eps_top_g    = epsilon_0 + Z_k(i,2)*kappa;
    eps_bottom_g = epsilon_0 + Z_k(i,1)*kappa;
    eps_top_12    = Te * eps_top_g;
    eps_bottom_12 = Te * eps_bottom_g;
    sigma12_top(:,:,i)    = Q_local * (eps_top_12    - dT_C*alpha_local);  % Pa
    sigma12_bottom(:,:,i) = Q_local * (eps_bottom_12 - dT_C*alpha_local);  % Pa
end

%% -------------------- Max-stress (local) --------------------
MSC_X = zeros(n,2);   % [tension_ratio, compression_ratio] for σ1
MSC_Y = zeros(n,2);   % [tension_ratio, compression_ratio] for σ2
MSC_S = zeros(n,2);   % [top, bottom] for |τ12|

for i = 1:n
    s1t = sigma12_top(1,1,i);   s1b = sigma12_bottom(1,1,i);
    s2t = sigma12_top(2,1,i);   s2b = sigma12_bottom(2,1,i);
    t12t= sigma12_top(3,1,i);   t12b= sigma12_bottom(3,1,i);

    r1_t = max(0, s1t)/Xt;  r1_t_b = max(0, s1b)/Xt;
    r1_c = max(0,-s1t)/abs(Xc); r1_c_b = max(0,-s1b)/abs(Xc);
    MSC_X(i,1) = max(r1_t, r1_t_b);
    MSC_X(i,2) = max(r1_c, r1_c_b);

    r2_t = max(0, s2t)/Yt;  r2_t_b = max(0, s2b)/Yt;
    r2_c = max(0,-s2t)/abs(Yc); r2_c_b = max(0,-s2b)/abs(Yc);
    MSC_Y(i,1) = max(r2_t, r2_t_b);
    MSC_Y(i,2) = max(r2_c, r2_c_b);

    MSC_S(i,1) = abs(t12t)/S12;
    MSC_S(i,2) = abs(t12b)/S12;

    if MSC_X(i,1) >= 1, fprintf('Ply %d fails fiber tension (σ1)\n', i); end
    if MSC_X(i,2) >= 1, fprintf('Ply %d fails fiber compression (σ1)\n', i); end
    if MSC_Y(i,1) >= 1, fprintf('Ply %d fails matrix tension (σ2)\n', i); end
    if MSC_Y(i,2) >= 1, fprintf('Ply %d fails matrix compression (σ2)\n', i); end
    if max(MSC_S(i,:)) >= 1, fprintf('Ply %d fails in shear (τ12)\n', i); end
end

%% -------------------- Displays --------------------
disp('ABD matrix (SI):'); disp(ABD);
disp('N resultants (incl. thermal), N/m:'), disp(N);
disp('M resultants (incl. thermal), N:'), disp(M);

% %% Checks if needed
% format long g
% disp('epsilon_0 ='); disp(epsilon_0)
% disp('kappa =');     disp(kappa)
% 
% res = ABD*[epsilon_0; kappa] + [Nt; Mt];
% fprintf('||res||_2 = %.3e\n', norm(res,2))
% fprintf('||kappa||_2 = %.3e\n', norm(kappa,2))
% 
% alpha_eff = epsilon_0 / dT_C;
% disp('alpha_eff (per °C) ='); disp(alpha_eff)
% 
% N_back = A*epsilon_0 + Nt;
% M_back = B*epsilon_0 + D*kappa + Mt;
% disp('N_back (should be ~0) ='); disp(N_back)
% disp('M_back (should be ~0) ='); disp(M_back)
% 
% 
% 
