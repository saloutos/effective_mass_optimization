%% Deriving EOMs for simplified finger contact model
% modifying IMF formulation for finger-manipuland collisions instead of
% leg-ground impacts

% Andrew SaLoutos
% 2/8/2021

% toy system consists of motor/finger (modeled as rack and pinion),
% manipuland (modeled as another mass), and a contact stiffness (modeled as
% a spring with spring contstant kc)

%% Define variables and system states

syms I_m m_f m_o   % inertias
syms k_c           % contact spring stiffness
syms r_m           % pinion radius (gear reduction)
syms th(t) x_o(t)  % state variables
syms tau_m         % motor input torque

dth = diff(th(t),t);
dx_o = diff(x_o(t),t);

%% Energy terms

T = 0.5*I_m*dth^2 + 0.5*m_f*(r_m*dth)^2 + 0.5*m_o*dx_o^2;
V = 0.5*k_c*(x_o - r_m*th)^2;

%% Lagrangian to get eoms
L = T-V;

% for motor angle
eq_th = diff(diff(L,dth),t) - diff(L,th) == tau_m;

% for object position
eq_x_o = diff(diff(L,dx_o),t) - diff(L,x_o) == 0;

%% Sort into mass matrix, stiffness matrix

M = [I_m + m_f*r_m^2, 0; 0, m_o];
K = [k_c*r_m^2, -k_c*r_m; -k_c*r_m, k_c];

%% Natural freq with finger mass, rotor, contact stiffness?

kc = 1e5;
rm = 1/2;
mo_nom = 1;
Im_nom = 1e-4;
mf_nom = 0.1;

If_nom = Im_nom + mf_nom*rm^2;

Im = logspace(-8,0);
mf = logspace(-4,0);
mo = logspace(-3,2);

If = logspace(-4,0);

[X1, Y1] = meshgrid(Im,mf);
[X2, Y2] = meshgrid(Im,mo);
[X3, Y3] = meshgrid(mf,mo);
[X4, Y4] = meshgrid(If,mo);

meff1 = ((X1 + Y1*rm^2)*mo_nom)./((X1 + Y1*rm^2) + mo_nom*rm^2);
wn1 = sqrt( kc ./ meff1 ); % in Hz
If_tot = (X1 + Y1*rm^2); % in kgm2

meff2 = ((X2 + mf_nom*rm^2).*Y2)./((X2 + mf_nom*rm^2) + Y2*rm^2);
wn2 = sqrt( kc ./ meff2 );

meff3 = ((Im_nom + X3*rm^2).*Y3)./((Im_nom + X3*rm^2) + Y3*rm^2);
wn3 = sqrt( kc ./ meff3 );

meff4 = (X4.*Y4)./(X4 + Y4*rm^2);
wn4 = sqrt( kc ./ meff4 );

figure;
subplot(1,3,1);
surf(X1,Y1,wn1);
set(gca,'xscale','log'); set(gca,'yscale','log'); 
set(gca,'colorscale','log');
view(0,90);
xlabel('Inertia (kgm^{2})'); ylabel('Finger mass (kg)');
xlim([1e-8, 1e0]); ylim([1e-4, 1e0]);
cb1 = colorbar; cb1.Label.String = 'Natural Freq (Hz)';

subplot(1,3,2);
surf(X2,Y2,wn2);
set(gca,'xscale','log'); set(gca,'yscale','log'); 
set(gca,'colorscale','log');
view(0,90);
xlabel('Inertia (kgm^{2})'); ylabel('Object mass (kg)');
xlim([1e-8, 1e0]); ylim([1e-3, 1e2]);
cb2 = colorbar; cb2.Label.String = 'Natural Freq (Hz)';

subplot(1,3,3);
surf(X3,Y3,wn3);
set(gca,'xscale','log'); set(gca,'yscale','log'); 
set(gca,'colorscale','log');
view(0,90);
xlabel('Finger mass (kg)'); ylabel('Object mass (kg)');
xlim([1e-4,1e0]); ylim([1e-3, 1e2]);
cb3 = colorbar; cb3.Label.String = 'Natural Freq (Hz)';

figure;
subplot(1,3,1);
surf(X1,Y1,If_tot);
set(gca,'xscale','log'); set(gca,'yscale','log');
set(gca,'colorscale','log');
view(0,90);
xlabel('Inertia (kgm^{2})'); ylabel('Finger mass (kg)');
xlim([1e-8, 1e0]); ylim([1e-4, 1e0]);
cb4 = colorbar; cb4.Label.String = 'Total Finger Inertia (kgm^{2})';

subplot(1,3,2);
surf(X4,Y4,meff4);
set(gca,'xscale','log'); set(gca,'yscale','log'); 
set(gca,'colorscale','log');
view(0,90);
xlabel('Total Finger Inertia (kgm^{2})'); ylabel('Object mass (kg)');
xlim([1e-4, 1]); ylim([1e-3, 1e2]);
cb5 = colorbar; cb5.Label.String = 'Effective  mass (kg)';

subplot(1,3,3);
surf(X4,Y4,wn4);
set(gca,'xscale','log'); set(gca,'yscale','log'); 
set(gca,'colorscale','log');
view(0,90);
xlabel('Total Finger Inertia (kgm^{2})'); ylabel('Object mass (kg)');
xlim([1e-4, 1]); ylim([1e-3, 1e2]);
cb6 = colorbar; cb6.Label.String = 'Natural Freq (Hz)';

figure;
surf(X4,Y4,wn4);
set(gca,'xscale','log'); set(gca,'yscale','log'); 
set(gca,'colorscale','log');
view(0,90);
xlabel('Total Finger Inertia (kgm^{2})'); ylabel('Object mass (kg)');
xlim([1e-4, 1]); ylim([1e-3, 1e2]);
cb6 = colorbar; cb6.Label.String = 'Natural Freq (Hz)';

%% Inform a design...

wn_des = 1000; % in Hz
kc_spec = 100000; % in N/m
rm_spec = 0.5; % N = 2

meff_necessary = kc_spec/(wn_des^2) % maximum effective mass

m_f_tot_spec = logspace(-4,1); % range of possible finger lumped masses
m_o_spec = logspace(-2,1); % range of desired object masses

[X5,Y5] = meshgrid(m_f_tot_spec,m_o_spec);

meff5 = (X5.*Y5)./(X5 + Y5*rm^2);
wn5 = sqrt( kc ./ meff5 );

figure;
surf(X5,Y5,wn5);
set(gca,'xscale','log'); set(gca,'yscale','log'); 
set(gca,'colorscale','log');
view(0,90);
xlabel('Lumped Finger Mass (kg)'); ylabel('Object mass (kg)');
xlim([1e-4, 1]); ylim([1e-2, 1e1]);
cb7 = colorbar; cb7.Label.String = 'Natural Freq in Collision (Hz)';




% Note: m_f_tot is defined as Im/(rm^2) + mf
meff_test = (X5.*Y5)./(X5 + Y5);

meff_good = meff_test;
meff_good(meff_test>meff_necessary) = nan;

figure;
subplot(1,2,1);
surf(X5,Y5,meff_good);
set(gca,'xscale','log'); set(gca,'yscale','log'); 
set(gca,'colorscale','log');
xlim([1e-4, 1]); ylim([1e-2, 1]);
view(0,90);
xlabel('Lumped Finger Mass (kg)'); ylabel('Object Mass (kg)');
cb8 = colorbar; cb8.Label.String = 'Effective Mass (kg)';
title_string = sprintf('Feasible Lumped Finger Masses for \\omega_{n} = %.0f Hz, m_{eff} \\leq %.4f kg', wn_des, meff_necessary);
title(title_string);

% now determine finger mass?
full_meff_cols = any(isnan(meff_good));
ind_last_full_col = find(full_meff_cols,1,'first')-1;
m_f_tot_necessary = max(meff_good(:,ind_last_full_col))

% If_necessary = 0.01; % from control at 250Hz, kc=1000, N=2 (r=0.5)

Im_spec = logspace(-6,-2);
mf_spec = logspace(-4,0);

[A,B] = meshgrid(Im_spec,mf_spec);

m_f_tot_test = (A/(rm_spec^2)) + B;
% If_test = A + B*rm_spec^2;
m_f_tot_good = m_f_tot_test;
m_f_tot_good(m_f_tot_test>m_f_tot_necessary) = nan;

subplot(1,2,2);
surf(A,B,m_f_tot_good);
set(gca,'xscale','log'); set(gca,'yscale','log'); 
set(gca,'colorscale','log');
xlim([1e-6, 1e-2]); ylim([1e-4, 1e0]);
view(0,90);
xlabel('Rotor Inertia (kgm^{2})'); ylabel('Finger Mass (kg)');
cb9 = colorbar; cb9.Label.String = 'Total Finger Mass (kg)';
title_string = sprintf('Feasible Combinatations for m_{finger,total} \\leq %.4f kg', m_f_tot_necessary);
title(title_string);


%% Other extensions...
% add damping?
% set initial conditions and simulate?
% simulate modified IMF with 2DOF finger? Lambda for a cube in 2D is pretty
% simple...
