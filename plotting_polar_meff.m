% polar plots of effective mass ellipsoids

% Andrew SaLoutos
% 5/30/2021

%% import arm functions and parameters

% UR3 functions work for arms without rotor inertias
% TODO: re-derive LIMS and large manipulator dynamic equations with rotor inertias

% addpath(genpath('UR3_arm_functions')) % for UR3 arm, these have normal names (no prefixes)
addpath(genpath('LIMS_rotor_functions'))

% UR5 planar arm
m1 = 8.058;             
m2 = 2.846;
m3 = 3.035;             
m_motor = 0.0; % not using motor mass
I1 = 0.08333*m1*(3*0.06^2 + 0.425^2);            
I2 = 0.08333*m2*(3*0.06^2 + 0.3922^2);
I3 = 0.08333*m3*(3*0.06^2 + 0.0997^2); % loosest approximation, but probably okay    
I_motor = 0.0; % not using motor inertia
Ir = 2.07e-5;           
N = 101;
Ir1 = Ir*N^2;
Ir2 = Ir*N^2;
Ir3 = Ir*N^2;
l_OA = 0.425;
l_AB = 0.3922;          
l_BC = 0.0997;
l_O_m1 = l_OA/2;         
l_A_m2 = l_AB/2;
l_B_m3 = l_BC/2;        
g = 9.81;
% parameters
% p_UR5   = [m1 m2 m3 m_motor I1 I2 I3 I_motor Ir N l_O_m1 l_A_m2 l_B_m3 l_OA l_AB l_BC g]'; 
p_UR5   = [m1 m2 m3 I1 I2 I3 Ir1 Ir2 Ir3 l_O_m1 l_A_m2 l_B_m3 l_OA l_AB l_BC g]'; % needs special dynamics
p_UR5_nr   = [m1 m2 m3 I1 I2 I3 0.0 0.0 0.0 l_O_m1 l_A_m2 l_B_m3 l_OA l_AB l_BC g]'; % needs special dynamics

% LIMS parameters (including reflected rotor inertia)
m1 = 1.28;          
m2 = 0.64;
m3 = 0.32;
m_motor = 0.0;
I1 = 18796/(1000*1000);
I2 = 6895/(1000*1000);
I3 = 694/(1000*1000);
I_motor = 0.0;
Ir = 0; N = 0;
Ir1 = 193444/(1000*1000);
Ir2 = 108732/(1000*1000);
Ir3 = 33896/(1000*1000);
l_OA = 0.336;
l_AB = 0.31;
l_BC = 0.205; % could reduce this
l_O_m1 = l_OA/2;
l_A_m2 = l_AB/2;
l_B_m3 = l_BC/2;
g = 9.81;
p_LIMS   = [m1 m2 m3 I1 I2 I3 Ir1 Ir2 Ir3 l_O_m1 l_A_m2 l_B_m3 l_OA l_AB l_BC g]'; % needs special dynamics
p_LIMS_nr = [m1 m2 m3 I1 I2 I3 0.0 0.0 0.0 l_O_m1 l_A_m2 l_B_m3 l_OA l_AB l_BC g]';

% Human parameters
m1 = 2.5;          
m2 = 1.45;
m3 = 0.53;
m_motor = 0.0;
I1 = 14845/(1000*1000);
I2 = 9355/(1000*1000);
I3 = 1370/(1000*1000);
I_motor = 0.0;
Ir = 0; N = 0;
l_OA = 0.366;
l_AB = 0.319; %(36.6/29.8)*0.417 - 0.193;
l_BC = 0.193;
l_O_m1 = l_OA/2;
l_A_m2 = l_AB/2;
l_B_m3 = l_BC/2;
% p_human  = [m1 m2 m3 m_motor I1 I2 I3 I_motor Ir N l_O_m1 l_A_m2 l_B_m3 l_OA l_AB l_BC g]';
p_human   = [m1 m2 m3 I1 I2 I3 0.0 0.0 0.0 l_O_m1 l_A_m2 l_B_m3 l_OA l_AB l_BC g]'; % needs special dynamics

% MIT Manipulator parameters
m1 = 6.68; 
m1_CB = 0; % when upper arm is q=1.15 rad from vertical, or 0.4208 from horizontal
m2 = 0.79;
m3 = 0.82;  
m_motor = 0.0;
I1 = 299652/(1000*1000);
I2 = 22500/(1000*1000);
I3 = 5774/(1000*1000);
I_motor = 0.0;
Ir = 0.0; N = 0;
Ir1 = (400/36)*0.022;
Ir2 = (400/36)*0.022;
Ir3 = (400/36)*0.0023;
l_OA = 0.45;
l_AB = 0.35;
l_BC = 0.16;
l_O_m1 = 0.17;
l_A_m2 = 0.11;
l_B_m3 = 0.07;
p_MIT   = [m1 m2 m3 I1 I2 I3 Ir1 Ir2 Ir3 l_O_m1 l_A_m2 l_B_m3 l_OA l_AB l_BC g]'; % needs special dynamics
p_MIT_nr = [m1 m2 m3 I1 I2 I3 0.0 0.0 0.0 l_O_m1 l_A_m2 l_B_m3 l_OA l_AB l_BC g]';
p_MIT_CB   = [m1_CB m2 m3 I1 I2 I3 Ir1 Ir2 Ir3 l_O_m1 l_A_m2 l_B_m3 l_OA l_AB l_BC g]'; % needs special dynamics



%% for a given configuration, calculate the effective mass of each arm

q = [-pi/2, pi/2, 0.0]'; % use this configuration for all arms

LLv_inv_UR5 = LLv_op_inv_arm_rotors(q,p_UR5);
LLv_inv_UR5_nr = LLv_op_inv_arm_rotors(q,p_UR5_nr);
LLv_inv_LIMS = LLv_op_inv_arm_rotors(q,p_LIMS);
LLv_inv_LIMS_nr = LLv_op_inv_arm_rotors(q,p_LIMS_nr);
LLv_inv_human = LLv_op_inv_arm_rotors(q,p_human);
LLv_inv_MIT = LLv_op_inv_arm_rotors(q,p_MIT);
LLv_inv_MIT_nr = LLv_op_inv_arm_rotors(q,p_MIT_nr);
LLv_inv_MIT_CB = LLv_op_inv_arm_rotors(q,p_MIT_CB);

thetas = linspace(0,2*pi,1000);

meff_UR5 = zeros(size(thetas));
meff_UR5_nr = zeros(size(thetas));
meff_LIMS = zeros(size(thetas));
meff_LIMS_nr = zeros(size(thetas));
meff_human = zeros(size(thetas));
meff_MIT = zeros(size(thetas));
meff_MIT_nr = zeros(size(thetas));
meff_MIT_CB = zeros(size(thetas));

for ii=1:length(thetas)
    
    u = [cos(thetas(ii)); sin(thetas(ii))];
    
    meff_UR5(ii) = 1/(u'*LLv_inv_UR5*u);
    meff_UR5_nr(ii) = 1/(u'*LLv_inv_UR5_nr*u);
    meff_LIMS(ii) = 1/(u'*LLv_inv_LIMS*u);
    meff_LIMS_nr(ii) = 1/(u'*LLv_inv_LIMS_nr*u);
    meff_human(ii) = 1/(u'*LLv_inv_human*u);
    meff_MIT(ii) = 1/(u'*LLv_inv_MIT*u);    
    meff_MIT_CB(ii) = 1/(u'*LLv_inv_MIT_CB*u);
    meff_MIT_nr(ii) = 1/(u'*LLv_inv_MIT_nr*u); 
    
end


%% use polarplot

figure; subplot(1,2,1);
polarplot(thetas, meff_UR5,'LineWidth',1.25); hold on;
polarplot(thetas, meff_MIT,'LineWidth',1.25);
polarplot(thetas, meff_LIMS,'LineWidth',1.25);
polarplot(thetas, meff_human,'LineWidth',1.25);
lgd = legend({'UR5', 'MIT', 'LIMS', 'Human'},'FontSize',12);
rticks([2.5,5,7.5,10]);
% title('Effective Mass (kg)');

ax = gca;
rlabel = ax.RAxis.Label;
rlabel.String = 'Effective Mass (kg)';
rlabel.Position = [101, 6.125, 0];

subplot(1,2,2);
polarplot(thetas, meff_MIT,'LineWidth',1.25, 'Color', [0.8500, 0.3250, 0.0980]); hold on;
polarplot(thetas, meff_LIMS,'LineWidth',1.25, 'Color', [0.9290, 0.6940, 0.1250]);
polarplot(thetas, meff_human,'LineWidth',1.25, 'Color', [0.4940, 0.1840, 0.5560]);
rlim([0,5]);
rticks([1,2,3,4,5]);
ax = gca;
rlabel = ax.RAxis.Label;
rlabel.String = 'Effective Mass (kg)';
rlabel.Position = [101, 0.6125*5, 0];

lgd.Position = [0.465, 0.15, 0.1, 0.1];

sgtitle(['Effective mass for different manipulators at q = [-90' char(176) ', 90' char(176) ', 0' char(176) ']']);

%% second plot

figure(7); clf;
polarplot(thetas, meff_MIT,'LineWidth',1.25); hold on;
polarplot(thetas, meff_MIT_CB,'LineWidth',1.25);
polarplot(thetas, meff_MIT_nr,'LineWidth',1.25);
rlim([0,6]);
rticks([1,2,3,4,5]);
ax = gca;
rlabel = ax.RAxis.Label;
rlabel.String = 'Effective Mass (kg)';
rlabel.Position = [100, 0.6125*5, 0];
lgd = legend({'With rotor inertia', 'With perfect counterbalance', 'No rotor inertia'}, 'Location', 'southoutside','FontSize',10);

title(['Effective mass for MIT Manipulator with N = 20 at q = [-90' char(176) ', 90' char(176) ', 0' char(176) ']']);







