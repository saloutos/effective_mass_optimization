% finding global minimum and maximum effective masses for arm model

% Andrew SaLoutos
% 4/23/2021

% make sure functions are on path
% addpath(genpath('arm_functions'))

% addpath(genpath('mc_arm_functions'))
addpath(genpath('UR3_arm_functions'))

% % define arm parameters
% m1 = 1;                 m2 = 1;
% m3 = 0.5;               m_motor = 0.5;
% I1 = 0.02;              I2 = 0.02;
% I3 = 0.01;              I_motor = 0.000625; % assuming radius r=0.05m
% Ir = 6.25e-6;           N = 6;
% l_O_m1 = 0.25;          l_A_m2 = 0.25;
% l_B_m3 = 0.25;          l_OA = 0.5;
% l_AB = 0.5;             l_BC = 0.5;
% g = 9.81; % do I want gravity to start?
% % arrange in vector
% p   = [m1 m2 m3 m_motor I1 I2 I3 I_motor Ir N l_O_m1 l_A_m2 l_B_m3 l_OA l_AB l_BC g]';

% load('mc_arm_functions\param_MC.mat');
% p = p_MC;

load('UR3_arm_functions\param_UR3.mat');
p = p_UR3;

% number of points
num_pts = 100;

% link angles
link1_angles = deg2rad( linspace(-90,90,num_pts) );% spaced to avoid singularities
link2_angles = deg2rad( linspace(-180,180,2*num_pts) );
link3_angles = deg2rad( linspace(-90,90,num_pts) );

n = length(link1_angles)*length(link2_angles)*length(link3_angles);

meff_min = zeros(3,n);
meff_max = zeros(3,n);
counter = 1;

for ii = 1:length(link1_angles)
    link1_temp = link1_angles(ii);
    for jj = 1:length(link2_angles)
        link2_temp = link2_angles(jj);
        for kk = 1:length(link3_angles)
            link3_temp = link3_angles(kk);
            
            % set angles vector
            qi = [link1_temp, link2_temp, link3_temp]';
            
            % get lambda
            LLv_inv = LLv_arm_op_inv(qi,p);
            [V,D] = eig(LLv_inv);
            
            if D(1,1)>D(2,2)
                meff_min(1,counter) = 1/D(1,1);
                meff_min(2:3,counter) = V(:,1);
                meff_max(1,counter) = 1/D(2,2);
                meff_max(2:3,counter) = V(:,2);
            else
                meff_min(1,counter) = 1/D(2,2);
                meff_min(2:3,counter) = V(:,2);
                meff_max(1,counter) = 1/D(1,1);
                meff_max(2:3,counter) = V(:,1);
            end
            
            
            counter = counter + 1;
        end
    end
    fprintf('Iteration %d / %d\n', counter, n);
end

% display true minimums and maximums
min(meff_min(1,:))
max(meff_max(1,:))

% could even run for different values of n in linspace, show meff-min
% converging and meff-max blowing up?



