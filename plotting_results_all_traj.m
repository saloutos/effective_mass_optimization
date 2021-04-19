% plotting interesting quantities across all trajectories

% Andrew SaLoutos
% 4/19/2021

% make sure functions are on path
addpath(genpath('arm_functions'))
addpath(genpath('block_functions'))

% import trajectory data
clear all

% load data
% load('multi_traj_data_linear.mat');
load('multi_traj_data_sinusoid.mat');

% each should contain TO_data_plain, TO_data_meff, TO_data_link3, and
% TO_data_meff_link3

% define arm parameters
m1 = 1;                 m2 = 1;
m3 = 0.5;               m_motor = 0.5;
I1 = 0.02;              I2 = 0.02;
I3 = 0.01;              I_motor = 0.000625; % assuming radius r=0.05m
Ir = 6.25e-6;           N = 6;
l_O_m1 = 0.25;          l_A_m2 = 0.25;
l_B_m3 = 0.25;          l_OA = 0.5;
l_AB = 0.5;             l_BC = 0.5;
g = 9.81; % do I want gravity to start?
% arrange in vector
p   = [m1 m2 m3 m_motor I1 I2 I3 I_motor Ir N l_O_m1 l_A_m2 l_B_m3 l_OA l_AB l_BC g]';


% iterate through each trajectory and compile quantities
N = size(TO_data_plain,2);

% quantities of interest?
solve_times = zeros(4,N);
avg_p_error = zeros(4,N);
avg_v_error = zeros(4,N);
avg_meff = zeros(4,N);

for ii=1:4 % for each cost function
    
    if (ii==1)
        TO_data = TO_data_plain;
    elseif (ii==2)
        TO_data = TO_data_meff;
    elseif (ii==3)
        TO_data = TO_data_link3;
    elseif (ii==4)
        TO_data = TO_data_meff_link3;
    end    
    
    for jj=1:N % for each trajectory
        
       solve_times(ii,jj) = TO_data(jj).solve_time;
       
       p_error = TO_data(jj).data.p(3,:);
       avg_p_error(ii,jj) = mean(p_error);
        
       v_error = TO_data(jj).data.v(3,:);
       avg_v_error(ii,jj) = mean(v_error);
       
       meff = TO_data(jj).data.meff(1,2:end); % first meff calc is garbage
       avg_meff(ii,jj) = mean(meff);
       
       % TODO: calculate effective momentum as well?
       
    end    
end

% calculate ratios for effective mass? each cost function over the plain TO
% can use avg(meff1)/avg(meff2) or avg( meff1[i]/meff2[i] )
% these two won't be equal

meff_ratios = zeros(3,N);

for ii=1:N
    % get meff of each cost function
    meff_base = TO_data_plain(ii).data.meff(1,2:end);
    meff_m = TO_data_meff(ii).data.meff(1,2:end);
    meff_l3 = TO_data_link3(ii).data.meff(1,2:end);
    meff_m_l3 = TO_data_meff_link3(ii).data.meff(1,2:end);
    % ratio at each 
    ratio_m = meff_m./meff_base;
    ratio_l3 = meff_l3./meff_base;
    ratio_m_l3 = meff_m_l3./meff_base;
    % save data
    meff_ratios(:,ii) = [mean(ratio_m); mean(ratio_l3); mean(ratio_m_l3)];
end


% plot quantities
figure(3); clf; hold on;
plot(solve_times(1,:),'o-'); plot(solve_times(2,:),'o-');
plot(solve_times(3,:),'o-'); plot(solve_times(4,:),'o-');
plot(ones(1,N),'r--');
xlabel('Traj #'); ylabel('Solve Time');
legend('TO plain', 'TO w/ m_{eff}', 'TO w/ link3', 'TO w/ both');

figure(4); clf;
subplot(2,2,1); hold on;
plot(meff_ratios(1,:));
plot(meff_ratios(2,:));
plot(meff_ratios(3,:));
plot(ones(1,N),'r--');
xlabel('Traj #'); ylabel('m_{eff} Ratio');
legend('TO w/ m_{eff}', 'TO w/ link3', 'TO w/ both');
title('Ratio of m_{eff} to Plain TO');

subplot(2,2,2); hold on;
plot(avg_meff(1,:)); plot(avg_meff(2,:));
plot(avg_meff(3,:)); plot(avg_meff(4,:));
xlabel('Traj #'); ylabel('Mean m_{eff}'); 
legend('TO plain', 'TO w/ m_{eff}', 'TO w/ link3', 'TO w/ both');
title('Average m_{eff}');

subplot(2,2,3); hold on;
plot(avg_p_error(1,:)); plot(avg_p_error(2,:));
plot(avg_p_error(3,:)); plot(avg_p_error(4,:));
xlabel('Traj #'); ylabel('Mean p_{err}'); 
legend('TO plain', 'TO w/ m_{eff}', 'TO w/ link3', 'TO w/ both');
title('Average Position Tracking Error');

subplot(2,2,4); hold on;
plot(avg_v_error(1,:)); plot(avg_v_error(2,:));
plot(avg_v_error(3,:)); plot(avg_v_error(4,:));
xlabel('Traj #'); ylabel('Mean v_{err}'); 
legend('TO plain', 'TO w/ m_{eff}', 'TO w/ link3', 'TO w/ both');
title('Average Velocity Tracking Error');