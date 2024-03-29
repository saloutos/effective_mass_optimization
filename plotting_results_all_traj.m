% plotting interesting quantities across all trajectories

% Andrew SaLoutos
% 4/19/2021

% make sure functions are on path
% addpath(genpath('arm_functions'))
addpath(genpath('block_functions'))

% addpath(genpath('mc_arm_functions'))
addpath(genpath('UR3_arm_functions'))

% import trajectory data
clear all

% load data
% load('multi_traj_data_linear.mat');
% load('multi_traj_data_linear_MC_2.mat');
% load('multi_traj_data_linear_UR3_2.mat');

% load('multi_traj_data_sinusoid.mat');
% load('multi_traj_data_sinusoid_MC_2.mat');
load('multi_traj_data_sinusoid_UR3_2.mat');

% each should contain TO_data_plain, TO_data_meff, TO_data_link3, and
% TO_data_meff_link3

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

% MC arm
% m1 = 0.195;             m2 = 0.262;
% m3 = 0.053;             m_motor = 0.527;
% I1 = 0.001170;          I2 = 0.001186;
% I3 = 0.000096;          I_motor = 0.000508;
% Ir = 0.000064;          N = 6;
% l_O_m1 = 0.092;         l_A_m2 = 0.201;
% l_B_m3 = 0.038;         l_OA = 0.2085;
% l_AB = 0.265;           l_BC = 0.1225;
% g = 9.81; % do I want gravity to start?
% % parameters
% p   = [m1 m2 m3 m_motor I1 I2 I3 I_motor Ir N l_O_m1 l_A_m2 l_B_m3 l_OA l_AB l_BC g]'; 

% % UR3 planar arm
m1 = 3.4445;            m2 = 1.437;
m3 = 1.9360;            m_motor = 0.0; % not using motor mass
I1 = 0.0219;            I2 = 0.0075;
I3 = 0.0064;            I_motor = 0.0; % not using motor inertia
Ir = 2.07e-5;           N = 101;
l_O_m1 = 0.113;         l_A_m2 = 0.163;
l_B_m3 = 0.0470;        l_OA = 0.24355;
l_AB = 0.2132;          l_BC = 0.08535;
g = 9.81; % do I want gravity to start?
% parameters
p   = [m1 m2 m3 m_motor I1 I2 I3 I_motor Ir N l_O_m1 l_A_m2 l_B_m3 l_OA l_AB l_BC g]'; 

% define object parameters
mo = 20.0; %1.0;
LLo = mo*eye(2);

% iterate through each trajectory and compile quantities
N = size(TO_data_plain,2);

% quantities of interest?
solve_times = zeros(4,N);
avg_p_error = zeros(4,N);
avg_v_error = zeros(4,N);
avg_meff = zeros(4,N);
avg_exit_vel = zeros(4,N);

exit_vel_data = [];
exit_vel_data_m = [];
exit_vel_data_l3 = [];
exit_vel_data_m_l3 = [];

bad_results = zeros(4,N);

summary_data = zeros(4,6);
summary_data_err = zeros(4,6);

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
        
        if ischar(TO_data(jj).data)
            disp(ii)
            disp(jj)
            bad_results(ii,jj) = 1;
            continue
        end
        
        solve_times(ii,jj) = TO_data(jj).solve_time;

        p_error = TO_data(jj).data.p(3,:);
        avg_p_error(ii,jj) = mean(p_error);

        v_error = TO_data(jj).data.v(3,:);
        avg_v_error(ii,jj) = mean(v_error);

        meff = TO_data(jj).data.meff(1,2:end); % first meff calc is garbage
        avg_meff(ii,jj) = mean(meff);
        
        % exit velocities
        exit_vels = zeros(3,length(v_error));
        for kk=1:length(v_error)
            v_temp = TO_data(jj).data.v(1:2,kk);       
            LLv_temp = reshape(TO_data(jj).data.LLv(:,kk),2,2);
            
            LLv_inv_temp = inv(LLv_temp);
            if (norm(v_temp)>1e-6)
                u_temp = v_temp/norm(v_temp);
                mf_temp = 1/(u_temp'*LLv_inv_temp*u_temp);
            else 
                mf_temp = 0;
            end
            ex_vel_temp = 2*(mf_temp/(mf_temp+mo))*v_temp;  
            
%             ex_vel_temp = 2*inv(LLv_temp+LLo)*LLv_temp*v_temp;
            exit_vels(:,kk) = [ex_vel_temp; norm(ex_vel_temp)];
            
        end
        avg_exit_vel(ii,jj) = mean(exit_vels(3,:));
        
        
        if (ii==1)
            exit_vel_data(:,:,jj) = exit_vels;
        elseif (ii==2)
            exit_vel_data_m(:,:,jj) = exit_vels;
        elseif (ii==3)
            exit_vel_data_l3(:,:,jj) = exit_vels;
        elseif (ii==4)
            exit_vel_data_m_l3(:,:,jj) = exit_vels;
        end
        
        % TODO: calculate effective momentum as well?
       
    end    
end

avg_exit_vel

% calculate ratios for effective mass? each cost function over the plain TO
% can use avg(meff1)/avg(meff2) or avg( meff1[i]/meff2[i] )
% these two won't be equal

meff_ratios = zeros(3,N);
exit_vel_ratios = zeros(3,N);

for ii=1:N
    % get meff of each cost function
    if bad_results(1,ii)==0
        meff_base = TO_data_plain(ii).data.meff(1,2:end);
        exit_vel_base = exit_vel_data(3,2:end,ii); 
    end
    if bad_results(2,ii)==0
        meff_m = TO_data_meff(ii).data.meff(1,2:end);
        exit_vel_m = exit_vel_data_m(3,2:end,ii);
    end
    if bad_results(3,ii)==0
        meff_l3 = TO_data_link3(ii).data.meff(1,2:end);
        exit_vel_l3 = exit_vel_data_l3(3,2:end,ii);
    end
    if bad_results(4,ii)==0
        meff_m_l3 = TO_data_meff_link3(ii).data.meff(1,2:end);
        exit_vel_m_l3 = exit_vel_data_m_l3(3,2:end,ii);
    end
    % ratio at each 
    if bad_results(1,ii)==0 && bad_results(2,ii)==0
        ratio_m = meff_m./meff_base;
        ev_ratio_m = exit_vel_m./exit_vel_base;
    else
        ratio_m = nan;
        ev_ratio_m = nan;
    end
    if bad_results(1,ii)==0 && bad_results(3,ii)==0
        ratio_l3 = meff_l3./meff_base;
        ev_ratio_l3 = exit_vel_l3./exit_vel_base;
    else
        ratio_l3 = nan;
        ev_ratio_l3 = nan;
    end
    if bad_results(1,ii)==0 && bad_results(4,ii)==0
        ratio_m_l3 = meff_m_l3./meff_base;
        ev_ratio_m_l3 = exit_vel_m_l3./exit_vel_base;
    else
        ratio_m_l3 = nan;
        ev_ratio_m_l3 = nan;
    end
    % save data
    meff_ratios(:,ii) = [mean(ratio_m); mean(ratio_l3); mean(ratio_m_l3)];
    exit_vel_ratios(:,ii) = [mean(ev_ratio_m); mean(ev_ratio_l3); mean(ev_ratio_m_l3)];
%     meffs = [meff_base; meff_m; meff_l3; meff_m_l3];
    
end

% set bad results to NaN
solve_times(solve_times==0) = NaN;
avg_meff(avg_meff==0) = NaN;
avg_p_error(avg_p_error==0) = NaN;
avg_v_error(avg_v_error==0) = NaN;
avg_exit_vel(avg_exit_vel==0) = NaN;

% get summary data
summary_data(:,1) = nanmean(solve_times,2);
summary_data(:,2) = min(solve_times,[],2);
summary_data(:,3) = max(solve_times,[],2);
summary_data(:,4) = nanmean(avg_meff,2);
summary_data(:,5) = min(avg_meff,[],2);
summary_data(:,6) = max(avg_meff,[],2);
summary_data(:,7) = nanmean(avg_exit_vel,2);
summary_data(:,8) = min(avg_exit_vel,[],2);
summary_data(:,9) = max(avg_exit_vel,[],2);

summary_data_err(:,1) = nanmean(avg_p_error,2);
summary_data_err(:,2) = min(avg_p_error,[],2);
summary_data_err(:,3) = max(avg_p_error,[],2);
summary_data_err(:,4) = nanmean(avg_v_error,2);
summary_data_err(:,5) = min(avg_v_error,[],2);
summary_data_err(:,6) = max(avg_v_error,[],2);

% display data for comparison
% summary_data
% summary_data_err

%% plot quantities
figure(3); clf; hold on;
plot(solve_times(1,:),'o-'); plot(solve_times(2,:),'o-');
plot(solve_times(3,:),'o-'); plot(solve_times(4,:),'o-');
plot(ones(1,N),'r--');
xlabel('Traj #'); ylabel('Solve Time (s)');
legend({'TO plain', 'TO w/ m_{eff}', 'TO w/ link3', 'TO w/ both'},'Location','northwest');
title('Optimization Solve Times for Random Trajectories');
% ylim([0,5]);

figure(4); clf;
subplot(2,3,1); hold on;
plot(meff_ratios(1,:));
plot(meff_ratios(2,:));
plot(meff_ratios(3,:));
plot(ones(1,N),'r--');
xlabel('Traj #'); ylabel('m_{eff} Ratio');
legend('TO w/ m_{eff}', 'TO w/ link3', 'TO w/ both');
title('Ratio of m_{eff} to Plain TO');
ylim([0,3]);

subplot(2,3,2); hold on;
plot(avg_meff(1,:)); plot(avg_meff(2,:));
plot(avg_meff(3,:)); plot(avg_meff(4,:));
xlabel('Traj #'); ylabel('Mean m_{eff} (kg)'); 
legend({'TO plain', 'TO w/ m_{eff}', 'TO w/ link3', 'TO w/ both'},'Location','northwest');
title('Average m_{eff}');

subplot(2,3,4); hold on;
plot(exit_vel_ratios(1,:));
plot(exit_vel_ratios(2,:));
plot(exit_vel_ratios(3,:));
plot(ones(1,N),'r--');
xlabel('Traj #'); ylabel('||v_{exit}|| Ratio');
legend('TO w/ m_{eff}', 'TO w/ link3', 'TO w/ both');
title('Ratio of ||v_{exit}|| to Plain TO');
ylim([0,3]);

subplot(2,3,5); hold on;
plot(avg_exit_vel(1,:)); plot(avg_exit_vel(2,:));
plot(avg_exit_vel(3,:)); plot(avg_exit_vel(4,:));
xlabel('Traj #'); ylabel('Mean ||v_{exit}|| (m/s)'); 
% legend('TO plain', 'TO w/ m_{eff}', 'TO w/ link3', 'TO w/ both');
title('Average ||v_{exit}||');

subplot(2,3,3); hold on;
plot(avg_p_error(1,:)); plot(avg_p_error(2,:));
plot(avg_p_error(3,:)); plot(avg_p_error(4,:));
xlabel('Traj #'); ylabel('Mean p_{err} (m)'); 
% legend('TO plain', 'TO w/ m_{eff}', 'TO w/ link3', 'TO w/ both');
title('Average Position Tracking Error');

subplot(2,3,6); hold on;
plot(avg_v_error(1,:)); plot(avg_v_error(2,:));
plot(avg_v_error(3,:)); plot(avg_v_error(4,:));
xlabel('Traj #'); ylabel('Mean v_{err} (m/s)'); 
% legend('TO plain', 'TO w/ m_{eff}', 'TO w/ link3', 'TO w/ both');
title('Average Velocity Tracking Error');

sgtitle('Comparing TO Performance Across Random Trajectories');

%% pseudocolor plots across trials and cost functions
figure(5);
subplot(3,1,1);
sz = size(meff_ratios);
C = [meff_ratios, zeros(sz(1),1); zeros(1,sz(2)+1)];
pcolor(C);
colorbar;
axis ij;
yticks([1.5,2.5,3.5]);
yticklabels({'meff','link3','both'});
xticks(5.5:5:30.5);
xticklabels({'5','10','15','20'}); %,'25','30'});
title('Meff ratios');

subplot(3,1,2);
sz = size(exit_vel_ratios);
C = [exit_vel_ratios, zeros(sz(1),1); zeros(1,sz(2)+1)];
pcolor(C);
colorbar;
axis ij;
yticks([1.5,2.5,3.5]);
yticklabels({'meff','link3','both'});
xticks(5.5:5:30.5);
xticklabels({'5','10','15','20'}); %,'25','30'});
title('Exit vel ratios');

% subplot(5,1,2);
% sz = size(avg_meff);
% C = [avg_meff, zeros(sz(1),1); zeros(1,sz(2)+1)];
% pcolor(C);
% colorbar;
% axis ij;
% yticks([1.5,2.5,3.5,4.5]);
% yticklabels({'plain','meff','link3','both'});
% xticks(5.5:5:30.5);
% xticklabels({'5','10','15','20','25','30'});
% title('Avg meff');
%
% subplot(5,1,3);
% sz = size(avg_p_error);
% C = [avg_p_error, zeros(sz(1),1); zeros(1,sz(2)+1)];
% pcolor(C);
% colorbar;
% axis ij;
% yticks([1.5,2.5,3.5,4.5]);
% yticklabels({'plain','meff','link3','both'});
% xticks(5.5:5:30.5);
% xticklabels({'5','10','15','20','25','30'});
% title('Avg p error');
% 
% subplot(5,1,4);
% sz = size(avg_v_error);
% C = [avg_v_error, zeros(sz(1),1); zeros(1,sz(2)+1)];
% pcolor(C);
% colorbar;
% axis ij;
% yticks([1.5,2.5,3.5,4.5]);
% yticklabels({'plain','meff','link3','both'});
% xticks(5.5:5:30.5);
% xticklabels({'5','10','15','20','25','30'});
% title('Avg v error');

%figure(6);
subplot(3,1,3);
sz = size(solve_times);
C = [solve_times, zeros(sz(1),1); zeros(1,sz(2)+1)];
pcolor(C);
colorbar;
axis ij;
yticks([1.5,2.5,3.5,4.5]);
yticklabels({'plain','meff','link3','both'});
xticks(5.5:5:30.5);
xticklabels({'5','10','15','20'}); %,'25','30'});
title('Solve Times');

% outliers make this difficult to accurately represent results









