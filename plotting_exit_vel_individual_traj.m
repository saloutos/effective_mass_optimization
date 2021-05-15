% plotting exit vel metrics for each TO output

% Andrew SaLoutos
% 5/5/2021

% make sure functions are on path
addpath(genpath('arm_functions'))
addpath(genpath('block_functions'))
addpath(genpath('arm_floating_functions'));

% import trajectory data
clear all

% load data
% load('multi_traj_data_linear_ikp.mat');
load('multi_traj_data_sinusoid_ikp.mat');

% make sure video filenames agree
filename1 = 'Figures/figure_sinusoid_traj_ikp_exit_vel_metrics.avi';



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

% define object inertia
mo = 1;
LLo = mo*eye(2);

% define base inertia
mb = 2; Ib = 0.5*mb*0.1^2;
pb = [p; mb; Ib];

% iterate through each trajectory and plot
N = size(TO_data_plain,2);
m = length(TO_data_plain(1).time);
exit_vel_data = zeros(16,m,N);

% struct to save frames for video
F1(N) = struct('cdata',[],'colormap',[]);

for ii=1:N
    if ischar(TO_data_plain(ii).data) || ischar(TO_data_meff(ii).data) || ischar(TO_data_link3(ii).data) || ischar(TO_data_meff_link3(ii).data)
        % bad TO outputs, just don't plot?
        figure(1); clf;
        plt_title = sprintf('Exit velocity metrics for trajectory #%d', ii);
        sgtitle(plt_title);
        % save figure as frame
        F1(ii) = getframe(gcf);
        continue
    end

    time_vec = TO_data_plain(ii).time;
    des_pts = TO_data_plain(ii).pts;
    des_vels = TO_data_plain(ii).vels;
    
    X = TO_data_plain(ii).data;
    X_m = TO_data_meff(ii).data;
    X_l3 = TO_data_link3(ii).data;
    X_m_l3 = TO_data_meff_link3(ii).data;

    act_pts = X.p;
    act_pts_m = X_m.p;
    act_pts_l3 = X_l3.p;
    act_pts_m_l3 = X_m_l3.p;
    
    act_vels = X.v(1:2,:);
    act_vels_m = X_m.v(1:2,:);
    act_vels_l3 = X_l3.v(1:2,:);
    act_vels_m_l3 = X_m_l3.v(1:2,:);

    q = X.q;
    q_m = X_m.q;
    q_l3 = X_l3.q;
    q_m_l3 = X_m_l3.q;  
    
    for jj=1:m % iterate through points in trajectory
        
        q_temp = [zeros(3,1); q(:,jj)];
        q_m_temp = [zeros(3,1); q_m(:,jj)];
        q_l3_temp = [zeros(3,1); q_l3(:,jj)];
        q_m_l3_temp = [zeros(3,1); q_m_l3(:,jj)];
        
        ev_temp = act_vels(:,jj)/norm(act_vels(:,jj));
        ev_m_temp = act_vels_m(:,jj)/norm(act_vels_m(:,jj));
        ev_l3_temp = act_vels_l3(:,jj)/norm(act_vels_l3(:,jj));
        ev_m_l3_temp = act_vels_m_l3(:,jj)/norm(act_vels_m_l3(:,jj));
        
        [phi_vf, Pvf, phi_vo, Pvo, phi_vol, Pvol] = eval_exit_vel_metrics(q_temp,pb,LLo);
        [phi_vf_m, Pvf_m, phi_vo_m, Pvo_m, phi_vol_m, Pvol_m] = eval_exit_vel_metrics(q_m_temp,pb,LLo);
        [phi_vf_l3, Pvf_l3, phi_vo_l3, Pvo_l3, phi_vol_l3, Pvol_l3] = eval_exit_vel_metrics(q_l3_temp,pb,LLo);
        [phi_vf_m_l3, Pvf_m_l3, phi_vo_m_l3, Pvo_m_l3, phi_vol_m_l3, Pvol_m_l3] = eval_exit_vel_metrics(q_m_l3_temp,pb,LLo);
        
        % get directional value, only for Pvol
        eta_vol = ev_temp'*Pvol*ev_temp;
        eta_vol_m = ev_m_temp'*Pvol_m*ev_m_temp;
        eta_vol_l3 = ev_l3_temp'*Pvol_l3*ev_l3_temp;
        eta_vol_m_l3 = ev_m_l3_temp'*Pvol_m_l3*ev_m_l3_temp;
        
        % get maximum and minimum eta_vol
        [V,D] = eig(Pvol);
        ax_lengths = [sqrt(D(1,1)),sqrt(D(2,2))];
        eta_max = max(ax_lengths);
        eta_min = min(ax_lengths);
        [V,D] = eig(Pvol_m);
        ax_lengths = [sqrt(D(1,1)),sqrt(D(2,2))];
        eta_max_m = max(ax_lengths);
        eta_min_m = min(ax_lengths);
        [V,D] = eig(Pvol_l3);
        ax_lengths = [sqrt(D(1,1)),sqrt(D(2,2))];
        eta_max_l3 = max(ax_lengths);
        eta_min_l3 = min(ax_lengths);
        [V,D] = eig(Pvol_m_l3);
        ax_lengths = [sqrt(D(1,1)),sqrt(D(2,2))];
        eta_max_m_l3 = max(ax_lengths);
        eta_min_m_l3 = min(ax_lengths);
        
        exit_vel_data(:,jj,ii) = [phi_vol; phi_vol_m; phi_vol_l3; phi_vol_m_l3; ...
                                  eta_vol; eta_vol_m; eta_vol_l3; eta_vol_m_l3; ...
                                  eta_max; eta_max_m; eta_max_l3; eta_max_m_l3; ...
                                  eta_min; eta_min_m; eta_min_l3; eta_min_m_l3];
                              
                      
    end 
    
    % plot
    data = exit_vel_data(:,:,ii);
    time_vec = TO_data_plain(ii).time;
    
    figure(1); clf;     
    
    subplot(1,4,1); hold on;
    plot(time_vec, data(1,:)); plot(time_vec, data(2,:));
    plot(time_vec, data(3,:)); plot(time_vec, data(4,:));
    xlabel('Time'); ylabel('Metric'); ylim([0, 1.0]);
    legend('plain','meff','link3','both');
    title('\phi_{vol}');
    
    subplot(1,4,2); hold on;
    plot(time_vec, data(5,:)); plot(time_vec, data(6,:));
    plot(time_vec, data(7,:)); plot(time_vec, data(8,:));
    xlabel('Time'); ylabel('Metric'); ylim([0, 1.0]);
%     legend('plain','meff','link3','both');
    title('\eta_{vol}');
    
    subplot(1,4,3); hold on;
    plot(time_vec, data(9,:)); plot(time_vec, data(10,:));
    plot(time_vec, data(11,:)); plot(time_vec, data(12,:));
    xlabel('Time'); ylabel('Metric'); ylim([0, 1.0]);
%     legend('plain','meff','link3','both');
    title('\eta_{max}');
    
    subplot(1,4,4); hold on;
    plot(time_vec, data(13,:)); plot(time_vec, data(14,:));
    plot(time_vec, data(15,:)); plot(time_vec, data(16,:));
    xlabel('Time'); ylabel('Metric'); ylim([0, 1.0]);
%     legend('plain','meff','link3','both');
    title('\eta_{min}');
    
    plt_title = sprintf('Exit velocity metrics for trajectory #%d', ii);
    sgtitle(plt_title);
    
    % save figure as frame
    F1(ii) = getframe(gcf);
    
%     pause;
    
    
end

% %% plot?
% figure(1);
% for ii=1:N
%     
%     data = exit_vel_data(:,:,ii);
%     time_vec = TO_data_plain(ii).time;
%     
%     clf;     
%     
%     subplot(1,4,1); hold on;
%     plot(time_vec, data(1,:)); plot(time_vec, data(2,:));
%     plot(time_vec, data(3,:)); plot(time_vec, data(4,:));
%     xlabel('Time'); ylabel('Metric'); ylim([0, 1.0]);
%     legend('plain','meff','link3','both');
%     title('\phi_{vol}');
%     
%     subplot(1,4,2); hold on;
%     plot(time_vec, data(5,:)); plot(time_vec, data(6,:));
%     plot(time_vec, data(7,:)); plot(time_vec, data(8,:));
%     xlabel('Time'); ylabel('Metric'); ylim([0, 1.0]);
% %     legend('plain','meff','link3','both');
%     title('\eta_{vol}');
%     
%     subplot(1,4,3); hold on;
%     plot(time_vec, data(9,:)); plot(time_vec, data(10,:));
%     plot(time_vec, data(11,:)); plot(time_vec, data(12,:));
%     xlabel('Time'); ylabel('Metric'); ylim([0, 1.0]);
% %     legend('plain','meff','link3','both');
%     title('\eta_{max}');
%     
%     subplot(1,4,4); hold on;
%     plot(time_vec, data(13,:)); plot(time_vec, data(14,:));
%     plot(time_vec, data(15,:)); plot(time_vec, data(16,:));
%     xlabel('Time'); ylabel('Metric'); ylim([0, 1.0]);
% %     legend('plain','meff','link3','both');
%     title('\eta_{min}');
%     
%     plt_title = sprintf('Exit velocity metrics for trajectory #%d', ii);
%     sgtitle(plt_title);
%     
%     pause;
% end

% write videos here
disp('Writing videos...');

v1 = VideoWriter(filename1);
v1.FrameRate = 30; % optional: control frame rate
open(v1)
for ii=1:N % for each trajectory
    for jj=1:60 % write 60 frames (two seconds) to video
        writeVideo(v1,F1(ii));
    end
end
close(v1)

disp('...done!');












