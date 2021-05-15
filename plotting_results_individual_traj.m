% plotting each TO output

% Andrew SaLoutos
% 4/17/2021

% make sure functions are on path
% addpath(genpath('arm_functions'))
addpath(genpath('block_functions'))

addpath(genpath('mc_arm_functions'))
% addpath(genpath('UR3_arm_functions'))

% import trajectory data
clear all

% load data
% load('multi_traj_data_linear.mat');
load('multi_traj_data_linear_MC.mat');
% load('multi_traj_data_linear_UR3.mat');

% load('multi_traj_data_sinusoid.mat');
% load('multi_traj_data_sinusoid_MC.mat');
% load('multi_traj_data_sinusoid_UR3.mat');

% make sure video filenames agree
filename1 = 'Figures/figure1_linear_traj_MC.avi';
filename2 = 'Figures/figure2_linear_traj_MC.avi';

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
m1 = 0.195;             m2 = 0.262;
m3 = 0.053;             m_motor = 0.527;
I1 = 0.001170;          I2 = 0.001186;
I3 = 0.000096;          I_motor = 0.000508;
Ir = 0.000064;          N = 6;
l_O_m1 = 0.092;         l_A_m2 = 0.201;
l_B_m3 = 0.038;         l_OA = 0.2085;
l_AB = 0.265;           l_BC = 0.1225;
g = 9.81; % do I want gravity to start?
% parameters
p   = [m1 m2 m3 m_motor I1 I2 I3 I_motor Ir N l_O_m1 l_A_m2 l_B_m3 l_OA l_AB l_BC g]'; 

% % UR3 planar arm
% m1 = 3.4445;            m2 = 1.437;
% m3 = 1.9360;            m_motor = 0.0; % not using motor mass
% I1 = 0.0219;            I2 = 0.0075;
% I3 = 0.0064;            I_motor = 0.0; % not using motor inertia
% Ir = 2.07e-5;           N = 101;
% l_O_m1 = 0.113;         l_A_m2 = 0.163;
% l_B_m3 = 0.0470;        l_OA = 0.24355;
% l_AB = 0.2132;          l_BC = 0.08535;
% g = 9.81; % do I want gravity to start?
% % parameters
% p   = [m1 m2 m3 m_motor I1 I2 I3 I_motor Ir N l_O_m1 l_A_m2 l_B_m3 l_OA l_AB l_BC g]'; 

% define object parameters
mo = 1;
LLo = mo*eye(2);



% iterate through each trajectory and plot
N = size(TO_data_plain,2);

% struct to save frames for video
F1(N) = struct('cdata',[],'colormap',[]);
F2(N) = struct('cdata',[],'colormap',[]);

for ii=1:N
    if ischar(TO_data_plain(ii).data) || ischar(TO_data_meff(ii).data) || ischar(TO_data_link3(ii).data) || ischar(TO_data_meff_link3(ii).data)
        % bad TO outputs, just don't plot?
        figure(1); clf;
        plt_title = sprintf('Optimized trajectory and effective mass for reference trajectory #%d', ii);
        sgtitle(plt_title);
        % save figure as frame
        F1(ii) = getframe(gcf);
        figure(2); clf;
        plt_title = sprintf('Tracking errors for reference trajectory #%d', ii);
        sgtitle(plt_title);
        % save figure as frame
        F2(ii) = getframe(gcf);
        continue
    end

    time_vec = TO_data_plain(ii).time;
    des_pts = TO_data_plain(ii).pts;
    
    X = TO_data_plain(ii).data;
    X_m = TO_data_meff(ii).data;
    X_l3 = TO_data_link3(ii).data;
    X_m_l3 = TO_data_meff_link3(ii).data;

    act_pts = X.p;
    act_pts_m = X_m.p;
    act_pts_l3 = X_l3.p;
    act_pts_m_l3 = X_m_l3.p;

    meff = X.meff;
    meff_m = X_m.meff;
    meff_l3 = X_l3.meff;
    meff_m_l3 = X_m_l3.meff;

    verr = X.v(3,:);
    verr_m = X_m.v(3,:);
    verr_l3 = X_l3.v(3,:);
    verr_m_l3 = X_m_l3.v(3,:);

    q0 = X.q(:,1);
    qf = X.q(:,end);
    q0_m = X_m.q(:,1);
    qf_m = X_m.q(:,end);
    q0_l3 = X_l3.q(:,1);
    qf_l3 = X_l3.q(:,end);
    q0_m_l3 = X_m_l3.q(:,1);
    qf_m_l3 = X_m_l3.q(:,end);

    LLv = X.LLv;
    LLv_m = X_m.LLv;
    LLv_l3 = X_l3.LLv;
    LLv_m_l3 = X_m_l3.LLv;    
    
    % get exit velocities [vx; vy; ||v||];
    exit_vel = zeros(3,length(time_vec));
    exit_vel_m = zeros(3,length(time_vec));
    exit_vel_l3 = zeros(3,length(time_vec));
    exit_vel_m_l3 = zeros(3,length(time_vec));
    
    for jj=1:length(time_vec)
        
        v_temp = X.v(1:2,jj);
        LLv_temp = reshape(LLv(:,jj),2,2);
        ex_vel_temp = 2*inv(LLv_temp+LLo)*LLv_temp*v_temp;
        exit_vel(:,jj) = [ex_vel_temp; norm(ex_vel_temp)];
        
        v_temp = X_m.v(1:2,jj);
        LLv_temp = reshape(LLv_m(:,jj),2,2);
        ex_vel_temp = 2*inv(LLv_temp+LLo)*LLv_temp*v_temp;
        exit_vel_m(:,jj) = [ex_vel_temp; norm(ex_vel_temp)];
        
        v_temp = X_l3.v(1:2,jj);
        LLv_temp = reshape(LLv_l3(:,jj),2,2);
        ex_vel_temp = 2*inv(LLv_temp+LLo)*LLv_temp*v_temp;
        exit_vel_l3(:,jj) = [ex_vel_temp; norm(ex_vel_temp)];
        
        v_temp = X_m_l3.v(1:2,jj);
        LLv_temp = reshape(LLv_m_l3(:,jj),2,2);
        ex_vel_temp = 2*inv(LLv_temp+LLo)*LLv_temp*v_temp;
        exit_vel_m_l3(:,jj) = [ex_vel_temp; norm(ex_vel_temp)];
        
    end
        
    % plot desired and actual trajectories
    figure(1); clf; 
    subplot(2,2,[1 3]); hold on;

    % trajectories
    plot(des_pts(1,:),des_pts(2,:),'g','LineWidth',2.0);
    plot(act_pts(1,:),act_pts(2,:),'Color',[0.1, 0.55, 0.1],'LineWidth',1.25);
    plot(act_pts_m(1,:),act_pts_m(2,:),'Color',[1.0, 0.27, 0],'LineWidth',1.25);
    plot(act_pts_l3(1,:),act_pts_l3(2,:),'Color',[0.68, 0.85, 0.9],'LineWidth',1.25);
    plot(act_pts_m_l3(1,:),act_pts_m_l3(2,:),'b');
    % arena
%     x_a = [0.4, 1.2, 1.2, 0.4, 0.4];
%     y_a = [-0.4, -0.4, 0.4, 0.4, -0.4];
    
    x_a = [0.2, 0.5, 0.5, 0.2, 0.2]; % for MC, UR3 plots
    y_a = [-0.15, -0.15, 0.15, 0.15, -0.15];
    
    plot(x_a,y_a,'r--');
    % initial and final joint configurations
    % lighter plots are without meff in optimization
    kp0_m_l3 = keypoints_arm(q0_m_l3,p);
    rA0_m_l3 = kp0_m_l3(:,1);
    rB0_m_l3 = kp0_m_l3(:,2);
    rC0_m_l3 = kp0_m_l3(:,3);
    plot([0 rA0_m_l3(1)],[0 rA0_m_l3(2)],'Color',[0.3,0.3,0.3],'LineWidth',2);
    plot([rA0_m_l3(1) rB0_m_l3(1)],[rA0_m_l3(2) rB0_m_l3(2)],'Color',[0.3,0.3,0.3],'LineWidth',2);
    plot([rB0_m_l3(1) rC0_m_l3(1)],[rB0_m_l3(2) rC0_m_l3(2)],'Color',[0.3,0.3,0.3],'LineWidth',2);
    plot([0 rA0_m_l3(1) rB0_m_l3(1)],[0 rA0_m_l3(2) rB0_m_l3(2)],'o','Color',[0.3,0.3,0.3],'MarkerSize',5,'MarkerFaceColor',[0.3,0.3,0.3]);
    kpf_m_l3 = keypoints_arm(qf_m_l3,p);
    rAf_m_l3 = kpf_m_l3(:,1);
    rBf_m_l3 = kpf_m_l3(:,2);
    rCf_m_l3 = kpf_m_l3(:,3);
    plot([0 rAf_m_l3(1)],[0 rAf_m_l3(2)],'Color',[0.45,0.45,0.45],'LineWidth',2);
    plot([rAf_m_l3(1) rBf_m_l3(1)],[rAf_m_l3(2) rBf_m_l3(2)],'Color',[0.45,0.45,0.45],'LineWidth',2);
    plot([rBf_m_l3(1) rCf_m_l3(1)],[rBf_m_l3(2) rCf_m_l3(2)],'Color',[0.45,0.45,0.45],'LineWidth',2);
    plot([0 rAf_m_l3(1) rBf_m_l3(1)],[0 rAf_m_l3(2) rBf_m_l3(2)],'o','Color',[0.45,0.45,0.45],'MarkerSize',5,'MarkerFaceColor',[0.45,0.45,0.45]);

    xlabel('X'); ylabel('Y'); legend({'Des','TO plain','TO w/ m_{eff}', 'TO w/ link3', 'TO w/ both'}, 'FontSize',6);
    axis equal; 
%     xlim([-0.5,1.5]); ylim([-1.0,1.0]);
%     xlim([-0.1, 1.5]); ylim([-0.7,0.7]);
    xlim([-0.1, 0.6]); ylim([-0.3,0.3]);

    
    % plot effective mass over time
    meff(1,1) = meff(1,2);
    meff_m(1,1) = meff_m(1,2); % shouldn't be any velocity at t=0 (just numerical noise from optimization), set to meff at next time-step for nice visuals 
    meff_l3(1,1) = meff_l3(1,2);
    meff_m_l3(1,1) = meff_m_l3(1,2);
    subplot(2,2,2); hold on;
    plot(time_vec,meff(1,:),'LineWidth',1.25); % add 'o-'
    plot(time_vec,meff_m(1,:),'LineWidth',1.25);
    plot(time_vec,meff_l3(1,:),'LineWidth',1.25);
    plot(time_vec,meff_m_l3(1,:),'LineWidth',1.25);
%     plot(time_vec,meff(2,:),'LineWidth',1.25);
%     plot(time_vec,meff(3,:),'LineWidth',1.25);
    xlabel('Time'); ylabel('m_{eff}'); 
%     legend('Actual','Min','Max');
    legend({'TO plain', 'TO w/ m_{eff}', 'TO w/ link3', 'TO w/ both'},'FontSize',6);
    ylim([0, 1.5]);
    title('Effective Mass');
%     ylim([0, 30]); % for UR3
    
    % plot norm of exit velocity over time
    subplot(2,2,4); hold on;
    plot(time_vec,exit_vel(3,:),'LineWidth',1.25);
    plot(time_vec,exit_vel_m(3,:),'LineWidth',1.25);
    plot(time_vec,exit_vel_l3(3,:),'LineWidth',1.25);
    plot(time_vec,exit_vel_m_l3(3,:),'LineWidth',1.25);
    xlabel('Time'); ylabel('||v_{exit}||');
    title('Exit Velocity');
%     legend({'TO plain', 'TO w/ m_{eff}', 'TO w/ link3', 'TO w/ both'},'FontSize',6);



    plt_title = sprintf('Optimized trajectory with effective mass and exit velocity for reference trajectory #%d', ii);
    sgtitle(plt_title);

    % save figure as frame
    F1(ii) = getframe(gcf);

    % also plot velocity and position tracking error over time
    figure(2); clf; % TODO: eventually combine with other figure?
    subplot(1,2,1); hold on;
    plot(time_vec, act_pts(3,:),'LineWidth',1.25);
    plot(time_vec, act_pts_m(3,:),'LineWidth',1.25);
    plot(time_vec, act_pts_l3(3,:),'LineWidth',1.25);
    plot(time_vec, act_pts_m_l3(3,:),'LineWidth',1.25);
    xlabel('Time'); ylabel('Position Error');
%     legend('TO plain', 'TO w/ m_{eff}', 'TO w/ link3', 'TO w/ both');
    subplot(1,2,2); hold on;
    plot(time_vec, verr, 'LineWidth',1.25);
    plot(time_vec, verr_m, 'LineWidth',1.25);
    plot(time_vec, verr_l3, 'LineWidth',1.25);
    plot(time_vec, verr_m_l3,'LineWidth',1.25);
    xlabel('Time'); ylabel('Velocity Error');
    legend('TO plain', 'TO w/ m_{eff}', 'TO w/ link3', 'TO w/ both');

    plt_title = sprintf('Tracking errors for reference trajectory #%d', ii);
    sgtitle(plt_title);

    % save figure as frame
    F2(ii) = getframe(gcf);
    
    
    % stop between each trajectory?
%     pause(0.5);
    pause;
    % TODO: make this a timed pause and record the animation as a video?
    
end

% % write videos here
% disp('Writing videos...');
% 
% v1 = VideoWriter(filename1);
% v1.FrameRate = 30; % optional: control frame rate
% open(v1)
% for ii=1:N % for each trajectory
%     for jj=1:60 % write 60 frames (two seconds) to video
%         writeVideo(v1,F1(ii));
%     end
% end
% close(v1)
% v2 = VideoWriter(filename2);
% v2.FrameRate = 30; % optional: control frame rate
% open(v2)
% for ii=1:N % for each trajectory
%     for jj=1:60 % write 60 frames (two seconds) to video
%         writeVideo(v2,F2(ii));
%     end
% end
% close(v2)
% disp('...done!');


