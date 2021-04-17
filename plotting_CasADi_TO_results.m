% plotting each TO output

% Andrew SaLoutos
% 4/17/2021

% make sure functions are on path
addpath(genpath('arm_functions'))
addpath(genpath('block_functions'))

% import trajectory data
clear all
load('multi_traj_data_unconstrained.mat');
load('multi_traj_data_no_meff.mat');
load('random_linear_traj.mat');

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


% iterate through each trajectory and plot
N = size(TO_data,2);
for ii=1:N
%     if (ii==6) || (ii==13) %could also use (TO_data.data == 'Bad')
%        % bad TO outputs
%        continue
%     end
    
    time_vec = TO_data(ii).time;
    des_pts = TO_data(ii).pts;
    
    
    X = TO_data(ii).data;
    X_nm = TO_data_no_meff(ii).data;
    act_pts = X.p;
    act_pts_nm = X_nm.p;
    meff = X.meff;
    meff_nm = X_nm.meff;
    q0 = X.q(:,1);
    qf = X.q(:,end);
    q0_nm = X_nm.q(:,1);
    qf_nm = X_nm.q(:,end);
    
    % plot desired and actual trajectories
    figure(1); clf; 
    subplot(1,2,1); hold on;
    
    % trajectories
    plot(des_pts(1,:),des_pts(2,:),'g','LineWidth',2.0);
    plot(act_pts(1,:),act_pts(2,:),'Color',[0.1, 0.55, 0.1],'LineWidth',1.25);
    plot(act_pts_nm(1,:),act_pts_nm(2,:),'Color',[1.0, 0.27, 0],'LineWidth',1.25);
    % arena
    x_a = [0.4, 1.2, 1.2, 0.4, 0.4];
    y_a = [-0.4, -0.4, 0.4, 0.4, -0.4];
    plot(x_a,y_a,'r--');
    % initial and final joint configurations
    % lighter plots are without meff in optimization
    kp0_nm = keypoints_arm(q0_nm,p);
    rA0_nm = kp0_nm(:,1);
    rB0_nm = kp0_nm(:,2);
    rC0_nm = kp0_nm(:,3);
    plot([0 rA0_nm(1)],[0 rA0_nm(2)],'Color',[0.25,0.25,0.25],'LineWidth',2);
    plot([rA0_nm(1) rB0_nm(1)],[rA0_nm(2) rB0_nm(2)],'Color',[0.25,0.25,0.25],'LineWidth',2);
    plot([rB0_nm(1) rC0_nm(1)],[rB0_nm(2) rC0_nm(2)],'Color',[0.25,0.25,0.25],'LineWidth',2);
    plot([0 rA0_nm(1) rB0_nm(1)],[0 rA0_nm(2) rB0_nm(2)],'o','Color',[0.25,0.25,0.25],'MarkerSize',5,'MarkerFaceColor',[0.25,0.25,0.25]);
    kpf_nm = keypoints_arm(qf_nm,p);
    rAf_nm = kpf_nm(:,1);
    rBf_nm = kpf_nm(:,2);
    rCf_nm = kpf_nm(:,3);
    plot([0 rAf_nm(1)],[0 rAf_nm(2)],'Color',[0.4,0.4,0.4],'LineWidth',2);
    plot([rAf_nm(1) rBf_nm(1)],[rAf_nm(2) rBf_nm(2)],'Color',[0.4,0.4,0.4],'LineWidth',2);
    plot([rBf_nm(1) rCf_nm(1)],[rBf_nm(2) rCf_nm(2)],'Color',[0.4,0.4,0.4],'LineWidth',2);
    plot([0 rAf_nm(1) rBf_nm(1)],[0 rAf_nm(2) rBf_nm(2)],'o','Color',[0.4,0.4,0.4],'MarkerSize',5,'MarkerFaceColor',[0.4,0.4,0.4]);
    kp0 = keypoints_arm(q0,p);
    rA0 = kp0(:,1);
    rB0 = kp0(:,2);
    rC0 = kp0(:,3);
    plot([0 rA0(1)],[0 rA0(2)],'Color',[0.4,0.4,0.4],'LineWidth',2);
    plot([rA0(1) rB0(1)],[rA0(2) rB0(2)],'Color',[0.4,0.4,0.4],'LineWidth',2);
    plot([rB0(1) rC0(1)],[rB0(2) rC0(2)],'Color',[0.4,0.4,0.4],'LineWidth',2);
    plot([0 rA0(1) rB0(1)],[0 rA0(2) rB0(2)],'o','Color',[0.4,0.4,0.4],'MarkerSize',5,'MarkerFaceColor',[0.4,0.4,0.4]);
    kpf = keypoints_arm(qf,p);
    rAf = kpf(:,1);
    rBf = kpf(:,2);
    rCf = kpf(:,3);
    plot([0 rAf(1)],[0 rAf(2)],'Color',[0.55,0.55,0.55],'LineWidth',2);
    plot([rAf(1) rBf(1)],[rAf(2) rBf(2)],'Color',[0.55,0.55,0.55],'LineWidth',2);
    plot([rBf(1) rCf(1)],[rBf(2) rCf(2)],'Color',[0.55,0.55,0.55],'LineWidth',2);
    plot([0 rAf(1) rBf(1)],[0 rAf(2) rBf(2)],'o','Color',[0.55,0.55,0.55],'MarkerSize',5,'MarkerFaceColor',[0.55,0.55,0.55]);
    
    xlabel('X'); ylabel('Y'); legend('Des','TO w/ m_{eff}','TO no m_{eff}');
    axis equal; 
    xlim([-0.5,1.5]); ylim([-1.0,1.0]);
    
    % plot effective mass over time
    subplot(1,2,2); hold on;
    plot(time_vec,meff(1,:),'LineWidth',1.25); % add 'o-'
    plot(time_vec,meff_nm(1,:),'LineWidth',1.25);
%     plot(time_vec,meff(2,:),'LineWidth',1.25);
%     plot(time_vec,meff(3,:),'LineWidth',1.25);
    xlabel('Time'); ylabel('m_{eff}'); 
%     legend('Actual','Min','Max');
    legend('TO w/ m_{eff}', 'TO no m_{eff}');
    
    plt_title = sprintf('Optimized trajectory and effective mass for reference trajectory #%d', ii);
    sgtitle(plt_title);
    
    % stop between each trajectory?
    pause;
    % TODO: make this a timed pause and record the animation as a video?
    
end


