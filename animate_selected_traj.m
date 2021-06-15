% plotting animations of trajectory optimization outputs

% Andrew SaLoutos
% 5/19/21

% make sure functions are on path
addpath(genpath('block_functions'))
addpath(genpath('mc_arm_functions'))

% import trajectory data
clear all

% start with linear data
load('multi_traj_data_linear_MC_2.mat');

% make sure video filenames agree
filename = 'Figures/animate_all_valid_traj_MC_2.avi';

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

% define object parameters
mo = 1;
LLo = mo*eye(2);

% select trajectories
sel_traj = [18,27,30]; %[1:7,9:23,25:30]; % entries that aren't NaN for either row; %, 27, 30]; %[9, 15, 20]; % linear
% minimum traj are [5, 27, 30]

% iterate through each trajectory and plot
% nf = length(sel_traj)*51; % 51 frames per trajectory?

% struct to save frames for video
F1(48*51) = struct('cdata',[],'colormap',[]);

% % joint limits of hardware
% q_ub = [ 1.880,  2.716,  1.016];
% q_lb = [-1.880, -2.716, -1.540];

% for loop through trajectories
counter = 0;
for ii=1:length(sel_traj)
    
    % for each trajectory
    time_vec = TO_data_plain(sel_traj(ii)).time;
    des_pts = TO_data_plain(sel_traj(ii)).pts;
        
    X = TO_data_plain(sel_traj(ii)).data;
    X_m_l3 = TO_data_meff_link3(sel_traj(ii)).data;
    
    act_pts = X.p;
    act_pts_m_l3 = X_m_l3.p;

    meff = X.meff;
    meff_m_l3 = X_m_l3.meff;

    verr = X.v(3,:);
    verr_m_l3 = X_m_l3.v(3,:);

    q0 = X.q(:,1);
    qf = X.q(:,end);
    q0_m_l3 = X_m_l3.q(:,1);
    qf_m_l3 = X_m_l3.q(:,end);

    q = X.q;
    q_m_l3 = X_m_l3.q;  
    
    
    dq = X.dq;
    dq_m_l3 = X_m_l3.dq;
    
    LLv = X.LLv;
    LLv_m_l3 = X_m_l3.LLv;    
    
    % get exit velocities [vx; vy; ||v||];
    exit_vel = zeros(3,length(time_vec));
    exit_vel_m_l3 = zeros(3,length(time_vec));
    
    for jj=1:length(time_vec)
        
        v_temp = X.v(1:2,jj);
        LLv_temp = reshape(LLv(:,jj),2,2);
        ex_vel_temp = 2*inv(LLv_temp+LLo)*LLv_temp*v_temp;
        exit_vel(:,jj) = [ex_vel_temp; norm(ex_vel_temp)];
        
        v_temp = X_m_l3.v(1:2,jj);
        LLv_temp = reshape(LLv_m_l3(:,jj),2,2);
        ex_vel_temp = 2*inv(LLv_temp+LLo)*LLv_temp*v_temp;
        exit_vel_m_l3(:,jj) = [ex_vel_temp; norm(ex_vel_temp)];
        
    end
    
%     joint_up_lims = zeros(3,length(time_vec));
%     joint_low_lims = zeros(3,length(time_vec));
%     
%     for jj=1:length(time_vec)
%         for kk=1:3
%             joint_up_lims(kk,jj) = (q(kk,jj)>q_ub(kk));
%             joint_low_lims(kk,jj) = (q(kk,jj)<q_lb(kk));
%         end
%     end
%     
%     joint_up_lims
%     joint_low_lims
    
    
%     % print trajectory data for hardware tests
%     fprintf('\n\n Plain Linear Trajectory %d \n\n', sel_traj(ii));
%     for jj=1:length(time_vec)
%         if (jj==1)
%             fprintf('{{ %.3ff, %.3ff, %.3ff, %.3ff, %.3ff, %.3ff, %.3ff},\n', time_vec(jj), q(1,jj), q(2,jj), q(3,jj), dq(1,jj), dq(2,jj), dq(3,jj));
%         elseif (jj==length(time_vec))
%             fprintf('{ %.3ff, %.3ff, %.3ff, %.3ff, %.3ff, %.3ff, %.3ff}};\n', time_vec(jj), q(1,jj), q(2,jj), q(3,jj), dq(1,jj), dq(2,jj), dq(3,jj));
%         else
%             fprintf('{ %.3ff, %.3ff, %.3ff, %.3ff, %.3ff, %.3ff, %.3ff},\n', time_vec(jj), q(1,jj), q(2,jj), q(3,jj), dq(1,jj), dq(2,jj), dq(3,jj));
%         end
%     end
%     fprintf('\n\n CA Linear Trajectory %d \n\n', sel_traj(ii));
%     for jj=1:length(time_vec)
%         if (jj==1)
%             fprintf('{{ %.3ff, %.3ff, %.3ff, %.3ff, %.3ff, %.3ff, %.3ff},\n', time_vec(jj), q_m_l3(1,jj), q_m_l3(2,jj), q_m_l3(3,jj), dq_m_l3(1,jj), dq_m_l3(2,jj), dq_m_l3(3,jj));
%         elseif (jj==length(time_vec))
%             fprintf('{ %.3ff, %.3ff, %.3ff, %.3ff, %.3ff, %.3ff, %.3ff}};\n', time_vec(jj), q_m_l3(1,jj), q_m_l3(2,jj), q_m_l3(3,jj), dq_m_l3(1,jj), dq_m_l3(2,jj), dq_m_l3(3,jj));
%         else
%             fprintf('{ %.3ff, %.3ff, %.3ff, %.3ff, %.3ff, %.3ff, %.3ff},\n', time_vec(jj), q_m_l3(1,jj), q_m_l3(2,jj), q_m_l3(3,jj), dq_m_l3(1,jj), dq_m_l3(2,jj), dq_m_l3(3,jj));
%         end
%     end
%     
    % iterate through actual trajetory and animate
    for jj=1:length(time_vec)
        % increment overall counter
        counter = counter+1;
    
        % plot desired and actual trajectories
        figure(1); clf; 
        subplot(2,2,[1 3]); hold on;

        % always plot trajectories
        plot(des_pts(1,:),des_pts(2,:),'g','LineWidth',2.0);
        plot(act_pts(1,:),act_pts(2,:),'Color',[0.0, 0.4470, 0.7410],'LineWidth',1.25);
        plot(act_pts_m_l3(1,:),act_pts_m_l3(2,:),'Color',[0.85, 0.3250, 0.0980],'LineWidth',1.25);
        % always plot arena    
        x_a = [0.2, 0.5, 0.5, 0.2, 0.2]; % for MC, UR3 plots
        y_a = [-0.15, -0.15, 0.15, 0.15, -0.15];
        plot(x_a,y_a,'r--');
        
        % initial and final joint configurations
        % darker plot is CA optimized
        kp_m_l3 = keypoints_arm(q_m_l3(:,jj),p);
        rA_m_l3 = kp_m_l3(:,1);
        rB_m_l3 = kp_m_l3(:,2);
        rC_m_l3 = kp_m_l3(:,3);
        plot([0 rA_m_l3(1)],[0 rA_m_l3(2)],'Color',[0.3,0.3,0.3],'LineWidth',2);
        plot([rA_m_l3(1) rB_m_l3(1)],[rA_m_l3(2) rB_m_l3(2)],'Color',[0.3,0.3,0.3],'LineWidth',2);
        plot([rB_m_l3(1) rC_m_l3(1)],[rB_m_l3(2) rC_m_l3(2)],'Color',[0.3,0.3,0.3],'LineWidth',2);
        plot([0 rA_m_l3(1) rB_m_l3(1)],[0 rA_m_l3(2) rB_m_l3(2)],'o','Color',[0.3,0.3,0.3],'MarkerSize',5,'MarkerFaceColor',[0.3,0.3,0.3]);
        % lighter plot is plain
        kp = keypoints_arm(q(:,jj),p);
        rA = kp(:,1);
        rB = kp(:,2);
        rC = kp(:,3);
        plot([0 rA(1)],[0 rA(2)],'Color',[0.45,0.45,0.45],'LineWidth',2);
        plot([rA(1) rB(1)],[rA(2) rB(2)],'Color',[0.45,0.45,0.45],'LineWidth',2);
        plot([rB(1) rC(1)],[rB(2) rC(2)],'Color',[0.45,0.45,0.45],'LineWidth',2);
        plot([0 rA(1) rB(1)],[0 rA(2) rB(2)],'o','Color',[0.45,0.45,0.45],'MarkerSize',5,'MarkerFaceColor',[0.45,0.45,0.45]);

        xlabel('X'); ylabel('Y'); legend({'Des','Plain','CA'}, 'FontSize',12);
        axis equal; 
        xlim([-0.1, 0.6]); ylim([-0.3,0.3]);


        % plot effective mass over time
        meff(1,1) = meff(1,2); % shouldn't be any velocity at t=0 (just numerical noise from optimization), set to meff at next time-step for nice visuals 
        meff_m_l3(1,1) = meff_m_l3(1,2);
        subplot(2,2,2); hold on;
        plot(time_vec,meff(1,:),'LineWidth',1.25,'Color',[0.0, 0.4470, 0.7410]); % add 'o-'
        plot(time_vec,meff_m_l3(1,:),'LineWidth',1.25,'Color',[0.85, 0.3250, 0.0980]);
        
        plot(time_vec(jj),meff(1,jj),'o','LineWidth',1.25,'Color',[0.0, 0.4470, 0.7410]); % plot for current time as well
        plot(time_vec(jj),meff_m_l3(1,jj),'o','LineWidth',1.25,'Color',[0.85, 0.3250, 0.0980]);
        
        xlabel('Time (s)'); ylabel('m_{eff}  (kg)'); 
        legend({'Plain', 'CA'},'FontSize',12,'Location','northwest');
        ylim([0, 1.0]);
        title('Effective Mass');

        % plot norm of exit velocity over time
        subplot(2,2,4); hold on;
        plot(time_vec,exit_vel(3,:),'LineWidth',1.25,'Color',[0.0, 0.4470, 0.7410]);
        plot(time_vec,exit_vel_m_l3(3,:),'LineWidth',1.25,'Color',[0.85, 0.3250, 0.0980]);
        
        plot(time_vec(jj),exit_vel(3,jj),'o','LineWidth',1.25,'Color',[0.0, 0.4470, 0.7410]);
        plot(time_vec(jj),exit_vel_m_l3(3,jj),'o','LineWidth',1.25,'Color',[0.85, 0.3250, 0.0980]);
        ylim([0, 0.4]);
        xlabel('Time (s)'); ylabel('||v_{exit}|| (m/s)');
        legend({'Plain', 'CA'},'FontSize',12,'Location','northwest');
        title('Exit Velocity, for Object With m=1.0kg');

        plt_title = sprintf('Optimized Trajectory with Effective Mass and Exit Velocity for Linear Reference Trajectory #%d', sel_traj(ii));
        sgtitle(plt_title);

        
%         t = annotation('textbox', [0.2, 0.1, 0.2, 0], 'string', 'Trajectories at 0.2x speed');
%         t.LineStyle = 'none';
%         t.FontSize = 14;
        
        % save figure as frame
        F1(counter) = getframe(gcf);

        % stop between each trajectory
%         pause(0.1);
        
    end % move to next trajectory
    
end

figure; plot(des_pts(1,:),des_pts(2,:))

% switch to sinusoid data
load('multi_traj_data_sinusoid_MC_2.mat');
sel_traj = [6,13,16]; %1:20; %[6, 12, 20]; % sinusoids
% minimum traj are [14,18,19]


for ii=1:length(sel_traj)
    
    % for each trajectory
    time_vec = TO_data_plain(sel_traj(ii)).time;
    des_pts = TO_data_plain(sel_traj(ii)).pts;
    
    X = TO_data_plain(sel_traj(ii)).data;
    X_m_l3 = TO_data_meff_link3(sel_traj(ii)).data;

    act_pts = X.p;
    act_pts_m_l3 = X_m_l3.p;

    meff = X.meff;
    meff_m_l3 = X_m_l3.meff;

    verr = X.v(3,:);
    verr_m_l3 = X_m_l3.v(3,:);

    q0 = X.q(:,1);
    qf = X.q(:,end);
    q0_m_l3 = X_m_l3.q(:,1);
    qf_m_l3 = X_m_l3.q(:,end);

    q = X.q;
    q_m_l3 = X_m_l3.q;  
    
    dq = X.dq;
    dq_m_l3 = X_m_l3.dq;
    
    LLv = X.LLv;
    LLv_m_l3 = X_m_l3.LLv;    
    
    % get exit velocities [vx; vy; ||v||];
    exit_vel = zeros(3,length(time_vec));
    exit_vel_m_l3 = zeros(3,length(time_vec));
    
    for jj=1:length(time_vec)
        
        v_temp = X.v(1:2,jj);
        LLv_temp = reshape(LLv(:,jj),2,2);
        ex_vel_temp = 2*inv(LLv_temp+LLo)*LLv_temp*v_temp;
        exit_vel(:,jj) = [ex_vel_temp; norm(ex_vel_temp)];
        
        v_temp = X_m_l3.v(1:2,jj);
        LLv_temp = reshape(LLv_m_l3(:,jj),2,2);
        ex_vel_temp = 2*inv(LLv_temp+LLo)*LLv_temp*v_temp;
        exit_vel_m_l3(:,jj) = [ex_vel_temp; norm(ex_vel_temp)];
        
    end
    
%     joint_up_lims = zeros(3,length(time_vec));
%     joint_low_lims = zeros(3,length(time_vec));
%     
%     for jj=1:length(time_vec)
%         for kk=1:3
%             joint_up_lims(kk,jj) = (q(kk,jj)>q_ub(kk));
%             joint_low_lims(kk,jj) = (q(kk,jj)<q_lb(kk));
%         end
%     end
%     
%     joint_up_lims
%     joint_low_lims
%     
%     % print trajectory data for hardware tests
%     fprintf('\n\n Plain Sinusoid Trajectory %d \n\n', sel_traj(ii));
%     for jj=1:length(time_vec)
%         if (jj==1)
%             fprintf('{{ %.3ff, %.3ff, %.3ff, %.3ff, %.3ff, %.3ff, %.3ff},\n', time_vec(jj), q(1,jj), q(2,jj), q(3,jj), dq(1,jj), dq(2,jj), dq(3,jj));
%         elseif (jj==length(time_vec))
%             fprintf('{ %.3ff, %.3ff, %.3ff, %.3ff, %.3ff, %.3ff, %.3ff}};\n', time_vec(jj), q(1,jj), q(2,jj), q(3,jj), dq(1,jj), dq(2,jj), dq(3,jj));
%         else
%             fprintf('{ %.3ff, %.3ff, %.3ff, %.3ff, %.3ff, %.3ff, %.3ff},\n', time_vec(jj), q(1,jj), q(2,jj), q(3,jj), dq(1,jj), dq(2,jj), dq(3,jj));
%         end
%     end
%     fprintf('\n\n CA Sinusoid Trajectory %d \n\n', sel_traj(ii));
%     for jj=1:length(time_vec)
%         if (jj==1)
%             fprintf('{{ %.3ff, %.3ff, %.3ff, %.3ff, %.3ff, %.3ff, %.3ff},\n', time_vec(jj), q_m_l3(1,jj), q_m_l3(2,jj), q_m_l3(3,jj), dq_m_l3(1,jj), dq_m_l3(2,jj), dq_m_l3(3,jj));
%         elseif (jj==length(time_vec))
%             fprintf('{ %.3ff, %.3ff, %.3ff, %.3ff, %.3ff, %.3ff, %.3ff}};\n', time_vec(jj), q_m_l3(1,jj), q_m_l3(2,jj), q_m_l3(3,jj), dq_m_l3(1,jj), dq_m_l3(2,jj), dq_m_l3(3,jj));
%         else
%             fprintf('{ %.3ff, %.3ff, %.3ff, %.3ff, %.3ff, %.3ff, %.3ff},\n', time_vec(jj), q_m_l3(1,jj), q_m_l3(2,jj), q_m_l3(3,jj), dq_m_l3(1,jj), dq_m_l3(2,jj), dq_m_l3(3,jj));
%         end
%     end
    
    % iterate through actual trajetory and animate
    for jj=1:length(time_vec)
        % increment overall counter
        counter = counter+1;
    
        % plot desired and actual trajectories
        figure(1); clf; 
        subplot(2,2,[1 3]); hold on;

        % always plot trajectories
        plot(des_pts(1,:),des_pts(2,:),'g','LineWidth',2.0);
        plot(act_pts(1,:),act_pts(2,:),'Color',[0.0, 0.4470, 0.7410],'LineWidth',1.25);
        plot(act_pts_m_l3(1,:),act_pts_m_l3(2,:),'Color',[0.85, 0.3250, 0.0980],'LineWidth',1.25);
        % always plot arena    
        x_a = [0.2, 0.5, 0.5, 0.2, 0.2]; % for MC, UR3 plots
        y_a = [-0.15, -0.15, 0.15, 0.15, -0.15];
        plot(x_a,y_a,'r--');
        
        % initial and final joint configurations
        % darker plot is CA optimized
        kp_m_l3 = keypoints_arm(q_m_l3(:,jj),p);
        rA_m_l3 = kp_m_l3(:,1);
        rB_m_l3 = kp_m_l3(:,2);
        rC_m_l3 = kp_m_l3(:,3);
        plot([0 rA_m_l3(1)],[0 rA_m_l3(2)],'Color',[0.3,0.3,0.3],'LineWidth',2);
        plot([rA_m_l3(1) rB_m_l3(1)],[rA_m_l3(2) rB_m_l3(2)],'Color',[0.3,0.3,0.3],'LineWidth',2);
        plot([rB_m_l3(1) rC_m_l3(1)],[rB_m_l3(2) rC_m_l3(2)],'Color',[0.3,0.3,0.3],'LineWidth',2);
        plot([0 rA_m_l3(1) rB_m_l3(1)],[0 rA_m_l3(2) rB_m_l3(2)],'o','Color',[0.3,0.3,0.3],'MarkerSize',5,'MarkerFaceColor',[0.3,0.3,0.3]);
        % lighter plot is plain
        kp = keypoints_arm(q(:,jj),p);
        rA = kp(:,1);
        rB = kp(:,2);
        rC = kp(:,3);
        plot([0 rA(1)],[0 rA(2)],'Color',[0.45,0.45,0.45],'LineWidth',2);
        plot([rA(1) rB(1)],[rA(2) rB(2)],'Color',[0.45,0.45,0.45],'LineWidth',2);
        plot([rB(1) rC(1)],[rB(2) rC(2)],'Color',[0.45,0.45,0.45],'LineWidth',2);
        plot([0 rA(1) rB(1)],[0 rA(2) rB(2)],'o','Color',[0.45,0.45,0.45],'MarkerSize',5,'MarkerFaceColor',[0.45,0.45,0.45]);

        xlabel('X'); ylabel('Y'); legend({'Des','Plain','CA'}, 'FontSize',12);
        axis equal; 
        xlim([-0.1, 0.6]); ylim([-0.3,0.3]);


        % plot effective mass over time
        meff(1,1) = meff(1,2); % shouldn't be any velocity at t=0 (just numerical noise from optimization), set to meff at next time-step for nice visuals 
        meff_m_l3(1,1) = meff_m_l3(1,2);
        subplot(2,2,2); hold on;
        plot(time_vec,meff(1,:),'LineWidth',1.25,'Color',[0.0, 0.4470, 0.7410]); % add 'o-'
        plot(time_vec,meff_m_l3(1,:),'LineWidth',1.25,'Color',[0.85, 0.3250, 0.0980]);
        
        plot(time_vec(jj),meff(1,jj),'o','LineWidth',1.25,'Color',[0.0, 0.4470, 0.7410]); % plot for current time as well
        plot(time_vec(jj),meff_m_l3(1,jj),'o','LineWidth',1.25,'Color',[0.85, 0.3250, 0.0980]);
        
        xlabel('Time (s)'); ylabel('m_{eff} (kg)'); 
        legend({'Plain', 'CA'},'FontSize',12,'Location','northwest');
        ylim([0, 1.0]);
        title('Effective Mass');

        % plot norm of exit velocity over time
        subplot(2,2,4); hold on;
        plot(time_vec,exit_vel(3,:),'LineWidth',1.25,'Color',[0.0, 0.4470, 0.7410]);
        plot(time_vec,exit_vel_m_l3(3,:),'LineWidth',1.25,'Color',[0.85, 0.3250, 0.0980]);
        
        plot(time_vec(jj),exit_vel(3,jj),'o','LineWidth',1.25,'Color',[0.0, 0.4470, 0.7410]);
        plot(time_vec(jj),exit_vel_m_l3(3,jj),'o','LineWidth',1.25,'Color',[0.85, 0.3250, 0.0980]);
        ylim([0, 0.4]);
        xlabel('Time (s)'); ylabel('||v_{exit}|| (m/s)');
        legend({'Plain', 'CA'},'FontSize',12,'Location','northwest');
        title('Exit Velocity, for Object With m=1.0kg');

        plt_title = sprintf('Optimized Trajectory with Effective Mass and Exit Velocity for Sinusoid Reference Trajectory #%d', sel_traj(ii));
        sgtitle(plt_title);

        % save figure as frame
        F1(counter) = getframe(gcf);

        % stop between each trajectory
%         pause(0.1);
        
    end % move to next trajectory
    
end

% figure; plot(des_pts(1,:),des_pts(2,:))

%% write videos here
disp('Writing videos...');

v1 = VideoWriter(filename);
v1.FrameRate = 30; % optional: control frame rate
open(v1)
for ii=1:(48*51) % for each timestep in each trajectory
    for jj=1:3 % play trajectories at 0.2x speed...151 frames total at 30fps
        writeVideo(v1,F1(ii));
    end
end
close(v1)
disp('...done!');


