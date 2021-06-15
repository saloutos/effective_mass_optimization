% plotting from all real trajectories

% Andrew SaLoutos
% 5/29/2021

%% arm functions and parameters

% clear all
addpath(genpath('mc_arm_functions'));

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

%% load desired data

% Linear #18
load('L18_real_data.mat'); % hardware data
load('linear18_video_data.mat'); % video data
obj_locs = [0.41, 0.35, 0.27; -0.07, 0.00, 0.10]; 

% Sinusoid #6
% load('S6_real_data.mat'); % hardware data
% load('sinusoid6_video_data.mat'); % video data
% obj_locs = [0.33, 0.37, 0.40; 0.05, 0.00, -0.10];

%% effective mass at collision time

% plain_col_time1 = mean( [plain_real(1).col_time, plain_real(2).col_time, plain_real(3).col_time, plain_real(4).col_time, plain_real(5).col_time] )-0.5;
% plain_col_time2 = mean( [plain_real(6).col_time, plain_real(7).col_time, plain_real(8).col_time, plain_real(9).col_time, plain_real(10).col_time] )-0.5;
% plain_col_time3 = mean( [plain_real(11).col_time, plain_real(12).col_time, plain_real(13).col_time, plain_real(13).col_time, plain_real(15).col_time] )-0.5;
% 
% CA_col_time1 = mean( [CA_real(1).col_time, CA_real(2).col_time, CA_real(3).col_time, CA_real(4).col_time, CA_real(5).col_time] )-0.5;
% CA_col_time2 = mean( [CA_real(6).col_time, CA_real(7).col_time, CA_real(8).col_time, CA_real(9).col_time, CA_real(10).col_time] )-0.5;
% CA_col_time3 = mean( [CA_real(11).col_time, CA_real(12).col_time, CA_real(13).col_time, CA_real(13).col_time, CA_real(15).col_time] )-0.5;
% 
% plain_des_meff_cols = interp1(plain_des.t, plain_des.meff, [plain_col_time1,plain_col_time2,plain_col_time3]);
% CA_des_meff_cols = interp1(CA_des.t, CA_des.meff, [CA_col_time1,CA_col_time2,CA_col_time3]);
% meff_des_rats = CA_des_meff_cols./plain_des_meff_cols;
% 
% meff_vals = zeros(15,2);
% for ii=1:15
%     meff_vals(ii,:) = [plain_real(ii).meff_col, CA_real(ii).meff_col];
% end
% 
% plain_meff_cols = [mean(meff_vals(1:5,1)), mean(meff_vals(6:10,1)), mean(meff_vals(11:15,1))];
% CA_meff_cols = [mean(meff_vals(1:5,2)), mean(meff_vals(6:10,2)), mean(meff_vals(11:15,2))];
% meff_rats = CA_meff_cols./plain_meff_cols;

%% plot for cartesian trajectories
figure(1); clf;


subplot(2,3,1);

% desired trajectory
hold on;
plot(plain_des.p(1,:),plain_des.p(2,:),'--','Color',[0.7,0.7,0.7],'LineWidth',1.5);
plot(CA_des.p(1,:),CA_des.p(2,:),'--','Color',[0.3,0.3,0.3],'LineWidth',1.5);

% object locations...not sure if I should include these
% subplot(1,2,1); 
% plot(obj_locs(1,:),obj_locs(2,:),'og','MarkerSize',17,'LineWidth',3);
% subplot(1,2,2); plot(obj_locs(1,:),obj_locs(2,:),'og','MarkerSize',17,'LineWidth',3);

% arena
% x_a = [0.2, 0.5, 0.5, 0.2, 0.2]; % for MC, UR3 plots
% y_a = [-0.15, -0.15, 0.15, 0.15, -0.15];
% plot(x_a,y_a,'k--');

% trials
for ii=1:15
    pl_trial = plot(plain_real(ii).p(1,333:999),plain_real(ii).p(2,333:999),'Color',[1.0,0,0,0.2]);
    ca_trial = plot(CA_real(ii).p(1,333:999),CA_real(ii).p(2,333:999),'Color',[0,0,1.0,0.2]);
end

% set axes equal
axis equal; xlim([0.2, 0.5]); ylim([-0.175,0.175]); 

% titles and axes lables
title('End-Effector Position, Linear #18');
xlabel('X (m)'); ylabel('Y (m)');
legend([pl_trial, ca_trial], {'Plain','CA'});

% plot for joint trajectories
figure(2); clf;
% joint positions
subplot(1,2,1); hold on;
plot(plain_des.t, plain_des.q(1,:),'b--','LineWidth',1.5);
plot(plain_des.t, plain_des.q(2,:),'r--','LineWidth',1.5);
plot(plain_des.t, plain_des.q(3,:),'g--','LineWidth',1.5);
plot(CA_des.t, CA_des.q(1,:),'b','LineWidth',1.5);
plot(CA_des.t, CA_des.q(2,:),'r','LineWidth',1.5);
plot(CA_des.t, CA_des.q(3,:),'g','LineWidth',1.5);
% trials
for ii=1:15
    plot(plain_real(ii).t(333:999)-plain_real(ii).t(333),plain_real(ii).q(1,333:999),'--','Color',[0.0,0.0,1.0,0.2]);
    plot(plain_real(ii).t(333:999)-plain_real(ii).t(333),plain_real(ii).q(2,333:999),'--','Color',[1.0,0.0,0.0,0.2]);
    plot(plain_real(ii).t(333:999)-plain_real(ii).t(333),plain_real(ii).q(3,333:999),'--','Color',[0.0,1.0,0.0,0.2]);
    plot(CA_real(ii).t(333:999)-CA_real(ii).t(333),CA_real(ii).q(1,333:999),'Color',[0.0,0.0,1.0,0.2]);
    plot(CA_real(ii).t(333:999)-CA_real(ii).t(333),CA_real(ii).q(2,333:999),'Color',[1.0,0.0,0.0,0.2]);
    plot(CA_real(ii).t(333:999)-CA_real(ii).t(333),CA_real(ii).q(3,333:999),'Color',[0.0,1.0,0.0,0.2]);  
end
% titles, labels, limits
xlim([0,1]);
legend({'q_1, plain','q_2, plain','q_3, plain','q_1, CA','q_2, CA','q_3, CA'},'Location','southoutside','NumColumns',2);
title('Joint Positions for Linear Traj #18, Plain');
xlabel('Time (s)'); ylabel('q (rad)');

% joint velocities
subplot(1,2,2); hold on;
plot(plain_des.t, plain_des.dq(1,:),'b--','LineWidth',1.5);
plot(plain_des.t, plain_des.dq(2,:),'r--','LineWidth',1.5);
plot(plain_des.t, plain_des.dq(3,:),'g--','LineWidth',1.5);
plot(CA_des.t, CA_des.dq(1,:),'b','LineWidth',1.5);
plot(CA_des.t, CA_des.dq(2,:),'r','LineWidth',1.5);
plot(CA_des.t, CA_des.dq(3,:),'g','LineWidth',1.5);
% trials
for ii=1:15
    plot(plain_real(ii).t(333:999)-plain_real(ii).t(333),plain_real(ii).dq(1,333:999),'--','Color',[0.0,0.0,1.0,0.2]);
    plot(plain_real(ii).t(333:999)-plain_real(ii).t(333),plain_real(ii).dq(2,333:999),'--','Color',[1.0,0.0,0.0,0.2]);
    plot(plain_real(ii).t(333:999)-plain_real(ii).t(333),plain_real(ii).dq(3,333:999),'--','Color',[0.0,1.0,0.0,0.2]);
    plot(CA_real(ii).t(333:999)-CA_real(ii).t(333),CA_real(ii).dq(1,333:999),'Color',[0.0,0.0,1.0,0.2]);
    plot(CA_real(ii).t(333:999)-CA_real(ii).t(333),CA_real(ii).dq(2,333:999),'Color',[1.0,0.0,0.0,0.2]);
    plot(CA_real(ii).t(333:999)-CA_real(ii).t(333),CA_real(ii).dq(3,333:999),'Color',[0.0,1.0,0.0,0.2]);  
end
% titles, labels, limits
xlim([0,1]);
legend({'dq_1, plain','dq_2, plain','dq_3, plain','dq_1, CA','dq_2, CA','dq_3, CA'},'Location','southoutside','NumColumns',2);
title('Joint Velocities for Linear Traj #18, Plain');
xlabel('Time (s)'); ylabel('dq (rad/s)');



% plot for object displacements
% could filter displacements?

% figure(3); clf; hold on;

figure(1);
subplot(2,3,[2 3]); hold on;

pos1 = [1,2,3,4,5,16,17,18,19,20];
pos2 = [6,7,8,9,10,21,22,23,24,25];
pos3 = [11,12,13,14,15,26,27,28,29,30];

mean_displacements = zeros(6,60);

for ii=1:10 % for all trials at a position
    nn = pos1(ii);
    if (ii<6)
        plot(video_data(nn).time-0.5, video_data(nn).obj_displacements,'--','Color',[1.0,0,0,0.4]); % plain trials  
        mean_displacements(1,:) = mean_displacements(1,:) + 0.2*video_data(nn).obj_displacements(1:60);
    else 
        plot(video_data(nn).time-0.5, video_data(nn).obj_displacements,'Color',[1.0,0,0,0.4]); % CA trials
        mean_displacements(2,:) = mean_displacements(2,:) + 0.2*video_data(nn).obj_displacements(1:60);
    end
end

for ii=1:10 % for all trials at a position
    nn = pos2(ii);
    if (ii<6)
        plot(video_data(nn).time-0.5, video_data(nn).obj_displacements,'--','Color',[0,0.5,0,0.4]); % plain trials
        mean_displacements(3,:) = mean_displacements(3,:) + 0.2*video_data(nn).obj_displacements(1:60);
    else 
        plot(video_data(nn).time-0.5, video_data(nn).obj_displacements,'Color',[0,0.5,0,0.4]); % CA trials
        mean_displacements(4,:) = mean_displacements(4,:) + 0.2*video_data(nn).obj_displacements(1:60);
    end
end

for ii=1:10 % for all trials at a position
    nn = pos3(ii);
    if (ii<6)
        plot(video_data(nn).time-0.5, video_data(nn).obj_displacements,'--','Color',[0,0,1.0,0.4]); % plain trials
        mean_displacements(5,:) = mean_displacements(5,:) + 0.2*video_data(nn).obj_displacements(1:60);
    else 
        plot(video_data(nn).time-0.5, video_data(nn).obj_displacements,'Color',[0,0,1.0,0.4]); % CA trials
        mean_displacements(6,:) = mean_displacements(6,:) + 0.2*video_data(nn).obj_displacements(1:60);
    end
end

mean_time = linspace(0,2,61);
mean_time = mean_time(1:60)-0.5;
p1_plain = plot(mean_time, mean_displacements(1,:),'r--','LineWidth',1.5);
p1_ca = plot(mean_time, mean_displacements(2,:),'r','LineWidth',1.5);
p2_plain = plot(mean_time, mean_displacements(3,:),'--','Color',[0,0.5,0],'LineWidth',1.5);
p2_ca = plot(mean_time, mean_displacements(4,:),'Color',[0,0.5,0],'LineWidth',1.5);
p3_plain = plot(mean_time, mean_displacements(5,:),'b--','LineWidth',1.5);
p3_ca = plot(mean_time, mean_displacements(6,:),'b','LineWidth',1.5);


% limits, labels, title
xlim([0,1]); xlabel('Time (s)'); ylabel('Object Displacement (mm)');
legend([p1_plain, p1_ca, p2_plain, p2_ca, p3_plain, p3_ca],{'Pos. 1, plain', 'Pos. 1, CA', 'Pos. 2, plain', 'Pos. 2, CA','Pos. 3, plain','Pos. 3, CA'},'Location','eastoutside');
title('Object Displacements, Linear #18');


%% now do same thing for sinusoid trajectories
% Sinusoid #6
load('S6_real_data.mat'); % hardware data
load('sinusoid6_video_data.mat'); % video data
obj_locs = [0.33, 0.37, 0.40; 0.05, 0.00, -0.10];

subplot(2,3,4);

% desired trajectory
hold on;
plot(plain_des.p(1,:),plain_des.p(2,:),'--','Color',[0.7,0.7,0.7],'LineWidth',1.5);
plot(CA_des.p(1,:),CA_des.p(2,:),'--','Color',[0.3,0.3,0.3],'LineWidth',1.5);

% trials
for ii=1:15
    pl_trial = plot(plain_real(ii).p(1,333:999),plain_real(ii).p(2,333:999),'Color',[1.0,0,0,0.2]);
    ca_trial = plot(CA_real(ii).p(1,333:999),CA_real(ii).p(2,333:999),'Color',[0,0,1.0,0.2]);
end

% set axes equal
axis equal; xlim([0.25, 0.45]); ylim([-0.125,0.125]); 
% titles and axes lables
title('End-Effector Position, Sinusoid #6');
xlabel('X (m)'); ylabel('Y (m)');
legend([pl_trial, ca_trial],{'Plain','CA'});

% second subplot
subplot(2,3,[5 6]); hold on;

pos1 = [1,2,3,4,5,16,17,18,19,20];
pos2 = [6,7,8,9,10,21,22,23,24,25];
pos3 = [11,12,13,14,15,26,27,28,29,30];

mean_displacements = zeros(6,60);

for ii=1:10 % for all trials at a position
    nn = pos1(ii);
    if (ii<6)
        plot(video_data(nn).time-0.5, video_data(nn).obj_displacements,'--','Color',[1.0,0,0,0.4]); % plain trials  
        mean_displacements(1,:) = mean_displacements(1,:) + 0.2*video_data(nn).obj_displacements(1:60);
    else 
        plot(video_data(nn).time-0.5, video_data(nn).obj_displacements,'Color',[1.0,0,0,0.4]); % CA trials
        mean_displacements(2,:) = mean_displacements(2,:) + 0.2*video_data(nn).obj_displacements(1:60);
    end
end

for ii=1:10 % for all trials at a position
    nn = pos2(ii);
    if (ii<6)
        plot(video_data(nn).time-0.5, video_data(nn).obj_displacements,'--','Color',[0,0.5,0,0.4]); % plain trials
        mean_displacements(3,:) = mean_displacements(3,:) + 0.2*video_data(nn).obj_displacements(1:60);
    else 
        plot(video_data(nn).time-0.5, video_data(nn).obj_displacements,'Color',[0,0.5,0,0.4]); % CA trials
        mean_displacements(4,:) = mean_displacements(4,:) + 0.2*video_data(nn).obj_displacements(1:60);
    end
end

for ii=1:10 % for all trials at a position
    nn = pos3(ii);
    if (ii<6)
        plot(video_data(nn).time-0.5, video_data(nn).obj_displacements,'--','Color',[0,0,1.0,0.4]); % plain trials
        mean_displacements(5,:) = mean_displacements(5,:) + 0.2*video_data(nn).obj_displacements(1:60);
    else 
        plot(video_data(nn).time-0.5, video_data(nn).obj_displacements,'Color',[0,0,1.0,0.4]); % CA trials
        mean_displacements(6,:) = mean_displacements(6,:) + 0.2*video_data(nn).obj_displacements(1:60);
    end
end

mean_time = linspace(0,2,61);
mean_time = mean_time(1:60)-0.5;
p1_plain = plot(mean_time, mean_displacements(1,:),'r--','LineWidth',1.5);
p1_ca = plot(mean_time, mean_displacements(2,:),'r','LineWidth',1.5);
p2_plain = plot(mean_time, mean_displacements(3,:),'--','Color',[0,0.5,0],'LineWidth',1.5);
p2_ca = plot(mean_time, mean_displacements(4,:),'Color',[0,0.5,0],'LineWidth',1.5);
p3_plain = plot(mean_time, mean_displacements(5,:),'b--','LineWidth',1.5);
p3_ca = plot(mean_time, mean_displacements(6,:),'b','LineWidth',1.5);

% limits, labels, title
xlim([0,1]); xlabel('Time (s)'); ylabel('Object Displacement (mm)');
legend([p1_plain, p1_ca, p2_plain, p2_ca, p3_plain, p3_ca],{'Pos. 1, plain', 'Pos. 1, CA', 'Pos. 2, plain', 'Pos. 2, CA','Pos. 3, plain','Pos. 3, CA'},'Location','eastoutside');
title('Object Displacements, Sinusoid #6');



% title for whole figure?
% sgtitle('Experimental Results (All Trials)');


figure(3);
% joint positions
subplot(1,2,1); hold on;
plot(plain_des.t, plain_des.q(1,:),'b--','LineWidth',1.5);
plot(plain_des.t, plain_des.q(2,:),'r--','LineWidth',1.5);
plot(plain_des.t, plain_des.q(3,:),'g--','LineWidth',1.5);
plot(CA_des.t, CA_des.q(1,:),'b','LineWidth',1.5);
plot(CA_des.t, CA_des.q(2,:),'r','LineWidth',1.5);
plot(CA_des.t, CA_des.q(3,:),'g','LineWidth',1.5);
% trials
for ii=1:15
    plot(plain_real(ii).t(333:999)-plain_real(ii).t(333),plain_real(ii).q(1,333:999),'--','Color',[0.0,0.0,1.0,0.2]);
    plot(plain_real(ii).t(333:999)-plain_real(ii).t(333),plain_real(ii).q(2,333:999),'--','Color',[1.0,0.0,0.0,0.2]);
    plot(plain_real(ii).t(333:999)-plain_real(ii).t(333),plain_real(ii).q(3,333:999),'--','Color',[0.0,1.0,0.0,0.2]);
    plot(CA_real(ii).t(333:999)-CA_real(ii).t(333),CA_real(ii).q(1,333:999),'Color',[0.0,0.0,1.0,0.2]);
    plot(CA_real(ii).t(333:999)-CA_real(ii).t(333),CA_real(ii).q(2,333:999),'Color',[1.0,0.0,0.0,0.2]);
    plot(CA_real(ii).t(333:999)-CA_real(ii).t(333),CA_real(ii).q(3,333:999),'Color',[0.0,1.0,0.0,0.2]);  
end
% titles, labels, limits
xlim([0,1]);
legend({'q_1, plain','q_2, plain','q_3, plain','q_1, CA','q_2, CA','q_3, CA'},'Location','southoutside','NumColumns',2);
title('Joint Positions for Sinusoid Traj #6, Plain');
xlabel('Time (s)'); ylabel('q (rad)');

% joint velocities
subplot(1,2,2); hold on;
plot(plain_des.t, plain_des.dq(1,:),'b--','LineWidth',1.5);
plot(plain_des.t, plain_des.dq(2,:),'r--','LineWidth',1.5);
plot(plain_des.t, plain_des.dq(3,:),'g--','LineWidth',1.5);
plot(CA_des.t, CA_des.dq(1,:),'b','LineWidth',1.5);
plot(CA_des.t, CA_des.dq(2,:),'r','LineWidth',1.5);
plot(CA_des.t, CA_des.dq(3,:),'g','LineWidth',1.5);
% trials
for ii=1:15
    plot(plain_real(ii).t(333:999)-plain_real(ii).t(333),plain_real(ii).dq(1,333:999),'--','Color',[0.0,0.0,1.0,0.2]);
    plot(plain_real(ii).t(333:999)-plain_real(ii).t(333),plain_real(ii).dq(2,333:999),'--','Color',[1.0,0.0,0.0,0.2]);
    plot(plain_real(ii).t(333:999)-plain_real(ii).t(333),plain_real(ii).dq(3,333:999),'--','Color',[0.0,1.0,0.0,0.2]);
    plot(CA_real(ii).t(333:999)-CA_real(ii).t(333),CA_real(ii).dq(1,333:999),'Color',[0.0,0.0,1.0,0.2]);
    plot(CA_real(ii).t(333:999)-CA_real(ii).t(333),CA_real(ii).dq(2,333:999),'Color',[1.0,0.0,0.0,0.2]);
    plot(CA_real(ii).t(333:999)-CA_real(ii).t(333),CA_real(ii).dq(3,333:999),'Color',[0.0,1.0,0.0,0.2]);  
end
% titles, labels, limits
xlim([0,1]);
legend({'dq_1, plain','dq_2, plain','dq_3, plain','dq_1, CA','dq_2, CA','dq_3, CA'},'Location','southoutside','NumColumns',2);
title('Joint Velocities for Sinusoid Traj #6, Plain');
xlabel('Time (s)'); ylabel('dq (rad/s)');
