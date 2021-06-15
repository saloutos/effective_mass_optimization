% plotting real data vs TO data

% Andrew SaLoutos
% 5/21/2021

addpath(genpath('mc_arm_functions'));
clear all


%% import data and declare parameters

% data is [t, q1, q2, q3, dq1, dq2, dq3]

% [data_plain,plain_col_time] = process_raw_traj_data('Hardware Trajectories\Testing Data\linear_traj_18_plain_no_coll_1.txt');
% [data_CA,CA_col_time] = process_raw_traj_data('Hardware Trajectories\Testing Data\linear_traj_18_CA_no_coll_1.txt');

[data_plain,plain_col_time] = process_raw_traj_data('Hardware Trajectories\Testing Data\Linear18\L18_plain_pos2_trial2.txt');
[data_CA,CA_col_time] = process_raw_traj_data('Hardware Trajectories\Testing Data\Linear18\L18_CA_pos2_trial2.txt');

traj_num = 18;

t_start = 0.5/0.0015;
t_end = 1.5/0.0015;

% load desired data
load('multi_traj_data_linear_MC_2.mat');
% load('multi_traj_data_sinusoid_MC_2.mat');

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

%% pull out desired data

% for plain trajectory
time_vec = TO_data_plain(traj_num).time;
qdes = TO_data_plain(traj_num).data.q;
dqdes = TO_data_plain(traj_num).data.dq;
pts_des = TO_data_plain(traj_num).data.p;

% from imported data
time_real = data_plain(:,1)';
qreal = data_plain(:,2:4)';
dqreal = data_plain(:,5:7)';

% get positions
pts_real = zeros(2,length(time_real));
for ii=1:length(time_real)
    pt_temp = position_tip(qreal(:,ii),p);
    pts_real(:,ii) = pt_temp(1:2);
end

% for CA trajectory
time_vec_CA = TO_data_meff_link3(traj_num).time;
qdes_CA = TO_data_meff_link3(traj_num).data.q;
dqdes_CA = TO_data_meff_link3(traj_num).data.dq;
pts_des_CA = TO_data_meff_link3(traj_num).data.p;

% from imported data
time_real_CA = data_CA(:,1)';
qreal_CA = data_CA(:,2:4)';
dqreal_CA = data_CA(:,5:7)';

% get positions
pts_real_CA = zeros(2,length(time_real_CA));
for ii=1:length(time_real_CA)
    pt_temp = position_tip(qreal_CA(:,ii),p);
    pts_real_CA(:,ii) = pt_temp(1:2);
end


%% generate plots
figure(1); clf; 
subplot(1,2,1); hold on;
plot(time_vec,qdes(1,:),'--b');
plot(time_real(333:1000)-time_real(333),qreal(1,333:1000),'b');
plot(time_vec,qdes(2,:),'--r');
plot(time_real(333:1000)-time_real(333),qreal(2,333:1000),'r');
plot(time_vec,qdes(3,:),'--g');
plot(time_real(333:1000)-time_real(333),qreal(3,333:1000),'g');
xlabel('Time (s)'); ylabel('Position (rad)');
xlim([0,1]);
legend({'q_1, des','q_1, real','q_2, des','q_2, real','q_3, des','q_3, real'},'Location','southoutside','NumColumns',3);
title('Joint Positions for Linear Traj #18, Plain');

subplot(1,2,2); hold on;
plot(time_vec_CA,qdes_CA(1,:),'--b');
plot(time_real_CA(333:1000)-time_real_CA(333),qreal_CA(1,333:1000),'b');
plot(time_vec_CA,qdes_CA(2,:),'--r');
plot(time_real_CA(333:1000)-time_real_CA(333),qreal_CA(2,333:1000),'r');
plot(time_vec_CA,qdes_CA(3,:),'--g');
plot(time_real_CA(333:1000)-time_real_CA(333),qreal_CA(3,333:1000),'g');
xlabel('Time (s)'); ylabel('Position (rad)');
xlim([0,1]);
legend({'q_1, des','q_1, real','q_2, des','q_2, real','q_3, des','q_3, real'},'Location','southoutside','NumColumns',3);
title('Joint Positions for Linear Traj #18, CA');

figure(2); clf; hold on;
pd_pl = plot(pts_des(1,:),pts_des(2,:),'k');
pd_CA = plot(pts_des_CA(1,:),pts_des_CA(2,:),'--k');
pr_pl = plot(pts_real(1,:),pts_real(2,:),'r');
pr_CA = plot(pts_real_CA(1,:),pts_real_CA(2,:),'g');
x_a = [0.2, 0.5, 0.5, 0.2, 0.2]; % for MC, UR3 plots
y_a = [-0.15, -0.15, 0.15, 0.15, -0.15];

plot(x_a,y_a,'r--');
xlabel('X'); ylabel('Y'); axis equal;
xlim([0.1, 0.6]); ylim([-0.25,0.25]); 
legend([pd_pl, pr_pl, pd_CA, pr_CA], {'Des, Plain','Real, Plain','Des, CA','Real, CA'});

title('Cartesian Tracking for Linear Traj #18');

%% error statistics

time_real_int = time_real(333:999)-time_real(333);
qdes_full = [interp1(time_vec,qdes(1,:),time_real_int); interp1(time_vec,qdes(2,:),time_real_int); interp1(time_vec,qdes(3,:),time_real_int)];
dqdes_full = [interp1(time_vec,dqdes(1,:),time_real_int); interp1(time_vec,dqdes(2,:),time_real_int); interp1(time_vec,dqdes(3,:),time_real_int)];

time_real_int_CA = time_real_CA(333:999)-time_real_CA(333);
qdes_full_CA = [interp1(time_vec_CA,qdes_CA(1,:),time_real_int_CA); interp1(time_vec_CA,qdes_CA(2,:),time_real_int_CA); interp1(time_vec_CA,qdes_CA(3,:),time_real_int_CA)];
dqdes_full_CA = [interp1(time_vec_CA,dqdes_CA(1,:),time_real_int_CA); interp1(time_vec_CA,dqdes_CA(2,:),time_real_int_CA); interp1(time_vec_CA,dqdes_CA(3,:),time_real_int_CA)];

qerr = abs( qreal(:,333:999) - qdes_full );
qerr_CA = abs( qreal_CA(:,333:999) - qdes_full_CA );
qerr_avg = mean(qerr,2)
qerr_CA_avg = mean(qerr_CA(:,1:634),2)

dqerr = abs( dqreal(:,333:999) - dqdes_full );
dqerr_CA = abs( dqreal_CA(:,333:999) - dqdes_full_CA );
dqerr_avg = mean(dqerr,2)
dqerr_CA_avg = mean(dqerr_CA(:,1:634),2)

% simulate collision checking

qlim = [0.05, 0.05, 0.05];
dqlim = 2.0;

col_det_plain = zeros(5,size(qerr,2));
q1_check_plain = (qerr(1,:)>qlim(1));
q2_check_plain = (qerr(2,:)>qlim(2));
q3_check_plain = (qerr(3,:)>qlim(3));

dq3_check_plain = (dqerr(3,:)>dqlim);

for ii=5:size(qerr,2)
    if sum(q1_check_plain(ii-4:ii))==5
        col_det_plain(1,ii) = 1;
    end
    if sum(q2_check_plain(ii-4:ii))==5
        col_det_plain(2,ii) = 1;
    end
    if sum(q3_check_plain(ii-4:ii))==5
        col_det_plain(3,ii) = 1;
    end

    if (col_det_plain(1,ii)==1)&&(col_det_plain(2,ii)==1)&&(col_det_plain(3,ii)==1)
        col_det_plain(4,ii) = 1;
    end   
end
 
for ii=5:size(dqerr,2)
    if  dq3_check_plain(ii)==1
        col_det_plain(5,ii:end) = 5;
        break
    end
end

% 
% figure(3); clf;
% subplot(4,2,1); hold on;
% plot(qerr(1,:));
% plot(col_det_plain(1,:));
% subplot(4,2,3); hold on;
% plot(qerr(2,:));
% plot(col_det_plain(2,:));
% subplot(4,2,5); hold on;
% plot(qerr(3,:));
% plot(col_det_plain(3,:));
% subplot(4,2,7);
% plot(col_det_plain(4,:));
% 
% subplot(4,2,2);
% plot(dqerr(1,:));
% subplot(4,2,4);
% plot(dqerr(2,:));
% subplot(4,2,6); 



col_det_CA = zeros(5,size(qerr_CA,2));
q1_check_CA = (qerr_CA(1,:)>qlim(1));
q2_check_CA = (qerr_CA(2,:)>qlim(2));
q3_check_CA = (qerr_CA(3,:)>qlim(3));

dq3_check_CA = (dqerr_CA(3,:)>dqlim);

for ii=5:size(qerr_CA,2)
    if sum(q1_check_CA(ii-4:ii))==5
        col_det_CA(1,ii) = 1;
    end
    if sum(q2_check_CA(ii-4:ii))==5
        col_det_CA(2,ii) = 1;
    end
    if sum(q3_check_CA(ii-4:ii))==5
        col_det_CA(3,ii) = 1;
    end

    if (col_det_CA(1,ii)==1)&&(col_det_CA(2,ii)==1)&&(col_det_CA(3,ii)==1)
        col_det_CA(4,ii) = 1;
    end   
end
    
for ii=5:size(dqerr_CA,2)
    if  dq3_check_CA(ii)==1
        col_det_CA(5,ii:end) = 5;
        break
    end
end


% figure(4); clf;
% subplot(4,2,1); hold on;
% plot(qerr_CA(1,:));
% plot(col_det_CA(1,:));
% subplot(4,2,3); hold on;
% plot(qerr_CA(2,:));
% plot(col_det_CA(2,:));
% subplot(4,2,5); hold on;
% plot(qerr_CA(3,:));
% plot(col_det_CA(3,:));
% subplot(4,2,7);
% plot(col_det_CA(4,:));


% subplot(4,2,2);
% plot(dqerr_CA(1,:));
% subplot(4,2,4);
% plot(dqerr_CA(2,:));
% subplot(4,2,6); 


figure(6); clf; 
subplot(1,2,1); hold on;
plot(dqerr(3,:));
plot(col_det_plain(5,:));

subplot(1,2,2); hold on;
plot(dqerr_CA(3,:));
plot(col_det_CA(5,:));

        
        
        
        
        
        
        