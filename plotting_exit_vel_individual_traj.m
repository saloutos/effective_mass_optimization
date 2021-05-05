% plotting each TO output

% Andrew SaLoutos
% 4/17/2021

% make sure functions are on path
addpath(genpath('arm_functions'))
addpath(genpath('block_functions'))
addpath(genpath('arm_floating_functions'));

% import trajectory data
clear all

% load data
load('multi_traj_data_linear_ikp.mat');
% load('multi_traj_data_sinusoid_ikp.mat');

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

% define base parameters
mb = 0; Ib = 0.5*mb*0.2^2;
mb2 = 1; Ib2 = 0.5*mb2*0.2^2;

pb = [p; mb; Ib];
pb2 = [p; mb2; Ib2];

% define object inertia
mo = 1;

LLo = mo*eye(2);
LLo2 = mo*2*eye(2);

% iterate through each trajectory and plot
N = size(TO_data_plain,2);
m = length(TO_data_plain(1).time);
exit_vel_data = zeros(12,m,N);

exit_vel_data2 = zeros(12,m,N);


for ii=1:N
    if ischar(TO_data_plain(ii).data) || ischar(TO_data_meff(ii).data) || ischar(TO_data_link3(ii).data) || ischar(TO_data_meff_link3(ii).data)
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
        
        [phi_vf, Pvf, nu_vf, phi_vo, Pvo] = eval_exit_vel_metrics(q_temp,pb,LLo,ev_temp);
        [phi_vf_m, Pvf_m, nu_vf_m, phi_vo_m, Pvo_m] = eval_exit_vel_metrics(q_m_temp,pb,LLo,ev_m_temp);
        [phi_vf_l3, Pvf_l3, nu_vf_l3, phi_vo_l3, Pvo_l3] = eval_exit_vel_metrics(q_l3_temp,pb,LLo,ev_l3_temp);
        [phi_vf_m_l3, Pvf_m_l3, nu_vf_m_l3, phi_vo_m_l3, Pvo_m_l3] = eval_exit_vel_metrics(q_m_l3_temp,pb,LLo,ev_m_l3_temp);
        
        exit_vel_data(:,jj,ii) = [phi_vf; phi_vf_m; phi_vf_l3; phi_vf_m_l3; ...
                                  nu_vf; nu_vf_m; nu_vf_l3; nu_vf_m_l3; ...
                                  phi_vo; phi_vo_m; phi_vo_l3; phi_vo_m_l3];
          
        % using second set of inertial parameters...  
        [phi_vf, Pvf, nu_vf, phi_vo, Pvo] = eval_exit_vel_metrics(q_temp,pb2,LLo,ev_temp);
        [phi_vf_m, Pvf_m, nu_vf_m, phi_vo_m, Pvo_m] = eval_exit_vel_metrics(q_m_temp,pb2,LLo,ev_m_temp);
        [phi_vf_l3, Pvf_l3, nu_vf_l3, phi_vo_l3, Pvo_l3] = eval_exit_vel_metrics(q_l3_temp,pb2,LLo,ev_l3_temp);
        [phi_vf_m_l3, Pvf_m_l3, nu_vf_m_l3, phi_vo_m_l3, Pvo_m_l3] = eval_exit_vel_metrics(q_m_l3_temp,pb2,LLo,ev_m_l3_temp);
        
        exit_vel_data2(:,jj,ii) = [phi_vf; phi_vf_m; phi_vf_l3; phi_vf_m_l3; ...
                                   nu_vf; nu_vf_m; nu_vf_l3; nu_vf_m_l3; ...
                                   phi_vo; phi_vo_m; phi_vo_l3; phi_vo_m_l3];
                              
                      
    end   
end

%% plot?
figure(1);
for ii=1:N
    
    data = exit_vel_data(:,:,ii);
    time_vec = TO_data_plain(ii).time;
    
    clf; 
    
    data2 = exit_vel_data2(:,:,ii);
    
    subplot(1,2,1); hold on;
    plot(time_vec, data(1,:)); plot(time_vec, data(2,:));
    plot(time_vec, data(3,:)); plot(time_vec, data(4,:));
    xlabel('Time'); ylabel('Metric');
    legend('plain','meff','link3','both');
    title('\phi_{vf}');
    
    subplot(1,2,2); hold on;
    plot(time_vec, data2(1,:)); plot(time_vec, data2(2,:));
    plot(time_vec, data2(3,:)); plot(time_vec, data2(4,:));
    xlabel('Time'); ylabel('Metric');
    legend('plain','meff','link3','both');
    title('\phi_{vf},2');
    
    
%     subplot(1,3,1); hold on;
%     plot(time_vec, data(1,:)); plot(time_vec, data(2,:));
%     plot(time_vec, data(3,:)); plot(time_vec, data(4,:));
%     xlabel('Time'); ylabel('Metric');
%     legend('plain','meff','link3','both');
%     title('\phi_{vf}');
%     
%     subplot(1,3,2); hold on;
%     plot(time_vec, data(5,:)); plot(time_vec, data(6,:));
%     plot(time_vec, data(7,:)); plot(time_vec, data(8,:));
%     xlabel('Time'); ylabel('Metric');
%     legend('plain','meff','link3','both');
%     title('\eta_{vf}');
%     
%     subplot(1,3,3); hold on;
%     plot(time_vec, data(9,:)); plot(time_vec, data(10,:));
%     plot(time_vec, data(11,:)); plot(time_vec, data(12,:));
%     xlabel('Time'); ylabel('Metric');
%     legend('plain','meff','link3','both');
%     title('\phi_{vo}');

    plt_title = sprintf('Exit velocity metrics for trajectory #%d', ii);
    sgtitle(plt_title);
    
    pause;
end



