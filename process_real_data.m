% process real data from text files into structs for plotting

% Andrew SaLoutos
% 5/21/2021

clear all
addpath(genpath('mc_arm_functions'));

%% import data and declare parameters

% saved data is [t, q1, q2, q3, dq1, dq2, dq3]

t_start = 0.5/0.0015;
t_end = 1.5/0.0015;

% load desired data

% load('multi_traj_data_linear_MC_2.mat'); % switch filenames below as well
% traj_num = 18;

load('multi_traj_data_sinusoid_MC_2.mat'); % switch filenames below as well
traj_num = 6;

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

%% iterate through plain trajectory data

plain_real = struct('pos',[],'trial',[],'col_time',[],'meff_col',[],'t',[],'q',[],'dq',[],'p',[],'v',[]);
plain_real = repmat(plain_real,1,15);

for ii=1:15
    posloc = floor((ii-1)/5)+1;
    trialnum = mod((ii-1),5)+1;
    next_filename = sprintf('Hardware Trajectories\\Testing Data\\Sinusoid6\\S6_plain_pos%d_trial%d.txt',posloc,trialnum);
    [data_plain,plain_col_time] = process_raw_traj_data(next_filename);
    
    % from imported data
    time_real = data_plain(:,1)';
    qreal = data_plain(:,2:4)';
    dqreal = data_plain(:,5:7)';

    % get positions and velocities
    pts_real = zeros(2,length(time_real));
    vels_real = zeros(2,length(time_real));
    for jj=1:length(time_real)
        pt_temp = position_tip(qreal(:,jj),p);
        pts_real(:,jj) = pt_temp(1:2);
        J_temp = jacobian_v_tip(qreal(:,jj),p);
        vels_real(:,jj) = J_temp*dqreal(:,jj);
    end
    
    % get effective mass at collision
    [~,cii] = min(abs(time_real-plain_col_time));
    
    vel_col = mean( vels_real(1:2,(cii-10):cii), 2 );
    u = vel_col/norm(vel_col);
    
    q_col = mean( qreal(:,(cii-10):cii), 2);
    dq_col = mean( dqreal(:,(cii-10):cii), 2);
    
    zcol = [q_col; dq_col];
    
    LLv_inv = LLv_arm_op_inv(zcol,p);
    meff_col = 1/(u'*LLv_inv*u);
    
    % store data
    plain_real(ii).pos = posloc;
    plain_real(ii).trial = trialnum;
    plain_real(ii).col_time = plain_col_time;
    plain_real(ii).meff_col = meff_col;
    plain_real(ii).t = time_real;
    plain_real(ii).q = qreal;
    plain_real(ii).dq = dqreal;
    plain_real(ii).p = pts_real;
    plain_real(ii).v = vels_real;
    
end

% for plain trajectory
plain_des = struct('t',[],'q',[],'dq',[],'p',[],'v',[],'meff',[]);
plain_des.t = TO_data_plain(traj_num).time;
plain_des.q = TO_data_plain(traj_num).data.q;
plain_des.dq = TO_data_plain(traj_num).data.dq;
plain_des.p = TO_data_plain(traj_num).data.p;
plain_des.v = TO_data_plain(traj_num).data.v;
plain_des.meff = TO_data_plain(traj_num).data.meff(1,:);

%% iterate through CA trajectory data

CA_real = struct('pos',[],'trial',[],'col_time',[],'meff_col',[],'t',[],'q',[],'dq',[],'p',[],'v',[]);
CA_real = repmat(CA_real,1,15);

for ii=1:15
    posloc = floor((ii-1)/5)+1;
    trialnum = mod((ii-1),5)+1;
    next_filename = sprintf('Hardware Trajectories\\Testing Data\\Sinusoid6\\S6_CA_pos%d_trial%d.txt',posloc,trialnum);
    [data_CA,CA_col_time] = process_raw_traj_data(next_filename);
    
    % from imported data
    time_real = data_CA(:,1)';
    qreal = data_CA(:,2:4)';
    dqreal = data_CA(:,5:7)';

    % get positions and velocities
    pts_real = zeros(2,length(time_real));
    vels_real = zeros(2,length(time_real));
    for jj=1:length(time_real)
        pt_temp = position_tip(qreal(:,jj),p);
        pts_real(:,jj) = pt_temp(1:2);
        J_temp = jacobian_v_tip(qreal(:,jj),p);
        vels_real(:,jj) = J_temp*dqreal(:,jj);
    end
    
    % get effective mass at collision
    [~,cii] = min(abs(time_real-CA_col_time));
    
    vel_col = mean( vels_real(1:2,(cii-10):cii), 2 );
    u = vel_col/norm(vel_col);
    
    q_col = mean( qreal(:,(cii-10):cii), 2);
    dq_col = mean( dqreal(:,(cii-10):cii), 2);
    
    zcol = [q_col; dq_col];
    
    LLv_inv = LLv_arm_op_inv(zcol,p);
    meff_col = 1/(u'*LLv_inv*u);
    
    % store data
    CA_real(ii).pos = posloc;
    CA_real(ii).trial = trialnum;
    CA_real(ii).col_time = CA_col_time;
    CA_real(ii).meff_col = meff_col;
    CA_real(ii).t = time_real;
    CA_real(ii).q = qreal;
    CA_real(ii).dq = dqreal;
    CA_real(ii).p = pts_real;
    CA_real(ii).v = vels_real;
    
end

% for CA trajectory
CA_des = struct('t',[],'q',[],'dq',[],'p',[],'v',[],'meff',[]);
CA_des.t = TO_data_meff_link3(traj_num).time;
CA_des.q = TO_data_meff_link3(traj_num).data.q;
CA_des.dq = TO_data_meff_link3(traj_num).data.dq;
CA_des.p = TO_data_meff_link3(traj_num).data.p;
CA_des.v = TO_data_meff_link3(traj_num).data.v;
CA_des.meff = TO_data_meff_link3(traj_num).data.meff(1,:);


        
        
        
        
        
        
        