% Andrew SaLoutos
% 3/25/21

% trajectory optimization for effective mass, using direct transcription formulation
% also using fmincon solver (better candidates?)

%% fixed parameters    
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

%% define trajectory 
num_pts = 21;
Tf = 1;
time_vec = linspace(0,Tf,num_pts); % include these in optimization variables?
dt = time_vec(2)-time_vec(1);
pts_x = linspace(1.0,0.8,num_pts); %ones(1,num_pts);
pts_y = linspace(-0.4,0.4,num_pts);

% thetas = linspace(0,2*pi,num_pts+1); 
% thetas = thetas(1:num_pts);
% center = [0.9, -0.4]; %[0.9;0.3];
% radius = 0.2;
% pts_x = center(1) + radius*cos(thetas);
% pts_y = center(2) + radius*sin(thetas);

pts_x = linspace(0.6,1.0,num_pts);
pts_y = 0.2*sin(((1/0.4)*2*pi)*pts_x);

% desired end-effector trajectory
pts = [pts_x;pts_y];

vels = diff(pts,1,2)/dt; % populate velocity vectors?
vels = [zeros(2,1), vels]; % pad with zeros for first point in trajectory

%% run optimization

% Cost weights:       meff, meff_dir, e3_dir, u_mag, ee_pos, eef_pos, dq_mag, vee_mag, vee_dir
alpha_vec =         [  5.0,      0.0,    0.0,   0.5, 1000.0,     0.0,     0.0,    1.0,     1.0]; 
                    % [10 (or 0), 0 (or 10), 0, 0.5, 1000, 0, 0, 1, 1]

% actual optimization
[TO_data, opt_time] = meff_optimization(pts, Tf, p, alpha_vec);

% display opt time
opt_time_avg = opt_time/num_pts;
fprintf('\nTotal optimization time: %.1f ms\n\n',opt_time*1000.0);
fprintf('Optimization time per point: %.1f ms\n\n',opt_time_avg*1000.0);

% extract solution data
z_data_TO = TO_data(2:10,:);
q_data_TO = z_data_TO(1:3,:);
dq_data_TO = z_data_TO(4:6,:);
u_data_TO = z_data_TO(7:9,:);
meff_data_TO = TO_data(11:15,:);
peff_data_TO = TO_data(16:18,:);
ee_data_TO = TO_data(19:21,:);

%% simulate arm with TO solution as reference
sim_dt = 0.001;
sim_data = simulate_arm(sim_dt,time_vec,z_data_TO,p);
% sim data is [t(i); qff(i); dqff(i); q(i); dq(i); ddq(i); uff(i); u(i); ee_pos(i)]

% extract simulation data
time_vec_sim = sim_data(1,:);
qff_data_sim = sim_data(2:4,:);
dqff_data_sim = sim_data(5:7,:);
q_data_sim = sim_data(8:10,:);
dq_data_sim = sim_data(11:13,:);
ddq_data_sim = sim_data(14:16,:);
uff_data_sim = sim_data(17:19,:);
u_data_sim = sim_data(20:22,:);
meff_data_sim = sim_data(23:27,:);
peff_data_sim = sim_data(28:30,:);
ee_data_sim = sim_data(31:34,:);

%% generate plots and animations
figure(2); clf;
subplot(3,1,1); hold on;
plot(time_vec, q_data_TO); 
xlabel('Time'); ylabel('Joint Angle'); legend('q1','q2','q3');
subplot(3,1,2); hold on;
plot(time_vec, dq_data_TO);
xlabel('Time'); ylabel('Joint Velocity'); legend('dq1','dq2','dq3');
subplot(3,1,3); hold on;
plot(time_vec, u_data_TO);
xlabel('Time'); ylabel('Joint Torque'); legend('u1','u2','u3');

figure(3); clf; 
subplot(3,1,1); hold on;
plot(time_vec,meff_data_TO(1:3,:));
xlabel('Time'); ylabel('Effective Mass'); legend('Actual','Min','Max');
subplot(3,1,2); hold on;
plot(time_vec, peff_data_TO(3,:));
xlabel('Time'); ylabel('Effective Momentum');
subplot(3,1,3); hold on;
plot(time_vec, ee_data_TO(3,:));
xlabel('Time'); ylabel('Endpoint Position Error');

% animation of kinematics from optimization and dynamic simulation
% filename = 'TO_mod_cost_sine_curve_meff_direct.gif'; % save animation as a gif
figure(1);
for ii=1:(num_pts-1)
    t_i = time_vec(ii);
    z_i = [q_data_TO(:,ii); dq_data_TO(:,ii)];
    
    plot_arm_kinematics(t_i,z_i,p,meff_data_TO(:,ii),ee_data_TO,pts,vels(:,ii));
%     if ii==1
%         gif(filename,'DelayTime',dt);
%     else
%         gif;
%     end
    pause(dt);
end

% plots from simulation 
figure(5); clf; 
subplot(2,1,1); hold on;
plot(time_vec, q_data_TO); % ff values
plot(time_vec_sim, q_data_sim, '--');
xlabel('Time'); ylabel('Joint Angle'); legend('q1,ff','q2,ff','q3,ff','q1','q2','q3');
subplot(2,1,2); hold on;
plot(time_vec, dq_data_TO); % ff values
plot(time_vec_sim, dq_data_sim, '--'); % actual values
xlabel('Time'); ylabel('Joint Velocity'); legend('dq1,ff','dq2,ff','dq3,ff','dq1','dq2','dq3');

figure(6); clf;
subplot(3,1,1); hold on;
plot(time_vec, u_data_TO(1,:)); % ff
plot(time_vec_sim, u_data_sim(1,:), '--'); % actual values
xlabel('Time'); ylabel('Joint Torque'); legend('u1,ff','u1');
subplot(3,1,2); hold on;
plot(time_vec, u_data_TO(2,:)); % ff
plot(time_vec_sim, u_data_sim(2,:), '--'); % actual values
xlabel('Time'); ylabel('Joint Torque'); legend('u2,ff','u2');
subplot(3,1,3); hold on;
plot(time_vec, u_data_TO(3,:)); % ff
plot(time_vec_sim, u_data_sim(3,:), '--'); % actual values
xlabel('Time'); ylabel('Joint Torque'); legend('u3,ff','u3');

figure(7); clf;
subplot(3,1,1); hold on;
plot(time_vec,meff_data_TO(1:3,:));
plot(time_vec_sim,meff_data_sim(1,:));
xlabel('Time'); ylabel('Effective Mass'); legend('TO','Min','Max','Sim');
subplot(3,1,2); hold on;
plot(time_vec, peff_data_TO(3,:));
plot(time_vec_sim, peff_data_sim(3,:));
xlabel('Time'); ylabel('Effective Momentum'); legend('TO','Sim');
subplot(3,1,3); hold on;
ee_error_vec_sim = ee_data_sim(1:2,:)-ee_data_sim(3:4,:);
ee_error_mag_sim = sqrt( ee_error_vec_sim(1,:).^2 + ee_error_vec_sim(2,:).^2 );
plot(time_vec_sim, ee_error_mag_sim);
xlabel('Time'); ylabel('Endpoint Position Error from TO');

figure(4);
skip = 20;
for ii=1:skip:size(sim_data,2)
    t_i = sim_data(1,ii);
    z_i = sim_data(8:13,ii);
    plot_arm_simulation(t_i,z_i,p,pts,ee_data_sim);
    pause(skip*sim_dt);    
end

%% helper functions 
function plot_arm_kinematics(t_i,z_i,p,meff_sol_i,ee_sol,pts,vels_i)

    % doesn't create a new figure/plot, need to call subplot/figure before and also add title after
    % plot arm posture in black, with ellipsoid, effective mass directions, and velocity vector

    %%% plot arm kinematics
    q_i = z_i(1:3);
    dq_i = z_i(4:6);
    
    %%% keypoints and lmdir for final pose
    kp = keypoints_arm(z_i,p);
    rA = kp(:,1);
    rB = kp(:,2);
    rC = kp(:,3);
    % calculate ellipse for operational space inertia
    LLv_inv = LLv_arm_op_inv(z_i,p); % just care about x and y for now
    [V,D] = eig(LLv_inv);
    % pull out lmdir
    lm_dir = meff_sol_i(4:5);
    % calculate endpoint velocity
    ve_temp = jacobian_tip(z_i,p)*dq_i;
    ve_i = ve_temp(1:2); % get actual tip velocity
    ve_des = vels_i;
    ev_i = ve_i/norm(ve_i);
    
    % using parametric equations of ellipse, calculate ellipse points relative to foot position
    th_ellipse = atan2(V(2,1),V(1,1)); % angle between first eigenvector and positive x axis
    gamma_ell = 0.1; % TODO: better way to implement scaling of the ellipse?
    l_x = gamma_ell*sqrt((1/D(1,1))); 
    l_y = gamma_ell*sqrt((1/D(2,2))); 
    jj = linspace(0, 2*pi, 100);
    % make this belted ellipsoid?
    x_ell = (l_x*cos(jj))*cos(th_ellipse) - (l_y*sin(jj))*sin(th_ellipse);
    y_ell = (l_x*cos(jj))*sin(th_ellipse) + (l_y*sin(jj))*cos(th_ellipse);


    %%% plot these values 
    figure(1);
    clf; hold on;
    
    cur_time = sprintf('t = %.3f', t_i);
    h_ti = text(-0.4,-0.8,cur_time);
    
    des_pts = plot(pts(1,:), pts(2,:),'g','LineWidth',2.0); % have this already plotted?
    act_pts = plot(ee_sol(1,:), ee_sol(2,:),'Color',[0.1, 0.55, 0.1]); 
    
    ell_color = [0.85, 0.33, 0.];
    h_OA = plot([0 rA(1)],[0 rA(2)],'k','LineWidth',2);
    h_AB = plot([rA(1) rB(1)],[rA(2) rB(2)],'k','LineWidth',2);
    h_BC = plot([rB(1) rC(1)],[rB(2) rC(2)],'k','LineWidth',2);
    h_dots = plot([0 rA(1) rB(1)],[0 rA(2) rB(2)],'ok','MarkerSize',5,'MarkerFaceColor','k');
    h_LL = plot([rC(1)+x_ell],[rC(2)+y_ell],'LineWidth',1.5,'Color',ell_color);

    h_lm1 = quiver([rC(1)],[rC(2)],[lm_dir(1)],[lm_dir(2)],0.25,'r','LineWidth',2);
    h_lm2 = quiver([rC(1)],[rC(2)],[-lm_dir(1)],[-lm_dir(2)],0.25,'r','LineWidth',2);
    h_v_i = quiver([rC(1)],[rC(2)],[ve_i(1)],[ve_i(2)],0.25,'b','LineWidth',2); % real velocity is in dark blue
    h_v_des = quiver([rC(1)],[rC(2)],[ve_des(1)],[ve_des(2)],0.25,'Color',[0.6, 0.8, 1.0],'LineWidth',2); % desired velocity is in light blue
    xlabel('X'); ylabel('Y'); axis equal; 
    xlim([-0.5,1.5]); ylim([-1.0,1.0]);

end

function plot_arm_simulation(t_i,z_i,p,pts,ee_data)

    % doesn't create a new figure/plot, need to call subplot/figure before and also add title after
    % plot arm posture in black, with ellipsoid, effective mass directions, and velocity vector

    %%% plot arm kinematics
    q_i = z_i(1:3);
    dq_i = z_i(4:6);
    
    %%% keypoints and lmdir for final pose
    kp = keypoints_arm(z_i,p);
    rA = kp(:,1);
    rB = kp(:,2);
    rC = kp(:,3);
    % calculate endpoint velocity
    ve_temp = jacobian_tip(z_i,p)*dq_i;
    ve_i = ve_temp(1:2); % get actual tip velocity
    ev_i = ve_i/norm(ve_i);

    %%% plot these values 
    figure(4);
    clf; hold on;
    
    cur_time = sprintf('t = %.3f', t_i);
    h_ti = text(-0.4,-0.8,cur_time);
    
    des_pts = plot(pts(1,:), pts(2,:),'g','LineWidth',2.0); % have this already plotted?
    act_pts = plot(ee_data(1,:), ee_data(2,:),'Color',[0.1, 0.55, 0.1]); 
    ff_pts = plot(ee_data(3,:), ee_data(4,:),'r');
    
    h_OA = plot([0 rA(1)],[0 rA(2)],'k','LineWidth',2);
    h_AB = plot([rA(1) rB(1)],[rA(2) rB(2)],'k','LineWidth',2);
    h_BC = plot([rB(1) rC(1)],[rB(2) rC(2)],'k','LineWidth',2);
    h_dots = plot([0 rA(1) rB(1)],[0 rA(2) rB(2)],'ok','MarkerSize',5,'MarkerFaceColor','k');

    h_v_i = quiver([rC(1)],[rC(2)],[ve_i(1)],[ve_i(2)],0.25,'b','LineWidth',2); % real velocity is in dark blue
    xlabel('X'); ylabel('Y'); axis equal; 
    xlim([-0.5,1.5]); ylim([-1.0,1.0]);

end

function plot_arm_object_simulation(t_i,z_i,p,pts,ee_data)

    % doesn't create a new figure/plot, need to call subplot/figure before and also add title after
    % plot arm posture in black, with ellipsoid, effective mass directions, and velocity vector

    %%% plot arm kinematics
    q_i = z_i(1:3);
    dq_i = z_i(4:6);
    
    %%% keypoints and lmdir for final pose
    kp = keypoints_arm(z_i,p);
    rA = kp(:,1);
    rB = kp(:,2);
    rC = kp(:,3);
    % calculate endpoint velocity
    ve_temp = jacobian_tip(z_i,p)*dq_i;
    ve_i = ve_temp(1:2); % get actual tip velocity
    ev_i = ve_i/norm(ve_i);

    %%% plot these values 
    figure(4);
    clf; hold on;
    
    cur_time = sprintf('t = %.3f', t_i);
    h_ti = text(-0.4,-0.8,cur_time);
    
    des_pts = plot(pts(1,:), pts(2,:),'g','LineWidth',2.0); % have this already plotted?
    act_pts = plot(ee_data(1,:), ee_data(2,:),'Color',[0.1, 0.55, 0.1]); 
    ff_pts = plot(ee_data(3,:), ee_data(4,:),'r');
    
    h_OA = plot([0 rA(1)],[0 rA(2)],'k','LineWidth',2);
    h_AB = plot([rA(1) rB(1)],[rA(2) rB(2)],'k','LineWidth',2);
    h_BC = plot([rB(1) rC(1)],[rB(2) rC(2)],'k','LineWidth',2);
    h_dots = plot([0 rA(1) rB(1)],[0 rA(2) rB(2)],'ok','MarkerSize',5,'MarkerFaceColor','k');

    h_v_i = quiver([rC(1)],[rC(2)],[ve_i(1)],[ve_i(2)],0.25,'b','LineWidth',2); % real velocity is in dark blue
    xlabel('X'); ylabel('Y'); axis equal; 
    xlim([-0.5,1.5]); ylim([-1.0,1.0]);

end
