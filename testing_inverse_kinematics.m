% testing inverse kinematics solutions given a trajectory


addpath(genpath('mc_arm_functions')) % for MC arm

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




%% import randomly generated trajectories

% traj_lib = load('random_linear_traj_MC_UR3.mat');
traj_lib = load('random_sinusoid_traj_MC_UR3.mat');

num_traj = size(traj_lib.trajectories,2);

%% define trajectory
num_pts = 51;
Tf = 1;
time_vec = linspace(0,Tf,num_pts); % include these in optimization variables?
dt = time_vec(2)-time_vec(1);

lin_prob = [8,24];
sin_prob = [11,12,15,19];

for jj=1:length(sin_prob)

    traj_ind = sin_prob(jj);    

    % if using linear trajectories
%     traj_sample = traj_lib.trajectories(:,traj_ind);
%     
%     pts_x = linspace(traj_sample(1),traj_sample(2),num_pts);
%     pts_y = linspace(traj_sample(3),traj_sample(4),num_pts);
%     
%     % desired end-effector trajectory
%     pts = [pts_x;pts_y];

    % if using sinusoidal trajectories
    pts = traj_lib.pts(:,:,traj_ind); % returns [pts_x;pts_y] of that trajectory

    % approximate velocities
    vels = diff(pts,1,2)/dt; % populate velocity vectors?
    vels = [zeros(2,1), vels]; % pad with zeros for first point in trajectory


    %% ik
    tic;
    z0 = zeros(6,num_pts);
    for ii=1:num_pts
        pt_i = pts(:,ii);
        vel_i = vels(:,ii);
%         z0(1:3,ii) = inverse_kinematics_init(p,pt_i); % get feasible ik solution for each point in trajectory...replace this with fmincon?
        z0(1:3,ii) = inverse_kinematics_perp(p,pt_i,vel_i);

        % given solution, get approximate joint velocities using jacobian pseudoinverse
        J_i = jacobian_tip(z0(1:3,ii),p);
        z0(4:6,ii) = pinv(J_i)*[vel_i;0];
    end
    ik_time = toc;

    % check ik with plot
    figure(1);
    pltt = sprintf('Traj #%d',traj_ind);
    for ii=1:num_pts
        t_i = time_vec(ii);
        plot_arm_kinematics_only(t_i,z0(:,ii),p,pts);
        title(pltt);
        pause(dt);
    end
   
end



%% helper functions

function q_ik = inverse_kinematics_init(p,ee_i)
    % pull out necessary parameters
    l1 = p(14);
    l2 = p(15);
    l3 = p(16);

    % pull out desired end-effector position
    % note: will have to do interpolation eventually here
    x_ee = ee_i(1);
    y_ee = ee_i(2);
    
    % set q1 initially
    % calculate q2 and q3
%     q1_ik = 0.0; %atan2(y_ee,x_ee)   
%     x_A = l1*cos(q1_ik);
%     y_A = l1*sin(q1_ik);
% 
%     % alternatively, could set th3 = q1+q2+q3 so that the third link is
%     % perpendicular to the desired velocity (to start)
%     % could also warm start with a first solve to find the configuration
%     % with the third link close to parallel to the desired velocity
% 
%     x_ee = x_ee-x_A % distances taken from second joint
%     y_ee = y_ee-y_A
% 
%     q3_ik = acos( (x_ee.^2 + y_ee.^2 - l2^2 - l3^2) / (2 * l2 * l3) ) % this can be +/-, leave as just + for now
%     q2_ik = atan2(y_ee,x_ee) - atan2( (l3*sin(q3_ik)), (l2+(l3*cos(q3_ik))) ) - q1_ik

    q3_ik = atan2(y_ee,x_ee);

    x_B = x_ee - l3*cos(q3_ik);
    y_B = y_ee - l3*sin(q3_ik);

    % ik from origin to [x_B,y_B]
    q2_ik = acos( (x_B.^2 + y_B.^2 - l1^2 - l2^2) / (2 * l1 * l2) ); % this can be +/-, leave as just + for now
    q1_ik = atan2(y_B,x_B) - atan2( (l2*sin(q2_ik)), (l1+(l2*cos(q2_ik))) );

    q3_ik = q3_ik - q2_ik - q1_ik;


    % outputs
    q_ik = [q1_ik;q2_ik;q3_ik];

end

function q_ik = inverse_kinematics_perp(p,ee_i,ve_i)
    % pull out necessary parameters
    l1 = p(14);
    l2 = p(15);
    l3 = p(16);

    % pull out desired end-effector position
    % note: will have to do interpolation eventually here
    x_ee = ee_i(1);
    y_ee = ee_i(2);

    vx_ee = ve_i(1);
    vy_ee = ve_i(2);
    
    qv = atan2(vy_ee,vx_ee); % angle of velocity vector
    
    % two candidate link 3 angles
    q3_1 = qv + pi/2;
    q3_2 = qv - pi/2;
    
    x_B1 = x_ee - l3*cos(q3_1);
    y_B1 = y_ee - l3*sin(q3_1);
    
    x_B2 = x_ee - l3*cos(q3_2);
    y_B2 = y_ee - l3*sin(q3_2);
    
    d1 = sqrt(x_B1^2 + y_B1^2);
    d2 = sqrt(x_B2^2 + y_B2^2);
    
    % select a direction
    if (d1<(l1+l2))||(d2<(l1+l2))
        if (d1<d2)
            q3_ik = q3_1;
            x_B = x_B1;
            y_B = y_B1;
        else
            q3_ik = q3_2;
            x_B = x_B2;
            y_B = y_B2;
        end
        
        q3_ik = q3_2;
        x_B = x_B2;
        y_B = y_B2;
        
        % ik from origin to [x_B,y_B]
        q2_ik = acos( (x_B.^2 + y_B.^2 - l1^2 - l2^2) / (2 * l1 * l2) ); % this can be +/-, leave as just + for now
        q1_ik = atan2(y_B,x_B) - atan2( (l2*sin(q2_ik)), (l1+(l2*cos(q2_ik))) );
        
        
        q3_ik = q3_ik - q2_ik - q1_ik;
        
    else
        % have link 3 aligned between end-effector location and origin
        
        % set q1 initially
        % calculate q2 and q3
        q3_ik = atan2(y_ee,x_ee);
        
        x_B = x_ee - l3*cos(q3_ik);
        y_B = y_ee - l3*sin(q3_ik);

        % ik from origin to [x_B,y_B]
        q2_ik = acos( (x_B.^2 + y_B.^2 - l1^2 - l2^2) / (2 * l1 * l2) ); % this can be +/-, leave as just + for now
        q1_ik = atan2(y_B,x_B) - atan2( (l2*sin(q2_ik)), (l1+(l2*cos(q2_ik))) );
        
        q3_ik = q3_ik - q2_ik - q1_ik;
        
    end  

    % outputs
    q_ik = [q1_ik;q2_ik;q3_ik];

end

function plot_arm_kinematics_only(t_i,z_i,p,pts)

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
   
    %%% plot these values 
    figure(1);
    clf; hold on;
    
    cur_time = sprintf('t = %.3f', t_i);
    
%     h_ti = text(-0.4,-0.8,cur_time);
    h_ti = text(-0.18,-0.25,cur_time);

    des_pts = plot(pts(1,:), pts(2,:),'g','LineWidth',2.0); % have this already plotted?
    
    % plot arena as well
    % TODO: don't have these dimensions hardcoded
    % TODO: have text showing traj_ind? include in title?
%     x_a = [0.4, 1.2, 1.2, 0.4, 0.4];
%     y_a = [-0.4, -0.4, 0.4, 0.4, -0.4];
    
    x_a = [0.2, 0.5, 0.5, 0.2, 0.2]; % for MC, UR3 plots
    y_a = [-0.15, -0.15, 0.15, 0.15, -0.15];
    
    arena_pts = plot(x_a,y_a,'r--');

    h_OA = plot([0 rA(1)],[0 rA(2)],'k','LineWidth',2);
    h_AB = plot([rA(1) rB(1)],[rA(2) rB(2)],'k','LineWidth',2);
    h_BC = plot([rB(1) rC(1)],[rB(2) rC(2)],'k','LineWidth',2);
    h_dots = plot([0 rA(1) rB(1)],[0 rA(2) rB(2)],'ok','MarkerSize',5,'MarkerFaceColor','k');

    xlabel('X'); ylabel('Y'); axis equal; 
    
%     xlim([-0.5,1.5]); ylim([-1.0,1.0]);
    xlim([-0.2, 0.8]); ylim([-0.3,0.3]);


end