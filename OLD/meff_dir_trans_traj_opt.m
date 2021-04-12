% Andrew SaLoutos
% 3/9/21

% trajectory optimization for effective mass, using direct transcription formulation
% also using fmincon solver (better candidates?)

% fixed parameters    
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

% define trajectory points here
num_pts = 31;
Tf = 1;
time_vec = linspace(0,Tf,num_pts); % include these in optimization variables?
dt = time_vec(2)-time_vec(1);
pts_x = linspace(1.0,0.8,num_pts); %ones(1,num_pts);
pts_y = linspace(-0.4,0.4,num_pts);

thetas = linspace(0,2*pi,num_pts+1); 
thetas = thetas(1:num_pts);
center = [1.0, -0.1]; %[0.9;0.3];
radius = 0.2;

% pts_x = center(1) + radius*cos(thetas);
% pts_y = center(2) + radius*sin(thetas);

pts = [pts_x;pts_y];
vels = diff(pts,1,2)/dt; % populate velocity vectors?
vels = [zeros(2,1), vels]; % pad with zeros for first point in trajectory
% vels_x = zeros(1,num_pts);
% vels_y = ones(1,num_pts);

% symbolic variables for optimization
% at each time point: n inputs, n joint angles, n joint velocities
syms z [9 num_pts] real % [q_i[1:3]; dq_i[1:3]; u_i[1:3]] for each point, i
z = z(:); % vectorize


% initial values
z0 = zeros(9,num_pts);
for ii=1:num_pts
    ee_i = [pts_x(ii);pts_y(ii)];
    z0(1:3,ii) = inverse_kinematics_init(p,ee_i); % get feasible ik solution for each point in trajectory
    % given solution, get approximate joint velocities using jacobian pseudoinverse
    J_i = jacobian_tip(z0(1:3,ii),p);
    z0(4:6,ii) = pinv(J_i)*[vels(:,ii);0];
end
z0 = z0(:);

% fmincon inputs
cost_func = @(z)meff_cost(z,p,time_vec,pts,vels);
nonlcon = @(z)joint_constraints(z,p,time_vec,pts,vels);
A = []; % no other linear inequalities
b = [];
Aeq = []; % constraints on initial and final positions?
beq = [];
[ub,lb] = joint_and_torque_limits(z,p,time_vec); % torque and joint velocity limits
options = optimoptions('fmincon','display','none');

% run optimization
tic;
[z_sol, fval_sol, sol_hist] = fmincon_with_hist(cost_func,z0,A,b,Aeq,beq,lb,ub,nonlcon,options);
opt_time = toc;
opt_time_avg = opt_time/num_pts;
z_sol = reshape(z_sol,9,num_pts);

% display opt time
fprintf('\nTotal optimization time: %.1f ms\n\n',opt_time*1000.0);
fprintf('Optimization time per point in trajectory: %.1f ms\n\n',opt_time_avg*1000.0);

% plot solutions:
% 1. plot joint angles, velocities, and input torques; 
% 2. plot effective mass vs minimum and maximum effective mass; 
% 3. animate optimization solution (for now, just animate positions? could do full simulation later)

q_sol = z_sol(1:3,:);
dq_sol = z_sol(4:6,:);
u_sol = z_sol(7:9,:);
meff_sol = zeros(5,num_pts); % each column is [meff; meff_min; meff_max; lm_dir]
ee_sol = zeros(2,num_pts);
for ii=1:num_pts
    z_i = [q_sol(:,ii); dq_sol(:,ii)];
    % calculate endpoint velocity
    ve_temp = jacobian_tip(z_i,p)*dq_sol(:,ii);
    ve_i = ve_temp(1:2); % get actual tip velocity
    ev_i = ve_i/norm(ve_i);
    % calculate effective mass, minimum, maximum
    LLv_inv = LLv_arm_op_inv(z_i,p);
    meff_sol(1,ii) = 1/(ev_i'*LLv_inv*ev_i);
    [V,D] = eig(LLv_inv);
    if D(1,1)>D(2,2)
        meff_sol(2,ii) = 1/D(1,1);
        meff_sol(3,ii) = 1/D(2,2);
        meff_sol(4:5,ii) = V(:,1);
    else
        meff_sol(2,ii) = 1/D(2,2);
        meff_sol(3,ii) = 1/D(1,1);
        meff_sol(4:5,ii) = V(:,2);
    end
    % get endpoint position
    ee_temp = position_tip(z_i,p);
    ee_sol(:,ii) = ee_temp(1:2);
end

figure(2); clf;
subplot(3,1,1); hold on;
plot(time_vec, q_sol); 
xlabel('Time'); ylabel('Joint Angle'); legend('q1','q2','q3');
subplot(3,1,2); hold on;
plot(time_vec, dq_sol);
xlabel('Time'); ylabel('Joint Velocity'); legend('dq1','dq2','dq3');
subplot(3,1,3); hold on;
plot(time_vec, u_sol);
xlabel('Time'); ylabel('Joint Torque'); legend('u1','u2','u3');

figure(3); clf; hold on;
plot(time_vec,meff_sol(1:3,:));
xlabel('Time'); ylabel('Effective Mass'); legend('Actual','Min','Max');

% animation

% filename = 'TO_circle_2.gif'; % save animation as a gif
figure(1);
for ii=1:(num_pts-1)
    z_i = [q_sol(:,ii); dq_sol(:,ii)];
    plot_arm_kinematics(z_i,p,meff_sol(:,ii),ee_sol,pts,vels(:,ii));
%     if ii==1
%         gif(filename,'DelayTime',dt);
%     else
%         gif;
%     end
    pause(dt);
end

%%%%% helper functions %%%%%
function [c,ceq] = joint_constraints(z,p,time_vec,pts,vels)
    
    num_pts = length(time_vec);
    dt = time_vec(2)-time_vec(1);
    
    z = reshape(z,9,num_pts);
    q = z(1:3,:);
    dq = z(4:6,:);
    u = z(7:9,:);
    
    c = [];
    ceq = [];

    % could also calculate gradients of constraints here?
    
    % dynamics constraints... using helper functions -> A*ddq = b -> x[i+1] = x[i] + dt*dx/dt[i]
    for ii=2:num_pts
        
        z_i1 = [q(:,ii-1); dq(:,ii-1)];
        u_i1 = u(:,ii-1);
        A = A_arm(z_i1,p);
        b = b_arm(z_i1,u_i1,p);
        ddq_i1 = A\b;

        %q(:,ii) = q(:,ii-1) + dt*dq(:,ii-1);
        %dq(:,ii) = dq(:,ii-1) + dt*ddq_i
        
        ceq = [ceq; q(:,ii)-q(:,ii-1)-dt*dq(:,ii-1); dq(:,ii)-dq(:,ii-1)-dt*ddq_i1];
        
    end
        
    % end effector position constraint? at least for initial point, makes sense for all points until a better exploration objective is defined
    % end effector velocity constraint? not sure if this should be a constraint or a cost
    for ii=1:num_pts
        % get end effector position at time_vec(ii)
        ee_i = pts(:,ii);
        z_i = [q(:,ii); dq(:,ii)];
        ee_temp = position_tip(z_i,p);
        ceq_q_i = ee_temp(1:2)-ee_i;
        
        % get end effector velocity (in x,y) at time_vec(ii)
%         ve_des = [vels_x(ii);vels_y(ii)];
%         ve_temp = jacobian_tip(z_i,p)*dq(:,ii);
%         ve_i = ve_temp(1:2);
%         ve_i = ve_i/norm(ve_i);
%         ve_des = ve_des/norm(ve_des);
%         ceq_dq_i = ve_i-ve_des; % direction of end effector velocity must match 
        
        ceq = [ceq; ceq_q_i];% ceq_dq_i];
    end
        
end

function total_cost = meff_cost(z,p,time_vec,pts,vels)
    
    num_pts = length(time_vec);
    dt = time_vec(2)-time_vec(1);
    
    z = reshape(z,9,num_pts);
    q = z(1:3,:);
    dq = z(4:6,:);
    u = z(7:9,:);
    
    meff_cost_vec = zeros(num_pts,1);
    dq_cost_vec = zeros(num_pts,1);
    u_cost_vec = zeros(num_pts,1);
    v_cost_vec = zeros(num_pts,1);
    
    % running cost on meff along vdir
    % also have running cost on alignment with traj velocity
    % starting with zero velocity, so cost starts with second point in
    % trajectory
    for ii=2:num_pts
        z_i = [q(:,ii); dq(:,ii)];
        ve_temp = jacobian_tip(z_i,p)*dq(:,ii);
        ve_i = ve_temp(1:2); % get actual tip velocity
        ve_des = vels(:,ii);
        ev_i = ve_i/norm(ve_i);
        ev_des = ve_des/norm(ve_des);
%         ev_i = ve_des/norm(ve_des);
        
        LLv_inv = LLv_arm_op_inv(z_i,p);
%         meff_cost_vec(ii) = -(ev_des'*LLv_inv*ev_des); % negative seems nicer to work with than inverse
        meff_cost_vec(ii) = -(ev_i'*LLv_inv*ev_i); % negative seems nicer to work with than inverse
        
        v_cost_vec(ii) = -ev_i'*ev_des; % maximize alignment of the two vectors
        
    end
    
    % running cost on dq between positions?
    R1 = eye(3);
%     R1(3,3) = 0.1; % don't really care about joint 3 velocity
    R1(1,1) = 0.2;
    for ii=2:num_pts
        dq_i = dq(:,ii);
        dq_cost_vec(ii) = dq_i'*R1*dq_i;
    end
    
    % cost on inputs
    R2 = eye(3);
    R2(3,3) = 0.1; % don't really care about joint 3 torque
    for ii=1:num_pts
        u_i = u(:,ii);
        u_cost_vec(ii) = u_i'*R2*u_i;
    end
    
    % add costs together with weights
    a1 = 10.0*ones(1,num_pts); 
%     a1(1) = 3; % higher weight on effective mass at the start of the trajectory
    c1 = a1*meff_cost_vec;
    
    a2 = 0.1*ones(1,num_pts);
    c2 = a2*dq_cost_vec;
    
    a3 = 0.05*ones(1,num_pts);
    c3 = a3*u_cost_vec;
    
    a4 = 0.2*ones(1,num_pts);
    c4 = a4*v_cost_vec;

    % return total cost
    total_cost = c1 + c2 + c3 + c4;
    
end

function [ub,lb] = joint_and_torque_limits(z,p,time_vec)

    num_pts = length(time_vec);
    dt = time_vec(2)-time_vec(1);
    
    z = reshape(z,9,num_pts);
    q = z(1:3,:);
    dq = z(4:6,:);
    u = z(7:9,:);
    
    ub = [];
    lb = [];
    
    % define limits?
    q_ub = [pi;pi;pi];
    q_lb = [-pi;-pi;-pi];
    
    dq_ub = [30;30;30];
    dq_lb = [-30;-30;-30];
    
    u_ub = [20;20;20];
    u_lb = [-20;-20;-20];
    
    for ii=1:num_pts
        
        q_i = q(:,ii);
        dq_i = dq(:,ii);
        u_i = u(:,ii);
        
        ub = [ub; q_ub; dq_ub; u_ub]; % format needs to match vectorized optimization variables
        lb = [lb; q_lb; dq_lb; u_lb];        
        
    end

end

function [xsol, fval, history] = fmincon_with_hist(cost_func,x0,A,b,Aeq,beq,lb,ub,nonlcon_func,prev_options)
    history = [];
    options = optimoptions(prev_options,'OutputFcn', @save_output);

    [xsol, fval] = fmincon(cost_func,x0,A,b,Aeq,beq,lb,ub,nonlcon_func,options);    

    function stop = save_output(x,optimvalues,state)
        stop = false;
        new_opt_data = [optimvalues.fval, optimvalues.constrviolation, x']; % for this use, new_opt_data = [cost, q1, q2, q3]
        if isequal(state,'iter')
          history = [history; new_opt_data];
        end
    end

end

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
    q1_ik = 0; 
    x_A = l1*cos(q1_ik);
    y_A = l1*sin(q1_ik);

    % alternatively, could set th3 = q1+q2+q3 so that the third link is
    % perpendicular to the desired velocity (to start)
    % could also warm start with a first solve to find the configuration
    % with the third link close to parallel to the desired velocity

    x_ee = x_ee-x_A; % distances taken from second joint
    y_ee = y_ee-y_A;

    q3_ik = acos( (x_ee.^2 + y_ee.^2 - l2^2 - l3^2) / (2 * l2 * l3) ); % this can be +/-, leave as just + for now
    q2_ik = atan2(y_ee,x_ee) - atan2( (l3*sin(q3_ik)), (l2+(l3*cos(q3_ik))) );

    % outputs
    q_ik = [q1_ik;q2_ik;q3_ik];
end

function plot_arm_kinematics(z_i,p,meff_sol_i,ee_sol,pts,vels_i)

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
    
    des_pts = plot(pts(1,:), pts(2,:),'g','LineWidth',2); % have this already plotted?
    act_pts = plot(ee_sol(1,:), ee_sol(2,:),'Color',[0.1, 0.55, 0.1]); 
    
    ell_color = [0.85, 0.33, 0.];
    h_OA = plot([0 rA(1)],[0 rA(2)],'k','LineWidth',2);
    h_AB = plot([rA(1) rB(1)],[rA(2) rB(2)],'k','LineWidth',2);
    h_BC = plot([rB(1) rC(1)],[rB(2) rC(2)],'k','LineWidth',2);
    h_dots = plot([0 rA(1) rB(1)],[0 rA(2) rB(2)],'ok','MarkerSize',5,'MarkerFaceColor','k');
    h_LL = plot([rC(1)+x_ell],[rC(2)+y_ell],'LineWidth',1.5,'Color',ell_color);

    h_lm1 = quiver([rC(1)],[rC(2)],[lm_dir(1)],[lm_dir(2)],0.25,'r','LineWidth',2);
    h_lm2 = quiver([rC(1)],[rC(2)],[-lm_dir(1)],[-lm_dir(2)],0.25,'r','LineWidth',2);
    h_v_i = quiver([rC(1)],[rC(2)],[ve_i(1)],[ve_i(2)],0.25,'b','LineWidth',2);
    h_v_des = quiver([rC(1)],[rC(2)],[ve_des(1)],[ve_des(2)],0.25,'Color',[0.6, 0.8, 1.0],'LineWidth',2);
    xlabel('X'); ylabel('Y'); axis equal; 
    xlim([-0.5,1.5]); ylim([-1.0,1.0]);

end
