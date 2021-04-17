clear
% close all

%% Add Libraries
addpath(genpath('casadi'));
addpath(genpath('spatial_v2'));
addpath(genpath('arm_4_functions'))
addpath(genpath('arm_functions'))
addpath(genpath('block_functions'))
import casadi.* 

%% define arm parameters
% TODO: Use spatial v2?
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

rad1 = 0.05;
rad2 = 0.05;
rad3 = 0.05;
rad_mult = 1.2;

%% Build Robot Model
% TODO: Turn this into a function
% TODO: Use this model for forward dynamics constraint
% TODO: How to include rotor dynamics...
model.NB = 3;                                  % number of bodies
model.gravity = [0 0 -9.81];                   % gravity
model.parent  = zeros(1,model.NB);             % parent body indices
model.jtype   = repmat({'  '},model.NB,1);     % joint types
model.Xtree   = repmat({eye(3)},model.NB,1);   % coordinate transforms
model.I       = repmat({zeros(3)},model.NB,1); % spatial inertias
nb = 0; % current body index

model.appearance.base = ...
    {'colour',[0.9 0.9 0.9],'cyl', [0 0 0;0 0 -0.3], 3*rad1};
% Joint 1
nb = nb + 1;
model.parent(nb) = nb - 1;
model.jtype{nb}  = 'r';
model.Xtree{nb}  = eye(3);
model.I{nb}      = mcI(m1, [l_O_m1 0], I1);
cyl1 = [0 0 0;...
    l_OA 0 0];
cyl1b = [0 0 -rad_mult*rad1;...
    0 0 rad_mult*rad1];
model.appearance.body{nb} = {'colour',[1.0 0.1 0.1],...
    'cyl', cyl1, rad1, 'cyl', cyl1b, rad1};
% Joint 2
nb = nb + 1;
model.parent(nb) = nb - 1;
model.jtype{nb}  = 'r';
model.Xtree{nb}  = plnr(0, [l_OA 0]);
model.I{nb}      = mcI(m2, [l_A_m2 0], I2);
cyl2 = [0 0 0;...
    l_AB 0 0];
cyl2b = [0 0 -rad_mult*rad2;...
    0 0 rad_mult*rad2];
model.appearance.body{nb} = {'colour',[0.1 1.0 0.1],...
    'cyl', cyl2, rad2, 'cyl', cyl2b, rad2};
% Joint 3
nb = nb + 1;
model.parent(nb) = nb - 1;
model.jtype{nb}  = 'r';
model.Xtree{nb}  = plnr(0, [l_AB 0]);
model.I{nb}      = mcI(m3, [l_B_m3 0], I3);
cyl3 = [0 0 0;...
    l_BC 0 0];
cyl3b = [0 0 -rad_mult*rad3;...
    0 0 rad_mult*rad3];
model.appearance.body{nb} = {'colour',[0.1 0.1 1.0],...
    'cyl', cyl3, rad3, 'cyl', cyl3b, rad3};


%% define object parameters
xc = 0.8; yc = 0; radius = 0.05;        % cylinder geometry parameters
mo = 1; Io = 0.5*mo*radius*radius;    % inertial parameters
mu_o = 0.5;                           % friction coefficient with surface
% arrange in vector
p_obj = [mo Io radius mu_o g xc yc]'; % appended starting location of block

%% define trajectory
num_pts = 51;
Tf = 1;
time_vec = linspace(0,Tf,num_pts); % include these in optimization variables?
dt = time_vec(2)-time_vec(1);

% pts_x = linspace(1.0,0.8,num_pts); %ones(1,num_pts);
% pts_y = linspace(-0.4,0.4,num_pts);

% thetas = linspace(0,2*pi,num_pts+1);
% thetas = thetas(1:num_pts);
% center = [0.9, -0.4]; %[0.9;0.3];
% radius = 0.2;
% pts_x = center(1) + radius*cos(thetas);
% pts_y = center(2) + radius*sin(thetas);
% 
% pts_x = linspace(0.6,1.0,num_pts); % big sine wave
% pts_y = 0.2*sin(((1/0.4)*2*pi)*pts_x);

% pts_x = linspace(0.7,0.9,num_pts); % tiny sine wave
% pts_y = 0.1*sin(((1/0.2)*2*pi)*pts_x);

% import randomly generated linear trajectory
traj_lib = load('random_linear_traj.mat');
traj_ind = 13; %randi(size(traj_lib.trajectories,2));
traj_sample = traj_lib.trajectories(:,traj_ind);

pts_x = linspace(traj_sample(1),traj_sample(2),num_pts);
pts_y = linspace(traj_sample(3),traj_sample(4),num_pts);

% desired end-effector trajectory
pts = [pts_x;pts_y];

vels = diff(pts,1,2)/dt; % populate velocity vectors?
vels = [zeros(2,1), vels]; % pad with zeros for first point in trajectory

%% Optimization Variables

opti = casadi.Opti();
% Optimization variables
% q       (3xn) % TODO: don't hard code numberof DOFs
% dq      (3xn)
% u     (3xn)
X = opti.variable(3*3, num_pts);
opt_var.X   = X(:,:);
opt_var.q   = X(1:3,:);
opt_var.dq  = X(4:6,:);
opt_var.u   = X(7:9,:);

% TODO: create optimization parameters
% (replace things like p and alpha_vec with opt_param)

% Cost weights:        meff,      p,     v,   u,  dq
alpha_vec =         [  10.0, 1000.0, 100.0, 0.5, 1.0];

%% Cost Function
% TODO: turn this into a function of its own
meff_cost_vec = casadi.MX(zeros(num_pts,1));
p_cost_vec = casadi.MX(zeros(num_pts,1));
v_cost_vec = casadi.MX(zeros(num_pts,1));
u_cost_vec = casadi.MX(zeros(num_pts,1));
dq_cost_vec = casadi.MX(zeros(num_pts,1));

% running cost on meff along vdir
% also have running cost on alignment with traj velocity
% starting with zero velocity, so cost starts with second point in
% trajectory

% TODO: remove constraints on initial and final positions, add higher costs?

for ii=2:num_pts
    z_i = [opt_var.q(:,ii); opt_var.dq(:,ii)];
    
    ve_temp = jacobian_tip(z_i,p)*opt_var.dq(:,ii);
    ve_i = ve_temp(1:2); % get actual tip velocity
    ve_des = vels(:,ii);
    ev_i = ve_i/norm(ve_i);
    
    LLv_inv = LLv_arm_op_inv(z_i,p);
    %         meff_cost_vec(ii) = 1/(ev_i'*LLv_inv*ev_i); % negative seems nicer to work with than inverse
    meff_cost_vec(ii) = -(ev_i'*LLv_inv*ev_i); % negative seems nicer to work with than inverse
    %         meff_cost_vec(ii) = -((ev_i'*LLv_inv*ev_i)^2); % can square to get better gradients?
    
    % end-effector velocity cost
    Qvdes = eye(2);
    v_cost_vec(ii) = (ve_i-ve_des)'*Qvdes*(ve_i-ve_des);
    
end

% running cost on dq between positions?
R1 = eye(3);
for ii=2:num_pts
    dq_i = opt_var.dq(:,ii);
    dq_cost_vec(ii) = dq_i'*R1*dq_i;
end


R2 = eye(3);
Qpdes = eye(2);
for ii=1:num_pts
    % cost on inputs
    u_i = opt_var.u(:,ii);
    u_cost_vec(ii) = u_i'*R2*u_i;
    % end-effector position cost
    z_i = [opt_var.q(:,ii); opt_var.dq(:,ii)];
    p_temp = position_tip(z_i,p);
    p_diff = p_temp(1:2) - pts(:,ii);
    p_cost_vec(ii) = p_diff'*Qpdes*p_diff; 
    % higher cost on initial and final positions
    if (ii==1)||(ii==num_pts)
        p_cost_vec(ii) = 100*p_cost_vec(ii);
    end
end

% add costs together with weights

% cost on effective mass
c1 = alpha_vec(1)*sum(meff_cost_vec);

% cost on end-effector position (relative to trajectory)
c2 = alpha_vec(2)*sum(p_cost_vec);

% cost on direction of end-effector velocity (relative to trajectory)
c3 = alpha_vec(3)*sum(v_cost_vec);

% cost on inputs
c4 = alpha_vec(4)*sum(u_cost_vec);

% cost on magnitude of joint velocities
c5 = alpha_vec(5)*sum(dq_cost_vec);

% return total cost
total_cost = c1 + c2 + c3 + c4 + c5;

% set cost function
opti.minimize( total_cost );

%% Nonlinear Constraints
% dynamics constraints... using helper functions -> A*ddq = b -> x[i+1] = x[i] + dt*dx/dt[i]
for ii=2:num_pts
    
    z_i1 = [opt_var.q(:,ii-1); opt_var.dq(:,ii-1)];
    u_i1 = opt_var.u(:,ii-1);
    A = A_arm(z_i1,p);
    A_inv = A_inv_arm(z_i1,p);
    b = b_arm(z_i1,u_i1,p);
    %ddq_i1 = A\b;
    ddq_i1 = A_inv*b;
    
    % TODO: might also just be easier to use the spatial v2 package
    
    %q(:,ii) = q(:,ii-1) + dt*dq(:,ii-1);
    %dq(:,ii) = dq(:,ii-1) + dt*ddq_i
    opti.subject_to(opt_var.q(:,ii)-opt_var.q(:,ii-1)-dt*opt_var.dq(:,ii-1) == zeros(3,1))
    opti.subject_to(opt_var.dq(:,ii)-opt_var.dq(:,ii-1)-dt*ddq_i1 == zeros(3,1))
    
end

% TODO: look at adding a constraint on maximum effective momentum?

ee_0 = pts(:,1);
z_0 = [opt_var.q(:,1); opt_var.dq(:,1)];
ee_temp0 = position_tip(z_0,p);
%ceq_q_0 = [ee_temp0(1:2)-ee_0; opt_var.dq(:,1)];
% opti.subject_to(ee_temp0(1:2)-ee_0 == zeros(2,1))
opti.subject_to(opt_var.dq(:,1) == zeros(3,1))

ee_f = pts(:,end);
z_f = [opt_var.q(:,end); opt_var.dq(:,end)];
ee_tempf = position_tip(z_f,p);
%ceq_q_f = [ee_tempf(1:2)-ee_f]; %; dq(:,end)];
% opti.subject_to(ee_tempf(1:2)-ee_f == zeros(2,1))

%% Joint and Torque Limits

% define limits
q_ub = [pi;pi;pi];
q_lb = [-pi;-pi;-pi];

dq_ub = [30;30;30];
dq_lb = [-30;-30;-30];

u_ub = [20;20;20];
u_lb = [-20;-20;-20];

for ii=1:num_pts
    
    q_i = opt_var.q(:,ii);
    dq_i = opt_var.dq(:,ii);
    u_i = opt_var.u(:,ii);
    
    opti.subject_to(q_i <= q_ub);
    opti.subject_to(dq_i <= dq_ub);
    opti.subject_to(u_i <= u_ub);
    
    opti.subject_to(q_i >= q_lb);
    opti.subject_to(dq_i >= dq_lb);
    opti.subject_to(u_i >= u_lb);
    
end

%% Generate Initial Guess
% initial values
z0 = zeros(9,num_pts);

for ii=1:num_pts
    pt_i = pts(:,ii);
    vel_i = vels(:,ii);
    z0(1:3,ii) = inverse_kinematics_init(p,pt_i); % get feasible ik solution for each point in trajectory
    %     z0(1:3,ii) = inverse_kinematics_perp(p,pt_i,vel_i);
    
    % given solution, get approximate joint velocities using jacobian pseudoinverse
    J_i = jacobian_tip(z0(1:3,ii),p);
    z0(4:6,ii) = pinv(J_i)*[vel_i;0];
end
%z0 = z0(:);

opti.set_initial(X,z0);

%% Setup Casadi and Ipopt Options
p_opts = struct('expand',true); % this speeds up ~x10
s_opts = struct('max_iter',3000,...
    'max_cpu_time',40.0,...
    'tol', 1e-4,... % (1e-6), 1e-4 works well
    'acceptable_tol', 1e-4,... % (1e-4)
    'constr_viol_tol', 1e-3,... % (1e-6), 1e3 works well
    'acceptable_iter', 5,... % (15), % 5 works well
    'nlp_scaling_method','gradient-based',... {'gradient-based','none','equilibration-based'};
    'nlp_scaling_max_gradient',50,... % (100), % 50 works well
    'bound_relax_factor', 1e-6,... % (1e-8), % 1e-6 works well
    'fixed_variable_treatment','relax_bounds',... % {'make_parameter','make_constraint','relax_bounds'}; % relax bounds works well
    'bound_frac',5e-1,... % (1e-2), 5e-1 works well
    'bound_push',5e-1,... % (1e-2), 5e-1 works well
    'mu_strategy','adaptive',... % {'monotone','adaptive'}; % adaptive works very well
    'mu_oracle','probing',... % {'quality-function','probing','loqo'}; % probing works very well
    'fixed_mu_oracle','probing',... % {'average_compl','quality-function','probing','loqo'}; % probing decent
    'adaptive_mu_globalization','obj-constr-filter',... % {'obj-constr-filter','kkt-error','never-monotone-mode'};
    'mu_init',1e-1,... % [1e-1 1e-2 1]
    'alpha_for_y','bound-mult',... % {'primal','bound-mult','min','max','full','min-dual-infeas','safer-min-dual-infeas','primal-and-full'}; % primal or bound-mult seems best
    'alpha_for_y_tol',1e1,... % (1e1)
    'recalc_y','no',... % {'no','yes'};
    'max_soc',4,... % (4)
    'accept_every_trial_step','no',... % {'no','yes'}
    'linear_solver','mumps',... % {'ma27','mumps','ma57','ma77','ma86'} % ma57 seems to work well
    'linear_system_scaling','slack-based',... {'mc19','none','slack-based'}; % Slack-based
    'linear_scaling_on_demand','yes',... % {'yes','no'};
    'max_refinement_steps',10,... % (10)
    'min_refinement_steps',1,... % (1)
    'warm_start_init_point', 'no'); % (no)

s_opts.file_print_level = 0;
s_opts.print_level = 3;
s_opts.print_frequency_iter = 5;
s_opts.print_timing_statistics ='no';
opti.solver('ipopt',p_opts,s_opts);

%% Solve with Opti Stack
disp('***********************');
disp('Solving with Opti Stack');
disp('***********************');
sol = opti.solve();
% try
%     %sol = opti.solve_limited();
%     sol = opti.solve();
% catch
%     disp('Did not solve - try debug solution')
% end
% TODO: setup this logic for catching error to handle infeasible solution

%% Decompose Solution
X_star = sol.value(X);
q_star(1:3,:) = sol.value(opt_var.q);
dq_star(1:3,:) = sol.value(opt_var.dq);
u_star(1:3,:) = sol.value(opt_var.u);

% calculate effective mass, effective momentum
meff_star = zeros(5,num_pts); % each column is [meff; meff_min; meff_max; lm_dir]
peff_star = zeros(3,num_pts); % px, py, |p|
p_star = zeros(3,num_pts); % xact, yact, error
v_star = zeros(3,num_pts); % vxact, vyact, verror
for ii=1:num_pts
    z_i = [q_star(:,ii); dq_star(:,ii)];
    % calculate endpoint velocity
    ve_temp = jacobian_tip(z_i,p)*dq_star(:,ii);
    ve_i = ve_temp(1:2); % get actual tip velocity
    ev_i = ve_i/norm(ve_i);
    % calculate effective mass, minimum, maximum
    LLv_inv = LLv_arm_op_inv(z_i,p);
    meff_star(1,ii) = 1/(ev_i'*LLv_inv*ev_i);
    [V,D] = eig(LLv_inv);
    if D(1,1)>D(2,2)
        meff_star(2,ii) = 1/D(1,1);
        meff_star(3,ii) = 1/D(2,2);
        meff_star(4:5,ii) = V(:,1);
    else
        meff_star(2,ii) = 1/D(2,2);
        meff_star(3,ii) = 1/D(1,1);
        meff_star(4:5,ii) = V(:,2);
    end
    % calculate effective momentum
    peff_star(1:2,ii) = ve_i*meff_star(1,ii);
    peff_star(3,ii) = norm(ve_i)*meff_star(1,ii);
    % get endpoint position
    p_temp = position_tip(z_i,p);
    p_star(1:2,ii) = p_temp(1:2);
    p_error = p_star(1:2,ii) - pts(:,ii);
    p_star(3,ii) = norm(p_error);
    % fill v_sol
    v_star(1:2,ii) = ve_i;
    v_error = v_star(1:2,ii) - vels(:,ii);
    v_star(3,ii) = norm(v_error);
end

% TODO: save TO data in struct

%% Animate Solution
% showmotion(model, time_vec, q_star);

%% Initial Plots
% same plots as before for optimization results
% TODO: eventually add forward simulation with block back in?

figure(2); clf;
subplot(3,1,1); hold on;
plot(time_vec, q_star); 
xlabel('Time'); ylabel('Joint Angle'); legend('q1','q2','q3');
subplot(3,1,2); hold on;
plot(time_vec, dq_star);
xlabel('Time'); ylabel('Joint Velocity'); legend('dq1','dq2','dq3');
subplot(3,1,3); hold on;
plot(time_vec, u_star);
xlabel('Time'); ylabel('Joint Torque'); legend('u1','u2','u3');

figure(3); clf; 
subplot(2,1,1); hold on;
plot(time_vec,meff_star(1,:),'o-','LineWidth',1.25);
plot(time_vec,meff_star(2:3,:),'LineWidth',1.25);
xlabel('Time'); ylabel('Effective Mass'); legend('Actual','Min','Max');
subplot(2,1,2); hold on;
plot(time_vec, peff_star(3,:),'LineWidth',1.25);
xlabel('Time'); ylabel('Effective Momentum');

figure(4); clf;
subplot(2,1,1); hold on;
plot(time_vec, p_star(3,:));
xlabel('Time'); ylabel('Endpoint Position Error');
subplot(2,1,2); hold on;
plot(time_vec, v_star(3,:));
xlabel('Time'); ylabel('Endpoint Velocity Error');

% animation of kinematics from optimization and dynamic simulation
% filename = 'TO_mod_cost_sine_curve_meff_direct.gif'; % save animation as a gif
figure(1);
% for ii=1:num_pts
for ii=1 % just generate the plot
    t_i = time_vec(ii);
    z_i = [q_star(:,ii); dq_star(:,ii)];
    plot_arm_kinematics(t_i,z_i,p,meff_star(:,ii),p_star,pts,vels(:,ii));
%     if ii==1
%         gif(filename,'DelayTime',dt);
%     else
%         gif;
%     end
    pause(dt);
end

%% Functions

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
    q1_ik = 0.0;
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
    
    % plot arena as well
    % TODO: don't have these dimensions hardcoded
    % TODO: have text showing traj_ind? include in title?
    x_a = [0.4, 1.2, 1.2, 0.4, 0.4];
    y_a = [-0.4, -0.4, 0.4, 0.4, -0.4];
    arena_pts = plot(x_a,y_a,'r--');
    
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





