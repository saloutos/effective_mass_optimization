% Andrew SaLoutos
% 3/8/21

% optimizing effective mass across a trajectory of points
% derived from "meff_opt_single_point.m" script


% trajectory of points, velocity directions are along lines between points (coarse approximation)

% in this script, the optimization runs on all of the trajectory points simulataneously

num_pts = 20;

%%% straight vertical line
pts_x = ones(1,num_pts);
pts_y = linspace(-0.4,0.4,num_pts);
vels_x = zeros(1,num_pts);
vels_y = ones(1,num_pts);

%%% circle
center = [0.9, 0.3]; % [0.7, -0.1] produced some weird behavior
radius = 0.2;
thetas = linspace(0,2*pi,num_pts+1);
thetas = thetas(1:num_pts);
pts_x   = center(1)*ones(size(thetas)) + radius*cos(thetas);
pts_y   = center(2)*ones(size(thetas)) + radius*sin(thetas);
vels_x = zeros(1,num_pts);
velx_y = zeros(1,num_pts);
for ii=1:num_pts
    if ii<num_pts
        vels_x(ii) = pts_x(ii+1)-pts_x(ii);
        vels_y(ii) = pts_y(ii+1)-pts_y(ii);
    else % if ii==num_pts
        vels_x(ii) = pts_x(1)-pts_x(ii);
        vels_y(ii) = pts_y(1)-pts_y(ii);
    end
end
    

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

% symbolic variables for optimization
syms q [3 num_pts] real % [q1_i; q2_i; q3_i] for each point, i
q = q(:); % [q1_1;q2_1;q3_1;...q1_N;q2_N;q3_N]

% initial values
q0 = zeros(3,num_pts);
% for ii=1:num_pts
%     eei = [pts_x(ii);pts_y(ii)];
%     %%%% TODO: run warm-start here instead of just naive ik?
%     q0(:,ii) = inverse_kinematics_init(p,eei);
% end
q0 = q0(:);

% fmincon inputs
cost_func = @(q)meff_cost(q,p,q0,pts_x,pts_y,vels_x,vels_y);
nonlcon = @(q)joint_constraints(q,p,q0,pts_x,pts_y,vels_x,vels_y);
A = []; % No other constraints
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
options = optimoptions('fmincon','display','none');

% run optimization
tic;
[q_sol, fval_sol, sol_hist] = fmincon_with_hist(cost_func,q0,A,b,Aeq,beq,lb,ub,nonlcon,options);
opt_time = toc;
opt_time_avg = opt_time/num_pts;
q_sol = reshape(q_sol,3,num_pts);

% loop through and run optimization for each point
pt_vec = zeros(2,num_pts);
ev_vec = zeros(2,num_pts);
cost_vec = zeros(1,num_pts);
meff_vec = zeros(1,num_pts);
meffmin_vec = zeros(1,num_pts);
qdiff_vec = zeros(1,num_pts);
figure; % open new figure

% save animation as a gif
% filename = 'circle_traj_gif_2.gif';

for ii=1:num_pts
       
   % store point and velocity
   pt_vec(:,ii) = [pts_x(ii); pts_y(ii)];
   vel = [vels_x(ii); vels_y(ii)];
   ev_vec(:,ii) = vel/norm(vel);
   
   % get joint angles for this point from optimization solution
   if ii==1
       q_vals = [q_sol(:,ii), q_sol(:,ii)];
   else
       q_vals = [q_sol(:,ii-1), q_sol(:,ii)];
   end
   
   qdiff = q_vals(:,2)-q_vals(:,1);
   qdiff_vec(ii) = qdiff'*qdiff;
   
   % animation at each point
   clf; hold on; 
   plot(pts_x,pts_y,'Xg','MarkerSize',8);   
   [meff_temp, meffmin_temp] = plot_ik_sol_only(q_vals,p,ev_vec(:,ii),pt_vec(:,ii));
   
%    if ii==1
%        gif(filename,'DelayTime',0.25);
%    else
%        gif;
%    end
   
   pause(0.5);
   
   % save effective mass
   meff_vec(ii) = meff_temp;
   meffmin_vec(ii) = meffmin_temp;
   
end

% also plot effective mass
figure; subplot(2,1,1);
plot(meff_vec);
xlabel('Point #'); ylabel('Effective Mass');
subplot(2,1,2);
plot(qdiff_vec);
xlabel('Point #'); ylabel("qdiff'*qdiff");

% display opt time
fprintf('\nTotal optimization time: %.1f ms\n\n',opt_time*1000.0);
fprintf('Optimization time per point in trajectory: %.1f ms\n\n',opt_time_avg*1000.0);

%%%%% helper functions %%%%%
function [c,ceq] = joint_constraints(q,p,q0,pts_x,pts_y,vels_x,vels_y)
    
    num_pts = length(pts_x);
    q = reshape(q,3,num_pts); % columns are joint angles for each point
    
    c = [];
    ceq = [];

    % could also calculate gradients of constraints here?

    % dynamics? using helper functions -> A*ddq = b -> x[i+1] = x[i] + dt*dx/dt[i]
    
    % add constraint for maximum changes in joint angles? (essentially joint velocity limits)
        
    
    % add end effector position constraint for each point in trajectory
    for ii=1:num_pts
       eei = [pts_x(ii);pts_y(ii)];
       zi = [q(:,ii);zeros(3,1)];
       ee_temp = position_tip(zi,p);
       ceqi = ee_temp(1:2)-eei;
       ceq = [ceq; ceqi];        
    end
    
end

function total_cost = meff_cost(q,p,q0,pts_x,pts_y,vels_x,vels_y)
    
    num_pts = length(pts_x);
    q = reshape(q,3,num_pts); % columns are joint angles for each point
    
    meff_cost_vec = zeros(num_pts,1);
    dq_cost_vec = zeros(num_pts,1);
    
    % could also calculate gradient and hessians of costs here?
    
    % running cost on meff along vdir
    for ii=1:num_pts
        zi = [q(:,ii); zeros(3,1)];
        veli = [vels_x(ii);vels_y(ii)];
        evi = veli/norm(veli);
        
        LLv_inv = LLv_arm_op_inv(zi,p);
%         meff_cost_vec(ii) = 1/(evi'*LLv_inv*evi); 
        meff_cost_vec(ii) = -(evi'*LLv_inv*evi); % negative might work better?
    end
    
    % running cost on dq between positions?
    Q = eye(3);
    for ii=2:num_pts
        qdiff = q(:,ii)-q(:,ii-1);
        dq_cost_vec(ii) = qdiff'*Q*qdiff;
    end
    
    % add costs together with weights
    a1 = 3*ones(1,num_pts); 
%     a1(1) = 3; % higher weight on effective mass at the start of the trajectory
    c1 = a1*meff_cost_vec;
    
    a2 = 0.05*ones(1,num_pts); % constant weights on dq costs
%     a2(1) = 10; % higher weight on first joint velocity
    c2 = a2*dq_cost_vec;

    % return total cost
    total_cost = c1 + c2;
    
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

% inverse kinematics for initialization
function q_ik = inverse_kinematics_init(p,rE)
    % pull out necessary parameters
    l1 = p(14);
    l2 = p(15);
    l3 = p(16);

    % pull out desired end-effector position
    % note: will have to do interpolation eventually here
    x_ee = rE(1);
    y_ee = rE(2);

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

% not sure if this is the plot function I want
% for a single point, plot the solution
function [meff1f, meffminf] = plot_ik_sol_only(q_vals,p,vdir,rE)

    % doesn't create a new figure/plot, need to call subplot/figure before and also add title after

    % should take initial posture and final posture (maybe intermediate?)
    % plot initial posture and ellipsoid in gray
    % plot final posture in black, with ellipsoid, effective mass directions, and velocity vector

    % q_vals should be [q_prev, q_cur]; where each q is a 3x1 column vector of joint angles
    num_pts = size(q_vals,2);

    %%% plot arm kinematics

    %%% keypoints and lmdir for initial pose
    kp0 = keypoints_arm(q_vals(:,1),p);
    rA0 = kp0(:,1);
    rB0 = kp0(:,2);
    rC0 = kp0(:,3);
    % calculate ellipse for operational space inertia
    z_ik0 = [q_vals(:,1); zeros(3,1)];
    LLv_inv0 = LLv_arm_op_inv(z_ik0,p); % just care about x and y for now
    [V0,D0] = eig(LLv_inv0);
    meff10 = 1/(vdir'*LLv_inv0*vdir); % effective mass?
    % plot direction of lowest apparent inertia
    if D0(1,1) > D0(2,2)
        lm_dir0 = V0(:,1); % should already have unit magnitude
    else
        lm_dir0 = V0(:,2);         
    end
    % using parametric equations of ellipse, calculate ellipse points relative to foot position
    th_ellipse0 = atan2(V0(2,1),V0(1,1)); % angle between first eigenvector and positive x axis
    gamma_ell0 = 0.1; % TODO: better way to implement scaling of the ellipse?
    l_x0 = gamma_ell0*sqrt((1/D0(1,1))); 
    l_y0 = gamma_ell0*sqrt((1/D0(2,2))); 
    jj = linspace(0, 2*pi, 100);
    % make this belted ellipsoid?
    x_ell0 = (l_x0*cos(jj))*cos(th_ellipse0) - (l_y0*sin(jj))*sin(th_ellipse0);
    y_ell0 = (l_x0*cos(jj))*sin(th_ellipse0) + (l_y0*sin(jj))*cos(th_ellipse0);


    %%% keypoints and lmdir for final pose
    kpf = keypoints_arm(q_vals(:,num_pts),p);
    rAf = kpf(:,1);
    rBf = kpf(:,2);
    rCf = kpf(:,3);
    % calculate ellipse for operational space inertia
    z_ikf = [q_vals(:,num_pts); zeros(3,1)];
    LLv_invf = LLv_arm_op_inv(z_ikf,p); % just care about x and y for now
    [Vf,Df] = eig(LLv_invf);
    meff1f = 1/(vdir'*LLv_invf*vdir); % effective mass?
    % plot direction of lowest apparent inertia
    if Df(1,1) > Df(2,2)
        lm_dirf = Vf(:,1); % should already have unit magnitude
        meffminf = 1/Df(1,1);
    else
        lm_dirf = Vf(:,2);
        meffminf = 1/Df(2,2);
    end
    % using parametric equations of ellipse, calculate ellipse points relative to foot position
    th_ellipsef = atan2(Vf(2,1),Vf(1,1)); % angle between first eigenvector and positive x axis
    gamma_ellf = 0.1; % TODO: better way to implement scaling of the ellipse?
    l_xf = gamma_ellf*sqrt((1/Df(1,1))); 
    l_yf = gamma_ellf*sqrt((1/Df(2,2))); 
    jj = linspace(0, 2*pi, 100);
    % make this belted ellipsoid?
    x_ellf = (l_xf*cos(jj))*cos(th_ellipsef) - (l_yf*sin(jj))*sin(th_ellipsef);
    y_ellf = (l_xf*cos(jj))*sin(th_ellipsef) + (l_yf*sin(jj))*cos(th_ellipsef);


    %%% plot these values    
    des_pt = plot(rE(1), rE(2), 'Xg','MarkerSize',6);
    gray_color = [0.5, 0.5, 0.5];
    h_OA_0 = plot([0 rA0(1)],[0 rA0(2)],'LineWidth',2,'Color',gray_color);
    h_AB_0 = plot([rA0(1) rB0(1)],[rA0(2) rB0(2)],'LineWidth',2,'Color',gray_color);
    h_BC_0 = plot([rB0(1) rC0(1)],[rB0(2) rC0(2)],'LineWidth',2,'Color',gray_color);
    h_0_dots = plot([0 rA0(1) rB0(1)],[0 rA0(2) rB0(2)],'o','MarkerSize',5,'Color',gray_color,'MarkerFaceColor',gray_color);
    h_LL_0 = plot([rC0(1)+x_ell0],[rC0(2)+y_ell0],'Color',gray_color*0.7); % slightly darker?

    ell_color = [0.85, 0.33, 0.];
    h_OA_f = plot([0 rAf(1)],[0 rAf(2)],'k','LineWidth',2);
    h_AB_f = plot([rAf(1) rBf(1)],[rAf(2) rBf(2)],'k','LineWidth',2);
    h_BC_f = plot([rBf(1) rCf(1)],[rBf(2) rCf(2)],'k','LineWidth',2);
    h_f_dots = plot([0 rAf(1) rBf(1)],[0 rAf(2) rBf(2)],'ok','MarkerSize',5,'MarkerFaceColor','k');
    h_LL_f = plot([rCf(1)+x_ellf],[rCf(2)+y_ellf],'LineWidth',1.5,'Color',ell_color);

    h_lm1 = quiver([rE(1)],[rE(2)],[lm_dirf(1)],[lm_dirf(2)],0.25,'r','LineWidth',2);
    h_lm2 = quiver([rE(1)],[rE(2)],[-lm_dirf(1)],[-lm_dirf(2)],0.25,'r','LineWidth',2);
    h_v = quiver([rE(1)],[rE(2)],[vdir(1)],[vdir(2)],0.25,'b','LineWidth',2);
    xlabel('X'); ylabel('Y'); axis equal; 
    xlim([-0.5,1.5]); ylim([-1.0,1.0]);

    h_meff1f = text(0.8,0.9,'','FontSize',8);
    h_meffminf = text(0.8,0.75,'','FontSize',8);
    set(h_meff1f,'String',  sprintf('m_{eff,1} = %.2f',meff1f));
    set(h_meffminf,'String',  sprintf('m_{eff,min} = %.2f',meffminf));

end


