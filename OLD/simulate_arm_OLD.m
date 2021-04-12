function simulate_arm()

    %% Setup
    % fixed parameters    
    m1 = 1;                 m2 = 1;
    m3 = 1;                 m_motor = 0.5;
    I1 = 0.2;               I2 = 0.2;
    I3 = 0.2;               I_motor = 0.05;
    Ir = 0.01;              N = 6;
    l_O_m1 = 0.25;          l_A_m2 = 0.25;
    l_B_m3 = 0.25;          l_OA = 0.5;
    l_AB = 0.5;             l_BC = 0.5;
    g = 9.81; % do I want gravity to start?
    
    % arrange in vector
    p   = [m1 m2 m3 m_motor I1 I2 I3 I_motor Ir N l_O_m1 l_A_m2 l_B_m3 l_OA l_AB l_BC g]';
    
    % sim parameters...circular trajectory
    p_traj.thetas = linspace(0,2*pi,61);
    p_traj.thetas = p_traj.thetas(1:60);
    
    
    p_traj.center = [0.7, -0.1]; %[1.25,0];
    v_dir = [0.4;1]; % specify velocity direction for brute-force method
    q1_lim = [-pi/3, pi/3]; % limits for search on q1 (first link angle)
    
%     rE = [1.0; 0.1]; % [x;y];
%     vdir = [1.0; -1.0]; % [vx; vy];
    
    
    p_traj.r = 0.2;
    p_traj.x   = p_traj.center(1)*ones(size(p_traj.thetas)); % + p_traj.r*cos(p_traj.thetas);
    p_traj.y   = p_traj.center(2)*ones(size(p_traj.thetas)); % + p_traj.r*sin(p_traj.thetas);
    
    
    
    %% Dynamic simulation
    % not dynamic yet
    
    tspan = [0 10];
    z0 = [0,0,0,0,0,0];
%     opts = odeset('AbsTol',1e-8,'RelTol',1e-6);
%     sol = ode45(@dynamics,tspan,z0,opts,p,p_traj);
    % ...use simpler euler integration here?
    
    % solve inverse kinematics
%     q_ik = inverse_kinematics(0,0,p,p_traj);

    q_ik = inverse_kinematics_redund(0,0,p,p_traj,q1_lim); % for a single point, investigate redundancies
    
    % evaluate operational space inertia matrix here?
    
    %% Compute tip position over time
    % now just plotting ik solutions
    
    % get keypoints for plotting
%     kp = [zeros(2,1), keypoints_arm(q_ik(:,1),p)];
%     figure(1); clf; hold on; axis equal;
%     plot(kp(1,:), kp(2,:),'k','LineWidth',1.5);
%     plot(kp(1,:), kp(2,:),'or','MarkerSize',10,'MarkerFaceColor','r');
%     plot(p_traj.x, p_traj.y, 'Xg','MarkerSize',10,'LineWidth',2);
%     
%     xlabel('X'); ylabel('Y'); xlim([-0.5,2]); ylim([-1.5,1.5])
    
    %% Animate Solution
    
%     animate_ik_LL_v_sol(q_ik,p,p_traj);
    animate_ik_LL_vin_sol(q_ik,v_dir,p,p_traj);
%     figure(3); clf;
%     hold on    
%     % target position Q1.4
%     plot(.025 , -.125, 'r.','MarkerSize',6); 
%     % Target traj. Q 1.6
%     TH = 0:.1:2*pi;
%     plot( p_traj.x_0 + p_traj.r * cos(TH), ...
%           p_traj.y_0 + p_traj.r * sin(TH),'k--'); 
%     
%     animateSol(sol,p);
    
end

function q_ik = inverse_kinematics(t,z0,p,p_traj)
    % pull out necessary parameters
    l1 = p(14);
    l2 = p(15);
    l3 = p(16);
    
    % pull out desired end-effector position
    % note: will have to do interpolation eventually here
    x_ee = p_traj.x;
    y_ee = p_traj.y;
    
    % assume q1 = 0
    % calculate q2 and q3
    q1 = zeros(1,length(p_traj.thetas));
    x_A = l1;
    y_A = 0;
    
    x_ee = x_ee-x_A; % distances taken from second joint
    y_ee = y_ee-y_A;
    
    q3 = acos( (x_ee.^2 + y_ee.^2 - l2^2 - l3^2) / (2 * l2 * l3) ); % this can be +/-, leave as just + for now
    q2 = atan2(y_ee,x_ee) - atan2( (l3*sin(q3)), (l2+(l3*cos(q3))) );
    
    % outputs
    q_ik = [q1;q2;q3];
end

function q_ik = inverse_kinematics_redund(t,z0,p,p_traj,q1_lim)
    % pull out necessary parameters
    l1 = p(14);
    l2 = p(15);
    l3 = p(16);
    
    % pull out desired end-effector position
    % note: will have to do interpolation eventually here
    x_ee = p_traj.x;
    y_ee = p_traj.y;
    
    % assume q1 = 0
    % calculate q2 and q3
    q1 = linspace(q1_lim(1), q1_lim(2), length(p_traj.thetas));
    x_A = l1*cos(q1);
    y_A = l1*sin(q1);
    
    x_ee = x_ee-x_A; % distances taken from second joint
    y_ee = y_ee-y_A;
    
    q3 = acos( (x_ee.^2 + y_ee.^2 - l2^2 - l3^2) / (2 * l2 * l3) ); % this can be +/-, leave as just + for now
    q2 = atan2(y_ee,x_ee) - atan2( (l3*sin(q3)), (l2+(l3*cos(q3))) ) - q1;
    
    % outputs
    q_ik = [q1;q2;q3];
end

function tau = control_law(t,z,p,p_traj)
    % Controller gains, Update as necessary for Problem 1
    K_x = 40; % Spring stiffness X
    K_y = 40; % Spring stiffness Y
    D_x = 4;  % Damping X
    D_y = 4;  % Damping Y
    
    % Desired position of foot is a circle
    rEd = [p_traj.x_0 p_traj.y_0 0]' + ...
            p_traj.r*[cos(p_traj.omega*t) sin(p_traj.omega*t) 0]';
    % Compute desired velocity of foot
    vEd = p_traj.r*[-sin(p_traj.omega*t)*p_traj.omega    ...
                     cos(p_traj.omega*t)*p_traj.omega   0]';
    % Desired acceleration
    aEd = p_traj.r*[-cos(p_traj.omega*t)*p_traj.omega^2 ...
                    -sin(p_traj.omega*t)*p_traj.omega^2 0]';
    
    % Actual position and velocity 
    rE = position_foot(z,p);
    vE = velocity_foot(z,p);
    
    % Jacobian matrix \partial r_E / \partial q
    J  = jacobian_foot(z,p);
    dJ = jacobian_dot_foot(z,p);
    dq = z(3:4);
    
    % Compute virtual foce for Question 1.4 and 1.5
    f  = [K_x * (rEd(1) - rE(1) ) - D_x * (vE(1) - vEd(1) ) ;
          K_y * (rEd(2) - rE(2) ) - D_y * (vE(2) - vEd(2) ) ];
    
    %% Task-space compensation and feed forward for Question 1.8
    % Get joint space components of equations of motion
    Mass_Joint_Sp = A_leg(z,p);
    Grav_Joint_Sp = Grav_leg(z,p);
    Corr_Joint_Sp = Corr_leg(z,p);

    % Task-space mass matrix (Equaiton 51 in Khatib's paper)
    Jinv   = inv(J);
    Lambda = Jinv' * Mass_Joint_Sp * Jinv;
    
    % Coriolis force in task-space (Equation 51)
    mu     = Jinv' * Corr_Joint_Sp - Lambda * dJ * dq;
    
    % Gravity force in task-space (Equation 51)
    rho    = Jinv' * Grav_Joint_Sp; 
    
    % Add task-space acceleration force feedforward, coriolis, and gravity compensation 
    % (Equation 50)
    f(1:2) = f(1:2) + Lambda*aEd(1:2) + mu + rho;
    
    % Map to joint torques  
    tau = J' * f;
end

function dz = dynamics(t,z,p,p_traj)
    th1 = z(1);     th2 = z(2);
    dth1= z(3);     dth2= z(4);
    
    % Get mass matrix
    A = A_leg(z,p);
    
    % Compute Controls
    tau = control_law(t,z,p,p_traj);
    
    % Get b = Q - V(q,qd) - G(q)
    b = b_leg(z,tau,p);
  
    % Compute the contact force (used for problem 2)
    % Fc = contact_force(z,p);
    Fc = 0;
    
%     Tauc = joint_limit_torque(z,p);
    Tauc=0;
    
    % Compute the contribution of the contact force to the generalied force
    J = jacobian_foot(z,p);
    QFc=J'*[0 Fc]';
    
    % Contribution of constraint torque is more simple.
    QTauc= [0 Tauc]';
    
    % Solve for qdd.
    qdd = A\(b+QFc+QTauc);
    dz = 0*z;
    
    % Form dz
    dz(1:2) = z(3:4);
    dz(3:4) = qdd;
end

function animateSol(sol,p)
    % Prepare plot handles
    hold on
    h_OB = plot([0],[0],'LineWidth',2);
    h_AC = plot([0],[0],'LineWidth',2);
    h_BD = plot([0],[0],'LineWidth',2);
    h_CE = plot([0],[0],'LineWidth',2);
   
    
    xlabel('x'); ylabel('y');
    h_title = title('t=0.0s');
    
    axis equal
    axis([-.2 .2 -.3 .1]);

    %Step through and update animation
    for t = 0:.01:sol.x(end)
        % interpolate to get state at current time.
        z = interp1(sol.x',sol.y',t)';
        keypoints = keypoints_leg(z,p);

        rA = keypoints(:,1); % Vector to base of cart
        rB = keypoints(:,2);
        rC = keypoints(:,3); % Vector to tip of pendulum
        rD = keypoints(:,4);
        rE = keypoints(:,5);

        set(h_title,'String',  sprintf('t=%.2f',t) ); % update title
        
        set(h_OB,'XData',[0 rB(1)]);
        set(h_OB,'YData',[0 rB(2)]);
        
        set(h_AC,'XData',[rA(1) rC(1)]);
        set(h_AC,'YData',[rA(2) rC(2)]);
        
        set(h_BD,'XData',[rB(1) rD(1)]);
        set(h_BD,'YData',[rB(2) rD(2)]);
        
        set(h_CE,'XData',[rC(1) rE(1)]);
        set(h_CE,'YData',[rC(2) rE(2)]);

        pause(.01)
    end
end

% animate ik solutions with operational space inertia matrix and cartesian
% velocity (approximation)
function animate_ik_LL_vin_sol(q_ik,v_dir,p,p_traj)

    figure(1); clf; 
    
    subplot(1,2,1); hold on;
    % Prepare plot handles
    des_pts = plot(p_traj.x, p_traj.y, 'Xg','MarkerSize',6);
    h_OA = plot([0],[0],'LineWidth',2);
    h_AB = plot([0],[0],'LineWidth',2);
    h_BC = plot([0],[0],'LineWidth',2);
    h_LL = plot([0],[0]);
    h_lm1 = quiver([0],[0],[0],[0],0.25,'r','LineWidth',2);
    h_lm2 = quiver([0],[0],[0],[0],0.25,'r','LineWidth',2);
    h_v = quiver([0],[0],[0],[0],0.25,'b','LineWidth',2);
    xlabel('X'); ylabel('Y'); xlim([-0.5,2]); ylim([-1.5,1.5]);
    axis equal;
    h_title = title('Point 1');
    h_meff1 = text(-0.25,-0.5,'','FontSize',12);
    h_meffmin = text(-0.25,-0.75,'','FontSize',12);
    h_meffmax = text(-0.25,-1,'','FontSize',12);
    
    subplot(1,2,2); hold on;
    h_m1 = plot([0],[0],'LineWidth',1.5);
    h_mm = plot([0],[0],'LineWidth',1.5);
    h_mx = plot([0],[0],'LineWidth',1.5);
    legend('m_{eff,v}','m_{eff,min}','m_{eff,max}');
    %Step through points and update animation
    num_pts = length(p_traj.thetas);
    xlim([0,num_pts]); ylim([0,10]);
    
    m1_vec = [];
    mm_vec = [];
    mx_vec = [];
    pt_vec = [];
    
    h_bestm = text(5,7.5,'','FontSize',12);
    h_bestmi = text(5,6.5,'','FontSize',12);
    
    
    for ii=1:num_pts
        
        % get keypoints for plotting
        kp = keypoints_arm(q_ik(:,ii),p);
                
        rA = kp(:,1);
        rB = kp(:,2);
        rC = kp(:,3);
        
        % tip vel desired (direction based on next point in trajectory)
        pt_cur = [p_traj.x(ii), p_traj.y(ii)];
%         if ii<num_pts
%             pt_next = [p_traj.x(ii+1), p_traj.y(ii+1)];
%         else
%             pt_next = [p_traj.x(1), p_traj.y(1)];
%         end
%         rough_vel = pt_next - pt_cur;
%         v_norm = rough_vel'/norm(rough_vel); % unit vector along vel direction
        
        v_norm = v_dir/norm(v_dir); % make sure input velocity direction is a unit vector
        
        % calculate ellipse for operational space inertia
        z_ik = [q_ik(:,ii); zeros(3,1)];
        LLv_inv = LLv_arm_op_inv(z_ik,p); % just care about x and y for now
        [V,D] = eig(LLv_inv);
        meff1 = 1/(v_norm'*LLv_inv*v_norm); % effective mass?
        
        % plot direction of lowest apparent inertia
        if D(1,1) > D(2,2)
            lm_dir = V(:,1); % should already have unit magnitude
            meffmin = 1/D(1,1);
            meffmax = 1/D(2,2);
        else
            lm_dir = V(:,2);
            meffmin = 1/D(2,2);
            meffmax = 1/D(1,1);            
        end
        
        pt_vec = [pt_vec, ii];
        m1_vec = [m1_vec, meff1];
        mm_vec = [mm_vec, meffmin];
        mx_vec = [mx_vec, meffmax];
        
        % using parametric equations of ellipse, calculate ellipse points relative to foot position
        th_ellipse = atan2(V(2,1),V(1,1)); % angle between first eigenvector and positive x axis
        gamma_ell = 0.1; % TODO: better way to implement scaling of the ellipse?
        l_x = gamma_ell*sqrt((1/D(1,1))); 
        l_y = gamma_ell*sqrt((1/D(2,2))); 
        jj = linspace(0, 2*pi, 100);
        % make this belted ellipsoid?
        x_ell = (l_x*cos(jj))*cos(th_ellipse) - (l_y*sin(jj))*sin(th_ellipse);
        y_ell = (l_x*cos(jj))*sin(th_ellipse) + (l_y*sin(jj))*cos(th_ellipse);
        
        
        set(h_title,'String',  sprintf('Point %d',ii) ); % update title
        
        set(h_meff1,'String',  sprintf('m_{eff,1} = %.2f',meff1));
        set(h_meffmin,'String',  sprintf('m_{eff,min} = %.2f',meffmin));
        set(h_meffmax,'String',  sprintf('m_{eff,max} = %.2f',meffmax));
        
        set(h_OA,'XData',[0 rA(1)]);
        set(h_OA,'YData',[0 rA(2)]);
        
        set(h_AB,'XData',[rA(1) rB(1)]);
        set(h_AB,'YData',[rA(2) rB(2)]);
        
        set(h_BC,'XData',[rB(1) rC(1)]);
        set(h_BC,'YData',[rB(2) rC(2)]);

        set(h_LL,'XData',rC(1)+x_ell);
        set(h_LL,'YData',rC(2)+y_ell);
        
        set(h_v,'XData',pt_cur(1));
        set(h_v,'YData',pt_cur(2));
        set(h_v,'UData',v_norm(1));
        set(h_v,'VData',v_norm(2));  
        
        set(h_lm1,'XData',pt_cur(1));
        set(h_lm1,'YData',pt_cur(2));
        set(h_lm1,'UData',lm_dir(1));
        set(h_lm1,'VData',lm_dir(2));
        set(h_lm2,'XData',pt_cur(1));
        set(h_lm2,'YData',pt_cur(2));
        set(h_lm2,'UData',-lm_dir(1));
        set(h_lm2,'VData',-lm_dir(2));
        subplot(1,2,1); xlim([-0.5,2]); ylim([-1.5,1.5]);
        
        set(h_m1,'XData',pt_vec);
        set(h_m1,'YData',m1_vec);
        set(h_mm,'XData',pt_vec);
        set(h_mm,'YData',mm_vec);
        set(h_mx,'XData',pt_vec);
        set(h_mx,'YData',mx_vec);
        subplot(1,2,2); xlim([0,num_pts]); ylim([0,10]);
        
        [bestm, bestm_ind] = min(m1_vec);
        
        bestm_mmin = mm_vec(bestm_ind);
        
        set(h_bestm,'String',  sprintf('Best m_{eff,v} = %.2f, minumum m_{eff}: %.2f',bestm,bestm_mmin));
        set(h_bestmi,'String',  sprintf('Config. index: %d',bestm_ind));
        
        pause(.1)
    end
end

% animate ik solutions with operational space inertia matrix and cartesian
% velocity (approximation)
function animate_ik_LL_v_sol(q_ik,p,p_traj)

    figure(1); clf; 
    
    subplot(1,2,1); hold on;
    % Prepare plot handles
    des_pts = plot(p_traj.x, p_traj.y, 'Xg','MarkerSize',6);
    h_OA = plot([0],[0],'LineWidth',2);
    h_AB = plot([0],[0],'LineWidth',2);
    h_BC = plot([0],[0],'LineWidth',2);
    h_LL = plot([0],[0]);
    h_lm1 = quiver([0],[0],[0],[0],0.25,'r','LineWidth',2);
    h_lm2 = quiver([0],[0],[0],[0],0.25,'r','LineWidth',2);
    h_v = quiver([0],[0],[0],[0],0.25,'b','LineWidth',2);
    xlabel('X'); ylabel('Y'); xlim([-0.5,2]); ylim([-1.5,1.5]);
    axis equal;
    h_title = title('Point 1');
    h_meff1 = text(-0.25,-0.5,'','FontSize',12);
    h_meffmin = text(-0.25,-0.75,'','FontSize',12);
    h_meffmax = text(-0.25,-1,'','FontSize',12);
    
    subplot(1,2,2); hold on;
    h_m1 = plot([0],[0],'LineWidth',1.5);
    h_mm = plot([0],[0],'LineWidth',1.5);
    h_mx = plot([0],[0],'LineWidth',1.5);
    legend('m_{eff,v}','m_{eff,min}','m_{eff,max}');
    %Step through points and update animation
    num_pts = length(p_traj.thetas);
    xlim([0,num_pts]); ylim([0,10]);
    
    m1_vec = [];
    mm_vec = [];
    mx_vec = [];
    pt_vec = [];
    
    for ii=1:num_pts
        
        % get keypoints for plotting
        kp = keypoints_arm(q_ik(:,ii),p);
                
        rA = kp(:,1);
        rB = kp(:,2);
        rC = kp(:,3);
        
        % tip vel desired (direction based on next point in trajectory)
        pt_cur = [p_traj.x(ii), p_traj.y(ii)];
        if ii<num_pts
            pt_next = [p_traj.x(ii+1), p_traj.y(ii+1)];
        else
            pt_next = [p_traj.x(1), p_traj.y(1)];
        end
        rough_vel = pt_next - pt_cur;
        v_norm = rough_vel'/norm(rough_vel); % unit vector along vel direction
        
        % calculate ellipse for operational space inertia
        z_ik = [q_ik(:,ii); zeros(3,1)];
        LLv_inv = LLv_arm_op_inv(z_ik,p); % just care about x and y for now
        [V,D] = eig(LLv_inv);
        meff1 = 1/(v_norm'*LLv_inv*v_norm); % effective mass?
        
        % plot direction of lowest apparent inertia
        if D(1,1) > D(2,2)
            lm_dir = V(:,1); % should already have unit magnitude
            meffmin = 1/D(1,1);
            meffmax = 1/D(2,2);
        else
            lm_dir = V(:,2);
            meffmin = 1/D(2,2);
            meffmax = 1/D(1,1);            
        end
        
        pt_vec = [pt_vec, ii];
        m1_vec = [m1_vec, meff1];
        mm_vec = [mm_vec, meffmin];
        mx_vec = [mx_vec, meffmax];
        
        % using parametric equations of ellipse, calculate ellipse points relative to foot position
        th_ellipse = atan2(V(2,1),V(1,1)); % angle between first eigenvector and positive x axis
        gamma_ell = 0.1; % TODO: better way to implement scaling of the ellipse?
        l_x = gamma_ell*sqrt((1/D(1,1))); 
        l_y = gamma_ell*sqrt((1/D(2,2))); 
        jj = linspace(0, 2*pi, 100);
        % make this belted ellipsoid?
        x_ell = (l_x*cos(jj))*cos(th_ellipse) - (l_y*sin(jj))*sin(th_ellipse);
        y_ell = (l_x*cos(jj))*sin(th_ellipse) + (l_y*sin(jj))*cos(th_ellipse);
        
        
        set(h_title,'String',  sprintf('Point %d',ii) ); % update title
        
        set(h_meff1,'String',  sprintf('m_{eff,1} = %.2f',meff1));
        set(h_meffmin,'String',  sprintf('m_{eff,min} = %.2f',meffmin));
        set(h_meffmax,'String',  sprintf('m_{eff,max} = %.2f',meffmax));
        
        set(h_OA,'XData',[0 rA(1)]);
        set(h_OA,'YData',[0 rA(2)]);
        
        set(h_AB,'XData',[rA(1) rB(1)]);
        set(h_AB,'YData',[rA(2) rB(2)]);
        
        set(h_BC,'XData',[rB(1) rC(1)]);
        set(h_BC,'YData',[rB(2) rC(2)]);

        set(h_LL,'XData',rC(1)+x_ell);
        set(h_LL,'YData',rC(2)+y_ell);
        
        set(h_v,'XData',pt_cur(1));
        set(h_v,'YData',pt_cur(2));
        set(h_v,'UData',v_norm(1));
        set(h_v,'VData',v_norm(2));  
        
        set(h_lm1,'XData',pt_cur(1));
        set(h_lm1,'YData',pt_cur(2));
        set(h_lm1,'UData',lm_dir(1));
        set(h_lm1,'VData',lm_dir(2));
        set(h_lm2,'XData',pt_cur(1));
        set(h_lm2,'YData',pt_cur(2));
        set(h_lm2,'UData',-lm_dir(1));
        set(h_lm2,'VData',-lm_dir(2));
        subplot(1,2,1); xlim([-0.5,2]); ylim([-1.5,1.5]);
        
        set(h_m1,'XData',pt_vec);
        set(h_m1,'YData',m1_vec);
        set(h_mm,'XData',pt_vec);
        set(h_mm,'YData',mm_vec);
        set(h_mx,'XData',pt_vec);
        set(h_mx,'YData',mx_vec);
        subplot(1,2,2); xlim([0,num_pts]); ylim([0,10]);
        
        pause(.1)
    end
end

% Animate through ik solutions
function animate_ik_sol(q_ik,p,p_traj)

    figure(1); hold on;
    % Prepare plot handles
    des_pts = plot(p_traj.x, p_traj.y, 'Xg','MarkerSize',10,'LineWidth',2);
    h_OA = plot([0],[0],'LineWidth',2);
    h_AB = plot([0],[0],'LineWidth',2);
    h_BC = plot([0],[0],'LineWidth',2);
    xlabel('X'); ylabel('Y'); xlim([-0.5,2]); ylim([-1.5,1.5]);
    axis equal;
    h_title = title('Point 1');

    %Step through points and update animation
    for ii=1:length(p_traj.thetas)
        
        % get keypoints for plotting
        kp = keypoints_arm(q_ik(:,ii),p);
                
        rA = kp(:,1);
        rB = kp(:,2);
        rC = kp(:,3);

        set(h_title,'String',  sprintf('Point %d',ii) ); % update title
        
        set(h_OA,'XData',[0 rA(1)]);
        set(h_OA,'YData',[0 rA(2)]);
        
        set(h_AB,'XData',[rA(1) rB(1)]);
        set(h_AB,'YData',[rA(2) rB(2)]);
        
        set(h_BC,'XData',[rB(1) rC(1)]);
        set(h_BC,'YData',[rB(2) rC(2)]);

        pause(0.25)
    end
end