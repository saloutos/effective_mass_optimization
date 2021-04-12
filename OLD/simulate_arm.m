function sim_data = simulate_arm(sim_dt,tspan_TO,z_sol,p)

    % Andrew SaLoutos
    % 3/24/21

    % dynamic simulation for a 3DOF arm, following a joint position and
    % velocity trajectory with feed-forward inputs

    % INPUTS:
    % z_sol = [q_sol(i); dq_sol(i); u_sol(i)]; from trajectory optimization
    % p   = [m1 m2 m3 m_motor I1 I2 I3 I_motor Ir N l_O_m1 l_A_m2 l_B_m3 l_OA l_AB l_BC g]'; % fixed parameters for dynamic model
    
    %% Perform dynamic simulation
    % sort out different time vectors/values
    Tf = tspan_TO(end);
    num_TO_steps = length(tspan_TO);  
    num_sim_steps = (Tf/sim_dt)+1;
    tspan_sim = 0:sim_dt:Tf;
    
    % extract reference state trajectory and input trajectory
    z_TO = z_sol(1:6,:);
    u_TO = z_sol(7:9,:);    
    
    % interpolate state trajectory for feedforward
    % TODO: could look at doing more than linear interp, like cubic spline?
    z_TO_ff = interp1(tspan_TO',z_TO',tspan_sim')';
    % zero-order hold input trajectory
    u_TO_ff = interp1(tspan_TO',u_TO',tspan_sim','previous')';
    
    u_sim = zeros(3,length(tspan_sim));
    z_sim = zeros(6,length(tspan_sim)); %[q(i); dq(i)]
    dz_sim = zeros(6,length(tspan_sim)); %[dq(i); ddq(i)]
    
    % set initial states
    z_sim(:,1) = z_TO(:,1);
    
    % controller values
    Kp = 1000*eye(3);
    Kd = 20*eye(3);
    sig_q = 0.0015; % approximating position measurement noise
    sig_dq = sig_q*(1/sim_dt)*(1/50); % approximating velocity measurement noise
    mu_s = 0.05; % static friction
    nu_v = 0.02; % viscous friction
    
    % initialize data storage
    q_meas_sim = zeros(3,num_sim_steps);
    dq_meas_sim = zeros(3,num_sim_steps);
    ee_sim = zeros(2,num_sim_steps);
    ee_ff = zeros(2,num_sim_steps);
    meff_sim = zeros(5,num_sim_steps);
    peff_sim = zeros(3,num_sim_steps);
    % forward simulation
    %TODO: could switch to ode45 here? 
    for ii=1:num_sim_steps
        % at each time-step...
        % TODO: make each step here into a helper function
        
        % limits used in optimization...
%         q_ub = [pi;pi;pi];
%         q_lb = [-pi;-pi;-pi];
% 
%         dq_ub = [30;30;30];
%         dq_lb = [-30;-30;-30];
% 
        u_ub = [20;20;20];
        u_lb = [-20;-20;-20];      
        
        % get new state measurement
        q_meas_sim(:,ii) = z_sim(1:3,ii) + normrnd(0,sig_q,[3,1]);
        dq_meas_sim(:,ii) = z_sim(4:6,ii) + normrnd(0,sig_dq,[3,1]);
        
        % evaluate control law
        u_temp = u_TO_ff(:,ii) + Kp*(z_TO_ff(1:3,ii)-q_meas_sim(:,ii)) + Kd*(z_TO_ff(4:6,ii)-dq_meas_sim(:,ii));
        u_sim(:,ii) = min( max( u_temp, u_lb ), u_ub); % enfore input limits
        %TODO: add noise here? have measurement model for z_sim(:,ii) used in control law
        
        % evaluate friction
        Qfric = -mu_s*sign(z_sim(4:6,ii)) - nu_v*z_sim(4:6,ii); 
        
        % get ddq
        A = A_arm(z_sim(:,ii),p); % mass matrix
        b = b_arm(z_sim(:,ii),u_sim(:,ii),p); % remainder of eoms
        % TODO: add friction as generalized force? -> b = b + friction
        Fc = [0, 0]; % contact force?
        Tauc = 0; % joint limit torques?
        J = jacobian_tip(z_sim(:,ii),p);
        QFc = J'*[Fc, 0]'; % contact force contribution to generalized force
        QTauc = [0, 0, Tauc]'; % constraint torque contribution
        ddq_i = A\(b+Qfric); % with contact force and constraint: A\(b+Qfric+QFc+QTauc)
        
        % forward step
        dz_sim(:,ii) = [z_sim(4:6,ii); ddq_i]; % [dq(i); ddq(i)]
        % TODO: add noise here?
        if ii<num_sim_steps
            z_sim(:,ii+1) = z_sim(:,ii) + sim_dt*dz_sim(:,ii);
        end  
        
        % TODO: store other variables of interest, like cartesian tip position and effective mass? 
        
         % calculate endpoint velocity
        ve_temp = jacobian_tip(z_sim(:,ii),p)*z_sim(4:6,ii);
        ve_i = ve_temp(1:2); % get actual tip velocity
        ev_i = ve_i/norm(ve_i);
        % calculate effective mass, minimum, maximum
        LLv_inv = LLv_arm_op_inv(z_sim(:,ii),p);
        meff_sim(1,ii) = 1/(ev_i'*LLv_inv*ev_i);
        [V,D] = eig(LLv_inv);
        if D(1,1)>D(2,2)
            meff_sim(2,ii) = 1/D(1,1);
            meff_sim(3,ii) = 1/D(2,2);
            meff_sim(4:5,ii) = V(:,1);
        else
            meff_sim(2,ii) = 1/D(2,2);
            meff_sim(3,ii) = 1/D(1,1);
            meff_sim(4:5,ii) = V(:,2);
        end
        % calculate effective momentum
        peff_sim(1:2,ii) = ve_i*meff_sim(1,ii);
        peff_sim(3,ii) = norm(ve_i)*meff_sim(1,ii);
        
        ee_temp = position_tip(z_sim(:,ii),p);
        ee_ff_temp = position_tip(z_TO_ff(:,ii),p);
        
        ee_sim(:,ii) = ee_temp(1:2);
        ee_ff(:,ii) = ee_ff_temp(1:2);
        
    end
    
    %% Generate any static plots?
    
    % TODO: add plots of joint trajectory tracking?
    
    %% Return single matrix with [t(i); qff(i); dqff(i); q(i); dq(i); ddq(i); uff(i); u(i)] for each sim time-step
    
    sim_data = [tspan_sim; z_TO_ff; z_sim; dz_sim(4:6,:); u_TO_ff; u_sim; meff_sim; peff_sim; ee_sim; ee_ff];    
    
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

function Fc = contact_force(z,p)

    %% Fixed parameters for contact
    K_c = 1000;
    D_c = 20;
    yC  = -.125;
    
    rE = position_foot(z,p);
    vE = velocity_foot(z,p);
    C = rE(2)-yC; % C gives height above ground
    dC= vE(2);    % dC gives foot velocity
    
    if C < 0 % if foot is below ground
        Fc = -K_c * C - D_c * dC;
        if Fc < 0 % if foot force is pulling
            Fc = 0; % set it to zero
        end
    else
        Fc = 0;
    end
end

function Tauc = joint_limit_torque(z,p)
    %% Fixed parameters for rotational spring damper at joint limit contact
    Kappa_c = 10;
    Dampa_c = .2;
    
    q2_min = 10 * pi/ 180;
    C = z(2) - q2_min; % C gives distance away from constraint
    dC= z(4);
    
    if C < 0 % if constraint is violated
        Tauc = -Kappa_c * C - Dampa_c * dC;
        if Tauc < 0 % if foot force is pulling
            Tauc = 0; % set it to zero
        end
    else
        Tauc = 0;
    end
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