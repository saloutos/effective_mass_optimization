function [TO_data, opt_time] = meff_optimization(pts, Tf, p, alpha_vec)

    % Andrew SaLoutos
    % 3/26/21
    
    % function for running fmincon trajectory optimization for effective
    % mass along a specified trajectory
    
    % inputs: trajectory of pts, final time of trajectory, parameter vector, cost function weight vector
    
    % outputs: concatenate [time_vec, z_sol, meff_sol, peff_sol, ee_sol] -> 21 x num_pts matrix
    
    % parse inputs
    num_pts = size(pts,2);
    time_vec = linspace(0,Tf,num_pts);
    dt = time_vec(2);
    
    vels = diff(pts,1,2)/dt; % populate velocity vectors?
    vels = [zeros(2,1), vels]; % pad with zeros for first point in trajectory
        
    % symbolic variables for optimization
    % at each time point: n inputs, n joint angles, n joint velocities
    syms z [9 num_pts] real % [q_i[1:3]; dq_i[1:3]; u_i[1:3]] for each point, i
    z = z(:); % vectorize

    % fmincon inputs
    cost_func = @(z)meff_cost(z,p,time_vec,pts,vels,alpha_vec);
    nonlcon = @(z)joint_constraints(z,p,time_vec,pts,vels);
    A = []; % no other linear inequalities
    b = [];
    Aeq = []; % constraints on initial and final positions?
    beq = [];
    [ub,lb] = joint_and_torque_limits(z,p,time_vec); % torque and joint velocity limits
    options = optimoptions('fmincon','display','none');

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
    z0 = z0(:);

    % run optimization
    tic;
    [z_sol, fval_sol, sol_hist] = fmincon_with_hist(cost_func,z0,A,b,Aeq,beq,lb,ub,nonlcon,options);
    opt_time = toc;
    z_sol = reshape(z_sol,9,num_pts);

    % calculate desired output data from solution
    q_sol = z_sol(1:3,:);
    dq_sol = z_sol(4:6,:);
    u_sol = z_sol(7:9,:);
    meff_sol = zeros(5,num_pts); % each column is [meff; meff_min; meff_max; lm_dir]
    peff_sol = zeros(3,num_pts); % px, py, |p|
    ee_sol = zeros(3,num_pts); % xact, yact, error
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
        % calculate effective momentum
        peff_sol(1:2,ii) = ve_i*meff_sol(1,ii);
        peff_sol(3,ii) = norm(ve_i)*meff_sol(1,ii);
        % get endpoint position
        ee_temp = position_tip(z_i,p);
        ee_sol(1:2,ii) = ee_temp(1:2);
        ee_error = ee_sol(1:2,ii) - pts(:,ii);
        ee_sol(3,ii) = norm(ee_error);
    end
    
    % save outputs
    TO_data = [time_vec; z_sol; meff_sol; peff_sol; ee_sol]; % concatenate everything
    
end


%%% Helper functions %%%
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
    
    ee_0 = pts(:,1);
    z_0 = [q(:,1); dq(:,1)];
    ee_temp0 = position_tip(z_0,p);
    ceq_q_0 = [ee_temp0(1:2)-ee_0; dq(:,1)];
        
    ee_f = pts(:,end);
    z_f = [q(:,end); dq(:,end)];
    ee_tempf = position_tip(z_f,p);
    ceq_q_f = [ee_tempf(1:2)-ee_f]; %; dq(:,end)];
    
    ceq = [ceq; ceq_q_0; ceq_q_f]; % only constrain initial and final positions and velocities
    
%     for ii=1:num_pts
%         % get end effector position at time_vec(ii)
%         ee_i = pts(:,ii);
%         z_i = [q(:,ii); dq(:,ii)];
%         ee_temp = position_tip(z_i,p);
%         ceq_q_i = ee_temp(1:2)-ee_i;
%         
%         % get end effector velocity (in x,y) at time_vec(ii)
%         ve_des = [vels_x(ii);vels_y(ii)];
%         ve_temp = jacobian_tip(z_i,p)*dq(:,ii);
%         ve_i = ve_temp(1:2);
%         ve_i = ve_i/norm(ve_i);
%         ve_des = ve_des/norm(ve_des);
%         ceq_dq_i = ve_i-ve_des; % direction of end effector velocity must match 
%         
%         ceq = [ceq; ceq_q_i];% ceq_dq_i];
%     end
        
end

function total_cost = meff_cost(z,p,time_vec,pts,vels,alpha_vec)
    
    num_pts = length(time_vec);
    dt = time_vec(2)-time_vec(1);
    
    z = reshape(z,9,num_pts);
    q = z(1:3,:);
    dq = z(4:6,:);
    u = z(7:9,:);
    
    meff_cost_vec = zeros(num_pts,1);
    meff_dir_cost_vec = zeros(num_pts,1);
    e3_cost_vec = zeros(num_pts,1);
    u_cost_vec = zeros(num_pts,1);
    ee_cost_vec = zeros(num_pts,1);
    eef_cost_vec = zeros(num_pts,1);
    dq_cost_vec = zeros(num_pts,1);
    v_cost_vec = zeros(num_pts,1);
    vmag_cost_vec = zeros(num_pts,1);
    
    
    
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
%         meff_cost_vec(ii) = -((ev_i'*LLv_inv*ev_i)^2); % can square to get better gradients?
        
        [V,D] = eig(LLv_inv);
        if D(1,1) > D(2,2)
            meff_max_dir = V(:,2);
            meff_min_dir = V(:,1);
        else
            meff_max_dir = V(:,1);
            meff_min_dir = V(:,2);
        end
        meff_dir_cost_vec(ii) = (meff_max_dir'*ev_i)^2; % minimize dot product of maximum effective mass direction and vdir
        Qmeff = eye(2);
        %meff_dir_cost_vec(ii) = (meff_min_dir-ev_i)'*Qmeff*(meff_min_dir-ev_i);
        
        v_cost_vec(ii) = -ev_i'*ev_des; % maximize alignment of the two vectors
        Qvdes = eye(2);
        %v_cost_vec(ii) = (ve_i-ve_des)'*Qvdes*(ve_i-ve_des);
        
        
        vmag_cost_vec(ii) = ve_i'*ve_i; % cost on overall magnitude of velocity
        
        th_link_3 = q(1,ii)+q(2,ii)+q(3,ii); % angle between third link and world x-axis
        e_link_3 = [cos(th_link_3); sin(th_link_3)];
        % minimize dot product of e_link_3 and vdir
        e3_cost_vec(ii) = (e_link_3'*ev_des)^2;
%         e3_cost_vec(ii) = (e_link_3'*ev_i)^2;
        
    end
    
    % running cost on dq between positions?
    R1 = eye(3);
    for ii=2:num_pts
        dq_i = dq(:,ii);
        dq_cost_vec(ii) = dq_i'*R1*dq_i;
    end
    
    % cost on inputs
    R2 = eye(3);
    for ii=1:num_pts
        u_i = u(:,ii);
        u_cost_vec(ii) = u_i'*R2*u_i;
    end
    
    % cost on deviation from nominal trajectory
    Q1 = 1.0*eye(2); % seems really heavy?
    Q2 = 1.0*eye(2);
    zf = [q(:,end); dq(:,end)];
    eef = position_tip(zf,p);
    eef = eef(1:2);
    for ii=1:num_pts
        z_i = [q(:,ii); dq(:,ii)];
        ee_temp = position_tip(z_i,p);
        ee_diff = ee_temp(1:2) - pts(:,ii);
        eef_diff = eef - pts(:,ii);
        ee_cost_vec(ii) = ee_diff'*Q1*ee_diff;  
        eef_cost_vec(ii) = eef_diff'*Q2*eef_diff;
    end
    
    % add costs together with weights
    
    % cost on effective mass
    a1 = alpha_vec(1)*ones(1,num_pts); 
    c1 = a1*meff_cost_vec;
    
    % cost on direction of minimum effective mass at end-effector
    a2 = alpha_vec(2)*ones(1,num_pts);
    c2 = a2*meff_dir_cost_vec;
    
    % cost on orientation of last arm link (relative to trajectory velocity)
    a3 = alpha_vec(3)*ones(1,num_pts);
    c3 = a3*e3_cost_vec;
    
    % cost on inputs
    a4 = alpha_vec(4)*ones(1,num_pts);
    c4 = a4*u_cost_vec;
        
    % cost on end-effector position (relative to trajectory)
    a5 = alpha_vec(5)*ones(1,num_pts);
    c5 = a5*ee_cost_vec;
    
    % cost on distance to end of trajectory
    a6 = alpha_vec(6)*ones(1,num_pts);
    c6 = a6*eef_cost_vec;
        
    % cost on magnitude of joint velocities
    a7 = alpha_vec(7)*ones(1,num_pts);
    c7 = a7*dq_cost_vec;
    
    % cost on magnitude of end-effector velocity
    a8 = alpha_vec(8)*ones(1,num_pts);
    c8 = a8*vmag_cost_vec;
    
    % cost on direction of end-effector velocity (relative to trajectory)
    a9 = alpha_vec(9)*ones(1,num_pts);
    c9 = a9*v_cost_vec;
    
    % other costs to add...
    % final time? (if times are included as optimization variables?)
    
    % return total cost
    total_cost = c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8 + c9;
    
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

function q_ik = inverse_kinematics_perp(p,pt_i,vel_i)
    % pull out necessary parameters
    l1 = p(14);
    l2 = p(15);
    l3 = p(16);

    % pull out desired end-effector position
    % note: will have to do interpolation eventually here
    x_ee = pt_i(1);
    y_ee = pt_i(2);

    % get vector perpendicular to velocity
    ev_i = vel_i/norm(vel_i);
    e3_1 = [-ev_i(2); ev_i(1)]; % two candidates
    e3_2 = [ev_i(2); -ev_i(1)];
    
    x_B1 = x_ee + l3*e3_1(1);
    x_B2 = x_ee + l3*e3_2(1);
    y_B1 = y_ee + l3*e3_1(2);
    y_B2 = y_ee + l3*e3_2(2);
    
    
    
    % check which candidate is closer to origin?
    rB1 = sqrt(x_B1^2 + y_B1^2);
    rB2 = sqrt(x_B2^2 + y_B2^2);
    if rB1<=(l1+l2) || rB2<=(l1+l2)
        if rB1<rB2
            x_B = x_B1;
            y_B = y_B1;
            e3 = e3_1;
        else
            x_B = x_B2;
            y_B = y_B2;
            e3 = e3_2;
        end
        
        q2_ik = acos((x_B.^2+y_B.^2-l1^2-l2^2)/(2*l1*l2));
        q1_ik = atan2(y_B,x_B) - atan2(l2*sin(q2_ik),l1+l2*cos(q2_ik));
        q3_ik = atan2(e3(2),e3(1))-q2_ik-q1_ik;        
    else
        % return an error here?
        disp('Neither perpendicular direction is feasible');
        
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
    end
    
    % outputs
    q_ik = [q1_ik;q2_ik;q3_ik];
    
end

function total_cost = warm_start(x,p,vel_i)
    
    th_link_3 = x(1)+x(2)+x(3); % angle between third link and world x-axis
    e_link_3 = [cos(th_link_3); sin(th_link_3)];
    ev_i = vel_i/norm(vel_i);
    
    % minimize dot product of e_link_3 and vdir
    alpha = 1;
    c1 = alpha*(e_link_3'*ev_i)^2;
     
    total_cost = c1;
    
end

function [c,ceq] = warm_start_joint_constraints(x,p,pt_i)
    % could also calculate gradients of constraints here?
    z = [x;zeros(3,1)];
    c = [];
    ee_temp = position_tip(z,p);
    ceq = ee_temp(1:2)-pt_i; % just care about x,y position
end