% Andrew SaLoutos
% 3/3/21

% optimizing effective mass across many points?
% derived from "meff_opt_single_point.m" script

% randomly selected points and velocity normal
pts_x = [0.4, 0.6, 1.0, 1.0, 0.6, 0.4];
pts_y = [-0.7, -0.4, -0.25, 0.25, 0.4, 0.7];

vels_x = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0];
vels_y = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0];

num_pts = length(pts_x);

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

% symbolic variables for optimization
syms q1 q2 q3 real
q = [q1; q2; q3];

% loop through and plot each point
figure; sgtitle('Minimum m_{eff} at Various Poses');
spr = 2;
spc = num_pts/spr;
for ii=1:num_pts
   
   pt = [pts_x(ii); pts_y(ii)];
   vdir = [vels_x(ii); vels_y(ii)];
   vdir = vdir/norm(vdir);
   
   [hist_start, hist_warm] = meff_opt_single_point_func(q,p,vdir,pt);
   meff_final = hist_warm(end,1);
   
   subplot(spr,spc,ii); hold on;
   plot_ik_sol_only(hist_warm,p,vdir,pt);
   sp_title = sprintf('Pt: (%.2f, %.2f); V: (%.2f, %.2f)', pt(1), pt(2), vdir(1), vdir(2));
   title(sp_title);
   
end



% wrap previous script in function here for easy use
function [hist_start, hist_warm] = meff_opt_single_point_func(q,p,vdir,pt)

    % generate cost function (and gradients and hessians?)
    % want to maximize dot product of vdir and minimum effective mass

    % generate constraint functions? (and gradients and hessians?)
    % want to constraint joint angles to satisfy forward kinematics
    % nonlinear constrains need to be of the form ceq=0 or c<=0

    % run optimization
    q0 = inverse_kinematics(q,p,vdir,pt); % give a feasible initial solution
    cost_func1 = @(q)warm_start(q,p,vdir,pt);
    cost_func2 = @(q)meff_cost(q,p,vdir,pt);
    nonlcon = @(q)joint_constraints(q,p,vdir,pt);
    A = []; % No other constraints
    b = [];
    Aeq = [];
    beq = [];
    lb = [];
    ub = [];
    options = optimoptions('fmincon','display','none');
    % warm start opt
    tic;
    [q_start, fval_start, hist_start] = fmincon_with_hist(cost_func1,q0,A,b,Aeq,beq,lb,ub,nonlcon,options);
    opt1 = toc;
    % effective mass opt
    tic;
    [q_sol_warm, fval_warm, hist_warm] = fmincon_with_hist(cost_func2,q_start,A,b,Aeq,beq,lb,ub,nonlcon,options);
    opt2 = toc;
    
    % generate plot titles based on pt, vdir here
    
    
    fprintf('Displaying "warm start" IK solution...q = [%.2f, %.2f, %.2f]\n',q_start(1),q_start(2),q_start(3)');
%     plot_sols_and_optims(hist_start,p,vdir,pt,'Warm Start, Cost: Alignment of third link with v_{dir}');
    fprintf('Displaying minimum effective mass solution...q = [%.2f, %.2f, %.2f]\n',q_sol_warm(1),q_sol_warm(2),q_sol_warm(3)');
%     plot_sols_and_optims(hist_warm,p,vdir,pt,'Cost: m_{eff}');

    %%%%% helper functions %%%%%
    function [c,ceq] = joint_constraints(q,p,vdir,rE)

        % could also calculate gradients of constraints here?

        z = [q;zeros(3,1)];
        c = [];
        rE_temp = position_tip(z,p);
        ceq = rE_temp(1:2)-rE; % just care about x,y position
    end

    function total_cost = warm_start(q,p,vdir,rE)

        th_link_3 = q(1)+q(2)+q(3); % angle between third link and world x-axis
        e_link_3 = [cos(th_link_3); sin(th_link_3)];

        % minimize dot product of e_link_3 and vdir
        a1 = 1;
        c1 = a1*(e_link_3'*vdir)^2;
        
        % minimize joint angles too
        Q = eye(3);
        a2 = 0.1;
        c2 = a2*q'*Q*q;

        total_cost = c1 + c2;

    end

    function total_cost = meff_cost(q,p,vdir,rE)

        % could also calculate gradient and hessians of costs here?

        z = [q; zeros(3,1)];

        LLv_inv = LLv_arm_op_inv(z,p); % just care about x and y for now
        meff = 1/(vdir'*LLv_inv*vdir); % effective mass
        
        [V,D] = eig(LLv_inv);
        if D(1,1) > D(2,2) % direction of minimum effective mass
            meff_max_dir = V(:,2); 
        else
            meff_max_dir = V(:,1);
        end
        
        % minimize effective mass along vdir
        a1 = 1;
        c1 = a1*meff; %alpha*meff2; % these two are not the same!
        
        % minimize dot product of maximum effective mass direction and vdir
        a2 = 0; %1;
        c2 = a2*(meff_max_dir'*vdir)^2;
        
        % add cost on overall joint angles (either absolute values or relative to previous angles?)
        Q = eye(3);
        a3 = 0; %0.1;
        c3 = a3*q'*Q*q;

        % return total cost
        total_cost = c1 + c2 + c3;
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

    function q_ik = inverse_kinematics(q,p,vdir,rE)
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

    % does this function need to be modified here?
    function plot_sols_and_optims(sol_hist,p,vdir,rE,plt_title)

        % should take initial posture and final posture (maybe intermediate?)
        % plot initial posture and ellipsoid in gray
        % plot final posture in black, with ellipsoid, effective mass directions, and velocity vector
        % plot evolution of cost function and effective mass
        % no animation, so we can just plot everything right away

        q_iters = sol_hist(:,3:5)';
        fval_iters = sol_hist(:,1);
        constr_iters = sol_hist(:,2);
        num_pts = size(q_iters,2);
        v_norm = vdir/norm(vdir); % make sure input velocity direction is a unit vector

        figure; clf; 
        sgtitle(plt_title);

        %%% first subplot with arm kinematics

        %%% keypoints and lmdir for initial pose
        kp0 = keypoints_arm(q_iters(:,1),p);
        rA0 = kp0(:,1);
        rB0 = kp0(:,2);
        rC0 = kp0(:,3);
        % calculate ellipse for operational space inertia
        z_ik0 = [q_iters(:,1); zeros(3,1)];
        LLv_inv0 = LLv_arm_op_inv(z_ik0,p); % just care about x and y for now
        [V0,D0] = eig(LLv_inv0);
        meff10 = 1/(v_norm'*LLv_inv0*v_norm); % effective mass?
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
        kpf = keypoints_arm(q_iters(:,num_pts),p);
        rAf = kpf(:,1);
        rBf = kpf(:,2);
        rCf = kpf(:,3);
        % calculate ellipse for operational space inertia
        z_ikf = [q_iters(:,num_pts); zeros(3,1)];
        LLv_invf = LLv_arm_op_inv(z_ikf,p); % just care about x and y for now
        [Vf,Df] = eig(LLv_invf);
        meff1f = 1/(v_norm'*LLv_invf*v_norm); % effective mass?
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
        subplot(1,2,1); hold on;
        title('Arm Kinematics');

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
        h_v = quiver([rE(1)],[rE(2)],[v_norm(1)],[v_norm(2)],0.25,'b','LineWidth',2);
        xlabel('X'); ylabel('Y'); axis equal; 
        xlim([-0.5,1.5]); ylim([-1.5,0.5]);

        h_meff1f = text(1.,-1.,'','FontSize',10);
        h_meffminf = text(1.,-1.25,'','FontSize',10);
        set(h_meff1f,'String',  sprintf('m_{eff,1} = %.2f',meff1f));
        set(h_meffminf,'String',  sprintf('m_{eff,min} = %.2f',meffminf));

        %%% second subplot with cost function values and effective mass

        %%% calculate effective mass for each iteration
        m1_vec = zeros(1,num_pts);
        pt_vec = 1:num_pts;
        for ii=1:num_pts
            % calculate ellipse for operational space inertia
            z_ik = [q_iters(:,ii); zeros(3,1)];
            LLv_inv = LLv_arm_op_inv(z_ik,p); % just care about x and y for now
            [V,D] = eig(LLv_inv);
            meff1 = 1/(v_norm'*LLv_inv*v_norm); % effective mass?
            m1_vec(ii) = meff1;        
        end

        %%% plot these values
        subplot(1,2,2); 
        title('Optimization Values');
        yyaxis left;
        h_fval = plot(pt_vec,fval_iters,'LineWidth',1.5);
        ylabel('Cost Function');
        yyaxis right;
        h_m1 = plot(pt_vec,constr_iters,'LineWidth',1.5);
        ylabel('Constraint Violations');
        xlabel('Optimization Iteration');
        xlim([0,num_pts]);

    end

end

function plot_ik_sol_only(sol_hist,p,vdir,rE)

    % doesn't create a new figure/plot, need to call subplot/figure before and also add title after

    % should take initial posture and final posture (maybe intermediate?)
    % plot initial posture and ellipsoid in gray
    % plot final posture in black, with ellipsoid, effective mass directions, and velocity vector

    q_iters = sol_hist(:,3:5)';
    num_pts = size(q_iters,2);

    %%% plot arm kinematics

    %%% keypoints and lmdir for initial pose
    kp0 = keypoints_arm(q_iters(:,1),p);
    rA0 = kp0(:,1);
    rB0 = kp0(:,2);
    rC0 = kp0(:,3);
    % calculate ellipse for operational space inertia
    z_ik0 = [q_iters(:,1); zeros(3,1)];
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
    kpf = keypoints_arm(q_iters(:,num_pts),p);
    rAf = kpf(:,1);
    rBf = kpf(:,2);
    rCf = kpf(:,3);
    % calculate ellipse for operational space inertia
    z_ikf = [q_iters(:,num_pts); zeros(3,1)];
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








