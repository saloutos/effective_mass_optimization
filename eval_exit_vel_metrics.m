function [phi_vf, Pvf, phi_vo, Pvo, phi_vol, Pvol] = eval_exit_vel_metrics(q,p,LLo) 
    % evaluate both exit velocity metrics at a given configuration
    % can also return the directional metric for a given finger velocity

    % q is [xb; yb; thb; th1; th2; th3] (floating-base state)
    % p is robot parameters in column vector, [m1 m2 m3 m_motor I1 I2 I3 I_motor Ir N l_O_m1 l_A_m2 l_B_m3 l_OA l_AB l_BC g mb Ib]'
    % LLo is OSIM for object
    % ev is velocity direction [vx; vy]
    
    % evaluate H, J
    H = H_arm_floating(q,p);
    J = J_arm_floating(q,p);
    Hbb = H(1:3,1:3);
    Jb = J(:,1:3);

    % evaluate OSIMs
    LL_inv = J/H*J';
    LLf = inv(LL_inv);
    LLfl_inv = Jb/Hbb*Jb';
    LLfl = inv(LLfl_inv);

    % evaluate partial OSIM ratios
    LLpf = (LLf+LLo)\LLf;
    LLpfl = (LLfl+LLo)\LLfl;

    % evaluate metrics
    Pvf = 2*(LLpfl-LLpf);
    Pvo = (LLpfl/LLpf)-eye(2);
    Pvol = eye(2)-(LLpf/LLpfl);
    phi_vf = det(Pvf);
    phi_vo = det(Pvo);
    phi_vol = det(Pvol);

%     % along a specified finger velocity direction
%     vf = vf/norm(vf); % normalize
%     nu_vf = vf'*Pvf*vf;

end