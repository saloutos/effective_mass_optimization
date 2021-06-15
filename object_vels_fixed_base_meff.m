function [vop, vfp] = object_vels_fixed_base_meff(q,p,mo,vf)

    % evaluate H, J
    H = A_arm(q,p);
    J = jacobian_v_tip(q,p);

    % evaluate OSIMs
    LL_inv = J/H*J';
    LLf = inv(LL_inv);

    % velocity direction
    u = vf/norm(vf);
    
    % effective mass
    mf = 1/(u'*LL_inv*u);
    
    % get mass ratio
    mp = mf/(mf+mo);
    
    % object velocity
    vop = 2*mp*vf;
    
    % finger velocity
    vfp = (2*mp - 1)*vf;
    
%     % evaluate partial OSIM ratios
%     LLpf = (LLf+LLo)\LLf;
%     % calculate object velocities
%     vo = 2*LLpf*vf;
%     % calculate finger velocities
%     vfp = (2*LLpf-eye(2))*vf;
    
end