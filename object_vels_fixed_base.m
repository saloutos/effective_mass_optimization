function [vo, vfp] = object_vels_fixed_base(q,p,LLo,vf)

    % evaluate H, J
    H = A_arm(q,p);
    J = jacobian_v_tip(q,p);

    % evaluate OSIMs
    LL_inv = J/H*J';
    LLf = inv(LL_inv);

    % evaluate partial OSIM ratios
    LLpf = (LLf+LLo)\LLf;
    
    % calculate object velocities
    vo = 2*LLpf*vf;
    % calculate finger velocities
    vfp = (2*LLpf-eye(2))*vf;
    
end