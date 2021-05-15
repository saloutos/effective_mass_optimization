function [vo, vol, vfp, vfpl] = object_vels(q,p,LLo,vf)

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
    
    % calculate object velocities
    vo = 2*LLpf*vf;
    vol = 2*LLpfl*vf;
    % calculate finger velocities
    vfp = (2*LLpf-eye(2))*vf;
    vfpl = (2*LLpfl-eye(2))*vf;
    
end