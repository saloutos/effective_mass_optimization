% plotting surfaces of exit vel metrics wrt object and floating-base inertias

% Andrew SaLoutos
% 5/5/2021

%% add functions

addpath(genpath('arm_floating_functions'))

%% test parameters
m1 = 1;                 m2 = 1;
m3 = 0.5;               m_motor = 0.5;
I1 = 0.02;              I2 = 0.02;
I3 = 0.01;              I_motor = 0.000625; % assuming radius r=0.05m
Ir = 6.25e-6;           N = 6;
l_O_m1 = 0.25;          l_A_m2 = 0.25;
l_B_m3 = 0.25;          l_OA = 0.5;
l_AB = 0.5;             l_BC = 0.5;
g = 9.81; % do I want gravity to start?

mb = 1; Ib = 0.5*mb*0.1^2;

% arrange in vector
p   = [m1 m2 m3 m_motor I1 I2 I3 I_motor Ir N l_O_m1 l_A_m2 l_B_m3 l_OA l_AB l_BC g mb Ib]';

% define object inertia
mo = 1.0;
LLo = mo*eye(2);

%% test configuration
xb = 0; yb = 0; thb = 0;
th1 = 1.1; %-0.5;
th2 = -2.2; %1.0;
th3 = 1.1; %-0.5;

q = [xb; yb; thb; th1; th2; th3];

% along a specified finger velocity direction
vf = [ 0.0; 1.0];
evf = vf/norm(vf);


%% surface 1: mo, mb, eta_vol

mb_vec = linspace(0,20,100); %linspace(0,5,20);
mo_vec = linspace(0.1,20,100); % linspace(0.1,2,20);
[X,Y] = meshgrid(mb_vec,mo_vec);
Z1 = zeros(size(X));
Z2 = zeros(size(X));
Z3 = zeros(size(X));

for ii=1:length(mb_vec)
    mb_temp = mb_vec(ii);
    Ib_temp = 0.5*mb_temp*0.1^2;
    p(18) = mb_temp;
    p(19) = Ib_temp;
    for jj=1:length(mo_vec)
        mo_temp = mo_vec(jj);
        LLo_temp = mo_temp*eye(2);
        
        % get post-collision object velocities
        [vo, vol, vfp, vfpl] = object_vels(q,p,LLo_temp,vf);

        % eval exit vel metrics
        [phi_vf, Pvf, phi_vo, Pvo, phi_vol, Pvol] = eval_exit_vel_metrics(q,p,LLo_temp);

        
        evo = vo/norm(vo);
        evol = vol/norm(vol);

        eta_vf = evf'*Pvf*evf;
        eta_vo = evo'*Pvo*evo;
        eta_vol = evol'*Pvol*evol;
        
        Z1(jj,ii) = eta_vol; % fill Z matrix
        Z2(jj,ii) = phi_vol; 
        Z3(jj,ii) = phi_vf;
        
    end
end

% figure; surf(X,Y,Z1);
% xlabel('m_b'); ylabel('m_o'); zlabel('\eta_{vol}');

figure; surf(X,Y,Z2);
xlabel('m_b'); ylabel('m_o'); zlabel('\phi_{vol}');
str1 = sprintf('\\phi_{vol} at q = [%.1f, %.1f, %.1f] for varying m_o and m_b', q(4),q(5),q(6));
title(str1);

figure; surf(X,Y,Z3);
xlabel('m_b'); ylabel('m_o'); zlabel('\phi_{vf}');
str2 = sprintf('\\phi_{vf} at q = [%.1f, %.1f, %.1f] for varying m_o and m_b', q(4),q(5),q(6));
title(str2);

%% surface 2: mo, link3 angle?, link2 angle?, vf angle?, eta_vol

mb = 2; Ib = 0.5*mb*0.1^2;
p(18) = mb; p(19) = Ib;

mo = 1.0; LLo = mo*eye(2); 
link2_angles = deg2rad( linspace(-180,180,50) );
link3_angles = deg2rad( linspace(-90,90,50) );

[X,Y] = meshgrid(link2_angles, link3_angles);
Z1 = zeros(size(X));
Z2 = zeros(size(X));

for ii=1:length( link2_angles )
    l2_temp = link2_angles(ii);
    for jj=1:length(link3_angles)
        l3_temp = link3_angles(jj);
        
        q(5) = l2_temp;
        q(6) = l3_temp;
        
        % get post-collision object velocities
        [vo, vol, vfp, vfpl] = object_vels(q,p,LLo,vf);

        % eval exit vel metrics
        [phi_vf, Pvf, phi_vo, Pvo, phi_vol, Pvol] = eval_exit_vel_metrics(q,p,LLo);

        eta_vf = evf'*Pvf*evf;
        eta_vo = evf'*Pvo*evf;
        eta_vol = evf'*Pvol*evf;
        
        Z1(jj,ii) = eta_vol; % fill Z matrix
        Z2(jj,ii) = phi_vol; 
        
    end
end

figure; surf(X,Y,Z1);
xlabel('q_2'); ylabel('q_3'); zlabel('\eta_{vol}');
str3 = sprintf('\\eta_{vol} with m_o=%.1f, m_b=%.1f, and vf = [%.1f, %.1f] for varying link angles', mo,mb,vf(1),vf(2));
title(str3);


figure; surf(X,Y,Z2);
xlabel('q_2'); ylabel('q_3'); zlabel('\phi_{vol}');
str4 = sprintf('\\phi_{vol} with m_o=%.1f and m_b=%.1f for varying link angles', mo,mb);
title(str4);










