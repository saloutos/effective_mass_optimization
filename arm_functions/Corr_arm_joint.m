function Corr_Joint_Sp = Corr_arm_joint(in1,in2)
%CORR_ARM_JOINT
%    CORR_JOINT_SP = CORR_ARM_JOINT(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    15-Mar-2021 22:37:37

Ir = in2(9,:);
N = in2(10,:);
dth1 = in1(4,:);
dth2 = in1(5,:);
dth3 = in1(6,:);
l_AB = in2(15,:);
l_A_m2 = in2(12,:);
l_B_m3 = in2(13,:);
l_OA = in2(14,:);
m2 = in2(2,:);
m3 = in2(3,:);
m_motor = in2(4,:);
th1 = in1(1,:);
th2 = in1(2,:);
th3 = in1(3,:);
t2 = sin(th2);
t3 = sin(th3);
t4 = th2+th3;
t5 = dth1.^2;
t6 = dth2.^2;
t7 = dth3.^2;
t8 = sin(t4);
t9 = l_AB.*l_B_m3.*m3.*t3.*t7;
t10 = dth1.*dth3.*l_AB.*l_B_m3.*m3.*t3.*2.0;
t11 = dth2.*dth3.*l_AB.*l_B_m3.*m3.*t3.*2.0;
t12 = -t10;
t13 = -t11;
t14 = -t9;
Corr_Joint_Sp = [t12+t13+t14-Ir.*N.^2.*th1-l_AB.*l_OA.*m3.*t2.*t6-l_A_m2.*l_OA.*m2.*t2.*t6-l_B_m3.*l_OA.*m3.*t6.*t8-l_B_m3.*l_OA.*m3.*t7.*t8-l_AB.*l_OA.*m_motor.*t2.*t6-dth1.*dth2.*l_AB.*l_OA.*m3.*t2.*2.0-dth1.*dth2.*l_A_m2.*l_OA.*m2.*t2.*2.0-dth1.*dth2.*l_B_m3.*l_OA.*m3.*t8.*2.0-dth1.*dth3.*l_B_m3.*l_OA.*m3.*t8.*2.0-dth2.*dth3.*l_B_m3.*l_OA.*m3.*t8.*2.0-dth1.*dth2.*l_AB.*l_OA.*m_motor.*t2.*2.0;t12+t13+t14+l_AB.*l_OA.*m3.*t2.*t5+l_A_m2.*l_OA.*m2.*t2.*t5+l_B_m3.*l_OA.*m3.*t5.*t8+l_AB.*l_OA.*m_motor.*t2.*t5;l_B_m3.*m3.*(l_AB.*t3.*t5+l_AB.*t3.*t6+l_OA.*t5.*t8+dth1.*dth2.*l_AB.*t3.*2.0)];
