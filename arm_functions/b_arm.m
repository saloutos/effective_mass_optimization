function b = b_arm(in1,in2,in3)
%B_ARM
%    B = B_ARM(IN1,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    03-May-2021 15:46:19

dth1 = in1(4,:);
dth2 = in1(5,:);
dth3 = in1(6,:);
l_AB = in3(15,:);
l_A_m2 = in3(12,:);
l_B_m3 = in3(13,:);
l_OA = in3(14,:);
m2 = in3(2,:);
m3 = in3(3,:);
m_motor = in3(4,:);
tau1 = in2(1,:);
tau2 = in2(2,:);
tau3 = in2(3,:);
th2 = in1(2,:);
th3 = in1(3,:);
t2 = sin(th2);
t3 = sin(th3);
t4 = th2+th3;
t5 = dth1.^2;
t6 = sin(t4);
t7 = dth3.*l_AB.*l_B_m3.*m3.*t3;
t8 = dth1.*l_AB.*l_B_m3.*m3.*t3.*2.0;
t9 = dth2.*l_AB.*l_B_m3.*m3.*t3.*2.0;
t10 = dth1.*l_B_m3.*l_OA.*m3.*t6;
t11 = dth2.*l_B_m3.*l_OA.*m3.*t6;
t12 = dth3.*l_B_m3.*l_OA.*m3.*t6;
t15 = l_B_m3.*l_OA.*m3.*t5.*t6;
t13 = dth3.*t10;
t14 = t10.*2.0;
t17 = -t15;
t16 = -t13;
b = [tau1+dth2.*(t11+t12+t14+dth1.*l_AB.*l_OA.*m3.*t2.*2.0+dth2.*l_AB.*l_OA.*m3.*t2+dth1.*l_A_m2.*l_OA.*m2.*t2.*2.0+dth2.*l_A_m2.*l_OA.*m2.*t2+dth1.*l_AB.*l_OA.*m_motor.*t2.*2.0+dth2.*l_AB.*l_OA.*m_motor.*t2)+dth3.*(t7+t8+t9+t11+t12+t14);t16+t17+tau2+dth2.*(t10+dth1.*l_AB.*l_OA.*m3.*t2+dth1.*l_A_m2.*l_OA.*m2.*t2+dth1.*l_AB.*l_OA.*m_motor.*t2)-dth2.*t10+dth3.*(t7+t8+t9+t10)-l_AB.*l_OA.*m3.*t2.*t5-l_A_m2.*l_OA.*m2.*t2.*t5-l_AB.*l_OA.*m_motor.*t2.*t5-dth1.*dth2.*l_AB.*l_OA.*m3.*t2-dth1.*dth2.*l_A_m2.*l_OA.*m2.*t2-dth1.*dth2.*l_AB.*l_OA.*m_motor.*t2;t16+t17+tau3-dth1.*t7-dth2.*t7+dth3.*(t10+dth1.*l_AB.*l_B_m3.*m3.*t3+dth2.*l_AB.*l_B_m3.*m3.*t3)-l_AB.*l_B_m3.*m3.*t3.*t5-dth2.^2.*l_AB.*l_B_m3.*m3.*t3-dth1.*dth2.*l_AB.*l_B_m3.*m3.*t3.*2.0];
