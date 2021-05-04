function E = energy_arm(in1,in2)
%ENERGY_ARM
%    E = ENERGY_ARM(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    03-May-2021 15:46:19

I1 = in2(5,:);
I2 = in2(6,:);
I3 = in2(7,:);
I_motor = in2(8,:);
Ir = in2(9,:);
N = in2(10,:);
dth1 = in1(4,:);
dth2 = in1(5,:);
dth3 = in1(6,:);
l_AB = in2(15,:);
l_A_m2 = in2(12,:);
l_B_m3 = in2(13,:);
l_OA = in2(14,:);
l_O_m1 = in2(11,:);
m1 = in2(1,:);
m2 = in2(2,:);
m3 = in2(3,:);
m_motor = in2(4,:);
th2 = in1(2,:);
th3 = in1(3,:);
t2 = cos(th2);
t3 = cos(th3);
t4 = th2+th3;
t5 = N.^2;
t6 = dth1.^2;
t7 = dth2.^2;
t8 = dth3.^2;
t9 = l_AB.^2;
t10 = l_A_m2.^2;
t11 = l_B_m3.^2;
t12 = l_OA.^2;
t13 = cos(t4);
E = (I1.*t6)./2.0+(I2.*t6)./2.0+(I2.*t7)./2.0+(I3.*t6)./2.0+(I3.*t7)./2.0+(I3.*t8)./2.0+I_motor.*t6+(I_motor.*t7)./2.0+Ir.*t6+(Ir.*t7)./2.0+(l_O_m1.^2.*m1.*t6)./2.0+I2.*dth1.*dth2+I3.*dth1.*dth2+I3.*dth1.*dth3+I3.*dth2.*dth3+I_motor.*dth1.*dth2+Ir.*dth1.*dth2+(Ir.*t5.*t6)./2.0+(Ir.*t5.*t7)./2.0+(Ir.*t5.*t8)./2.0+(m2.*t6.*t10)./2.0+(m3.*t6.*t9)./2.0+(m2.*t7.*t10)./2.0+(m3.*t7.*t9)./2.0+(m2.*t6.*t12)./2.0+(m3.*t6.*t11)./2.0+(m3.*t6.*t12)./2.0+(m3.*t7.*t11)./2.0+(m3.*t8.*t11)./2.0+(m_motor.*t6.*t9)./2.0+(m_motor.*t7.*t9)./2.0+m_motor.*t6.*t12+dth1.*dth2.*m2.*t10+dth1.*dth2.*m3.*t9+dth1.*dth2.*m3.*t11+dth1.*dth3.*m3.*t11+dth2.*dth3.*m3.*t11+dth1.*dth2.*m_motor.*t9+Ir.*N.*dth1.*dth2+Ir.*N.*dth1.*dth3+Ir.*N.*dth2.*dth3+l_AB.*l_B_m3.*m3.*t3.*t6+l_AB.*l_B_m3.*m3.*t3.*t7+l_AB.*l_OA.*m3.*t2.*t6+l_A_m2.*l_OA.*m2.*t2.*t6+l_B_m3.*l_OA.*m3.*t6.*t13+l_AB.*l_OA.*m_motor.*t2.*t6+dth1.*dth2.*l_AB.*l_B_m3.*m3.*t3.*2.0+dth1.*dth3.*l_AB.*l_B_m3.*m3.*t3+dth2.*dth3.*l_AB.*l_B_m3.*m3.*t3+dth1.*dth2.*l_AB.*l_OA.*m3.*t2+dth1.*dth2.*l_A_m2.*l_OA.*m2.*t2+dth1.*dth2.*l_B_m3.*l_OA.*m3.*t13+dth1.*dth3.*l_B_m3.*l_OA.*m3.*t13+dth1.*dth2.*l_AB.*l_OA.*m_motor.*t2;
