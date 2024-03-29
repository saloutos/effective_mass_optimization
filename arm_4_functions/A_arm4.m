function A = A_arm4(in1,in2)
%A_ARM4
%    A = A_ARM4(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    03-May-2021 15:42:43

I1 = in2(6,:);
I2 = in2(7,:);
I3 = in2(8,:);
I4 = in2(9,:);
I_motor = in2(10,:);
Ir = in2(11,:);
N = in2(12,:);
l_AB = in2(18,:);
l_A_m2 = in2(14,:);
l_BC = in2(19,:);
l_B_m3 = in2(15,:);
l_C_m4 = in2(16,:);
l_OA = in2(17,:);
l_O_m1 = in2(13,:);
m1 = in2(1,:);
m2 = in2(2,:);
m3 = in2(3,:);
m4 = in2(4,:);
m_motor = in2(5,:);
th2 = in1(2,:);
th3 = in1(3,:);
th4 = in1(4,:);
t2 = cos(th2);
t3 = cos(th3);
t4 = cos(th4);
t5 = th2+th3;
t6 = th3+th4;
t7 = I_motor.*2.0;
t8 = Ir.*2.0;
t9 = N.^2;
t10 = l_AB.^2;
t11 = l_A_m2.^2;
t12 = l_BC.^2;
t13 = l_B_m3.^2;
t14 = l_C_m4.^2;
t15 = l_OA.^2;
t16 = Ir.*N;
t17 = cos(t5);
t18 = cos(t6);
t19 = t5+th4;
t21 = Ir.*t9;
t22 = m3.*t10;
t23 = m4.*t10;
t24 = m2.*t11;
t25 = m4.*t12;
t26 = m3.*t13;
t27 = m4.*t14;
t28 = m_motor.*t12;
t29 = l_AB.*l_BC.*m4.*t3;
t30 = l_AB.*l_B_m3.*m3.*t3;
t31 = l_BC.*l_C_m4.*m4.*t4;
t32 = l_AB.*l_OA.*m3.*t2;
t33 = l_AB.*l_OA.*m4.*t2;
t34 = l_A_m2.*l_OA.*m2.*t2;
t35 = l_AB.*l_BC.*m_motor.*t3;
t36 = m_motor.*t10.*2.0;
t41 = l_AB.*l_OA.*m_motor.*t2.*2.0;
t20 = cos(t19);
t37 = t29.*2.0;
t38 = t30.*2.0;
t39 = t31.*2.0;
t40 = t35.*2.0;
t42 = l_AB.*l_C_m4.*m4.*t18;
t43 = l_BC.*l_OA.*m4.*t17;
t44 = l_B_m3.*l_OA.*m3.*t17;
t45 = l_BC.*l_OA.*m_motor.*t17;
t48 = I4+t16+t27+t31;
t46 = t42.*2.0;
t47 = l_C_m4.*l_OA.*m4.*t20;
t49 = t42+t48;
t51 = I3+I4+I_motor+Ir+t16+t25+t26+t27+t28+t29+t30+t35+t39+t42;
t50 = t47+t49;
t52 = t43+t44+t45+t47+t51;
t53 = I2+I3+I4+t7+t8+t16+t22+t23+t24+t25+t26+t27+t28+t32+t33+t34+t36+t37+t38+t39+t40+t41+t43+t44+t45+t46+t47;
A = reshape([I1+I2+I3+I4+I_motor.*3.0+Ir.*3.0+t21+t22+t23+t24+t25+t26+t27+t28+t32.*2.0+t33.*2.0+t34.*2.0+t36+t37+t38+t39+t40+t43.*2.0+t44.*2.0+t45.*2.0+t46+t47.*2.0+m2.*t15+m3.*t15+m4.*t15+m_motor.*t15.*3.0+l_O_m1.^2.*m1+l_AB.*l_OA.*m_motor.*t2.*4.0,t53,t52,t50,t53,I2+I3+I4+t7+t8+t21+t22+t23+t24+t25+t26+t27+t28+t36+t37+t38+t39+t40+t46,t51,t49,t52,t51,I3+I4+I_motor+Ir+t21+t25+t26+t27+t28+t39,t48,t50,t49,t48,I4+t21+t27],[4,4]);
