function Corr_Joint_Sp = Corr_arm_joint4(in1,in2)
%CORR_ARM_JOINT4
%    CORR_JOINT_SP = CORR_ARM_JOINT4(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    03-May-2021 15:42:45

dth1 = in1(5,:);
dth2 = in1(6,:);
dth3 = in1(7,:);
dth4 = in1(8,:);
l_AB = in2(18,:);
l_A_m2 = in2(14,:);
l_BC = in2(19,:);
l_B_m3 = in2(15,:);
l_C_m4 = in2(16,:);
l_OA = in2(17,:);
m2 = in2(2,:);
m3 = in2(3,:);
m4 = in2(4,:);
m_motor = in2(5,:);
th2 = in1(2,:);
th3 = in1(3,:);
th4 = in1(4,:);
t2 = sin(th2);
t3 = sin(th3);
t4 = sin(th4);
t5 = th2+th3;
t6 = th3+th4;
t7 = dth1.^2;
t8 = dth2.^2;
t9 = dth3.^2;
t10 = dth4.^2;
t11 = sin(t5);
t12 = sin(t6);
t13 = t5+th4;
t15 = l_AB.*l_BC.*m4.*t3.*t9;
t16 = l_AB.*l_B_m3.*m3.*t3.*t9;
t17 = l_BC.*l_C_m4.*m4.*t4.*t10;
t18 = l_AB.*l_BC.*m_motor.*t3.*t9;
t19 = dth1.*dth3.*l_AB.*l_BC.*m4.*t3.*2.0;
t20 = dth2.*dth3.*l_AB.*l_BC.*m4.*t3.*2.0;
t21 = dth1.*dth3.*l_AB.*l_B_m3.*m3.*t3.*2.0;
t22 = dth2.*dth3.*l_AB.*l_B_m3.*m3.*t3.*2.0;
t23 = dth1.*dth4.*l_BC.*l_C_m4.*m4.*t4.*2.0;
t24 = dth2.*dth4.*l_BC.*l_C_m4.*m4.*t4.*2.0;
t25 = dth3.*dth4.*l_BC.*l_C_m4.*m4.*t4.*2.0;
t26 = dth1.*dth3.*l_AB.*l_BC.*m_motor.*t3.*2.0;
t27 = dth2.*dth3.*l_AB.*l_BC.*m_motor.*t3.*2.0;
t14 = sin(t13);
t28 = -t19;
t29 = -t20;
t30 = -t21;
t31 = -t22;
t32 = -t23;
t33 = -t24;
t34 = -t25;
t35 = -t26;
t36 = -t27;
t37 = l_AB.*l_C_m4.*m4.*t9.*t12;
t38 = l_AB.*l_C_m4.*m4.*t10.*t12;
t39 = l_BC.*l_OA.*m4.*t7.*t11;
t40 = l_B_m3.*l_OA.*m3.*t7.*t11;
t41 = l_BC.*l_OA.*m_motor.*t7.*t11;
t42 = dth1.*dth3.*l_AB.*l_C_m4.*m4.*t12.*2.0;
t43 = dth1.*dth4.*l_AB.*l_C_m4.*m4.*t12.*2.0;
t44 = dth2.*dth3.*l_AB.*l_C_m4.*m4.*t12.*2.0;
t45 = dth2.*dth4.*l_AB.*l_C_m4.*m4.*t12.*2.0;
t46 = dth3.*dth4.*l_AB.*l_C_m4.*m4.*t12.*2.0;
t47 = -t15;
t48 = -t16;
t49 = -t17;
t50 = -t18;
t51 = -t42;
t52 = -t43;
t53 = -t44;
t54 = -t45;
t55 = -t46;
t56 = l_C_m4.*l_OA.*m4.*t7.*t14;
t57 = -t37;
t58 = -t38;
Corr_Joint_Sp = [t28+t29+t30+t31+t32+t33+t34+t35+t36+t47+t48+t49+t50+t51+t52+t53+t54+t55+t57+t58-l_AB.*l_OA.*m3.*t2.*t8-l_AB.*l_OA.*m4.*t2.*t8-l_A_m2.*l_OA.*m2.*t2.*t8-l_BC.*l_OA.*m4.*t8.*t11-l_BC.*l_OA.*m4.*t9.*t11-l_B_m3.*l_OA.*m3.*t8.*t11-l_B_m3.*l_OA.*m3.*t9.*t11-l_C_m4.*l_OA.*m4.*t8.*t14-l_C_m4.*l_OA.*m4.*t9.*t14-l_C_m4.*l_OA.*m4.*t10.*t14-l_AB.*l_OA.*m_motor.*t2.*t8.*2.0-l_BC.*l_OA.*m_motor.*t8.*t11-l_BC.*l_OA.*m_motor.*t9.*t11-dth1.*dth2.*l_AB.*l_OA.*m3.*t2.*2.0-dth1.*dth2.*l_AB.*l_OA.*m4.*t2.*2.0-dth1.*dth2.*l_A_m2.*l_OA.*m2.*t2.*2.0-dth1.*dth2.*l_BC.*l_OA.*m4.*t11.*2.0-dth1.*dth3.*l_BC.*l_OA.*m4.*t11.*2.0-dth2.*dth3.*l_BC.*l_OA.*m4.*t11.*2.0-dth1.*dth2.*l_B_m3.*l_OA.*m3.*t11.*2.0-dth1.*dth3.*l_B_m3.*l_OA.*m3.*t11.*2.0-dth2.*dth3.*l_B_m3.*l_OA.*m3.*t11.*2.0-dth1.*dth2.*l_C_m4.*l_OA.*m4.*t14.*2.0-dth1.*dth3.*l_C_m4.*l_OA.*m4.*t14.*2.0-dth1.*dth4.*l_C_m4.*l_OA.*m4.*t14.*2.0-dth2.*dth3.*l_C_m4.*l_OA.*m4.*t14.*2.0-dth2.*dth4.*l_C_m4.*l_OA.*m4.*t14.*2.0-dth3.*dth4.*l_C_m4.*l_OA.*m4.*t14.*2.0-dth1.*dth2.*l_AB.*l_OA.*m_motor.*t2.*4.0-dth1.*dth2.*l_BC.*l_OA.*m_motor.*t11.*2.0-dth1.*dth3.*l_BC.*l_OA.*m_motor.*t11.*2.0-dth2.*dth3.*l_BC.*l_OA.*m_motor.*t11.*2.0;t28+t29+t30+t31+t32+t33+t34+t35+t36+t39+t40+t41+t47+t48+t49+t50+t51+t52+t53+t54+t55+t56+t57+t58+l_AB.*l_OA.*m3.*t2.*t7+l_AB.*l_OA.*m4.*t2.*t7+l_A_m2.*l_OA.*m2.*t2.*t7+l_AB.*l_OA.*m_motor.*t2.*t7.*2.0;t32+t33+t34+t39+t40+t41+t49+t56+l_AB.*l_BC.*m4.*t3.*t7+l_AB.*l_BC.*m4.*t3.*t8+l_AB.*l_B_m3.*m3.*t3.*t7+l_AB.*l_B_m3.*m3.*t3.*t8+l_AB.*l_C_m4.*m4.*t7.*t12+l_AB.*l_C_m4.*m4.*t8.*t12+l_AB.*l_BC.*m_motor.*t3.*t7+l_AB.*l_BC.*m_motor.*t3.*t8+dth1.*dth2.*l_AB.*l_BC.*m4.*t3.*2.0+dth1.*dth2.*l_AB.*l_B_m3.*m3.*t3.*2.0+dth1.*dth2.*l_AB.*l_C_m4.*m4.*t12.*2.0+dth1.*dth2.*l_AB.*l_BC.*m_motor.*t3.*2.0;l_C_m4.*m4.*(l_AB.*t7.*t12+l_AB.*t8.*t12+l_BC.*t4.*t7+l_BC.*t4.*t8+l_BC.*t4.*t9+l_OA.*t7.*t14+dth1.*dth2.*l_AB.*t12.*2.0+dth1.*dth2.*l_BC.*t4.*2.0+dth1.*dth3.*l_BC.*t4.*2.0+dth2.*dth3.*l_BC.*t4.*2.0)];
