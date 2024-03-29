function A = H_arm_floating(in1,in2)
%H_ARM_FLOATING
%    A = H_ARM_FLOATING(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    04-May-2021 16:15:48

I1 = in2(5,:);
I2 = in2(6,:);
I3 = in2(7,:);
I_motor = in2(8,:);
Ib = in2(19,:);
Ir = in2(9,:);
N = in2(10,:);
l_AB = in2(15,:);
l_A_m2 = in2(12,:);
l_B_m3 = in2(13,:);
l_OA = in2(14,:);
l_O_m1 = in2(11,:);
m1 = in2(1,:);
m2 = in2(2,:);
m3 = in2(3,:);
m_motor = in2(4,:);
mb = in2(18,:);
th1 = in1(4,:);
th2 = in1(5,:);
th3 = in1(6,:);
thb = in1(3,:);
t2 = th1+thb;
t3 = I_motor.*2.0;
t4 = N.^2;
t5 = l_AB.^2;
t6 = l_A_m2.^2;
t7 = l_B_m3.^2;
t8 = l_OA.^2;
t9 = l_O_m1.^2;
t10 = m_motor.*3.0;
t11 = Ir.*N;
t12 = cos(t2);
t13 = sin(t2);
t14 = t2+th2;
t18 = Ir.*t4;
t43 = m1+m2+m3+mb+t10;
t15 = cos(t14);
t16 = sin(t14);
t17 = t14+th3;
t19 = t12.^2;
t20 = t13.^2;
t21 = l_OA.*t12;
t22 = l_OA.*t13;
t29 = l_O_m1.*m1.*t12;
t33 = l_O_m1.*m1.*t13;
t23 = cos(t17);
t24 = sin(t17);
t25 = t21.*2.0;
t26 = t15.^2;
t27 = t22.*2.0;
t28 = t16.^2;
t30 = m_motor.*t21;
t31 = l_AB.*t15;
t32 = l_A_m2.*t15;
t34 = m_motor.*t22;
t35 = l_AB.*t16;
t36 = l_A_m2.*t16;
t52 = -t33;
t56 = t8.*t19.*2.0;
t57 = t9.*t19.*2.0;
t58 = t8.*t20.*2.0;
t59 = t9.*t20.*2.0;
t37 = m2.*t32;
t38 = m_motor.*t31;
t39 = l_B_m3.*t23;
t40 = m2.*t36;
t41 = m_motor.*t35;
t42 = l_B_m3.*t24;
t44 = t31.*2.0;
t45 = t32.*2.0;
t46 = t35.*2.0;
t47 = t36.*2.0;
t53 = -t34;
t61 = t21+t31;
t62 = t21+t32;
t63 = t22+t35;
t64 = t22+t36;
t95 = t56+t58;
t96 = t57+t59;
t48 = t39.*2.0;
t49 = t42.*2.0;
t50 = m3.*t39;
t51 = m3.*t42;
t54 = -t40;
t55 = -t41;
t65 = t61.^2;
t66 = t62.^2;
t67 = t63.^2;
t68 = t64.^2;
t69 = t31+t39;
t70 = t25+t44;
t71 = t25+t45;
t72 = t35+t42;
t73 = t27+t46;
t74 = t27+t47;
t85 = t39+t61;
t86 = t42+t63;
t87 = t44.*t61;
t88 = t45.*t62;
t91 = t46.*t63;
t92 = t47.*t64;
t106 = (m1.*t96)./2.0;
t107 = (m_motor.*t95)./2.0;
t60 = -t51;
t75 = t65.*2.0;
t76 = t66.*2.0;
t77 = t67.*2.0;
t78 = t68.*2.0;
t79 = t44+t48;
t80 = t46+t49;
t81 = (m2.*t71)./2.0;
t82 = (m_motor.*t70)./2.0;
t83 = (m2.*t74)./2.0;
t84 = (m_motor.*t73)./2.0;
t97 = t85.^2;
t99 = t86.^2;
t100 = t48.*t69;
t101 = t48+t70;
t102 = t49.*t72;
t103 = t49+t73;
t111 = t48.*t85;
t112 = t49.*t86;
t115 = t69.*t85.*2.0;
t116 = t72.*t86.*2.0;
t121 = t87+t91;
t122 = t88+t92;
t89 = -t83;
t90 = -t84;
t93 = (m3.*t79)./2.0;
t94 = (m3.*t80)./2.0;
t104 = t99.*2.0;
t105 = t97.*2.0;
t108 = (m3.*t103)./2.0;
t109 = (m3.*t101)./2.0;
t113 = t75+t77;
t114 = t76+t78;
t123 = (m2.*t122)./2.0;
t124 = (m_motor.*t121)./2.0;
t125 = t100+t102;
t130 = t111+t112;
t134 = t115+t116;
t98 = -t94;
t110 = -t108;
t117 = t37+t38+t93;
t118 = (m2.*t114)./2.0;
t119 = (m_motor.*t113)./2.0;
t126 = t104+t105;
t127 = (m3.*t125)./2.0;
t131 = (m3.*t130)./2.0;
t135 = (m3.*t134)./2.0;
t136 = t29+t30+t81+t82+t109;
t120 = t54+t55+t98;
t128 = (m3.*t126)./2.0;
t129 = I3+t11+t127;
t132 = I3+t131;
t137 = t52+t53+t89+t90+t110;
t138 = I2+I3+I_motor+t123+t124+t135;
t133 = t11+t132;
t139 = Ir+t11+t138;
t140 = I1+I2+I3+t3+t106+t107+t118+t119+t128;
A = reshape([t43,0.0,t137,t137,t120,t60,0.0,t43,t136,t136,t117,t50,t137,t136,I1+I2+I3+I_motor.*3.0+Ib+t106+t107+t118+t119+t128,t140,t138,t132,t137,t136,t140,Ir.*2.0+t18+t140,t139,t133,t120,t117,t138,t139,I2+I3+I_motor+Ir+t18+(m3.*(t69.^2.*2.0+t72.^2.*2.0))./2.0+(m2.*(t6.*t26.*2.0+t6.*t28.*2.0))./2.0+(m_motor.*(t5.*t26.*2.0+t5.*t28.*2.0))./2.0,t129,t60,t50,t132,t133,t129,I3+t18+(m3.*(t7.*t23.^2.*2.0+t7.*t24.^2.*2.0))./2.0],[6,6]);
