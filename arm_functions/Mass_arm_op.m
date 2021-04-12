function Mass_Op_Sp = Mass_arm_op(in1,in2)
%MASS_ARM_OP
%    MASS_OP_SP = MASS_ARM_OP(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    15-Mar-2021 22:37:37

I1 = in2(5,:);
I2 = in2(6,:);
I3 = in2(7,:);
I_motor = in2(8,:);
Ir = in2(9,:);
N = in2(10,:);
l_AB = in2(15,:);
l_A_m2 = in2(12,:);
l_BC = in2(16,:);
l_B_m3 = in2(13,:);
l_OA = in2(14,:);
l_O_m1 = in2(11,:);
m1 = in2(1,:);
m2 = in2(2,:);
m3 = in2(3,:);
m_motor = in2(4,:);
th1 = in1(1,:);
th2 = in1(2,:);
th3 = in1(3,:);
t2 = cos(th1);
t3 = cos(th2);
t4 = cos(th3);
t5 = sin(th1);
t6 = th1+th2;
t7 = th2+th3;
t8 = I_motor.*2.0;
t9 = Ir.*2.0;
t10 = N.^2;
t11 = l_AB.^2;
t12 = l_A_m2.^2;
t13 = l_B_m3.^2;
t14 = l_OA.^2;
t15 = l_O_m1.^2;
t16 = Ir.*N;
t17 = l_OA.*t2;
t18 = cos(t6);
t19 = cos(t7);
t20 = l_OA.*t5;
t21 = sin(t6);
t22 = t6+th3;
t25 = Ir.*t10;
t26 = m3.*t11;
t27 = m2.*t12;
t28 = m3.*t13;
t29 = m2.*t14;
t30 = m3.*t14;
t31 = m1.*t15;
t32 = m_motor.*t11;
t34 = l_AB.*l_B_m3.*m3.*t4;
t35 = l_AB.*l_OA.*m3.*t3;
t36 = l_A_m2.*l_OA.*m2.*t3;
t37 = l_AB.*l_OA.*m_motor.*t3;
t39 = m_motor.*t14.*2.0;
t23 = cos(t22);
t24 = sin(t22);
t33 = l_AB.*t18;
t38 = l_AB.*t21;
t40 = t34.*2.0;
t41 = t35.*2.0;
t42 = t36.*2.0;
t43 = t37.*2.0;
t44 = l_B_m3.*l_OA.*m3.*t19;
t47 = t18.*t20;
t48 = t17.*t21;
t67 = I3+t25+t28;
t70 = I3+t16+t28+t34;
t45 = t5.*t33;
t46 = t2.*t38;
t49 = t44.*2.0;
t50 = t20.*t33;
t51 = t17.*t38;
t52 = l_BC.*t5.*t23;
t53 = l_BC.*t2.*t24;
t54 = l_BC.*t20.*t23;
t55 = l_BC.*t17.*t24;
t57 = -t48;
t58 = t17+t33;
t59 = t20+t38;
t60 = l_BC.*t21.*t23;
t61 = l_BC.*t18.*t24;
t64 = l_BC.*t23.*t38;
t65 = l_BC.*t24.*t33;
t76 = t44+t70;
t81 = I2+I_motor+Ir+t26+t27+t32+t40+t67;
t94 = I2+I3+I_motor+Ir+t16+t26+t27+t28+t32+t35+t36+t37+t40+t44;
t56 = -t46;
t62 = -t51;
t63 = -t53;
t66 = -t55;
t68 = -t61;
t69 = -t65;
t72 = t47+t57;
t103 = I1+I2+I3+t8+t9+t26+t27+t28+t29+t30+t31+t32+t39+t40+t41+t42+t43+t49;
t71 = t45+t56;
t73 = t50+t62;
t75 = 1.0./t72;
t78 = t60+t68;
t85 = t54+t64+t66+t69;
t74 = 1.0./t71;
t77 = 1.0./t73;
t82 = t52+t63+t71;
t88 = t18.*t75.*t76;
t89 = t21.*t75.*t76;
t95 = t75.*t76.*t78;
t104 = t18.*t75.*t94;
t105 = t21.*t75.*t94;
t112 = t75.*t78.*t94;
t115 = t18.*t75.*t103;
t116 = t21.*t75.*t103;
t118 = t75.*t78.*t103;
t79 = t2.*t67.*t74;
t80 = t5.*t67.*t74;
t83 = t2.*t70.*t74;
t84 = t5.*t70.*t74;
t86 = t2.*t74.*t76;
t87 = t5.*t74.*t76;
t90 = t59.*t70.*t77;
t91 = t58.*t70.*t77;
t96 = t59.*t77.*t81;
t97 = t58.*t77.*t81;
t98 = t67.*t74.*t82;
t101 = t70.*t74.*t82;
t102 = t74.*t76.*t82;
t106 = t70.*t77.*t85;
t108 = t58.*t77.*t94;
t109 = t59.*t77.*t94;
t113 = t77.*t81.*t85;
t117 = t77.*t85.*t94;
t92 = -t90;
t93 = -t91;
t99 = -t96;
t100 = -t97;
t107 = -t106;
t110 = -t108;
t111 = -t109;
t114 = -t113;
t119 = -t117;
t120 = t79+t88+t93;
t121 = t80+t89+t92;
t122 = t83+t100+t104;
t123 = t84+t99+t105;
t124 = t95+t98+t107;
t125 = t86+t110+t115;
t126 = t87+t111+t116;
t127 = t101+t112+t114;
t128 = t102+t118+t119;
Mass_Op_Sp = reshape([t2.*t74.*t120+t18.*t75.*t125-t58.*t77.*t122,t2.*t74.*t121+t18.*t75.*t126-t58.*t77.*t123,-t2.*t74.*t124-t18.*t75.*t128+t58.*t77.*t127,t5.*t74.*t120+t21.*t75.*t125-t59.*t77.*t122,t5.*t74.*t121+t21.*t75.*t126-t59.*t77.*t123,-t5.*t74.*t124-t21.*t75.*t128+t59.*t77.*t127,-t74.*t82.*t120-t75.*t78.*t125+t77.*t85.*t122,-t74.*t82.*t121-t75.*t78.*t126+t77.*t85.*t123,t74.*t82.*t124+t75.*t78.*t128-t77.*t85.*t127],[3,3]);
