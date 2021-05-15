function dJ = jacobian_dot_tip(in1,in2)
%JACOBIAN_DOT_TIP
%    DJ = JACOBIAN_DOT_TIP(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    12-May-2021 19:41:07

dth1 = in1(4,:);
dth2 = in1(5,:);
dth3 = in1(6,:);
l_AB = in2(15,:);
l_BC = in2(16,:);
l_OA = in2(14,:);
th1 = in1(1,:);
th2 = in1(2,:);
th3 = in1(3,:);
t2 = th1+th2;
t3 = cos(t2);
t4 = sin(t2);
t5 = t2+th3;
t6 = cos(t5);
t7 = sin(t5);
t8 = l_AB.*t3;
t9 = l_AB.*t4;
t10 = l_BC.*t6;
t11 = l_BC.*t7;
t12 = dth3.*t10;
t13 = dth3.*t11;
t16 = t8+t10;
t17 = t9+t11;
t14 = -t12;
t15 = -t13;
t18 = dth2.*t16;
t19 = dth2.*t17;
t20 = -t18;
t21 = -t19;
dJ = reshape([t14+t20-dth1.*(t16+l_OA.*cos(th1)),t15+t21-dth1.*(t17+l_OA.*sin(th1)),0.0,t14+t20-dth1.*t16,t15+t21-dth1.*t17,0.0,t14-dth1.*t10-dth2.*t10,t15-dth1.*t11-dth2.*t11,0.0],[3,3]);
