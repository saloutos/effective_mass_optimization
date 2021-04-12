function dJ = jacobian_dot_tip4(in1,in2)
%JACOBIAN_DOT_TIP4
%    DJ = JACOBIAN_DOT_TIP4(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    29-Mar-2021 16:48:03

dth1 = in1(5,:);
dth2 = in1(6,:);
dth3 = in1(7,:);
dth4 = in1(8,:);
l_AB = in2(18,:);
l_BC = in2(19,:);
l_CD = in2(20,:);
l_OA = in2(17,:);
th1 = in1(1,:);
th2 = in1(2,:);
th3 = in1(3,:);
th4 = in1(4,:);
t2 = th1+th2;
t3 = cos(t2);
t4 = sin(t2);
t5 = t2+th3;
t6 = cos(t5);
t7 = sin(t5);
t8 = t5+th4;
t9 = l_AB.*t3;
t10 = l_AB.*t4;
t11 = cos(t8);
t12 = sin(t8);
t13 = l_BC.*t6;
t14 = l_BC.*t7;
t15 = l_CD.*t11;
t16 = l_CD.*t12;
t17 = dth4.*t15;
t18 = dth4.*t16;
t21 = t13+t15;
t22 = t14+t16;
t19 = -t17;
t20 = -t18;
t23 = dth3.*t21;
t24 = dth3.*t22;
t27 = t9+t21;
t28 = t10+t22;
t25 = -t23;
t26 = -t24;
t29 = dth2.*t28;
t30 = dth2.*t27;
t31 = -t29;
t32 = -t30;
dJ = reshape([t19+t25+t32-dth1.*(t27+l_OA.*cos(th1)),t20+t26+t31-dth1.*(t28+l_OA.*sin(th1)),0.0,t19+t25+t32-dth1.*t27,t20+t26+t31-dth1.*t28,0.0,t19+t25-dth1.*t21-dth2.*t21,t20+t26-dth1.*t22-dth2.*t22,0.0,t19-dth1.*t15-dth2.*t15-dth3.*t15,t20-dth1.*t16-dth2.*t16-dth3.*t16,0.0],[3,4]);