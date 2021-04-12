function J = jacobian_tip4(in1,in2)
%JACOBIAN_TIP4
%    J = JACOBIAN_TIP4(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    29-Mar-2021 16:48:03

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
t17 = -t10;
t15 = l_CD.*t11;
t16 = l_CD.*t12;
t18 = -t14;
t19 = -t16;
J = reshape([t17+t18+t19-l_OA.*sin(th1),t9+t13+t15+l_OA.*cos(th1),1.0,t17+t18+t19,t9+t13+t15,1.0,t18+t19,t13+t15,1.0,t19,t15,1.0],[3,4]);
