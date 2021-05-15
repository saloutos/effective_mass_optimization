function Jv = jacobian_v_tip(in1,in2)
%JACOBIAN_V_TIP
%    JV = JACOBIAN_V_TIP(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    12-May-2021 18:28:58

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
t12 = -t9;
t13 = -t11;
Jv = reshape([t12+t13-l_OA.*sin(th1),t8+t10+l_OA.*cos(th1),t12+t13,t8+t10,t13,t10],[2,3]);
