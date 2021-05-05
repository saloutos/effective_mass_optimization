function Jv = J_arm_floating(in1,in2)
%J_ARM_FLOATING
%    JV = J_ARM_FLOATING(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    04-May-2021 16:15:49

l_AB = in2(15,:);
l_BC = in2(16,:);
l_OA = in2(14,:);
th1 = in1(4,:);
th2 = in1(5,:);
th3 = in1(6,:);
thb = in1(3,:);
t2 = th1+thb;
t3 = cos(t2);
t4 = sin(t2);
t5 = t2+th2;
t6 = cos(t5);
t7 = sin(t5);
t8 = t5+th3;
t9 = l_OA.*t3;
t10 = l_OA.*t4;
t11 = cos(t8);
t12 = sin(t8);
t13 = l_AB.*t6;
t14 = l_AB.*t7;
t17 = -t10;
t15 = l_BC.*t11;
t16 = l_BC.*t12;
t18 = -t14;
t19 = -t16;
t20 = t9+t13+t15;
t21 = t17+t18+t19;
Jv = reshape([1.0,0.0,0.0,1.0,t21,t20,t21,t20,t18+t19,t13+t15,t19,t15],[2,6]);
