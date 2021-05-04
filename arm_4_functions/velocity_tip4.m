function dqD = velocity_tip4(in1,in2)
%VELOCITY_TIP4
%    DQD = VELOCITY_TIP4(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    03-May-2021 15:42:44

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
dqD = [-dth3.*(t14+t16)-dth1.*(t10+t14+t16+l_OA.*sin(th1))-dth4.*t16-dth2.*(t10+t14+t16);dth3.*(t13+t15)+dth1.*(t9+t13+t15+l_OA.*cos(th1))+dth4.*t15+dth2.*(t9+t13+t15);dth1+dth2+dth3+dth4];
