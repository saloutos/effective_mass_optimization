function keypoints = keypoints_arm4(in1,in2)
%KEYPOINTS_ARM4
%    KEYPOINTS = KEYPOINTS_ARM4(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    03-May-2021 15:42:45

l_AB = in2(18,:);
l_BC = in2(19,:);
l_CD = in2(20,:);
l_OA = in2(17,:);
th1 = in1(1,:);
th2 = in1(2,:);
th3 = in1(3,:);
th4 = in1(4,:);
t2 = cos(th1);
t3 = sin(th1);
t4 = th1+th2;
t5 = l_OA.*t2;
t6 = cos(t4);
t7 = l_OA.*t3;
t8 = sin(t4);
t9 = t4+th3;
t10 = cos(t9);
t11 = sin(t9);
t12 = t9+th4;
t13 = l_AB.*t6;
t14 = l_AB.*t8;
t15 = l_BC.*t10;
t16 = l_BC.*t11;
keypoints = reshape([t5,t7,t5+t13,t7+t14,t5+t13+t15,t7+t14+t16,t5+t13+t15+l_CD.*cos(t12),t7+t14+t16+l_CD.*sin(t12)],[2,4]);
