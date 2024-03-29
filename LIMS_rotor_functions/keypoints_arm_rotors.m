function keypoints = keypoints_arm_rotors(in1,in2)
%KEYPOINTS_ARM_ROTORS
%    KEYPOINTS = KEYPOINTS_ARM_ROTORS(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    30-May-2021 10:56:54

l_AB = in2(14,:);
l_BC = in2(15,:);
l_OA = in2(13,:);
th1 = in1(1,:);
th2 = in1(2,:);
th3 = in1(3,:);
t2 = cos(th1);
t3 = sin(th1);
t4 = th1+th2;
t5 = l_OA.*t2;
t6 = cos(t4);
t7 = l_OA.*t3;
t8 = sin(t4);
t9 = t4+th3;
t10 = l_AB.*t6;
t11 = l_AB.*t8;
keypoints = reshape([t5,t7,t5+t10,t7+t11,t5+t10+l_BC.*cos(t9),t7+t11+l_BC.*sin(t9)],[2,3]);
