function qC = position_tip_arm_rotors(in1,in2)
%POSITION_TIP_ARM_ROTORS
%    QC = POSITION_TIP_ARM_ROTORS(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    30-May-2021 10:56:53

l_AB = in2(14,:);
l_BC = in2(15,:);
l_OA = in2(13,:);
th1 = in1(1,:);
th2 = in1(2,:);
th3 = in1(3,:);
t2 = th1+th2;
t3 = t2+th3;
qC = [l_AB.*cos(t2)+l_BC.*cos(t3)+l_OA.*cos(th1);l_AB.*sin(t2)+l_BC.*sin(t3)+l_OA.*sin(th1);t3];
