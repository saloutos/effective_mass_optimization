function E = energy_arm_rotors(in1,in2)
%ENERGY_ARM_ROTORS
%    E = ENERGY_ARM_ROTORS(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    30-May-2021 10:56:53

I1 = in2(4,:);
I2 = in2(5,:);
I3 = in2(6,:);
Ir1 = in2(7,:);
Ir2 = in2(8,:);
Ir3 = in2(9,:);
dth1 = in1(4,:);
dth2 = in1(5,:);
dth3 = in1(6,:);
l_AB = in2(14,:);
l_A_m2 = in2(11,:);
l_B_m3 = in2(12,:);
l_OA = in2(13,:);
l_O_m1 = in2(10,:);
m1 = in2(1,:);
m2 = in2(2,:);
m3 = in2(3,:);
th1 = in1(1,:);
th2 = in1(2,:);
th3 = in1(3,:);
t2 = cos(th1);
t3 = sin(th1);
t4 = dth1+dth2;
t5 = th1+th2;
t6 = dth1.^2;
t7 = l_O_m1.^2;
t8 = dth3+t4;
t9 = l_OA.*t2;
t10 = cos(t5);
t11 = l_OA.*t3;
t12 = sin(t5);
t13 = t5+th3;
t16 = t4.^2;
t14 = cos(t13);
t15 = sin(t13);
t17 = l_AB.*t10;
t18 = l_AB.*t12;
t19 = t8.^2;
t20 = l_B_m3.*t14;
t21 = l_B_m3.*t15;
E = (I1.*t6)./2.0+(I2.*t16)./2.0+(I3.*t19)./2.0+(Ir1.*t6)./2.0+(Ir2.*t16)./2.0+(Ir3.*t19)./2.0+(m3.*((dth2.*(t17+t20)+dth3.*t20+dth1.*(t9+t17+t20)).^2+(dth2.*(t18+t21)+dth3.*t21+dth1.*(t11+t18+t21)).^2))./2.0+(m1.*(t2.^2.*t6.*t7+t3.^2.*t6.*t7))./2.0+(m2.*((dth1.*(t9+l_A_m2.*t10)+dth2.*l_A_m2.*t10).^2+(dth1.*(t11+l_A_m2.*t12)+dth2.*l_A_m2.*t12).^2))./2.0;