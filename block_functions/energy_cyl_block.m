function E = energy_cyl_block(in1,in2)
%ENERGY_CYL_BLOCK
%    E = ENERGY_CYL_BLOCK(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    31-Mar-2021 23:01:39

I = in2(2,:);
dth = in1(6,:);
dx = in1(4,:);
dy = in1(5,:);
m = in2(1,:);
E = (I.*dth.^2)./2.0+(m.*(dx.^2+dy.^2))./2.0;
