function A = A_cyl_block(in1,in2)
%A_CYL_BLOCK
%    A = A_CYL_BLOCK(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    31-Mar-2021 23:01:38

I = in2(2,:);
m = in2(1,:);
A = reshape([m,0.0,0.0,0.0,m,0.0,0.0,0.0,I],[3,3]);
