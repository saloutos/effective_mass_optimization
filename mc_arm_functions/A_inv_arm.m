function Mass_Joint_Sp_inv = A_inv_arm(in1,in2)
%A_INV_ARM
%    MASS_JOINT_SP_INV = A_INV_ARM(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    12-May-2021 18:28:57

I1 = in2(5,:);
I2 = in2(6,:);
I3 = in2(7,:);
I_motor = in2(8,:);
Ir = in2(9,:);
N = in2(10,:);
l_AB = in2(15,:);
l_A_m2 = in2(12,:);
l_B_m3 = in2(13,:);
l_OA = in2(14,:);
l_O_m1 = in2(11,:);
m1 = in2(1,:);
m2 = in2(2,:);
m3 = in2(3,:);
th2 = in1(2,:);
th3 = in1(3,:);
t2 = cos(th2);
t3 = cos(th3);
t4 = th2+th3;
t5 = Ir.^2;
t6 = N.^2;
t8 = l_AB.^2;
t9 = l_A_m2.^2;
t10 = l_B_m3.^2;
t11 = l_OA.^2;
t12 = l_O_m1.^2;
t13 = m2.^2;
t14 = m3.^2;
t15 = m3.^3;
t20 = I1.*9.0e+6;
t21 = I_motor.*9.0e+6;
t22 = Ir.*9.0e+6;
t26 = Ir.*N.*1.4e+7;
t27 = I2.*2.25e+9;
t29 = Ir.*N.*1.26e+8;
t35 = Ir.*N.*3.5e+9;
t36 = I1.*I3.*5.625e+10;
t37 = I3.*I_motor.*5.625e+10;
t38 = I3.*Ir.*5.625e+10;
t39 = I3.*1.00894978125e+11;
t43 = I1.*2.88957537e+11;
t44 = I_motor.*2.88957537e+11;
t45 = Ir.*2.88957537e+11;
t55 = I2.*I3.*7.03125e+12;
t57 = Ir.*N.*8.98979004e+11;
t60 = I3.*Ir.*N.*8.75e+10;
t64 = I3.*Ir.*N.*7.875e+11;
t67 = I1.*I2.*2.025e+13;
t68 = I2.*I_motor.*2.025e+13;
t69 = I2.*Ir.*2.025e+13;
t82 = I2.*Ir.*N.*6.3e+13;
t85 = I3.*Ir.*N.*1.09375e+13;
t102 = I1.*I3.*9.08054803125e+14;
t103 = I3.*I_motor.*9.08054803125e+14;
t104 = I3.*Ir.*9.08054803125e+14;
t119 = I3.*Ir.*N.*2.8250593875e+15;
t139 = I1.*I2.*I3.*6.328125e+16;
t140 = I2.*I3.*I_motor.*6.328125e+16;
t141 = I2.*I3.*Ir.*6.328125e+16;
t199 = I2.*I3.*Ir.*N.*1.96875e+17;
t7 = t6.^2;
t16 = t2.^2;
t17 = t3.^2;
t18 = cos(t4);
t23 = t11.*2.583e+6;
t24 = l_OA.*t2.*5.73426e+5;
t25 = l_OA.*t2.*5.160834e+6;
t28 = -t26;
t30 = t6.*t22;
t31 = m2.*t11.*9.0e+6;
t32 = m3.*t11.*9.0e+6;
t33 = m1.*t12.*9.0e+6;
t34 = l_OA.*t2.*1.433565e+8;
t40 = Ir.*t6.*1.96e+8;
t41 = l_AB.*l_OA.*m3.*t2.*9.0e+6;
t42 = l_A_m2.*l_OA.*m2.*t2.*9.0e+6;
t46 = m3.*t8.*2.25e+9;
t47 = m2.*t9.*2.25e+9;
t49 = l_AB.*l_OA.*m3.*t2.*8.1e+7;
t50 = l_A_m2.*l_OA.*m2.*t2.*8.1e+7;
t51 = I3.*l_OA.*t2.*3.5839125e+9;
t52 = l_AB.*l_B_m3.*m3.*t3.*2.25e+9;
t53 = l_AB.*l_OA.*m3.*t2.*2.25e+9;
t54 = l_A_m2.*l_OA.*m2.*t2.*2.25e+9;
t56 = I3.*t11.*1.614375e+10;
t58 = I3.*l_OA.*t2.*3.22552125e+10;
t59 = t11.*8.2930813119e+10;
t61 = m3.*t10.*1.00894978125e+11;
t62 = I3.*l_OA.*t2.*4.479890625e+11;
t63 = -t57;
t65 = -t60;
t66 = t6.*t38;
t70 = m2.*t11.*2.88957537e+11;
t71 = m3.*t11.*2.88957537e+11;
t72 = m1.*t12.*2.88957537e+11;
t73 = I1.*m3.*t10.*5.625e+10;
t74 = I3.*m2.*t11.*5.625e+10;
t75 = I3.*m3.*t11.*5.625e+10;
t76 = I3.*m1.*t12.*5.625e+10;
t77 = I_motor.*m3.*t10.*5.625e+10;
t78 = Ir.*m3.*t10.*5.625e+10;
t79 = Ir.*t6.*9.88163429e+11;
t83 = I2.*t11.*5.81175e+12;
t84 = l_OA.*m3.*t2.*t10.*3.5839125e+9;
t86 = I1.*l_AB.*l_B_m3.*m3.*t3.*5.625e+10;
t87 = I3.*l_AB.*l_OA.*m3.*t2.*5.625e+10;
t88 = I3.*l_A_m2.*l_OA.*m2.*t2.*5.625e+10;
t89 = I_motor.*l_AB.*l_B_m3.*m3.*t3.*5.625e+10;
t90 = Ir.*l_AB.*l_B_m3.*m3.*t3.*5.625e+10;
t92 = m3.*t10.*t11.*1.614375e+10;
t93 = l_OA.*m3.*t2.*t10.*3.22552125e+10;
t95 = Ir.*N.*m3.*t10.*8.75e+10;
t96 = -t82;
t97 = l_AB.*l_B_m3.*l_OA.*m3.*t2.*t3.*3.5839125e+9;
t98 = I3.*l_AB.*l_OA.*m3.*t2.*5.0625e+11;
t99 = I3.*l_A_m2.*l_OA.*m2.*t2.*5.0625e+11;
t100 = Ir.*N.*l_OA.*t2.*4.013982e+12;
t101 = I3.*Ir.*t6.*1.225e+12;
t105 = I3.*m3.*t8.*7.03125e+12;
t106 = I3.*m2.*t9.*7.03125e+12;
t107 = I2.*m3.*t10.*7.03125e+12;
t109 = l_AB.*l_B_m3.*m3.*t3.*t11.*1.614375e+10;
t110 = I2.*Ir.*t6.*6.925e+13;
t111 = l_OA.*m3.*t2.*t10.*4.479890625e+11;
t112 = l_AB.*l_B_m3.*l_OA.*m3.*t2.*t3.*3.22552125e+10;
t113 = I1.*Ir.*t6.*4.9e+13;
t114 = I_motor.*Ir.*t6.*4.9e+13;
t115 = Ir.*N.*l_AB.*l_B_m3.*m3.*t3.*8.75e+10;
t117 = Ir.*N.*m3.*t10.*7.875e+11;
t122 = I1.*m3.*t8.*2.025e+13;
t123 = I1.*m2.*t9.*2.025e+13;
t124 = I2.*m2.*t11.*2.025e+13;
t125 = I2.*m3.*t11.*2.025e+13;
t126 = I2.*m1.*t12.*2.025e+13;
t127 = I_motor.*m3.*t8.*2.025e+13;
t128 = I_motor.*m2.*t9.*2.025e+13;
t129 = Ir.*m3.*t8.*2.025e+13;
t130 = Ir.*m2.*t9.*2.025e+13;
t133 = I3.*t11.*2.60611728496875e+14;
t134 = I3.*l_AB.*l_OA.*m3.*t2.*7.03125e+12;
t135 = I3.*l_A_m2.*l_OA.*m2.*t2.*7.03125e+12;
t136 = m2.*m3.*t10.*t11.*5.625e+10;
t137 = m1.*m3.*t10.*t12.*5.625e+10;
t138 = Ir.*l_OA.*t2.*t6.*6.243972e+12;
t143 = Ir.*N.*l_AB.*l_B_m3.*m3.*t3.*7.875e+11;
t146 = Ir.*N.*m3.*t8.*6.3e+13;
t147 = Ir.*N.*m2.*t9.*6.3e+13;
t149 = m3.*t8.*t11.*5.81175e+12;
t150 = m2.*t9.*t11.*5.81175e+12;
t152 = I1.*l_AB.*l_B_m3.*m3.*t3.*2.025e+13;
t153 = I_motor.*l_AB.*l_B_m3.*m3.*t3.*2.025e+13;
t154 = Ir.*l_AB.*l_B_m3.*m3.*t3.*2.025e+13;
t155 = -t119;
t157 = Ir.*N.*m3.*t10.*1.09375e+13;
t158 = l_AB.*l_B_m3.*m2.*m3.*t3.*t11.*5.625e+10;
t159 = l_AB.*l_B_m3.*m1.*m3.*t3.*t12.*5.625e+10;
t160 = l_A_m2.*l_OA.*m2.*m3.*t2.*t10.*5.625e+10;
t161 = t10.*t11.*t14.*5.625e+10;
t162 = I3.*Ir.*t6.*3.105323215625e+15;
t165 = Ir.*N.*l_AB.*l_B_m3.*m3.*t3.*6.3e+13;
t166 = Ir.*N.*l_AB.*l_OA.*m3.*t2.*6.3e+13;
t167 = Ir.*N.*l_A_m2.*l_OA.*m2.*t2.*6.3e+13;
t168 = Ir.*t6.*t11.*1.4063e+13;
t169 = l_AB.*l_B_m3.*m3.*t3.*t11.*5.81175e+12;
t173 = l_A_m2.*l_OA.*m2.*m3.*t2.*t10.*5.0625e+11;
t174 = l_AB.*l_A_m2.*l_B_m3.*l_OA.*m2.*m3.*t2.*t3.*5.625e+10;
t175 = l_AB.*l_OA.*t2.*t10.*t14.*5.625e+10;
t176 = l_AB.*l_B_m3.*t3.*t11.*t14.*5.625e+10;
t177 = Ir.*m3.*t6.*t10.*1.225e+12;
t178 = I1.*m3.*t10.*9.08054803125e+14;
t179 = I3.*m2.*t11.*9.08054803125e+14;
t180 = I3.*m3.*t11.*9.08054803125e+14;
t181 = I3.*m1.*t12.*9.08054803125e+14;
t182 = I_motor.*m3.*t10.*9.08054803125e+14;
t183 = Ir.*m3.*t10.*9.08054803125e+14;
t184 = m2.*m3.*t9.*t10.*7.03125e+12;
t187 = Ir.*m3.*t6.*t8.*6.925e+13;
t188 = Ir.*m2.*t6.*t9.*6.925e+13;
t189 = Ir.*m2.*t6.*t11.*4.9e+13;
t190 = Ir.*m3.*t6.*t11.*4.9e+13;
t191 = Ir.*m1.*t6.*t12.*4.9e+13;
t194 = I2.*I3.*t11.*1.816171875e+16;
t198 = l_AB.*l_A_m2.*l_B_m3.*l_OA.*m2.*m3.*t2.*t3.*5.0625e+11;
t200 = l_AB.*l_OA.*t2.*t10.*t14.*5.0625e+11;
t201 = Ir.*N.*m3.*t10.*2.8250593875e+15;
t203 = m2.*m3.*t8.*t11.*2.025e+13;
t204 = m1.*m3.*t8.*t12.*2.025e+13;
t205 = m2.*m3.*t9.*t11.*2.025e+13;
t206 = m1.*m2.*t9.*t12.*2.025e+13;
t207 = Ir.*l_AB.*l_B_m3.*m3.*t3.*t6.*1.225e+12;
t209 = m3.*t10.*t11.*2.60611728496875e+14;
t210 = l_A_m2.*l_OA.*m2.*m3.*t2.*t10.*7.03125e+12;
t211 = t8.*t10.*t14.*7.03125e+12;
t216 = I1.*I3.*m3.*t8.*6.328125e+16;
t217 = I1.*I3.*m2.*t9.*6.328125e+16;
t218 = I1.*I2.*m3.*t10.*6.328125e+16;
t219 = I2.*I3.*m2.*t11.*6.328125e+16;
t220 = I2.*I3.*m3.*t11.*6.328125e+16;
t221 = I2.*I3.*m1.*t12.*6.328125e+16;
t222 = I3.*I_motor.*m3.*t8.*6.328125e+16;
t223 = I3.*I_motor.*m2.*t9.*6.328125e+16;
t224 = I2.*I_motor.*m3.*t10.*6.328125e+16;
t225 = I3.*Ir.*m3.*t8.*6.328125e+16;
t226 = I3.*Ir.*m2.*t9.*6.328125e+16;
t227 = I2.*Ir.*m3.*t10.*6.328125e+16;
t229 = -t199;
t230 = l_AB.*l_B_m3.*m2.*m3.*t3.*t11.*2.025e+13;
t231 = l_AB.*l_B_m3.*m1.*m3.*t3.*t12.*2.025e+13;
t234 = t8.*t11.*t14.*2.025e+13;
t235 = t9.*t11.*t13.*2.025e+13;
t236 = l_B_m3.*l_OA.*t2.*t3.*t8.*t14.*5.625e+10;
t239 = l_AB.*l_OA.*t2.*t10.*t14.*7.03125e+12;
t240 = I3.*Ir.*N.*l_OA.*t2.*1.254369375e+16;
t241 = Ir.*l_AB.*l_B_m3.*m3.*t3.*t6.*6.925e+13;
t244 = I1.*I3.*Ir.*t6.*1.53125e+17;
t245 = I3.*I_motor.*Ir.*t6.*1.53125e+17;
t247 = l_B_m3.*l_OA.*t2.*t3.*t8.*t14.*5.0625e+11;
t248 = I2.*I3.*Ir.*t6.*2.1640625e+17;
t251 = l_AB.*l_B_m3.*t3.*t11.*t14.*2.025e+13;
t254 = m2.*m3.*t10.*t11.*9.08054803125e+14;
t255 = m1.*m3.*t10.*t12.*9.08054803125e+14;
t257 = Ir.*l_AB.*l_OA.*m3.*t2.*t6.*9.8e+13;
t258 = Ir.*l_A_m2.*l_OA.*m2.*t2.*t6.*9.8e+13;
t260 = Ir.*m3.*t6.*t10.*3.105323215625e+15;
t262 = I3.*m3.*t8.*t11.*1.816171875e+16;
t263 = I3.*m2.*t9.*t11.*1.816171875e+16;
t264 = I2.*m3.*t10.*t11.*1.816171875e+16;
t268 = I3.*Ir.*N.*m3.*t8.*1.96875e+17;
t269 = I3.*Ir.*N.*m2.*t9.*1.96875e+17;
t270 = I2.*Ir.*N.*m3.*t10.*1.96875e+17;
t271 = I3.*Ir.*l_OA.*t2.*t6.*1.95124125e+16;
t274 = I3.*Ir.*t6.*t11.*4.3946875e+16;
t276 = t10.*t11.*t14.*9.08054803125e+14;
t282 = I1.*m2.*m3.*t9.*t10.*6.328125e+16;
t283 = I3.*m2.*m3.*t8.*t11.*6.328125e+16;
t284 = I3.*m1.*m3.*t8.*t12.*6.328125e+16;
t285 = I3.*m2.*m3.*t9.*t11.*6.328125e+16;
t286 = I3.*m1.*m2.*t9.*t12.*6.328125e+16;
t287 = I2.*m2.*m3.*t10.*t11.*6.328125e+16;
t288 = I2.*m1.*m3.*t10.*t12.*6.328125e+16;
t289 = I_motor.*m2.*m3.*t9.*t10.*6.328125e+16;
t290 = Ir.*m2.*m3.*t9.*t10.*6.328125e+16;
t291 = I3.*Ir.*N.*l_AB.*l_OA.*m3.*t2.*1.96875e+17;
t292 = I3.*Ir.*N.*l_A_m2.*l_OA.*m2.*t2.*1.96875e+17;
t299 = Ir.*N.*l_OA.*m3.*t2.*t10.*1.254369375e+16;
t305 = I1.*Ir.*m3.*t6.*t10.*1.53125e+17;
t306 = I3.*Ir.*m2.*t6.*t11.*1.53125e+17;
t307 = I3.*Ir.*m3.*t6.*t11.*1.53125e+17;
t308 = I3.*Ir.*m1.*t6.*t12.*1.53125e+17;
t309 = I_motor.*Ir.*m3.*t6.*t10.*1.53125e+17;
t312 = I1.*t8.*t10.*t14.*6.328125e+16;
t313 = I3.*t8.*t11.*t14.*6.328125e+16;
t314 = I3.*t9.*t11.*t13.*6.328125e+16;
t315 = I2.*t10.*t11.*t14.*6.328125e+16;
t316 = I_motor.*t8.*t10.*t14.*6.328125e+16;
t317 = Ir.*t8.*t10.*t14.*6.328125e+16;
t318 = I3.*Ir.*m3.*t6.*t8.*2.1640625e+17;
t319 = I3.*Ir.*m2.*t6.*t9.*2.1640625e+17;
t320 = I2.*Ir.*m3.*t6.*t10.*2.1640625e+17;
t333 = m2.*m3.*t9.*t10.*t11.*1.816171875e+16;
t334 = I3.*Ir.*l_AB.*l_OA.*m3.*t2.*t6.*3.0625e+17;
t335 = I3.*Ir.*l_A_m2.*l_OA.*m2.*t2.*t6.*3.0625e+17;
t336 = Ir.*N.*m2.*m3.*t9.*t10.*1.96875e+17;
t337 = Ir.*l_OA.*m3.*t2.*t6.*t10.*1.95124125e+16;
t339 = Ir.*m3.*t6.*t10.*t11.*4.3946875e+16;
t343 = m1.*m2.*m3.*t9.*t10.*t12.*6.328125e+16;
t344 = t8.*t10.*t11.*t14.*1.816171875e+16;
t345 = Ir.*N.*l_A_m2.*l_OA.*m2.*m3.*t2.*t10.*1.96875e+17;
t346 = Ir.*N.*t8.*t10.*t14.*1.96875e+17;
t348 = t8.*t10.*t11.*t15.*6.328125e+16;
t350 = Ir.*m2.*m3.*t6.*t10.*t11.*1.53125e+17;
t351 = Ir.*m1.*m3.*t6.*t10.*t12.*1.53125e+17;
t352 = m2.*t8.*t10.*t11.*t14.*6.328125e+16;
t353 = m1.*t8.*t10.*t12.*t14.*6.328125e+16;
t354 = m2.*t9.*t10.*t11.*t14.*6.328125e+16;
t355 = m3.*t9.*t10.*t11.*t13.*6.328125e+16;
t356 = Ir.*N.*l_AB.*l_OA.*t2.*t10.*t14.*1.96875e+17;
t357 = Ir.*m2.*m3.*t6.*t9.*t10.*2.1640625e+17;
t363 = Ir.*t6.*t10.*t11.*t14.*1.53125e+17;
t364 = Ir.*t6.*t8.*t10.*t14.*2.1640625e+17;
t366 = Ir.*l_A_m2.*l_OA.*m2.*m3.*t2.*t6.*t10.*3.0625e+17;
t375 = Ir.*l_AB.*l_OA.*t2.*t6.*t10.*t14.*3.0625e+17;
t19 = t18.^2;
t48 = -t40;
t80 = l_B_m3.*l_OA.*m3.*t18.*7.98159825e+8;
t81 = l_B_m3.*l_OA.*m3.*t18.*1.125e+9;
t91 = l_B_m3.*l_OA.*m3.*t18.*7.183438425e+9;
t116 = t11.*t16.*8.2204344369e+10;
t118 = -t95;
t120 = t6.*t78;
t121 = I2.*l_B_m3.*l_OA.*m3.*t18.*5.625e+10;
t131 = -t100;
t132 = -t101;
t142 = t5.*t7.*4.9e+13;
t144 = -t115;
t145 = I2.*l_B_m3.*l_OA.*m3.*t18.*5.0625e+11;
t151 = t6.*t90;
t163 = Ir.*N.*l_B_m3.*l_OA.*m3.*t18.*8.75e+10;
t164 = l_B_m3.*m3.*t2.*t11.*t18.*3.5839125e+9;
t171 = -t146;
t172 = -t147;
t185 = l_AB.*m3.*t11.*t16.*2.580417e+12;
t186 = l_A_m2.*m2.*t11.*t16.*2.580417e+12;
t195 = -t165;
t196 = -t166;
t197 = -t167;
t202 = l_B_m3.*l_OA.*m2.*m3.*t9.*t18.*5.625e+10;
t208 = -t177;
t214 = Ir.*N.*l_B_m3.*l_OA.*m3.*t18.*3.15e+13;
t215 = I3.*t11.*t16.*2.56888576153125e+14;
t228 = l_B_m3.*l_OA.*m2.*m3.*t9.*t18.*5.0625e+11;
t232 = l_B_m3.*l_OA.*t8.*t14.*t18.*5.625e+10;
t233 = -t201;
t238 = -t207;
t246 = l_B_m3.*l_OA.*t8.*t14.*t18.*5.0625e+11;
t250 = l_A_m2.*l_B_m3.*m2.*m3.*t2.*t11.*t18.*5.625e+10;
t252 = Ir.*l_B_m3.*l_OA.*m3.*t6.*t18.*1.225e+12;
t256 = -t240;
t259 = l_B_m3.*m3.*t2.*t11.*t18.*1.2902085e+12;
t261 = I3.*t5.*t7.*1.53125e+17;
t266 = l_AB.*l_B_m3.*t2.*t11.*t14.*t18.*5.625e+10;
t267 = l_AB.*l_OA.*t3.*t10.*t14.*t18.*5.625e+10;
t275 = l_AB.*l_A_m2.*m2.*m3.*t11.*t16.*4.05e+13;
t278 = Ir.*l_B_m3.*l_OA.*m3.*t6.*t18.*4.9e+13;
t279 = t17.*t211;
t280 = m3.*t10.*t11.*t16.*2.56888576153125e+14;
t281 = l_AB.*l_OA.*t3.*t10.*t14.*t18.*5.0625e+11;
t295 = -t268;
t296 = -t269;
t297 = -t270;
t301 = I3.*l_AB.*m3.*t11.*t16.*8.063803125e+15;
t302 = I3.*l_A_m2.*m2.*t11.*t16.*8.063803125e+15;
t303 = t8.*t10.*t14.*t17.*-7.03125e+12;
t311 = l_A_m2.*l_B_m3.*m2.*m3.*t2.*t11.*t18.*2.025e+13;
t321 = -t291;
t322 = -t292;
t323 = l_AB.*l_OA.*t3.*t10.*t14.*t18.*7.03125e+12;
t324 = -t299;
t325 = t16.*t234;
t326 = t16.*t235;
t330 = l_AB.*l_B_m3.*t2.*t11.*t14.*t18.*2.025e+13;
t332 = m3.*t5.*t7.*t10.*1.53125e+17;
t340 = t8.*t11.*t14.*t16.*-2.025e+13;
t341 = t9.*t11.*t13.*t16.*-2.025e+13;
t347 = -t336;
t349 = l_A_m2.*m2.*m3.*t10.*t11.*t16.*8.063803125e+15;
t358 = -t345;
t359 = -t346;
t360 = I3.*l_AB.*l_A_m2.*m2.*m3.*t11.*t16.*1.265625e+17;
t361 = l_AB.*t10.*t11.*t14.*t16.*8.063803125e+15;
t365 = -t356;
t367 = t17.*t312;
t368 = t16.*t313;
t369 = t16.*t314;
t370 = t17.*t316;
t371 = t17.*t317;
t376 = I1.*t8.*t10.*t14.*t17.*-6.328125e+16;
t377 = I3.*t8.*t11.*t14.*t16.*-6.328125e+16;
t378 = I3.*t9.*t11.*t13.*t16.*-6.328125e+16;
t379 = I_motor.*t8.*t10.*t14.*t17.*-6.328125e+16;
t380 = Ir.*t8.*t10.*t14.*t17.*-6.328125e+16;
t383 = t17.*t344;
t384 = t17.*t346;
t385 = t16.*t348;
t386 = t17.*t348;
t387 = l_AB.*t2.*t3.*t10.*t11.*t14.*t18.*8.063803125e+15;
t389 = t17.*t352;
t390 = t17.*t353;
t391 = t16.*t355;
t392 = t8.*t10.*t11.*t14.*t17.*-1.816171875e+16;
t393 = l_AB.*l_A_m2.*m2.*t10.*t11.*t14.*t16.*1.265625e+17;
t394 = Ir.*N.*l_AB.*l_OA.*t3.*t10.*t14.*t18.*1.96875e+17;
t395 = t8.*t10.*t11.*t15.*t16.*-6.328125e+16;
t396 = t8.*t10.*t11.*t15.*t17.*-6.328125e+16;
t397 = m2.*t8.*t10.*t11.*t14.*t17.*-6.328125e+16;
t398 = m1.*t8.*t10.*t12.*t14.*t17.*-6.328125e+16;
t399 = m3.*t9.*t10.*t11.*t13.*t16.*-6.328125e+16;
t400 = t17.*t364;
t404 = Ir.*t6.*t8.*t10.*t14.*t17.*-2.1640625e+17;
t406 = Ir.*l_AB.*l_OA.*t3.*t6.*t10.*t14.*t18.*3.0625e+17;
t411 = t2.*t3.*t8.*t10.*t11.*t15.*t18.*1.265625e+17;
t412 = l_AB.*l_A_m2.*m2.*t2.*t3.*t10.*t11.*t14.*t18.*1.265625e+17;
t94 = -t80;
t108 = -t91;
t148 = -t116;
t156 = -t121;
t170 = -t145;
t192 = -t163;
t193 = -t164;
t212 = -t185;
t213 = -t186;
t237 = -t202;
t242 = -t214;
t243 = -t215;
t249 = -t228;
t253 = -t232;
t265 = -t246;
t272 = -t250;
t273 = -t252;
t277 = -t259;
t293 = -t266;
t294 = -t267;
t298 = -t275;
t300 = t19.*t161;
t304 = -t280;
t310 = -t281;
t327 = t10.*t11.*t14.*t19.*-5.625e+10;
t328 = -t301;
t329 = -t302;
t331 = -t311;
t338 = -t323;
t342 = -t330;
t362 = -t349;
t372 = -t360;
t373 = t19.*t276;
t374 = -t361;
t381 = t10.*t11.*t14.*t19.*-9.08054803125e+14;
t382 = t19.*t315;
t388 = I2.*t10.*t11.*t14.*t19.*-6.328125e+16;
t401 = -t393;
t402 = t19.*t348;
t403 = t19.*t354;
t405 = t8.*t10.*t11.*t15.*t19.*-6.328125e+16;
t407 = t19.*t363;
t408 = m2.*t9.*t10.*t11.*t14.*t19.*-6.328125e+16;
t409 = -t406;
t410 = Ir.*t6.*t10.*t11.*t14.*t19.*-1.53125e+17;
t413 = t27+t34+t35+t39+t46+t47+t52+t53+t54+t55+t61+t62+t81+t85+t105+t106+t107+t111+t134+t135+t157+t184+t210+t211+t239+t303+t338+3.2106393e+7;
t414 = t25+t29+t48+t49+t50+t58+t64+t93+t98+t99+t108+t112+t117+t132+t143+t170+t173+t198+t200+t208+t238+t247+t249+t265+t273+t310;
t415 = t20+t21+t22+t23+t24+t28+t30+t31+t32+t33+t36+t37+t38+t41+t42+t51+t56+t65+t66+t73+t74+t75+t76+t77+t78+t84+t86+t87+t88+t89+t90+t92+t94+t97+t109+t118+t120+t136+t137+t144+t151+t156+t158+t159+t160+t161+t174+t175+t176+t192+t193+t236+t237+t253+t272+t293+t294+t327;
t416 = t43+t44+t45+t59+t63+t67+t68+t69+t70+t71+t72+t79+t83+t96+t102+t103+t104+t110+t113+t114+t122+t123+t124+t125+t126+t127+t128+t129+t130+t131+t133+t138+t139+t140+t141+t142+t148+t149+t150+t152+t153+t154+t155+t162+t168+t169+t171+t172+t178+t179+t180+t181+t182+t183+t187+t188+t189+t190+t191+t194+t195+t196+t197+t203+t204+t205+t206+t209+t212+t213+t216+t217+t218+t219+t220+t221+t222+t223+t224+t225+t226+t227+t229+t230+t231+t233+t234+t235+t241+t242+t243+t244+t245+t248+t251+t254+t255+t256+t257+t258+t260+t261+t262+t263+t264+t271+t274+t276+t277+t278+t282+t283+t284+t285+t286+t287+t288+t289+t290+t295+t296+t297+t298+t304+t305+t306+t307+t308+t309+t312+t313+t314+t315+t316+t317+t318+t319+t320+t321+t322+t324+t328+t329+t331+t332+t333+t334+t335+t337+t339+t340+t341+t342+t343+t344+t347+t348+t350+t351+t352+t353+t354+t355+t357+t358+t359+t362+t363+t364+t365+t366+t372+t374+t375+t376+t377+t378+t379+t380+t381+t384+t387+t388+t392+t394+t395+t396+t397+t398+t399+t401+t404+t405+t408+t409+t410+t411+t412;
t417 = 1.0./t416;
t418 = t413.*t417.*9.0e+3;
t420 = t414.*t417.*1.25e+5;
t421 = t415.*t417.*1.125e+6;
t419 = -t418;
t422 = -t421;
Mass_Joint_Sp_inv = reshape([t417.*(I2.*2.025e+10+I3.*9.08054803125e+11+Ir.*t6.*4.9e+10+m2.*t9.*2.025e+10+m3.*t8.*2.025e+10+m3.*t10.*9.08054803125e+11+I2.*I3.*6.328125e+13+I3.*Ir.*t6.*1.53125e+14+I3.*m2.*t9.*6.328125e+13+I3.*m3.*t8.*6.328125e+13+I2.*m3.*t10.*6.328125e+13+t8.*t10.*t14.*6.328125e+13+l_AB.*l_B_m3.*m3.*t3.*2.025e+10+m2.*m3.*t9.*t10.*6.328125e+13-t8.*t10.*t14.*t17.*6.328125e+13+Ir.*m3.*t6.*t10.*1.53125e+14+2.88957537e+8).*1.0e+3,t419,t420,t419,t417.*(I1.*2.5e+8+I2.*2.5e+8+I3.*1.1210553125e+10+I_motor.*2.5e+8+Ir.*2.5e+8+t11.*7.175e+7+I3.*t11.*2.2421875e+11+Ir.*t6.*2.5e+8+l_OA.*t2.*3.1857e+7+m2.*t9.*2.5e+8+m3.*t8.*2.5e+8+m1.*t12.*2.5e+8+m2.*t11.*2.5e+8+m3.*t10.*1.1210553125e+10+m3.*t11.*2.5e+8+I1.*I3.*7.8125e+11+I2.*I3.*7.8125e+11+I3.*I_motor.*7.8125e+11+I3.*Ir.*7.8125e+11+I3.*Ir.*t6.*7.8125e+11+I3.*l_OA.*t2.*9.9553125e+10+I1.*m3.*t10.*7.8125e+11+I3.*m2.*t9.*7.8125e+11+I3.*m3.*t8.*7.8125e+11+I2.*m3.*t10.*7.8125e+11+I3.*m1.*t12.*7.8125e+11+I3.*m2.*t11.*7.8125e+11+I3.*m3.*t11.*7.8125e+11+I_motor.*m3.*t10.*7.8125e+11+Ir.*m3.*t10.*7.8125e+11+m3.*t10.*t11.*2.2421875e+11+t8.*t10.*t14.*7.8125e+11+t10.*t11.*t14.*7.8125e+11+l_AB.*l_B_m3.*m3.*t3.*2.5e+8+l_AB.*l_OA.*m3.*t2.*5.0e+8+l_A_m2.*l_OA.*m2.*t2.*5.0e+8+l_B_m3.*l_OA.*m3.*t18.*2.5e+8+l_OA.*m3.*t2.*t10.*9.9553125e+10+m2.*m3.*t9.*t10.*7.8125e+11+m1.*m3.*t10.*t12.*7.8125e+11+m2.*m3.*t10.*t11.*7.8125e+11-t8.*t10.*t14.*t17.*7.8125e+11-t10.*t11.*t14.*t19.*7.8125e+11+Ir.*m3.*t6.*t10.*7.8125e+11+I3.*l_AB.*l_OA.*m3.*t2.*1.5625e+12+I3.*l_A_m2.*l_OA.*m2.*t2.*1.5625e+12+l_AB.*l_OA.*t2.*t10.*t14.*1.5625e+12-l_AB.*l_OA.*t3.*t10.*t14.*t18.*1.5625e+12+l_A_m2.*l_OA.*m2.*m3.*t2.*t10.*1.5625e+12+3.567377e+6).*8.1e+4,t422,t420,t422,t417.*(I1.*2.90577537e+11+I_motor.*2.90577537e+11+Ir.*2.90577537e+11+t11.*8.3395753119e+10+t67+t68+t69+t83+t96+t110+t113+t114+t122+t123+t124+t125+t126+t127+t128+t129+t130+t131+t138+t142+t148+t149+t150+t168+t171+t172+t187+t188+t189+t190+t191+t196+t197+t203+t204+t205+t206+t212+t213+t234+t235+t257+t258+t298+t340+t341+I3.*t11.*5.81175e+12+Ir.*t6.*9.93703429e+11+m1.*t12.*2.90577537e+11+m2.*t11.*2.90577537e+11+m3.*t11.*2.90577537e+11+I1.*I3.*2.025e+13+I3.*I_motor.*2.025e+13+I3.*Ir.*2.025e+13-Ir.*N.*9.04019004e+11-I3.*Ir.*N.*6.3e+13+I3.*Ir.*t6.*6.925e+13+I1.*m3.*t10.*2.025e+13+I3.*m1.*t12.*2.025e+13+I3.*m2.*t11.*2.025e+13+I3.*m3.*t11.*2.025e+13+I_motor.*m3.*t10.*2.025e+13+Ir.*m3.*t10.*2.025e+13+m3.*t10.*t11.*5.81175e+12+t10.*t11.*t14.*2.025e+13+m1.*m3.*t10.*t12.*2.025e+13+m2.*m3.*t10.*t11.*2.025e+13-t10.*t11.*t14.*t19.*2.025e+13-Ir.*N.*m3.*t10.*6.3e+13+Ir.*m3.*t6.*t10.*6.925e+13+I1.*l_AB.*l_B_m3.*m3.*t3.*4.05e+13+I_motor.*l_AB.*l_B_m3.*m3.*t3.*4.05e+13+Ir.*l_AB.*l_B_m3.*m3.*t3.*4.05e+13+l_AB.*l_B_m3.*m3.*t3.*t11.*1.16235e+13+l_AB.*l_B_m3.*t3.*t11.*t14.*4.05e+13-l_B_m3.*m3.*t2.*t11.*t18.*2.580417e+12-l_AB.*l_B_m3.*t2.*t11.*t14.*t18.*4.05e+13-Ir.*N.*l_AB.*l_B_m3.*m3.*t3.*1.26e+14-Ir.*N.*l_B_m3.*l_OA.*m3.*t18.*6.3e+13+Ir.*l_AB.*l_B_m3.*m3.*t3.*t6.*1.385e+14+Ir.*l_B_m3.*l_OA.*m3.*t6.*t18.*9.8e+13+l_AB.*l_B_m3.*m1.*m3.*t3.*t12.*4.05e+13+l_AB.*l_B_m3.*m2.*m3.*t3.*t11.*4.05e+13-l_A_m2.*l_B_m3.*m2.*m3.*t2.*t11.*t18.*4.05e+13).*3.125e+3],[3,3]);
