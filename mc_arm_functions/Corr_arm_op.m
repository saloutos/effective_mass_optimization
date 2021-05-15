function Corr_Op_Sp = Corr_arm_op(in1,in2)
%CORR_ARM_OP
%    CORR_OP_SP = CORR_ARM_OP(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    12-May-2021 18:29:43

I1 = in2(5,:);
I2 = in2(6,:);
I3 = in2(7,:);
I_motor = in2(8,:);
Ir = in2(9,:);
N = in2(10,:);
dth1 = in1(4,:);
dth2 = in1(5,:);
dth3 = in1(6,:);
l_AB = in2(15,:);
l_A_m2 = in2(12,:);
l_BC = in2(16,:);
l_B_m3 = in2(13,:);
l_OA = in2(14,:);
l_O_m1 = in2(11,:);
m1 = in2(1,:);
m2 = in2(2,:);
m3 = in2(3,:);
th1 = in1(1,:);
th2 = in1(2,:);
th3 = in1(3,:);
t2 = cos(th1);
t3 = cos(th2);
t4 = cos(th3);
t5 = sin(th1);
t6 = sin(th2);
t7 = sin(th3);
t8 = th1+th2;
t9 = th2+th3;
t10 = Ir.^2;
t11 = N.^2;
t13 = dth1.^2;
t14 = dth2.^2;
t15 = dth3.^2;
t16 = l_AB.^2;
t17 = l_A_m2.^2;
t18 = l_B_m3.^2;
t19 = l_OA.^2;
t20 = l_O_m1.^2;
t21 = m2.^2;
t22 = m3.^2;
t23 = m3.^3;
t74 = Ir.*N.*(1.4e+1./9.0);
t94 = I1.*9.0e+6;
t95 = I_motor.*9.0e+6;
t96 = Ir.*9.0e+6;
t107 = I1.*2.5e+8;
t108 = I2.*2.5e+8;
t109 = I_motor.*2.5e+8;
t110 = Ir.*2.5e+8;
t113 = Ir.*N.*1.4e+7;
t117 = I2.*2.25e+9;
t122 = Ir.*N.*1.26e+8;
t125 = I3.*1.1210553125e+10;
t127 = I2.*2.025e+10;
t133 = Ir.*N.*3.5e+9;
t134 = I1.*I3.*5.625e+10;
t135 = I3.*I_motor.*5.625e+10;
t136 = I3.*Ir.*5.625e+10;
t138 = I3.*1.00894978125e+11;
t139 = I1.*2.90577537e+11;
t140 = I_motor.*2.90577537e+11;
t141 = Ir.*2.90577537e+11;
t153 = I1.*2.88957537e+11;
t154 = I_motor.*2.88957537e+11;
t155 = Ir.*2.88957537e+11;
t156 = I3.*9.08054803125e+11;
t169 = I2.*I3.*7.03125e+12;
t172 = Ir.*N.*8.98979004e+11;
t176 = Ir.*N.*9.04019004e+11;
t181 = I3.*Ir.*N.*8.75e+10;
t187 = I1.*I3.*7.8125e+11;
t188 = I2.*I3.*7.8125e+11;
t189 = I3.*I_motor.*7.8125e+11;
t190 = I3.*Ir.*7.8125e+11;
t203 = I3.*Ir.*N.*7.875e+11;
t207 = I2.*I3.*6.328125e+13;
t208 = I1.*I2.*2.025e+13;
t209 = I1.*I3.*2.025e+13;
t210 = I2.*I_motor.*2.025e+13;
t211 = I3.*I_motor.*2.025e+13;
t212 = I2.*Ir.*2.025e+13;
t213 = I3.*Ir.*2.025e+13;
t228 = I2.*Ir.*N.*6.3e+13;
t229 = I3.*Ir.*N.*6.3e+13;
t235 = I3.*Ir.*N.*1.09375e+13;
t267 = I1.*I3.*9.08054803125e+14;
t268 = I3.*I_motor.*9.08054803125e+14;
t269 = I3.*Ir.*9.08054803125e+14;
t291 = I3.*Ir.*N.*2.8250593875e+15;
t320 = I1.*I2.*I3.*6.328125e+16;
t321 = I2.*I3.*I_motor.*6.328125e+16;
t322 = I2.*I3.*Ir.*6.328125e+16;
t395 = I2.*I3.*Ir.*N.*1.96875e+17;
t12 = t11.^2;
t24 = t3.^2;
t25 = t4.^2;
t26 = l_OA.*t2;
t27 = cos(t8);
t28 = cos(t9);
t29 = l_OA.*t5;
t30 = sin(t8);
t31 = sin(t9);
t32 = t8+th3;
t35 = Ir.*t11;
t36 = m3.*t16;
t37 = m2.*t17;
t38 = m3.*t18;
t39 = m2.*t19;
t40 = m3.*t19;
t41 = m1.*t20;
t44 = l_AB.*l_B_m3.*m3.*t4;
t45 = l_AB.*l_OA.*m3.*t3;
t46 = l_A_m2.*l_OA.*m2.*t3;
t48 = l_AB.*t7.*t13;
t49 = l_AB.*t7.*t14;
t50 = dth1.*dth2.*l_AB.*t7.*2.0;
t55 = dth3.*l_AB.*l_B_m3.*m3.*t7;
t56 = dth2.*l_AB.*l_OA.*m3.*t6;
t57 = dth2.*l_A_m2.*l_OA.*m2.*t6;
t70 = dth1.*l_AB.*l_B_m3.*m3.*t7.*2.0;
t71 = dth2.*l_AB.*l_B_m3.*m3.*t7.*2.0;
t72 = dth1.*l_AB.*l_OA.*m3.*t6.*2.0;
t73 = dth1.*l_A_m2.*l_OA.*m2.*t6.*2.0;
t77 = l_AB.*l_B_m3.*m3.*t7.*t15;
t78 = l_AB.*l_OA.*m3.*t6.*t13;
t79 = l_A_m2.*l_OA.*m2.*t6.*t13;
t99 = t19.*2.583e+6;
t103 = l_OA.*t3.*5.73426e+5;
t104 = l_OA.*t3.*5.160834e+6;
t116 = t19.*(2.87e+2./1.0e+3);
t120 = l_OA.*t3.*3.1857e+7;
t121 = -t113;
t123 = t19.*7.175e+7;
t132 = l_OA.*t3.*1.433565e+8;
t165 = I3.*l_OA.*t3.*3.5839125e+9;
t171 = I3.*t19.*1.614375e+10;
t175 = I3.*l_OA.*t3.*3.22552125e+10;
t177 = t19.*8.2930813119e+10;
t180 = t19.*8.3395753119e+10;
t182 = l_OA.*t3.*1.27428e-1;
t183 = l_OA.*t3.*6.3714e-2;
t192 = dth1.*l_OA.*t6.*1.27428e-1;
t193 = dth2.*l_OA.*t6.*6.3714e-2;
t194 = I3.*l_OA.*t3.*9.9553125e+10;
t196 = I3.*l_OA.*t3.*4.479890625e+11;
t200 = -t172;
t202 = -t176;
t204 = -t181;
t205 = I3.*t19.*2.2421875e+11;
t231 = I2.*t19.*5.81175e+12;
t232 = I3.*t19.*5.81175e+12;
t246 = l_OA.*t6.*t13.*6.3714e-2;
t250 = -t228;
t251 = -t229;
t265 = Ir.*N.*l_OA.*t3.*4.013982e+12;
t314 = I3.*t19.*2.60611728496875e+14;
t341 = -t291;
t348 = t18.*t19.*t22.*5.625e+10;
t370 = l_AB.*l_OA.*t3.*t18.*t22.*5.625e+10;
t371 = l_AB.*l_B_m3.*t4.*t19.*t22.*5.625e+10;
t390 = I2.*I3.*t19.*1.816171875e+16;
t397 = l_AB.*l_OA.*t3.*t18.*t22.*5.0625e+11;
t408 = t16.*t18.*t22.*7.8125e+11;
t409 = t18.*t19.*t22.*7.8125e+11;
t415 = t16.*t18.*t22.*7.03125e+12;
t437 = l_AB.*l_OA.*t3.*t18.*t22.*1.5625e+12;
t438 = -t395;
t441 = t16.*t18.*t22.*6.328125e+13;
t444 = t16.*t19.*t22.*2.025e+13;
t445 = t17.*t19.*t21.*2.025e+13;
t446 = t18.*t19.*t22.*2.025e+13;
t447 = l_B_m3.*l_OA.*t3.*t4.*t16.*t22.*5.625e+10;
t451 = l_AB.*l_OA.*t3.*t18.*t22.*7.03125e+12;
t453 = I3.*Ir.*N.*l_OA.*t3.*1.254369375e+16;
t460 = l_AB.*l_B_m3.*t4.*t19.*t22.*4.05e+13;
t464 = l_B_m3.*l_OA.*t3.*t4.*t16.*t22.*5.0625e+11;
t468 = l_AB.*l_B_m3.*t4.*t19.*t22.*2.025e+13;
t497 = t18.*t19.*t22.*9.08054803125e+14;
t538 = I1.*t16.*t18.*t22.*6.328125e+16;
t539 = I3.*t16.*t19.*t22.*6.328125e+16;
t540 = I3.*t17.*t19.*t21.*6.328125e+16;
t541 = I2.*t18.*t19.*t22.*6.328125e+16;
t542 = I_motor.*t16.*t18.*t22.*6.328125e+16;
t543 = Ir.*t16.*t18.*t22.*6.328125e+16;
t577 = t16.*t18.*t19.*t22.*1.816171875e+16;
t579 = Ir.*N.*t16.*t18.*t22.*1.96875e+17;
t583 = t16.*t18.*t19.*t23.*6.328125e+16;
t593 = Ir.*N.*l_AB.*l_OA.*t3.*t18.*t22.*1.96875e+17;
t33 = cos(t32);
t34 = sin(t32);
t42 = t28.^2;
t43 = l_AB.*t27;
t47 = l_AB.*t30;
t51 = t44.*2.0;
t52 = t45.*2.0;
t53 = t46.*2.0;
t62 = l_B_m3.*l_OA.*m3.*t28;
t65 = t27.*t29;
t66 = t26.*t30;
t75 = l_OA.*t13.*t31;
t80 = dth1.*t55.*2.0;
t81 = dth2.*t55.*2.0;
t82 = dth2.*l_B_m3.*l_OA.*m3.*t31;
t83 = dth3.*l_B_m3.*l_OA.*m3.*t31;
t97 = dth1.*l_B_m3.*l_OA.*m3.*t31.*2.0;
t100 = -t77;
t119 = t35.*(1.96e+2./8.1e+1);
t124 = I3+t38+3.2e-4;
t126 = t35.*9.0e+6;
t128 = t39.*9.0e+6;
t129 = t40.*9.0e+6;
t130 = t41.*9.0e+6;
t142 = t35.*1.96e+8;
t143 = t35.*2.5e+8;
t144 = t45.*9.0e+6;
t145 = t46.*9.0e+6;
t146 = t36.*2.5e+8;
t147 = t37.*2.5e+8;
t148 = t39.*2.5e+8;
t149 = t40.*2.5e+8;
t150 = t41.*2.5e+8;
t157 = t36.*2.25e+9;
t158 = t37.*2.25e+9;
t160 = t45.*8.1e+7;
t161 = t46.*8.1e+7;
t162 = t44.*2.5e+8;
t163 = t45.*5.0e+8;
t164 = t46.*5.0e+8;
t166 = t44.*2.25e+9;
t167 = t45.*2.25e+9;
t168 = t46.*2.25e+9;
t170 = t38.*1.1210553125e+10;
t173 = t36.*2.025e+10;
t174 = t37.*2.025e+10;
t186 = t35.*4.9e+10;
t191 = t38.*1.00894978125e+11;
t195 = t44.*2.025e+10;
t197 = t39.*2.90577537e+11;
t198 = t40.*2.90577537e+11;
t199 = t41.*2.90577537e+11;
t206 = I3.*t35.*5.625e+10;
t214 = t39.*2.88957537e+11;
t215 = t40.*2.88957537e+11;
t216 = t41.*2.88957537e+11;
t217 = I1.*t38.*5.625e+10;
t218 = I3.*t39.*5.625e+10;
t219 = I3.*t40.*5.625e+10;
t220 = I3.*t41.*5.625e+10;
t221 = I_motor.*t38.*5.625e+10;
t222 = Ir.*t38.*5.625e+10;
t223 = t38.*9.08054803125e+11;
t225 = t35.*9.88163429e+11;
t233 = l_OA.*t3.*t38.*3.5839125e+9;
t234 = t35.*9.93703429e+11;
t236 = I1.*t44.*5.625e+10;
t237 = I3.*t45.*5.625e+10;
t238 = I3.*t46.*5.625e+10;
t239 = I_motor.*t44.*5.625e+10;
t240 = Ir.*t44.*5.625e+10;
t244 = t19.*t38.*1.614375e+10;
t245 = l_OA.*t3.*t38.*3.22552125e+10;
t248 = I3.*t35.*7.8125e+11;
t249 = Ir.*N.*t38.*8.75e+10;
t252 = l_OA.*t3.*t44.*3.5839125e+9;
t253 = I3.*t45.*5.0625e+11;
t254 = I3.*t46.*5.0625e+11;
t255 = I3+t38+t44+1.6e-4;
t256 = I3.*t36.*7.8125e+11;
t257 = I3.*t37.*7.8125e+11;
t258 = I1.*t38.*7.8125e+11;
t259 = I2.*t38.*7.8125e+11;
t260 = I3.*t39.*7.8125e+11;
t261 = I3.*t40.*7.8125e+11;
t262 = I3.*t41.*7.8125e+11;
t263 = I_motor.*t38.*7.8125e+11;
t264 = Ir.*t38.*7.8125e+11;
t266 = I3.*t35.*1.225e+12;
t270 = I3.*t36.*7.03125e+12;
t271 = I3.*t37.*7.03125e+12;
t272 = I2.*t38.*7.03125e+12;
t274 = l_OA.*t3.*t38.*9.9553125e+10;
t275 = l_AB.*l_B_m3.*t4.*t40.*1.614375e+10;
t276 = I2.*t35.*6.925e+13;
t277 = I3.*t35.*6.925e+13;
t279 = l_OA.*t3.*t38.*4.479890625e+11;
t280 = l_OA.*t3.*t44.*3.22552125e+10;
t281 = I1.*t35.*4.9e+13;
t282 = I_motor.*t35.*4.9e+13;
t284 = Ir.*N.*t44.*8.75e+10;
t285 = t19.*t24.*8.2204344369e+10;
t286 = Ir.*N.*t38.*7.875e+11;
t288 = I3.*t45.*1.5625e+12;
t289 = I3.*t46.*1.5625e+12;
t290 = t19.*t38.*2.2421875e+11;
t292 = t35.*t38.*5.625e+10;
t293 = I3.*t36.*6.328125e+13;
t294 = I3.*t37.*6.328125e+13;
t295 = I2.*t38.*6.328125e+13;
t297 = I1.*t36.*2.025e+13;
t298 = I1.*t37.*2.025e+13;
t299 = I1.*t38.*2.025e+13;
t300 = I2.*t39.*2.025e+13;
t301 = I2.*t40.*2.025e+13;
t302 = I3.*t39.*2.025e+13;
t303 = I3.*t40.*2.025e+13;
t304 = I2.*t41.*2.025e+13;
t305 = I3.*t41.*2.025e+13;
t306 = I_motor.*t36.*2.025e+13;
t307 = I_motor.*t37.*2.025e+13;
t308 = I_motor.*t38.*2.025e+13;
t309 = Ir.*t36.*2.025e+13;
t310 = Ir.*t37.*2.025e+13;
t311 = Ir.*t38.*2.025e+13;
t312 = -t265;
t315 = I3.*t45.*7.03125e+12;
t316 = I3.*t46.*7.03125e+12;
t317 = t38.*t39.*5.625e+10;
t318 = t38.*t41.*5.625e+10;
t319 = l_OA.*t3.*t35.*6.243972e+12;
t323 = t10.*t12.*4.9e+13;
t324 = I1.*t44.*4.05e+13;
t325 = I_motor.*t44.*4.05e+13;
t326 = Ir.*t44.*4.05e+13;
t327 = Ir.*N.*t44.*7.875e+11;
t330 = Ir.*N.*t36.*6.3e+13;
t331 = Ir.*N.*t37.*6.3e+13;
t332 = Ir.*N.*t38.*6.3e+13;
t334 = t19.*t36.*5.81175e+12;
t335 = t19.*t37.*5.81175e+12;
t336 = t19.*t38.*5.81175e+12;
t337 = t35.*t44.*5.625e+10;
t338 = I1.*t44.*2.025e+13;
t339 = I_motor.*t44.*2.025e+13;
t340 = Ir.*t44.*2.025e+13;
t343 = Ir.*N.*t38.*1.09375e+13;
t344 = t39.*t44.*5.625e+10;
t345 = t41.*t44.*5.625e+10;
t346 = t38.*t46.*5.625e+10;
t347 = I3.*t35.*1.53125e+14;
t350 = I3.*t35.*3.105323215625e+15;
t352 = l_B_m3.*t3.*t28.*t40.*3.5839125e+9;
t353 = Ir.*N.*t44.*6.3e+13;
t354 = Ir.*N.*t45.*6.3e+13;
t355 = Ir.*N.*t46.*6.3e+13;
t356 = t19.*t35.*1.4063e+13;
t357 = l_AB.*l_B_m3.*t4.*t40.*5.81175e+12;
t358 = t35.*t38.*7.8125e+11;
t365 = t38.*t46.*5.0625e+11;
t366 = t44.*t46.*5.625e+10;
t367 = t37.*t38.*7.8125e+11;
t368 = t38.*t39.*7.8125e+11;
t369 = t38.*t41.*7.8125e+11;
t372 = t35.*t38.*1.225e+12;
t373 = I1.*t38.*9.08054803125e+14;
t374 = I3.*t39.*9.08054803125e+14;
t375 = I3.*t40.*9.08054803125e+14;
t376 = I3.*t41.*9.08054803125e+14;
t377 = I_motor.*t38.*9.08054803125e+14;
t378 = Ir.*t38.*9.08054803125e+14;
t379 = t37.*t38.*7.03125e+12;
t380 = l_AB.*t24.*t40.*2.580417e+12;
t381 = l_A_m2.*t24.*t39.*2.580417e+12;
t382 = t35.*t36.*6.925e+13;
t383 = t35.*t37.*6.925e+13;
t384 = t35.*t38.*6.925e+13;
t385 = t35.*t39.*4.9e+13;
t386 = t35.*t40.*4.9e+13;
t387 = t35.*t41.*4.9e+13;
t394 = t44.*t46.*5.0625e+11;
t396 = t38.*t46.*1.5625e+12;
t398 = Ir.*N.*t38.*2.8250593875e+15;
t399 = l_AB.*l_B_m3.*t4.*t40.*1.16235e+13;
t400 = t37.*t38.*6.328125e+13;
t402 = t36.*t39.*2.025e+13;
t403 = t36.*t41.*2.025e+13;
t404 = t37.*t40.*2.025e+13;
t405 = t37.*t41.*2.025e+13;
t406 = t38.*t39.*2.025e+13;
t407 = t38.*t41.*2.025e+13;
t410 = t35.*t44.*1.225e+12;
t411 = Ir.*N.*t44.*1.26e+14;
t413 = t19.*t38.*2.60611728496875e+14;
t414 = t38.*t46.*7.03125e+12;
t419 = I3.*t19.*t24.*2.56888576153125e+14;
t421 = I1.*I3.*t36.*6.328125e+16;
t422 = I1.*I3.*t37.*6.328125e+16;
t423 = I1.*I2.*t38.*6.328125e+16;
t424 = I2.*I3.*t39.*6.328125e+16;
t425 = I2.*I3.*t40.*6.328125e+16;
t426 = I2.*I3.*t41.*6.328125e+16;
t427 = I3.*I_motor.*t36.*6.328125e+16;
t428 = I3.*I_motor.*t37.*6.328125e+16;
t429 = I2.*I_motor.*t38.*6.328125e+16;
t430 = I3.*Ir.*t36.*6.328125e+16;
t431 = I3.*Ir.*t37.*6.328125e+16;
t432 = I2.*Ir.*t38.*6.328125e+16;
t434 = t39.*t44.*4.05e+13;
t435 = t41.*t44.*4.05e+13;
t439 = t39.*t44.*2.025e+13;
t440 = t41.*t44.*2.025e+13;
t442 = l_B_m3.*l_OA.*t16.*t22.*t28.*5.625e+10;
t452 = t35.*t38.*1.53125e+14;
t454 = t35.*t44.*6.925e+13;
t457 = I1.*I3.*t35.*1.53125e+17;
t458 = I3.*I_motor.*t35.*1.53125e+17;
t462 = l_B_m3.*t3.*t28.*t40.*2.580417e+12;
t463 = l_B_m3.*l_OA.*t16.*t22.*t28.*5.0625e+11;
t465 = I2.*I3.*t35.*2.1640625e+17;
t467 = l_A_m2.*l_B_m3.*m3.*t3.*t28.*t39.*5.625e+10;
t469 = t35.*t44.*1.385e+14;
t474 = t38.*t39.*9.08054803125e+14;
t475 = t38.*t41.*9.08054803125e+14;
t476 = -t453;
t477 = t35.*t45.*9.8e+13;
t478 = t35.*t46.*9.8e+13;
t479 = l_B_m3.*t3.*t28.*t40.*1.2902085e+12;
t480 = t35.*t38.*3.105323215625e+15;
t481 = I3.*t10.*t12.*1.53125e+17;
t483 = I3.*t19.*t36.*1.816171875e+16;
t484 = I3.*t19.*t37.*1.816171875e+16;
t485 = I2.*t19.*t38.*1.816171875e+16;
t487 = l_AB.*l_B_m3.*t3.*t19.*t22.*t28.*5.625e+10;
t488 = l_AB.*l_OA.*t4.*t18.*t22.*t28.*5.625e+10;
t489 = I3.*Ir.*N.*t36.*1.96875e+17;
t490 = I3.*Ir.*N.*t37.*1.96875e+17;
t491 = I2.*Ir.*N.*t38.*1.96875e+17;
t492 = I3.*l_OA.*t3.*t35.*1.95124125e+16;
t495 = I3.*t19.*t35.*4.3946875e+16;
t496 = l_AB.*l_A_m2.*m3.*t24.*t39.*4.05e+13;
t501 = t25.*t415;
t502 = t19.*t24.*t38.*2.56888576153125e+14;
t503 = l_AB.*l_OA.*t4.*t18.*t22.*t28.*5.0625e+11;
t504 = I1.*t37.*t38.*6.328125e+16;
t505 = I3.*t36.*t39.*6.328125e+16;
t506 = I3.*t36.*t41.*6.328125e+16;
t507 = I3.*t37.*t40.*6.328125e+16;
t508 = I3.*t37.*t41.*6.328125e+16;
t509 = I2.*t38.*t39.*6.328125e+16;
t510 = I2.*t38.*t41.*6.328125e+16;
t511 = I_motor.*t37.*t38.*6.328125e+16;
t512 = Ir.*t37.*t38.*6.328125e+16;
t513 = I3.*Ir.*N.*t45.*1.96875e+17;
t514 = I3.*Ir.*N.*t46.*1.96875e+17;
t522 = Ir.*N.*l_OA.*t3.*t38.*1.254369375e+16;
t523 = t25.*t408;
t525 = I3.*l_AB.*t24.*t40.*8.063803125e+15;
t526 = I3.*l_A_m2.*t24.*t39.*8.063803125e+15;
t527 = l_A_m2.*l_B_m3.*m3.*t3.*t28.*t39.*4.05e+13;
t528 = t16.*t18.*t22.*t25.*-7.03125e+12;
t529 = l_AB.*l_OA.*t4.*t18.*t22.*t28.*1.5625e+12;
t531 = I1.*t35.*t38.*1.53125e+17;
t532 = I3.*t35.*t39.*1.53125e+17;
t533 = I3.*t35.*t40.*1.53125e+17;
t534 = I3.*t35.*t41.*1.53125e+17;
t535 = I_motor.*t35.*t38.*1.53125e+17;
t537 = l_A_m2.*l_B_m3.*m3.*t3.*t28.*t39.*2.025e+13;
t544 = I3.*t35.*t36.*2.1640625e+17;
t545 = I3.*t35.*t37.*2.1640625e+17;
t546 = I2.*t35.*t38.*2.1640625e+17;
t549 = l_AB.*l_OA.*t4.*t18.*t22.*t28.*7.03125e+12;
t551 = t25.*t441;
t552 = t24.*t444;
t553 = t24.*t445;
t554 = t16.*t18.*t22.*t25.*-7.8125e+11;
t560 = l_AB.*l_B_m3.*t3.*t19.*t22.*t28.*2.025e+13;
t562 = t10.*t12.*t38.*1.53125e+17;
t563 = t19.*t37.*t38.*1.816171875e+16;
t564 = I3.*t35.*t45.*3.0625e+17;
t565 = I3.*t35.*t46.*3.0625e+17;
t566 = Ir.*N.*t37.*t38.*1.96875e+17;
t567 = l_OA.*t3.*t35.*t38.*1.95124125e+16;
t569 = t19.*t35.*t38.*4.3946875e+16;
t571 = t16.*t18.*t22.*t25.*-6.328125e+13;
t572 = t16.*t19.*t22.*t24.*-2.025e+13;
t573 = t17.*t19.*t21.*t24.*-2.025e+13;
t574 = l_AB.*l_B_m3.*t3.*t19.*t22.*t28.*4.05e+13;
t576 = t37.*t38.*t41.*6.328125e+16;
t578 = Ir.*N.*t38.*t46.*1.96875e+17;
t584 = l_A_m2.*t24.*t38.*t39.*8.063803125e+15;
t587 = t35.*t38.*t39.*1.53125e+17;
t588 = t35.*t38.*t41.*1.53125e+17;
t589 = t16.*t18.*t22.*t39.*6.328125e+16;
t590 = t16.*t18.*t22.*t41.*6.328125e+16;
t591 = t18.*t19.*t22.*t37.*6.328125e+16;
t592 = t17.*t19.*t21.*t38.*6.328125e+16;
t594 = t35.*t37.*t38.*2.1640625e+17;
t596 = -t579;
t598 = I3.*l_AB.*l_A_m2.*m3.*t24.*t39.*1.265625e+17;
t599 = l_AB.*t18.*t19.*t22.*t24.*8.063803125e+15;
t601 = t18.*t19.*t22.*t35.*1.53125e+17;
t602 = t16.*t18.*t22.*t35.*2.1640625e+17;
t603 = -t593;
t604 = t35.*t38.*t46.*3.0625e+17;
t605 = t25.*t538;
t606 = t24.*t539;
t607 = t24.*t540;
t608 = t25.*t542;
t609 = t25.*t543;
t613 = l_AB.*l_OA.*t3.*t18.*t22.*t35.*3.0625e+17;
t614 = I1.*t16.*t18.*t22.*t25.*-6.328125e+16;
t615 = I3.*t16.*t19.*t22.*t24.*-6.328125e+16;
t616 = I3.*t17.*t19.*t21.*t24.*-6.328125e+16;
t617 = I_motor.*t16.*t18.*t22.*t25.*-6.328125e+16;
t618 = Ir.*t16.*t18.*t22.*t25.*-6.328125e+16;
t621 = t25.*t577;
t622 = t25.*t579;
t623 = t24.*t583;
t624 = t25.*t583;
t625 = l_AB.*t3.*t4.*t18.*t19.*t22.*t28.*8.063803125e+15;
t630 = t16.*t18.*t19.*t22.*t25.*-1.816171875e+16;
t631 = l_AB.*l_A_m2.*t18.*t22.*t24.*t39.*1.265625e+17;
t632 = Ir.*N.*l_AB.*l_OA.*t4.*t18.*t22.*t28.*1.96875e+17;
t633 = t16.*t18.*t19.*t23.*t24.*-6.328125e+16;
t634 = t16.*t18.*t19.*t23.*t25.*-6.328125e+16;
t635 = t16.*t18.*t22.*t25.*t39.*-6.328125e+16;
t636 = t16.*t18.*t22.*t25.*t41.*-6.328125e+16;
t637 = t17.*t19.*t21.*t24.*t38.*-6.328125e+16;
t642 = t16.*t18.*t22.*t25.*t35.*-2.1640625e+17;
t644 = l_AB.*l_OA.*t4.*t18.*t22.*t28.*t35.*3.0625e+17;
t651 = t3.*t4.*t16.*t18.*t19.*t23.*t28.*1.265625e+17;
t652 = l_AB.*l_A_m2.*t3.*t4.*t18.*t22.*t28.*t39.*1.265625e+17;
t54 = l_BC.*t33;
t58 = l_BC.*t34;
t63 = t5.*t43;
t64 = t2.*t47;
t76 = t62.*2.0;
t84 = t29.*t43;
t85 = t26.*t47;
t90 = -t80;
t91 = -t81;
t93 = -t66;
t98 = l_B_m3.*m3.*t75;
t101 = t26+t43;
t102 = t29+t47;
t159 = -t142;
t201 = t62.*2.5e+8;
t226 = t62.*7.98159825e+8;
t227 = t62.*1.125e+9;
t241 = t62.*7.183438425e+9;
t287 = -t249;
t296 = I2.*t62.*5.625e+10;
t313 = -t266;
t328 = -t284;
t329 = I2.*t62.*5.0625e+11;
t333 = -t285;
t351 = Ir.*N.*t62.*8.75e+10;
t362 = -t330;
t363 = -t331;
t364 = -t332;
t389 = -t352;
t391 = -t353;
t392 = -t354;
t393 = -t355;
t401 = t37.*t62.*5.625e+10;
t412 = -t372;
t416 = -t380;
t417 = -t381;
t418 = Ir.*N.*t62.*3.15e+13;
t420 = Ir.*N.*t62.*6.3e+13;
t436 = t37.*t62.*5.0625e+11;
t443 = -t398;
t449 = -t410;
t450 = -t411;
t456 = -t419;
t470 = t35.*t62.*1.225e+12;
t472 = -t442;
t482 = -t462;
t486 = -t463;
t493 = -t467;
t498 = -t479;
t499 = t35.*t62.*4.9e+13;
t500 = t62+t255;
t515 = -t487;
t516 = -t488;
t517 = -t489;
t518 = -t490;
t519 = -t491;
t520 = t35.*t62.*9.8e+13;
t521 = -t496;
t524 = t42.*t348;
t530 = -t502;
t536 = -t503;
t547 = -t513;
t548 = -t514;
t550 = -t522;
t555 = t18.*t19.*t22.*t42.*-5.625e+10;
t556 = -t525;
t557 = -t526;
t558 = -t527;
t559 = -t529;
t561 = -t537;
t568 = -t549;
t570 = t42.*t409;
t575 = -t560;
t580 = -t566;
t581 = t42.*t446;
t582 = t18.*t19.*t22.*t42.*-7.8125e+11;
t585 = -t574;
t586 = t48+t49+t50+t75;
t595 = -t578;
t597 = t18.*t19.*t22.*t42.*-2.025e+13;
t600 = -t584;
t610 = -t598;
t611 = t42.*t497;
t612 = -t599;
t619 = t18.*t19.*t22.*t42.*-9.08054803125e+14;
t620 = t42.*t541;
t626 = I2.*t18.*t19.*t22.*t42.*-6.328125e+16;
t627 = t25.*t589;
t628 = t25.*t590;
t629 = t24.*t592;
t638 = t25.*t602;
t639 = -t631;
t640 = t42.*t583;
t641 = t42.*t591;
t643 = t16.*t18.*t19.*t23.*t42.*-6.328125e+16;
t646 = t42.*t601;
t647 = t18.*t19.*t22.*t37.*t42.*-6.328125e+16;
t648 = -t644;
t650 = t18.*t19.*t22.*t35.*t42.*-1.53125e+17;
t661 = t55+t70+t71+t82+t83+t97;
t666 = I2+I3+t36+t37+t38+t51+t119+1.4349508e-2;
t676 = I2+I3+t36+t37+t38+t45+t46+t51+t62+t74+t183+1.4349508e-2;
t683 = t56+t57+t72+t73+t82+t83+t97+t192+t193;
t724 = t127+t156+t173+t174+t186+t195+t207+t223+t293+t294+t295+t347+t400+t441+t452+t571+2.88957537e+8;
t59 = dth1.*t54;
t60 = dth2.*t54;
t61 = dth3.*t54;
t67 = dth1.*t58;
t68 = dth2.*t58;
t69 = dth3.*t58;
t86 = t5.*t54;
t87 = t2.*t58;
t88 = t29.*t54;
t89 = t26.*t58;
t92 = -t64;
t105 = t30.*t54;
t106 = t27.*t58;
t111 = -t85;
t114 = t47.*t54;
t115 = t43.*t58;
t151 = t43+t54;
t152 = t47+t58;
t224 = t54+t101;
t230 = t58+t102;
t243 = t65+t93;
t247 = -t226;
t273 = -t241;
t342 = -t296;
t359 = -t329;
t388 = -t351;
t448 = -t401;
t455 = -t418;
t459 = -t420;
t466 = -t436;
t494 = -t470;
t664 = dth3.*t661;
t674 = t78+t79+t90+t91+t98+t100+t246;
t684 = dth2.*t683;
t688 = I1+I2+I3+I_motor+Ir+t35+t36+t37+t38+t39+t40+t41+t51+t52+t53+t76+t116+t182+1.4349508e-2;
t748 = t117+t132+t133+t138+t157+t158+t166+t167+t168+t169+t191+t196+t227+t235+t270+t271+t272+t279+t315+t316+t343+t379+t414+t415+t451+t528+t568+3.2106393e+7;
t754 = t107+t108+t109+t110+t120+t123+t125+t143+t146+t147+t148+t149+t150+t162+t163+t164+t170+t187+t188+t189+t190+t194+t201+t205+t248+t256+t257+t258+t259+t260+t261+t262+t263+t264+t274+t288+t289+t290+t358+t367+t368+t369+t396+t408+t409+t437+t554+t559+t582+3.567377e+6;
t112 = -t87;
t118 = -t89;
t131 = -t106;
t137 = -t115;
t178 = dth1.*t151;
t179 = dth2.*t151;
t184 = dth1.*t152;
t185 = dth2.*t152;
t242 = t63+t92;
t278 = dth1.*t224;
t283 = dth1.*t230;
t349 = t84+t111;
t361 = 1.0./t243;
t461 = t59+t60+t61;
t471 = t67+t68+t69;
t698 = t664+t684;
t749 = t104+t122+t159+t160+t161+t175+t203+t245+t253+t254+t273+t280+t286+t313+t327+t359+t365+t394+t397+t412+t449+t464+t466+t486+t494+t536;
t757 = t94+t95+t96+t99+t103+t121+t126+t128+t129+t130+t134+t135+t136+t144+t145+t165+t171+t204+t206+t217+t218+t219+t220+t221+t222+t233+t236+t237+t238+t239+t240+t244+t247+t252+t275+t287+t292+t317+t318+t328+t337+t342+t344+t345+t346+t348+t366+t370+t371+t388+t389+t447+t448+t472+t493+t515+t516+t555;
t761 = t139+t140+t141+t180+t197+t198+t199+t202+t208+t209+t210+t211+t212+t213+t231+t232+t234+t250+t251+t276+t277+t281+t282+t297+t298+t299+t300+t301+t302+t303+t304+t305+t306+t307+t308+t309+t310+t311+t312+t319+t323+t324+t325+t326+t333+t334+t335+t336+t356+t362+t363+t364+t382+t383+t384+t385+t386+t387+t392+t393+t399+t402+t403+t404+t405+t406+t407+t416+t417+t434+t435+t444+t445+t446+t450+t459+t460+t469+t477+t478+t482+t520+t521+t558+t572+t573+t585+t597;
t762 = t153+t154+t155+t177+t200+t208+t210+t212+t214+t215+t216+t225+t231+t250+t267+t268+t269+t276+t281+t282+t297+t298+t300+t301+t304+t306+t307+t309+t310+t312+t314+t319+t320+t321+t322+t323+t333+t334+t335+t338+t339+t340+t341+t350+t356+t357+t362+t363+t373+t374+t375+t376+t377+t378+t382+t383+t385+t386+t387+t390+t391+t392+t393+t402+t403+t404+t405+t413+t416+t417+t421+t422+t423+t424+t425+t426+t427+t428+t429+t430+t431+t432+t438+t439+t440+t443+t444+t445+t454+t455+t456+t457+t458+t465+t468+t474+t475+t476+t477+t478+t480+t481+t483+t484+t485+t492+t495+t497+t498+t499+t504+t505+t506+t507+t508+t509+t510+t511+t512+t517+t518+t519+t521+t530+t531+t532+t533+t534+t535+t538+t539+t540+t541+t542+t543+t544+t545+t546+t547+t548+t550+t556+t557+t561+t562+t563+t564+t565+t567+t569+t572+t573+t575+t576+t577+t580+t583+t587+t588+t589+t590+t591+t592+t594+t595+t596+t600+t601+t602+t603+t604+t610+t612+t613+t614+t615+t616+t617+t618+t619+t622+t625+t626+t630+t632+t633+t634+t635+t636+t637+t639+t642+t643+t647+t648+t650+t651+t652;
t360 = 1.0./t242;
t433 = 1.0./t349;
t473 = t105+t131;
t653 = t86+t112+t242;
t654 = t69+t184+t185;
t655 = t61+t178+t179;
t658 = t61+t179+t278;
t659 = t69+t185+t283;
t660 = t88+t114+t118+t137;
t665 = t27.*t361.*t500;
t667 = t30.*t361.*t500;
t686 = t27.*t361.*t676;
t687 = t30.*t361.*t676;
t696 = t27.*t361.*t688;
t697 = t30.*t361.*t688;
t763 = 1.0./t762;
t645 = t5.*t124.*t360;
t649 = t2.*t124.*t360;
t656 = t2.*t255.*t360;
t657 = t5.*t255.*t360;
t662 = t2.*t360.*t500;
t663 = t5.*t360.*t500;
t668 = t102.*t255.*t433;
t669 = t101.*t255.*t433;
t672 = t361.*t473.*t500;
t673 = t124.*t360.*t653;
t675 = t255.*t360.*t653;
t677 = t360.*t500.*t653;
t678 = t101.*t433.*t666;
t679 = t102.*t433.*t666;
t682 = t255.*t433.*t660;
t689 = t101.*t433.*t676;
t690 = t102.*t433.*t676;
t693 = t433.*t660.*t666;
t695 = t361.*t473.*t676;
t701 = t433.*t660.*t676;
t703 = t361.*t473.*t688;
t764 = t724.*t763.*1.0e+3;
t767 = t748.*t763.*9.0e+3;
t771 = t151.*t748.*t763.*-9.0e+3;
t772 = t152.*t748.*t763.*-9.0e+3;
t773 = t749.*t763.*1.25e+5;
t780 = t754.*t763.*8.1e+4;
t784 = t152.*t754.*t763.*-8.1e+4;
t785 = t151.*t754.*t763.*-8.1e+4;
t786 = t757.*t763.*1.125e+6;
t792 = t151.*t757.*t763.*-1.125e+6;
t793 = t152.*t757.*t763.*-1.125e+6;
t794 = t761.*t763.*3.125e+3;
t670 = -t668;
t671 = -t669;
t680 = -t678;
t681 = -t679;
t685 = -t682;
t691 = -t689;
t692 = -t690;
t694 = -t693;
t702 = -t701;
t765 = t224.*t764;
t766 = t230.*t764;
t768 = -t767;
t769 = t151.*t767;
t770 = t152.*t767;
t774 = t230.*t767;
t775 = t224.*t767;
t776 = t54.*t773;
t777 = t58.*t773;
t778 = t224.*t773;
t779 = t230.*t773;
t781 = -t780;
t782 = t152.*t780;
t783 = t151.*t780;
t787 = -t786;
t788 = t54.*t786;
t789 = t58.*t786;
t790 = t151.*t786;
t791 = t152.*t786;
t795 = t54.*t794;
t796 = t58.*t794;
t699 = t645+t667+t670;
t700 = t649+t665+t671;
t710 = t672+t673+t685;
t711 = t656+t680+t686;
t712 = t657+t681+t687;
t725 = t662+t691+t696;
t726 = t663+t692+t697;
t731 = t675+t694+t695;
t742 = t677+t702+t703;
t797 = t764+t768+t773;
t798 = t765+t771+t776;
t799 = t766+t772+t777;
t800 = t767+t781+t786;
t801 = t775+t785+t788;
t802 = t774+t784+t789;
t803 = t773+t787+t794;
t804 = t778+t792+t795;
t805 = t779+t793+t796;
t704 = t2.*t360.*t700;
t705 = t5.*t360.*t700;
t706 = t2.*t360.*t699;
t707 = t5.*t360.*t699;
t708 = t360.*t653.*t700;
t709 = t360.*t653.*t699;
t713 = t2.*t360.*t710;
t714 = t5.*t360.*t710;
t715 = t101.*t433.*t711;
t716 = t102.*t433.*t711;
t717 = t101.*t433.*t712;
t718 = t102.*t433.*t712;
t723 = t360.*t653.*t710;
t727 = t433.*t660.*t711;
t728 = t433.*t660.*t712;
t732 = t27.*t361.*t726;
t733 = t30.*t361.*t726;
t734 = t27.*t361.*t725;
t735 = t30.*t361.*t725;
t736 = t361.*t473.*t726;
t737 = t361.*t473.*t725;
t738 = t101.*t433.*t731;
t739 = t102.*t433.*t731;
t743 = t433.*t660.*t731;
t745 = t27.*t361.*t742;
t746 = t30.*t361.*t742;
t747 = t361.*t473.*t742;
t719 = -t715;
t720 = -t716;
t721 = -t717;
t722 = -t718;
t729 = -t727;
t730 = -t728;
t740 = -t738;
t741 = -t739;
t744 = -t743;
t750 = t706+t721+t732;
t751 = t707+t722+t733;
t752 = t704+t719+t734;
t753 = t705+t720+t735;
t755 = t709+t730+t736;
t756 = t708+t729+t737;
t758 = t714+t741+t746;
t759 = t713+t740+t745;
t760 = t723+t744+t747;
Corr_Op_Sp = [t674.*(-t750.*t801+t752.*t802+t759.*t800)+t698.*(-t750.*t798+t752.*t799+t759.*t797)+dth3.*(t461.*t752+t471.*t753)+dth2.*(t654.*t753+t655.*t752)+dth1.*(t658.*t752+t659.*t753)-l_B_m3.*m3.*t586.*(-t750.*t804+t752.*t805+t759.*t803);t674.*(-t751.*t801+t753.*t802+t758.*t800)+t698.*(-t751.*t798+t753.*t799+t758.*t797)+dth3.*(t461.*t750+t471.*t751)+dth2.*(t654.*t751+t655.*t750)+dth1.*(t658.*t750+t659.*t751)-l_B_m3.*m3.*t586.*(-t751.*t804+t753.*t805+t758.*t803);-t674.*(-t755.*t801+t756.*t802+t760.*t800)-t698.*(-t755.*t798+t756.*t799+t760.*t797)-dth3.*(t461.*t759+t471.*t758)-dth2.*(t654.*t758+t655.*t759)-dth1.*(t658.*t759+t659.*t758)+l_B_m3.*m3.*t586.*(-t755.*t804+t756.*t805+t760.*t803)];
