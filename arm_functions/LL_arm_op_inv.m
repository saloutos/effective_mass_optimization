function Mass_Op_Sp_inv = LL_arm_op_inv(in1,in2)
%LL_ARM_OP_INV
%    MASS_OP_SP_INV = LL_ARM_OP_INV(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    03-May-2021 15:46:37

I1 = in2(5,:);
I2 = in2(6,:);
I3 = in2(7,:);
I_motor = in2(8,:);
Ir = in2(9,:);
N = in2(10,:);
l_AB = in2(15,:);
l_A_m2 = in2(12,:);
l_BC = in2(16,:);
l_B_m3 = in2(13,:);
l_OA = in2(14,:);
l_O_m1 = in2(11,:);
m1 = in2(1,:);
m2 = in2(2,:);
m3 = in2(3,:);
m_motor = in2(4,:);
th1 = in1(1,:);
th2 = in1(2,:);
th3 = in1(3,:);
t2 = cos(th1);
t3 = cos(th2);
t4 = cos(th3);
t5 = sin(th1);
t6 = th1+th2;
t7 = th2+th3;
t8 = I_motor.^2;
t9 = Ir.^2;
t10 = Ir.^3;
t11 = N.^2;
t12 = N.^3;
t15 = l_AB.^2;
t16 = l_A_m2.^2;
t17 = l_B_m3.^2;
t18 = l_OA.^2;
t19 = l_O_m1.^2;
t20 = m2.^2;
t21 = m3.^2;
t22 = m3.^3;
t23 = m_motor.^2;
t24 = I1.*I2;
t25 = I1.*I3;
t26 = I2.*I3;
t27 = I1.*I_motor;
t28 = I2.*I_motor;
t29 = I3.*I_motor;
t30 = I1.*Ir;
t31 = I2.*Ir;
t32 = I3.*Ir;
t42 = I_motor.*Ir.*N;
t51 = I_motor.*Ir.*2.0;
t13 = t11.^2;
t14 = t11.^3;
t33 = t3.^2;
t34 = t4.^2;
t35 = I3.*t24;
t36 = I_motor.*t25;
t37 = I_motor.*t26;
t38 = Ir.*t25;
t39 = Ir.*t26;
t40 = N.*t30;
t41 = N.*t32;
t43 = l_OA.*t2;
t44 = cos(t6);
t45 = cos(t7);
t46 = l_OA.*t5;
t47 = sin(t6);
t48 = t6+th3;
t49 = t29.*2.0;
t50 = t32.*2.0;
t55 = N.*t31.*2.0;
t57 = t42.*2.0;
t58 = I3.*t8;
t59 = I3.*t9;
t60 = N.*t9;
t71 = Ir.*N.*t29.*4.0;
t72 = t11.*t30;
t73 = t11.*t31;
t74 = t11.*t32;
t75 = I_motor.*Ir.*t11;
t76 = I1.*m3.*t15;
t77 = I3.*m3.*t15;
t78 = I1.*m2.*t16;
t79 = I3.*m2.*t16;
t80 = I1.*m3.*t17;
t81 = I2.*m3.*t17;
t82 = I2.*m2.*t18;
t83 = I2.*m3.*t18;
t84 = I3.*m2.*t18;
t85 = I3.*m3.*t18;
t86 = I2.*m1.*t19;
t87 = I3.*m1.*t19;
t88 = I1.*m_motor.*t15;
t89 = I3.*m_motor.*t15;
t90 = I_motor.*m3.*t15;
t91 = I_motor.*m2.*t16;
t92 = I_motor.*m3.*t17;
t93 = I_motor.*m2.*t18;
t94 = I_motor.*m3.*t18;
t95 = I_motor.*m1.*t19;
t96 = Ir.*m3.*t15;
t97 = Ir.*m2.*t16;
t98 = Ir.*m3.*t17;
t99 = Ir.*m2.*t18;
t100 = Ir.*m3.*t18;
t101 = Ir.*m1.*t19;
t102 = I_motor.*m_motor.*t15;
t103 = Ir.*m_motor.*t15;
t104 = m3.*t15.*t25;
t105 = m2.*t16.*t25;
t106 = m3.*t17.*t24;
t107 = m2.*t18.*t26;
t108 = m3.*t18.*t26;
t109 = m1.*t19.*t26;
t110 = m_motor.*t15.*t25;
t111 = m3.*t15.*t29;
t112 = m2.*t16.*t29;
t113 = m3.*t17.*t27;
t114 = m3.*t17.*t28;
t115 = m2.*t18.*t29;
t116 = m3.*t18.*t29;
t117 = m1.*t19.*t29;
t118 = m3.*t15.*t32;
t119 = m2.*t16.*t32;
t120 = m3.*t17.*t30;
t121 = m3.*t17.*t31;
t122 = m2.*t18.*t32;
t123 = m3.*t18.*t32;
t124 = m1.*t19.*t32;
t125 = m_motor.*t15.*t29;
t126 = m_motor.*t15.*t32;
t138 = t11.*t51;
t140 = t9.*t11;
t141 = t9.*t12;
t144 = I1.*l_AB.*l_B_m3.*m3.*t4;
t145 = I3.*l_AB.*l_OA.*m3.*t3;
t146 = I3.*l_A_m2.*l_OA.*m2.*t3;
t147 = I_motor.*l_AB.*l_B_m3.*m3.*t4;
t148 = I3.*l_AB.*l_OA.*m_motor.*t3;
t149 = Ir.*l_AB.*l_B_m3.*m3.*t4;
t151 = I2.*m_motor.*t18.*2.0;
t152 = I3.*m_motor.*t18.*2.0;
t154 = I_motor.*m_motor.*t18.*2.0;
t155 = Ir.*m_motor.*t18.*2.0;
t156 = Ir.*t11.*t24;
t157 = Ir.*t11.*t27;
t158 = Ir.*t11.*t28;
t162 = Ir.*N.*l_AB.*l_OA.*m3.*t3;
t163 = Ir.*N.*l_A_m2.*l_OA.*m2.*t3;
t164 = Ir.*N.*l_AB.*l_OA.*m_motor.*t3;
t165 = m_motor.*t18.*t26.*2.0;
t167 = m3.*t17.*t51;
t174 = m3.*t8.*t17;
t175 = m3.*t9.*t17;
t187 = Ir.*t11.*t29.*4.0;
t191 = Ir.*t8.*t11;
t210 = m3.*t17.*t42.*4.0;
t224 = m2.*m3.*t16.*t17;
t225 = m2.*m3.*t15.*t18;
t226 = m1.*m3.*t15.*t19;
t227 = m2.*m3.*t16.*t18;
t228 = m1.*m2.*t16.*t19;
t229 = m2.*m3.*t17.*t18;
t230 = m1.*m3.*t17.*t19;
t231 = m3.*m_motor.*t15.*t17;
t232 = m2.*m_motor.*t15.*t18;
t233 = m1.*m_motor.*t15.*t19;
t271 = Ir.*l_AB.*l_OA.*m3.*t3.*t11;
t272 = Ir.*l_A_m2.*l_OA.*m2.*t3.*t11;
t273 = Ir.*l_AB.*l_OA.*m_motor.*t3.*t11;
t290 = l_AB.*l_B_m3.*m2.*m3.*t4.*t18;
t291 = l_AB.*l_B_m3.*m1.*m3.*t4.*t19;
t292 = l_A_m2.*l_OA.*m2.*m3.*t3.*t17;
t293 = l_AB.*l_OA.*m3.*m_motor.*t3.*t17;
t294 = t15.*t17.*t21;
t295 = t15.*t18.*t21;
t296 = t16.*t18.*t20;
t297 = t17.*t18.*t21;
t298 = m3.*m_motor.*t15.*t18.*3.0;
t299 = m2.*m_motor.*t16.*t18.*2.0;
t300 = m3.*m_motor.*t17.*t18.*2.0;
t333 = l_AB.*l_B_m3.*m3.*t4.*t42.*-2.0;
t346 = l_AB.*l_A_m2.*l_B_m3.*l_OA.*m2.*m3.*t3.*t4;
t347 = l_AB.*l_OA.*t3.*t17.*t21;
t348 = l_AB.*l_B_m3.*t4.*t18.*t21;
t353 = l_AB.*l_B_m3.*m3.*m_motor.*t4.*t18.*2.0;
t354 = l_AB.*l_B_m3.*m3.*m_motor.*t4.*t18.*4.0;
t355 = t15.*t18.*t23.*2.0;
t386 = l_B_m3.*l_OA.*m3.*m_motor.*t3.*t4.*t15;
t425 = l_B_m3.*l_OA.*t3.*t4.*t15.*t21;
t426 = t15.*t17.*t18.*t22;
t52 = cos(t48);
t53 = sin(t48);
t54 = Ir.*t49;
t56 = t41.*2.0;
t61 = t45.^2;
t62 = -t55;
t63 = -t41;
t65 = -t57;
t66 = t60.*2.0;
t67 = l_AB.*t44;
t68 = l_AB.*t47;
t69 = N.*t38.*2.0;
t70 = N.*t39.*2.0;
t127 = N.*t98;
t128 = N.*t99;
t129 = N.*t100;
t130 = N.*t101;
t134 = -t71;
t135 = t73.*2.0;
t136 = t11.*t50;
t137 = N.*t59.*4.0;
t139 = t75.*3.0;
t142 = t9.*t13;
t143 = t10.*t14;
t150 = t92.*2.0;
t153 = t98.*2.0;
t161 = N.*t149;
t166 = m_motor.*t18.*t49;
t168 = m_motor.*t18.*t50;
t169 = N.*t96.*2.0;
t170 = N.*t97.*2.0;
t172 = N.*t103.*2.0;
t173 = N.*t155;
t176 = -t74;
t178 = t140.*2.0;
t179 = t144.*2.0;
t180 = t145.*2.0;
t181 = t146.*2.0;
t182 = t147.*2.0;
t183 = t148.*2.0;
t184 = t149.*2.0;
t185 = t11.*t38.*2.0;
t186 = t11.*t39.*2.0;
t188 = I2.*t140;
t190 = I_motor.*t140;
t194 = t162.*2.0;
t195 = t163.*2.0;
t196 = t164.*2.0;
t204 = m3.*t17.*t40.*2.0;
t205 = m3.*t17.*t55;
t211 = m_motor.*t18.*t41.*4.0;
t212 = t11.*t96;
t213 = t11.*t97;
t214 = t11.*t98;
t215 = t11.*t99;
t216 = t11.*t100;
t217 = t11.*t101;
t218 = t11.*t103;
t219 = I2.*l_B_m3.*l_OA.*m3.*t45;
t220 = I_motor.*l_B_m3.*l_OA.*m3.*t45;
t221 = Ir.*l_B_m3.*l_OA.*m3.*t45;
t222 = -t140;
t223 = -t141;
t234 = I2.*t141.*2.0;
t236 = t13.*t59.*3.0;
t237 = t11.*t59.*7.0;
t238 = t12.*t59.*6.0;
t239 = I_motor.*t141.*2.0;
t247 = m3.*t17.*t78;
t248 = m2.*t18.*t77;
t249 = m1.*t19.*t77;
t250 = m3.*t18.*t79;
t251 = m1.*t19.*t79;
t252 = m2.*t18.*t81;
t253 = m1.*t19.*t81;
t254 = m_motor.*t17.*t76;
t255 = m3.*t17.*t91;
t256 = m_motor.*t15.*t84;
t257 = m_motor.*t15.*t87;
t258 = m2.*t18.*t92;
t259 = m1.*t19.*t92;
t260 = m3.*t17.*t97;
t261 = m2.*t18.*t98;
t262 = m1.*t19.*t98;
t263 = m_motor.*t17.*t90;
t264 = m_motor.*t17.*t96;
t265 = l_AB.*l_B_m3.*m3.*t4.*t40.*2.0;
t268 = l_AB.*l_B_m3.*m3.*t4.*t57;
t270 = t11.*t149;
t274 = m3.*t15.*t41.*-2.0;
t275 = m2.*t16.*t41.*-2.0;
t277 = N.*t121.*-2.0;
t278 = m2.*t18.*t41.*-2.0;
t279 = m3.*t18.*t41.*-2.0;
t280 = m1.*t19.*t41.*-2.0;
t281 = m_motor.*t15.*t41.*-2.0;
t282 = -t210;
t287 = m3.*t17.*t60.*4.0;
t289 = t11.*t155;
t301 = m3.*t15.*t72;
t302 = m2.*t16.*t72;
t303 = m2.*t18.*t73;
t304 = m3.*t18.*t73;
t305 = m1.*t19.*t73;
t306 = m_motor.*t15.*t72;
t307 = m3.*t15.*t75;
t308 = m2.*t16.*t75;
t309 = m2.*t18.*t75;
t310 = m3.*t18.*t75;
t311 = m1.*t19.*t75;
t312 = m_motor.*t15.*t75;
t317 = I1.*t294;
t318 = I3.*t295;
t319 = I3.*t296;
t320 = I2.*t297;
t321 = I_motor.*t294;
t322 = I_motor.*t297;
t323 = Ir.*t294;
t324 = Ir.*t297;
t325 = m_motor.*t18.*t77.*3.0;
t326 = m_motor.*t18.*t79.*2.0;
t327 = m_motor.*t18.*t81.*2.0;
t331 = l_AB.*l_OA.*m3.*t3.*t41.*-2.0;
t332 = l_A_m2.*l_OA.*m2.*t3.*t41.*-2.0;
t334 = l_AB.*l_OA.*m_motor.*t3.*t41.*-2.0;
t338 = t271.*2.0;
t339 = t272.*2.0;
t340 = t273.*2.0;
t349 = t290.*2.0;
t350 = t291.*2.0;
t351 = t292.*2.0;
t352 = t293.*2.0;
t358 = m3.*t17.*t72.*2.0;
t364 = m3.*t17.*t75.*4.0;
t366 = m_motor.*t18.*t74.*4.0;
t367 = m_motor.*t18.*t138;
t369 = m3.*t15.*t140;
t370 = m2.*t16.*t140;
t374 = m_motor.*t15.*t140;
t375 = I3.*t355;
t383 = l_AB.*l_B_m3.*m3.*t4.*t60.*-2.0;
t384 = l_B_m3.*l_OA.*m2.*m3.*t16.*t45;
t385 = l_B_m3.*l_OA.*m3.*m_motor.*t15.*t45;
t387 = t347.*2.0;
t388 = t348.*2.0;
t389 = l_AB.*l_B_m3.*m3.*t4.*t72.*2.0;
t392 = l_AB.*l_B_m3.*m3.*t4.*t138;
t396 = m3.*t15.*t141.*2.0;
t398 = m2.*t16.*t141.*2.0;
t401 = m3.*t17.*t140.*7.0;
t402 = m3.*t17.*t141.*6.0;
t403 = m_motor.*t15.*t141.*2.0;
t406 = m1.*t19.*t224;
t407 = m_motor.*t17.*t225;
t408 = m_motor.*t17.*t226;
t409 = l_AB.*l_A_m2.*m2.*m3.*t18.*t33.*2.0;
t410 = l_AB.*l_A_m2.*m2.*m_motor.*t18.*t33.*2.0;
t424 = l_B_m3.*l_OA.*t15.*t21.*t45;
t433 = l_AB.*l_B_m3.*m3.*t4.*t140.*6.0;
t435 = l_AB.*l_OA.*m3.*t3.*t141.*2.0;
t437 = l_A_m2.*l_OA.*m2.*t3.*t141.*2.0;
t439 = l_AB.*l_OA.*m_motor.*t3.*t141.*2.0;
t445 = m3.*m_motor.*t15.*t18.*t33.*2.0;
t446 = m_motor.*t18.*t224.*2.0;
t459 = l_AB.*l_A_m2.*m3.*t33.*t84.*2.0;
t460 = l_AB.*l_A_m2.*m_motor.*t33.*t84.*2.0;
t461 = l_A_m2.*l_B_m3.*m2.*m3.*t3.*t18.*t45;
t462 = l_AB.*l_B_m3.*m3.*m_motor.*t3.*t18.*t45;
t465 = Ir.*t11.*t295;
t466 = Ir.*t11.*t296;
t475 = l_AB.*l_B_m3.*m3.*t4.*t141.*8.0;
t480 = t34.*t294;
t481 = t33.*t295;
t482 = t33.*t296;
t483 = t15.*t18.*t23.*t33;
t484 = m2.*t18.*t294;
t485 = m1.*t19.*t294;
t486 = m2.*t16.*t297;
t487 = m3.*t17.*t296;
t489 = Ir.*N.*t347.*-2.0;
t490 = Ir.*N.*t348.*-2.0;
t497 = m_motor.*t18.*t33.*t77.*2.0;
t498 = l_AB.*l_B_m3.*t3.*t18.*t21.*t45;
t499 = l_AB.*l_OA.*t4.*t17.*t21.*t45;
t507 = l_B_m3.*l_OA.*m3.*t45.*t141.*4.0;
t515 = Ir.*t11.*t355;
t517 = m3.*t17.*t355;
t518 = m_motor.*t18.*t294.*3.0;
t545 = l_AB.*l_A_m2.*m_motor.*t33.*t229.*2.0;
t551 = t33.*t426;
t552 = t34.*t426;
t561 = l_AB.*l_A_m2.*m2.*t33.*t297.*2.0;
t575 = m_motor.*t18.*t33.*t294.*2.0;
t608 = l_AB.*l_A_m2.*m2.*t3.*t4.*t45.*t297.*2.0;
t610 = t3.*t4.*t45.*t426.*2.0;
t611 = m_motor.*t3.*t4.*t18.*t45.*t294.*2.0;
t64 = -t56;
t131 = -t66;
t132 = -t69;
t133 = -t70;
t159 = l_BC.*t52;
t160 = l_BC.*t53;
t171 = t127.*2.0;
t177 = -t137;
t189 = I1.*t142;
t192 = t161.*2.0;
t193 = t161.*4.0;
t197 = -t169;
t198 = -t170;
t199 = -t127;
t201 = -t172;
t202 = m3.*t15.*t56;
t203 = m2.*t16.*t56;
t206 = m2.*t18.*t56;
t207 = m3.*t18.*t56;
t208 = m1.*t19.*t56;
t209 = m_motor.*t15.*t56;
t235 = I2.*t142.*2.0;
t240 = I_motor.*t142.*3.0;
t241 = -t161;
t244 = -t194;
t245 = -t195;
t246 = -t196;
t266 = l_AB.*l_OA.*m3.*t3.*t56;
t267 = l_A_m2.*l_OA.*m2.*t3.*t56;
t269 = l_AB.*l_OA.*m_motor.*t3.*t56;
t276 = -t204;
t283 = -t211;
t284 = t212.*2.0;
t285 = t213.*2.0;
t286 = t11.*t153;
t288 = t218.*2.0;
t313 = N.*t221;
t314 = -t234;
t315 = -t238;
t316 = -t239;
t328 = m_motor.*t18.*t150;
t329 = m_motor.*t18.*t153;
t330 = -t265;
t335 = t11.*t184;
t336 = l_AB.*l_B_m3.*m3.*t4.*t66;
t337 = t270.*4.0;
t341 = -t214;
t342 = -t287;
t343 = -t219;
t344 = -t220;
t345 = -t221;
t356 = m3.*t15.*t136;
t357 = m2.*t16.*t136;
t359 = m3.*t17.*t135;
t360 = m2.*t18.*t136;
t361 = m3.*t18.*t136;
t362 = m1.*t19.*t136;
t363 = m_motor.*t15.*t136;
t365 = m_motor.*t18.*t135;
t371 = m2.*t18.*t142;
t372 = m3.*t18.*t142;
t373 = m1.*t19.*t142;
t376 = t11.*t221;
t377 = m3.*t17.*t170;
t380 = m_motor.*t17.*t169;
t381 = m_motor.*t18.*t127.*4.0;
t382 = -t270;
t390 = l_AB.*l_OA.*m3.*t3.*t136;
t391 = l_A_m2.*l_OA.*m2.*t3.*t136;
t393 = l_AB.*l_OA.*m_motor.*t3.*t136;
t397 = m3.*t15.*t142.*2.0;
t399 = m2.*t16.*t142.*2.0;
t400 = m3.*t17.*t142.*3.0;
t404 = m_motor.*t15.*t142.*2.0;
t405 = m_motor.*t18.*t142.*2.0;
t412 = l_AB.*l_B_m3.*m3.*t4.*t128.*2.0;
t413 = l_AB.*l_B_m3.*m3.*t4.*t130.*2.0;
t416 = l_AB.*l_B_m3.*m_motor.*t4.*t129.*4.0;
t417 = N.*t323.*2.0;
t418 = N.*t324.*2.0;
t419 = N.*t260.*-2.0;
t420 = m2.*t18.*t127.*-2.0;
t421 = m1.*t19.*t127.*-2.0;
t422 = N.*t264.*-2.0;
t427 = m2.*t18.*t212;
t428 = m1.*t19.*t212;
t429 = m3.*t18.*t213;
t430 = m1.*t19.*t213;
t431 = m_motor.*t15.*t215;
t432 = m_motor.*t15.*t217;
t434 = l_AB.*l_B_m3.*m3.*t4.*t142.*4.0;
t436 = l_AB.*l_OA.*m3.*t3.*t142.*2.0;
t438 = l_A_m2.*l_OA.*m2.*t3.*t142.*2.0;
t440 = l_AB.*l_OA.*m_motor.*t3.*t142.*2.0;
t441 = -t396;
t442 = -t398;
t443 = -t402;
t444 = -t403;
t447 = -t409;
t448 = -t410;
t449 = Ir.*N.*t387;
t450 = Ir.*N.*t388;
t454 = l_A_m2.*l_OA.*m2.*t3.*t127.*-2.0;
t455 = l_AB.*l_OA.*m_motor.*t3.*t127.*-2.0;
t463 = -t384;
t464 = -t385;
t472 = m_motor.*t18.*t212.*3.0;
t474 = m_motor.*t18.*t214.*4.0;
t476 = -t435;
t477 = -t437;
t478 = -t439;
t488 = -t445;
t491 = t34.*t317;
t492 = t33.*t318;
t493 = t33.*t319;
t494 = t34.*t321;
t495 = I3.*t483;
t496 = t34.*t323;
t500 = -t459;
t501 = -t460;
t502 = t461.*2.0;
t503 = t462.*2.0;
t504 = -t424;
t505 = l_B_m3.*l_OA.*m3.*t45.*t178;
t506 = l_B_m3.*l_OA.*m3.*t45.*t142.*2.0;
t508 = l_AB.*l_B_m3.*m3.*t4.*t215.*2.0;
t509 = l_AB.*l_B_m3.*m3.*t4.*t217.*2.0;
t512 = l_AB.*l_B_m3.*m_motor.*t4.*t216.*4.0;
t513 = t11.*t323.*2.0;
t514 = t11.*t324.*2.0;
t516 = -t475;
t519 = -t497;
t520 = t498.*2.0;
t521 = t499.*2.0;
t522 = -t461;
t524 = -t462;
t526 = Ir.*t11.*t387;
t527 = Ir.*t11.*t388;
t528 = -t507;
t529 = t61.*t297;
t530 = -t480;
t531 = -t481;
t532 = -t482;
t533 = -t483;
t534 = l_AB.*l_A_m2.*m3.*t33.*t215.*2.0;
t535 = l_AB.*l_A_m2.*m_motor.*t33.*t215.*2.0;
t536 = t61.*t320;
t537 = t61.*t322;
t538 = t61.*t324;
t546 = -t498;
t548 = -t499;
t553 = l_A_m2.*l_B_m3.*m3.*t3.*t45.*t128.*2.0;
t554 = l_AB.*l_B_m3.*m_motor.*t3.*t45.*t129.*2.0;
t557 = m2.*t18.*t480;
t558 = m1.*t19.*t480;
t559 = m3.*t17.*t482;
t560 = m3.*t17.*t483;
t562 = -t545;
t568 = t33.*t465;
t569 = t33.*t466;
t570 = Ir.*t11.*t483;
t576 = m_motor.*t18.*t480.*2.0;
t577 = -t561;
t578 = l_A_m2.*l_B_m3.*m3.*t3.*t45.*t215.*2.0;
t579 = l_AB.*l_B_m3.*m_motor.*t3.*t45.*t216.*2.0;
t581 = t61.*t426;
t582 = -t551;
t583 = -t552;
t585 = m_motor.*t18.*t33.*t212.*-2.0;
t586 = t61.*t486;
t587 = m_motor.*t18.*t61.*t294;
t592 = -t575;
t602 = Ir.*t11.*t498.*-2.0;
t603 = Ir.*t11.*t499.*-2.0;
t200 = -t171;
t242 = -t192;
t243 = -t193;
t368 = t313.*2.0;
t378 = m2.*t18.*t171;
t379 = m1.*t19.*t171;
t394 = -t313;
t411 = t376.*2.0;
t414 = l_A_m2.*l_OA.*m2.*t3.*t171;
t415 = l_AB.*l_OA.*m_motor.*t3.*t171;
t423 = -t381;
t451 = t11.*t345;
t452 = -t412;
t453 = -t413;
t456 = -t416;
t457 = -t417;
t458 = -t418;
t467 = t67+t159;
t468 = m3.*t17.*t285;
t469 = m2.*t18.*t286;
t470 = m1.*t19.*t286;
t471 = m_motor.*t17.*t284;
t473 = m_motor.*t18.*t285;
t479 = t68+t160;
t510 = l_A_m2.*l_OA.*m2.*t3.*t286;
t511 = l_AB.*l_OA.*m_motor.*t3.*t286;
t523 = -t502;
t525 = -t503;
t539 = -t491;
t540 = -t492;
t541 = -t493;
t542 = -t494;
t543 = -t495;
t544 = -t496;
t547 = -t520;
t549 = -t521;
t550 = t34.*t417;
t555 = -t534;
t556 = -t535;
t565 = -t529;
t566 = Ir.*N.*t520;
t567 = Ir.*N.*t521;
t571 = m_motor.*t18.*t33.*t284;
t572 = -t536;
t573 = -t537;
t574 = -t538;
t580 = t61.*t418;
t584 = t11.*t496.*2.0;
t588 = m2.*t18.*t530;
t589 = m1.*t19.*t530;
t590 = m3.*t17.*t532;
t591 = m3.*t17.*t533;
t593 = -t576;
t594 = Ir.*t11.*t520;
t595 = Ir.*t11.*t521;
t596 = -t578;
t597 = -t579;
t599 = -t568;
t600 = -t569;
t601 = Ir.*t11.*t533;
t604 = -t581;
t605 = t61.*t514;
t606 = -t586;
t607 = -t587;
t609 = t11.*t538.*-2.0;
t395 = -t368;
t563 = t43+t467;
t564 = t46+t479;
t598 = -t584;
t612 = t26+t29+t32+t64+t73+t75+t77+t79+t81+t89+t92+t98+t136+t142+t200+t212+t213+t218+t224+t231+t242+t286+t294+t335+t530;
t613 = t41+t127+t140+t145+t146+t148+t161+t162+t163+t164+t176+t223+t292+t293+t313+t341+t343+t344+t345+t346+t347+t382+t386+t425+t451+t463+t464+t504+t548;
t614 = t26+t29+t32+t63+t73+t74+t75+t77+t79+t81+t89+t92+t98+t141+t145+t146+t148+t199+t212+t213+t214+t218+t224+t231+t242+t271+t272+t273+t292+t293+t294+t335+t347+t376+t394+t530+t548;
t616 = t25+t29+t32+t40+t42+t60+t63+t74+t80+t84+t85+t87+t92+t98+t128+t129+t130+t141+t144+t145+t146+t147+t148+t149+t152+t162+t163+t164+t173+t199+t214+t222+t229+t230+t241+t270+t290+t291+t292+t293+t297+t300+t343+t344+t345+t346+t347+t348+t353+t386+t425+t463+t464+t504+t522+t524+t546+t548+t565;
t615 = t25+t26+t49+t50+t64+t72+t73+t77+t79+t80+t81+t84+t85+t87+t89+t136+t138+t140+t142+t150+t152+t153+t180+t181+t183+t200+t212+t213+t215+t216+t217+t218+t224+t229+t230+t231+t242+t286+t289+t294+t297+t300+t335+t338+t339+t340+t351+t352+t387+t395+t411+t530+t549+t565;
t617 = t8+t9+t24+t25+t27+t28+t29+t30+t31+t32+t51+t62+t64+t65+t72+t76+t78+t80+t82+t83+t84+t85+t86+t87+t88+t90+t91+t92+t93+t94+t95+t96+t97+t98+t99+t100+t101+t102+t103+t131+t135+t136+t139+t142+t151+t152+t154+t155+t178+t179+t182+t184+t197+t198+t200+t201+t215+t216+t217+t225+t226+t227+t228+t229+t230+t232+t233+t243+t244+t245+t246+t284+t285+t286+t288+t289+t295+t296+t297+t298+t299+t300+t337+t338+t339+t340+t349+t350+t354+t355+t388+t395+t411+t447+t448+t488+t523+t525+t531+t532+t533+t547+t565;
t618 = t35+t36+t37+t38+t39+t54+t58+t59+t104+t105+t106+t107+t108+t109+t110+t111+t112+t113+t114+t115+t116+t117+t118+t119+t120+t121+t122+t123+t124+t125+t126+t132+t133+t134+t143+t156+t157+t158+t165+t166+t167+t168+t174+t175+t177+t185+t186+t187+t188+t189+t190+t191+t235+t236+t237+t240+t247+t248+t249+t250+t251+t252+t253+t254+t255+t256+t257+t258+t259+t260+t261+t262+t263+t264+t274+t275+t276+t277+t278+t279+t280+t281+t282+t283+t301+t302+t303+t304+t305+t306+t307+t308+t309+t310+t311+t312+t314+t315+t316+t317+t318+t319+t320+t321+t322+t323+t324+t325+t326+t327+t328+t329+t330+t331+t332+t333+t334+t342+t356+t357+t358+t359+t360+t361+t362+t363+t364+t365+t366+t367+t369+t370+t371+t372+t373+t374+t375+t383+t389+t390+t391+t392+t393+t397+t399+t400+t401+t404+t405+t406+t407+t408+t419+t420+t421+t422+t423+t426+t427+t428+t429+t430+t431+t432+t433+t434+t436+t438+t440+t441+t442+t443+t444+t446+t452+t453+t454+t455+t456+t457+t458+t465+t466+t468+t469+t470+t471+t472+t473+t474+t476+t477+t478+t484+t485+t486+t487+t489+t490+t500+t501+t505+t506+t508+t509+t510+t511+t512+t513+t514+t515+t516+t517+t518+t519+t526+t527+t528+t539+t540+t541+t542+t543+t544+t550+t553+t554+t555+t556+t562+t566+t567+t572+t573+t574+t577+t580+t582+t583+t585+t588+t589+t590+t591+t592+t593+t596+t597+t598+t599+t600+t601+t602+t603+t604+t606+t607+t608+t609+t610+t611;
t619 = 1.0./t618;
t620 = t612.*t619;
t623 = t613.*t619;
t626 = t614.*t619;
t636 = t615.*t619;
t642 = t616.*t619;
t650 = t617.*t619;
t621 = t563.*t620;
t622 = t564.*t620;
t624 = t159.*t623;
t625 = t160.*t623;
t627 = -t626;
t628 = t563.*t623;
t629 = t564.*t623;
t630 = t479.*t626;
t631 = t467.*t626;
t634 = t563.*t626;
t635 = t564.*t626;
t637 = -t636;
t638 = t467.*t636;
t639 = t479.*t636;
t643 = -t642;
t644 = t159.*t642;
t645 = t160.*t642;
t646 = t467.*t642;
t647 = t479.*t642;
t651 = t159.*t650;
t652 = t160.*t650;
t632 = t479.*t627;
t633 = t467.*t627;
t640 = t467.*t637;
t641 = t479.*t637;
t648 = t467.*t643;
t649 = t479.*t643;
t653 = t620+t623+t627;
t656 = t626+t637+t642;
t659 = t623+t643+t650;
t654 = t621+t624+t633;
t655 = t622+t625+t632;
t657 = t634+t640+t644;
t658 = t635+t641+t645;
t660 = t628+t648+t651;
t661 = t629+t649+t652;
Mass_Op_Sp_inv = reshape([t160.*t661-t479.*t658+t564.*t655,-t160.*t660+t479.*t657-t564.*t654,-t160.*t659+t479.*t656-t564.*t653,-t159.*t661+t467.*t658-t563.*t655,t159.*t660-t467.*t657+t563.*t654,t159.*t659-t467.*t656+t563.*t653,-t622-t625-t629+t630+t647-t652+t658,t638+t654+t660+t159.*t643+t563.*t627,t620+t623.*2.0-t626.*2.0+t636-t642.*2.0+t650],[3,3]);
