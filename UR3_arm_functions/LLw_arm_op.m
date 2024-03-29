function Mass_Rot_Op_Sp = LLw_arm_op(in1,in2)
%LLW_ARM_OP
%    MASS_ROT_OP_SP = LLW_ARM_OP(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    12-May-2021 19:41:44

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
m_motor = in2(4,:);
th2 = in1(2,:);
th3 = in1(3,:);
t2 = cos(th2);
t3 = cos(th3);
t4 = th2+th3;
t5 = I_motor.^2;
t6 = Ir.^2;
t7 = Ir.^3;
t8 = N.^2;
t9 = N.^3;
t12 = l_AB.^2;
t13 = l_A_m2.^2;
t14 = l_B_m3.^2;
t15 = l_OA.^2;
t16 = l_O_m1.^2;
t17 = m2.^2;
t18 = m3.^2;
t19 = m3.^3;
t20 = m_motor.^2;
t21 = I1.*I3;
t22 = I2.*I3;
t23 = I3.*I_motor;
t24 = I3.*Ir;
t10 = t8.^2;
t11 = t8.^3;
t25 = t2.^2;
t26 = t3.^2;
t27 = I2.*t21;
t28 = I_motor.*t21;
t29 = I_motor.*t22;
t30 = Ir.*t21;
t31 = Ir.*t22;
t32 = N.*t24;
t33 = cos(t4);
t34 = Ir.*t23.*2.0;
t36 = I3.*t5;
t37 = I3.*t6;
t43 = Ir.*N.*t23.*4.0;
t44 = I1.*Ir.*t8;
t45 = I2.*Ir.*t8;
t46 = t8.*t24;
t47 = I_motor.*Ir.*t8;
t48 = I3.*m3.*t12;
t49 = I3.*m2.*t13;
t50 = I1.*m3.*t14;
t51 = I2.*m3.*t14;
t52 = I3.*m2.*t15;
t53 = I3.*m3.*t15;
t54 = I3.*m1.*t16;
t55 = I3.*m_motor.*t12;
t56 = I_motor.*m3.*t14;
t57 = Ir.*m3.*t14;
t58 = m3.*t12.*t21;
t59 = m2.*t13.*t21;
t61 = m2.*t15.*t22;
t62 = m3.*t15.*t22;
t63 = m1.*t16.*t22;
t64 = m_motor.*t12.*t21;
t65 = m3.*t12.*t23;
t66 = m2.*t13.*t23;
t69 = m2.*t15.*t23;
t70 = m3.*t15.*t23;
t71 = m1.*t16.*t23;
t72 = m3.*t12.*t24;
t73 = m2.*t13.*t24;
t76 = m2.*t15.*t24;
t77 = m3.*t15.*t24;
t78 = m1.*t16.*t24;
t79 = m_motor.*t12.*t23;
t80 = m_motor.*t12.*t24;
t87 = t6.*t8;
t88 = t6.*t9;
t91 = I3.*l_AB.*l_OA.*m3.*t2;
t92 = I3.*l_A_m2.*l_OA.*m2.*t2;
t93 = I3.*l_AB.*l_OA.*m_motor.*t2;
t94 = I3.*m_motor.*t15.*2.0;
t98 = Ir.*N.*l_AB.*l_B_m3.*m3.*t3;
t99 = Ir.*N.*l_AB.*l_OA.*m3.*t2;
t100 = Ir.*N.*l_A_m2.*l_OA.*m2.*t2;
t101 = Ir.*N.*l_AB.*l_OA.*m_motor.*t2;
t102 = m_motor.*t15.*t22.*2.0;
t103 = m_motor.*t15.*t23.*2.0;
t105 = m_motor.*t15.*t24.*2.0;
t107 = m3.*t5.*t14;
t108 = m3.*t6.*t14;
t112 = Ir.*t8.*t23.*4.0;
t116 = Ir.*t5.*t8;
t130 = Ir.*m3.*t8.*t12;
t131 = Ir.*m2.*t8.*t13;
t133 = Ir.*m2.*t8.*t15;
t134 = Ir.*m3.*t8.*t15;
t135 = Ir.*m1.*t8.*t16;
t136 = Ir.*m_motor.*t8.*t12;
t140 = m2.*m3.*t13.*t14;
t141 = m2.*m3.*t14.*t15;
t142 = m1.*m3.*t14.*t16;
t143 = m3.*m_motor.*t12.*t14;
t175 = Ir.*l_AB.*l_B_m3.*m3.*t3.*t8;
t188 = Ir.*m_motor.*t8.*t15.*2.0;
t189 = l_A_m2.*l_OA.*m2.*m3.*t2.*t14;
t190 = l_AB.*l_OA.*m3.*m_motor.*t2.*t14;
t191 = t12.*t14.*t18;
t192 = t14.*t15.*t18;
t193 = m3.*m_motor.*t14.*t15.*2.0;
t211 = I3.*t12.*t15.*t18;
t212 = I3.*t13.*t15.*t17;
t229 = N.*l_AB.*l_B_m3.*m3.*t3.*t6.*2.0;
t230 = Ir.*l_AB.*l_OA.*m3.*t2.*t8.*2.0;
t231 = Ir.*l_A_m2.*l_OA.*m2.*t2.*t8.*2.0;
t232 = Ir.*l_AB.*l_OA.*m_motor.*t2.*t8.*2.0;
t237 = l_AB.*l_A_m2.*l_B_m3.*l_OA.*m2.*m3.*t2.*t3;
t238 = l_AB.*l_OA.*t2.*t14.*t18;
t258 = I3.*t12.*t15.*t20.*2.0;
t268 = l_B_m3.*l_OA.*m3.*m_motor.*t2.*t3.*t12;
t302 = l_B_m3.*l_OA.*t2.*t3.*t12.*t18;
t303 = t12.*t14.*t15.*t19;
t324 = Ir.*N.*l_AB.*l_B_m3.*t3.*t15.*t18.*2.0;
t336 = Ir.*t8.*t12.*t15.*t18;
t337 = Ir.*t8.*t13.*t15.*t17;
t353 = m3.*t13.*t14.*t15.*t17;
t377 = Ir.*t8.*t12.*t15.*t20.*2.0;
t379 = m3.*t12.*t14.*t15.*t20.*2.0;
t383 = Ir.*l_AB.*l_B_m3.*t3.*t8.*t15.*t18.*2.0;
t35 = t32.*2.0;
t38 = t33.^2;
t39 = -t32;
t41 = N.*t30.*2.0;
t42 = N.*t31.*2.0;
t60 = I2.*t50;
t67 = I_motor.*t50;
t68 = I_motor.*t51;
t74 = Ir.*t50;
t75 = Ir.*t51;
t81 = N.*t57;
t84 = -t43;
t85 = t46.*2.0;
t86 = N.*t37.*4.0;
t89 = t6.*t10;
t90 = t7.*t11;
t95 = I2.*t44;
t96 = I_motor.*t44;
t97 = I_motor.*t45;
t104 = Ir.*t56.*2.0;
t110 = t8.*t30.*2.0;
t111 = t8.*t31.*2.0;
t113 = I2.*t87;
t115 = I_motor.*t87;
t117 = t98.*2.0;
t128 = Ir.*N.*t56.*4.0;
t129 = m_motor.*t15.*t32.*4.0;
t132 = t8.*t57;
t137 = I2.*l_B_m3.*l_OA.*m3.*t33;
t138 = I_motor.*l_B_m3.*l_OA.*m3.*t33;
t139 = Ir.*l_B_m3.*l_OA.*m3.*t33;
t144 = I2.*t88.*2.0;
t146 = t10.*t37.*3.0;
t147 = t8.*t37.*7.0;
t148 = t9.*t37.*6.0;
t149 = I_motor.*t88.*2.0;
t152 = m2.*t13.*t50;
t153 = m2.*t15.*t48;
t154 = m1.*t16.*t48;
t155 = m3.*t15.*t49;
t156 = m1.*t16.*t49;
t157 = m2.*t15.*t51;
t158 = m1.*t16.*t51;
t159 = m_motor.*t12.*t50;
t160 = m2.*t13.*t56;
t161 = m_motor.*t12.*t52;
t162 = m_motor.*t12.*t54;
t163 = m2.*t15.*t56;
t164 = m1.*t16.*t56;
t165 = m2.*t13.*t57;
t166 = m2.*t15.*t57;
t167 = m1.*t16.*t57;
t168 = m_motor.*t12.*t56;
t169 = m_motor.*t12.*t57;
t176 = m3.*t12.*t32.*-2.0;
t177 = m2.*t13.*t32.*-2.0;
t180 = m2.*t15.*t32.*-2.0;
t181 = m3.*t15.*t32.*-2.0;
t182 = m1.*t16.*t32.*-2.0;
t183 = m_motor.*t12.*t32.*-2.0;
t187 = N.*t108.*4.0;
t194 = m3.*t12.*t44;
t195 = m2.*t13.*t44;
t196 = m2.*t15.*t45;
t197 = m3.*t15.*t45;
t198 = m1.*t16.*t45;
t199 = m_motor.*t12.*t44;
t200 = m3.*t12.*t47;
t201 = m2.*t13.*t47;
t202 = m2.*t15.*t47;
t203 = m3.*t15.*t47;
t204 = m1.*t16.*t47;
t205 = m_motor.*t12.*t47;
t210 = I1.*t191;
t213 = I2.*t192;
t214 = I_motor.*t191;
t215 = I_motor.*t192;
t216 = Ir.*t191;
t217 = Ir.*t192;
t218 = m_motor.*t15.*t48.*3.0;
t219 = m_motor.*t15.*t49.*2.0;
t220 = m_motor.*t15.*t51.*2.0;
t221 = m_motor.*t15.*t56.*2.0;
t222 = m_motor.*t15.*t57.*2.0;
t223 = I1.*t98.*-2.0;
t224 = l_AB.*l_OA.*m3.*t2.*t32.*-2.0;
t225 = l_A_m2.*l_OA.*m2.*t2.*t32.*-2.0;
t226 = I_motor.*t98.*-2.0;
t227 = l_AB.*l_OA.*m_motor.*t2.*t32.*-2.0;
t228 = t175.*2.0;
t241 = m3.*t14.*t44.*2.0;
t242 = m3.*t14.*t45.*2.0;
t247 = m3.*t14.*t47.*4.0;
t248 = m_motor.*t15.*t45.*2.0;
t249 = m_motor.*t15.*t46.*4.0;
t250 = m_motor.*t15.*t47.*2.0;
t252 = m3.*t12.*t87;
t253 = m2.*t13.*t87;
t257 = m_motor.*t12.*t87;
t265 = -t229;
t266 = l_B_m3.*l_OA.*m2.*m3.*t13.*t33;
t267 = l_B_m3.*l_OA.*m3.*m_motor.*t12.*t33;
t269 = l_AB.*l_B_m3.*m3.*t3.*t44.*2.0;
t272 = l_AB.*l_B_m3.*m3.*t3.*t47.*2.0;
t275 = m3.*t12.*t88.*2.0;
t277 = m2.*t13.*t88.*2.0;
t280 = m3.*t14.*t87.*7.0;
t281 = m3.*t14.*t88.*6.0;
t282 = m_motor.*t12.*t88.*2.0;
t285 = m1.*t16.*t140;
t286 = m_motor.*t12.*t141;
t287 = m_motor.*t12.*t142;
t293 = m_motor.*t15.*t98.*4.0;
t301 = l_B_m3.*l_OA.*t12.*t18.*t33;
t304 = m2.*t15.*t130;
t305 = m1.*t16.*t130;
t306 = m3.*t15.*t131;
t307 = m1.*t16.*t131;
t308 = m_motor.*t12.*t133;
t309 = m_motor.*t12.*t135;
t310 = l_AB.*l_B_m3.*m3.*t3.*t87.*6.0;
t312 = l_AB.*l_OA.*m3.*t2.*t88.*2.0;
t314 = l_A_m2.*l_OA.*m2.*t2.*t88.*2.0;
t316 = l_AB.*l_OA.*m_motor.*t2.*t88.*2.0;
t322 = m_motor.*t15.*t140.*2.0;
t323 = Ir.*N.*t238.*2.0;
t325 = m2.*t15.*t98.*-2.0;
t326 = m1.*t16.*t98.*-2.0;
t332 = l_AB.*l_A_m2.*m3.*t25.*t52.*2.0;
t333 = l_AB.*l_A_m2.*m_motor.*t25.*t52.*2.0;
t342 = m_motor.*t15.*t130.*3.0;
t343 = m_motor.*t15.*t131.*2.0;
t345 = l_AB.*l_B_m3.*m3.*t3.*t88.*8.0;
t349 = t26.*t191;
t350 = m2.*t15.*t191;
t351 = m1.*t16.*t191;
t352 = m2.*t13.*t192;
t355 = -t324;
t357 = t25.*t211;
t358 = t25.*t212;
t360 = I3.*t12.*t15.*t20.*t25;
t362 = m_motor.*t15.*t25.*t48.*2.0;
t363 = l_AB.*l_OA.*t3.*t14.*t18.*t33;
t367 = l_B_m3.*l_OA.*m3.*t33.*t87.*2.0;
t369 = l_B_m3.*l_OA.*m3.*t33.*t88.*4.0;
t370 = l_AB.*l_B_m3.*m3.*t3.*t133.*2.0;
t371 = l_AB.*l_B_m3.*m3.*t3.*t135.*2.0;
t374 = l_AB.*l_B_m3.*m_motor.*t3.*t134.*4.0;
t380 = m_motor.*t15.*t191.*3.0;
t382 = Ir.*t8.*t238.*2.0;
t387 = l_AB.*l_A_m2.*m3.*t25.*t133.*2.0;
t388 = l_AB.*l_A_m2.*m_motor.*t25.*t133.*2.0;
t398 = l_AB.*l_A_m2.*m_motor.*t25.*t141.*2.0;
t401 = t25.*t303;
t402 = t26.*t303;
t403 = Ir.*N.*l_A_m2.*l_B_m3.*m2.*m3.*t2.*t15.*t33.*2.0;
t404 = Ir.*N.*l_AB.*l_B_m3.*m3.*m_motor.*t2.*t15.*t33.*2.0;
t409 = t25.*t353;
t410 = m3.*t12.*t14.*t15.*t20.*t25;
t411 = l_AB.*l_A_m2.*m2.*t25.*t192.*2.0;
t414 = Ir.*N.*l_AB.*l_B_m3.*t2.*t15.*t18.*t33.*2.0;
t416 = t25.*t336;
t417 = t25.*t337;
t418 = Ir.*t8.*t12.*t15.*t20.*t25;
t419 = m_motor.*t15.*t25.*t130.*2.0;
t423 = m_motor.*t15.*t25.*t191.*2.0;
t426 = l_A_m2.*l_B_m3.*m3.*t2.*t33.*t133.*2.0;
t427 = l_AB.*l_B_m3.*m_motor.*t2.*t33.*t134.*2.0;
t442 = Ir.*l_AB.*l_B_m3.*t2.*t8.*t15.*t18.*t33.*2.0;
t456 = l_AB.*l_A_m2.*m2.*t2.*t3.*t33.*t192.*2.0;
t458 = t2.*t3.*t33.*t303.*2.0;
t459 = m_motor.*t2.*t3.*t15.*t33.*t191.*2.0;
t40 = -t35;
t82 = -t41;
t83 = -t42;
t106 = t81.*2.0;
t109 = -t86;
t114 = I1.*t89;
t118 = -t81;
t120 = m3.*t12.*t35;
t121 = m2.*t13.*t35;
t122 = N.*t74.*2.0;
t123 = N.*t75.*2.0;
t124 = m2.*t15.*t35;
t125 = m3.*t15.*t35;
t126 = m1.*t16.*t35;
t127 = m_motor.*t12.*t35;
t145 = I2.*t89.*2.0;
t150 = I_motor.*t89.*3.0;
t151 = -t117;
t170 = I1.*t117;
t171 = l_AB.*l_OA.*m3.*t2.*t35;
t172 = l_A_m2.*l_OA.*m2.*t2.*t35;
t173 = I_motor.*t117;
t174 = l_AB.*l_OA.*m_motor.*t2.*t35;
t184 = -t128;
t185 = -t129;
t186 = t132.*2.0;
t206 = N.*t139;
t207 = -t144;
t208 = -t148;
t209 = -t149;
t233 = -t187;
t234 = -t137;
t235 = -t138;
t236 = -t139;
t239 = m3.*t12.*t85;
t240 = m2.*t13.*t85;
t243 = m2.*t15.*t85;
t244 = m3.*t15.*t85;
t245 = m1.*t16.*t85;
t246 = m_motor.*t12.*t85;
t254 = m2.*t15.*t89;
t255 = m3.*t15.*t89;
t256 = m1.*t16.*t89;
t259 = t8.*t139;
t264 = m_motor.*t15.*t81.*4.0;
t270 = l_AB.*l_OA.*m3.*t2.*t85;
t271 = l_A_m2.*l_OA.*m2.*t2.*t85;
t273 = l_AB.*l_OA.*m_motor.*t2.*t85;
t276 = m3.*t12.*t89.*2.0;
t278 = m2.*t13.*t89.*2.0;
t279 = m3.*t14.*t89.*3.0;
t283 = m_motor.*t12.*t89.*2.0;
t284 = m_motor.*t15.*t89.*2.0;
t289 = m2.*t15.*t117;
t290 = m1.*t16.*t117;
t294 = N.*t216.*2.0;
t295 = N.*t217.*2.0;
t296 = m2.*t13.*t81.*-2.0;
t297 = m2.*t15.*t81.*-2.0;
t298 = m1.*t16.*t81.*-2.0;
t299 = m_motor.*t12.*t81.*-2.0;
t311 = l_AB.*l_B_m3.*m3.*t3.*t89.*4.0;
t313 = l_AB.*l_OA.*m3.*t2.*t89.*2.0;
t315 = l_A_m2.*l_OA.*m2.*t2.*t89.*2.0;
t317 = l_AB.*l_OA.*m_motor.*t2.*t89.*2.0;
t318 = -t275;
t319 = -t277;
t320 = -t281;
t321 = -t282;
t327 = l_A_m2.*l_OA.*m2.*t2.*t81.*-2.0;
t328 = l_AB.*l_OA.*m_motor.*t2.*t81.*-2.0;
t329 = -t293;
t334 = -t266;
t335 = -t267;
t344 = m_motor.*t15.*t132.*4.0;
t346 = -t312;
t347 = -t314;
t348 = -t316;
t354 = -t323;
t356 = t26.*t210;
t359 = t26.*t214;
t361 = t26.*t216;
t364 = -t332;
t365 = -t333;
t366 = -t301;
t368 = l_B_m3.*l_OA.*m3.*t33.*t89.*2.0;
t375 = t8.*t216.*2.0;
t376 = t8.*t217.*2.0;
t378 = -t345;
t381 = -t362;
t384 = -t369;
t385 = t38.*t192;
t386 = -t349;
t389 = t38.*t213;
t390 = t38.*t215;
t391 = t38.*t217;
t393 = -t357;
t394 = -t358;
t396 = -t360;
t399 = -t363;
t405 = -t387;
t406 = -t388;
t407 = m2.*t15.*t349;
t408 = m1.*t16.*t349;
t412 = -t398;
t415 = Ir.*N.*t363.*2.0;
t424 = m_motor.*t15.*t349.*2.0;
t425 = -t411;
t429 = t38.*t303;
t430 = -t401;
t431 = -t402;
t433 = -t419;
t434 = t38.*t352;
t435 = m_motor.*t15.*t38.*t191;
t438 = -t409;
t439 = -t410;
t440 = -t423;
t443 = Ir.*t8.*t363.*2.0;
t444 = -t426;
t445 = -t427;
t447 = -t416;
t448 = -t417;
t449 = -t418;
t450 = -t442;
t119 = -t106;
t178 = -t122;
t179 = -t123;
t251 = t206.*2.0;
t260 = m2.*t13.*t106;
t261 = m2.*t15.*t106;
t262 = m1.*t16.*t106;
t263 = m_motor.*t12.*t106;
t288 = t259.*2.0;
t291 = l_A_m2.*l_OA.*m2.*t2.*t106;
t292 = l_AB.*l_OA.*m_motor.*t2.*t106;
t300 = -t264;
t330 = -t294;
t331 = -t295;
t338 = m2.*t13.*t186;
t339 = m2.*t15.*t186;
t340 = m1.*t16.*t186;
t341 = m_motor.*t12.*t186;
t372 = l_A_m2.*l_OA.*m2.*t2.*t186;
t373 = l_AB.*l_OA.*m_motor.*t2.*t186;
t392 = -t356;
t395 = -t359;
t397 = -t361;
t400 = t26.*t294;
t413 = -t385;
t420 = -t389;
t421 = -t390;
t422 = -t391;
t428 = t38.*t295;
t432 = t8.*t361.*2.0;
t436 = m2.*t15.*t386;
t437 = m1.*t16.*t386;
t441 = -t424;
t451 = -t443;
t452 = -t429;
t453 = t38.*t376;
t454 = -t434;
t455 = -t435;
t457 = t8.*t391.*-2.0;
t274 = -t251;
t446 = -t432;
t460 = t27+t28+t29+t30+t31+t34+t36+t37+t58+t59+t60+t61+t62+t63+t64+t65+t66+t67+t68+t69+t70+t71+t72+t73+t74+t75+t76+t77+t78+t79+t80+t82+t83+t84+t90+t95+t96+t97+t102+t103+t104+t105+t107+t108+t109+t110+t111+t112+t113+t114+t115+t116+t145+t146+t147+t150+t152+t153+t154+t155+t156+t157+t158+t159+t160+t161+t162+t163+t164+t165+t166+t167+t168+t169+t176+t177+t178+t179+t180+t181+t182+t183+t184+t185+t194+t195+t196+t197+t198+t199+t200+t201+t202+t203+t204+t205+t207+t208+t209+t210+t211+t212+t213+t214+t215+t216+t217+t218+t219+t220+t221+t222+t223+t224+t225+t226+t227+t233+t239+t240+t241+t242+t243+t244+t245+t246+t247+t248+t249+t250+t252+t253+t254+t255+t256+t257+t258+t265+t269+t270+t271+t272+t273+t276+t278+t279+t280+t283+t284+t285+t286+t287+t296+t297+t298+t299+t300+t303+t304+t305+t306+t307+t308+t309+t310+t311+t313+t315+t317+t318+t319+t320+t321+t322+t325+t326+t327+t328+t329+t330+t331+t336+t337+t338+t339+t340+t341+t342+t343+t344+t346+t347+t348+t350+t351+t352+t353+t354+t355+t364+t365+t367+t368+t370+t371+t372+t373+t374+t375+t376+t377+t378+t379+t380+t381+t382+t383+t384+t392+t393+t394+t395+t396+t397+t400+t403+t404+t405+t406+t412+t414+t415+t420+t421+t422+t425+t428+t430+t431+t433+t436+t437+t438+t439+t440+t441+t444+t445+t446+t447+t448+t449+t450+t451+t452+t454+t455+t456+t457+t458+t459;
t461 = 1.0./t460;
Mass_Rot_Op_Sp = 1.0./(t461.*(t21+t22+t23.*2.0+t24.*2.0+t40+t44+t45+t47.*2.0+t48+t49+t50+t51+t52+t53+t54+t55+t56.*2.0+t57.*2.0+t85+t87+t89+t91.*2.0+t92.*2.0+t93.*2.0+t94+t119+t130+t131+t133+t134+t135+t136+t140+t141+t142+t143+t151+t186+t188+t189.*2.0+t190.*2.0+t191+t192+t193+t228+t230+t231+t232+t238.*2.0+t274+t288-t363.*2.0+t386+t413)-t461.*(t22+t23+t24+t39+t45+t46+t47+t48+t49+t51+t55+t56+t57+t88+t91+t92+t93+t118+t130+t131+t132+t136+t140+t143+t151+t189+t190+t191-t206+t228+t238+t259+t386+t399+Ir.*l_AB.*l_OA.*m3.*t2.*t8+Ir.*l_A_m2.*l_OA.*m2.*t2.*t8+Ir.*l_AB.*l_OA.*m_motor.*t2.*t8).*2.0+t461.*(t22+t23+t24+t40+t45+t47+t48+t49+t51+t55+t56+t57+t85+t89+t119+t130+t131+t136+t140+t143+t151+t186+t191+t228+t386)-t461.*(t21+t23+t24+t39+t46+t50+t52+t53+t54+t56+t57-t87+t88+t91+t92+t93+t94-t98+t99+t100+t101+t118+t132+t141+t142+t175+t189+t190+t192+t193+t234+t235+t236+t237+t238+t268+t302+t334+t335+t366+t399+t413+N.*t6+I1.*Ir.*N+I_motor.*Ir.*N+Ir.*N.*m1.*t16+Ir.*N.*m2.*t15+Ir.*N.*m3.*t15+Ir.*N.*m_motor.*t15.*2.0+I1.*l_AB.*l_B_m3.*m3.*t3+I_motor.*l_AB.*l_B_m3.*m3.*t3+Ir.*l_AB.*l_B_m3.*m3.*t3+l_AB.*l_B_m3.*t3.*t15.*t18-l_AB.*l_B_m3.*t2.*t15.*t18.*t33+l_AB.*l_B_m3.*m1.*m3.*t3.*t16+l_AB.*l_B_m3.*m2.*m3.*t3.*t15+l_AB.*l_B_m3.*m3.*m_motor.*t3.*t15.*2.0-l_A_m2.*l_B_m3.*m2.*m3.*t2.*t15.*t33-l_AB.*l_B_m3.*m3.*m_motor.*t2.*t15.*t33).*2.0+t461.*(t32-t46+t81+t87-t88+t91+t92+t93+t98+t99+t100+t101-t132-t175+t189+t190+t206+t234+t235+t236+t237+t238+t268+t302+t334+t335+t366+t399+t8.*t236).*2.0+t461.*(t5+t6+t21+t23+t24+t40+t44+t45.*2.0+t47.*3.0+t50+t52+t53+t54+t56+t57+t85+t87.*2.0+t89+t94-t98.*4.0-t99.*2.0-t100.*2.0-t101.*2.0+t119+t130.*2.0+t131.*2.0+t133+t134+t135+t136.*2.0+t141+t142+t175.*4.0+t186+t188+t192+t193+t230+t231+t232+t274+t288+t413-N.*t6.*2.0+I1.*I2+I1.*I_motor+I2.*I_motor+I1.*Ir+I2.*Ir+I_motor.*Ir.*2.0-I2.*Ir.*N.*2.0-I_motor.*Ir.*N.*2.0+I1.*m2.*t13+I1.*m3.*t12+I2.*m1.*t16+I2.*m2.*t15+I2.*m3.*t15+I1.*m_motor.*t12+I_motor.*m2.*t13+I_motor.*m3.*t12+I2.*m_motor.*t15.*2.0+I_motor.*m1.*t16+I_motor.*m2.*t15+I_motor.*m3.*t15+Ir.*m2.*t13+Ir.*m3.*t12+Ir.*m1.*t16+Ir.*m2.*t15+Ir.*m3.*t15+I_motor.*m_motor.*t12+I_motor.*m_motor.*t15.*2.0+Ir.*m_motor.*t12+Ir.*m_motor.*t15.*2.0+t12.*t15.*t18+t13.*t15.*t17+t12.*t15.*t20.*2.0+m1.*m2.*t13.*t16+m1.*m3.*t12.*t16+m2.*m3.*t12.*t15+m2.*m3.*t13.*t15+m1.*m_motor.*t12.*t16+m2.*m_motor.*t12.*t15+m2.*m_motor.*t13.*t15.*2.0+m3.*m_motor.*t12.*t15.*3.0-t12.*t15.*t18.*t25-t13.*t15.*t17.*t25-t12.*t15.*t20.*t25-Ir.*N.*m2.*t13.*2.0-Ir.*N.*m3.*t12.*2.0-Ir.*N.*m_motor.*t12.*2.0+I1.*l_AB.*l_B_m3.*m3.*t3.*2.0+I_motor.*l_AB.*l_B_m3.*m3.*t3.*2.0+Ir.*l_AB.*l_B_m3.*m3.*t3.*2.0+l_AB.*l_B_m3.*t3.*t15.*t18.*2.0-m3.*m_motor.*t12.*t15.*t25.*2.0-l_AB.*l_B_m3.*t2.*t15.*t18.*t33.*2.0-l_AB.*l_A_m2.*m2.*m3.*t15.*t25.*2.0+l_AB.*l_B_m3.*m1.*m3.*t3.*t16.*2.0+l_AB.*l_B_m3.*m2.*m3.*t3.*t15.*2.0-l_AB.*l_A_m2.*m2.*m_motor.*t15.*t25.*2.0+l_AB.*l_B_m3.*m3.*m_motor.*t3.*t15.*4.0-l_A_m2.*l_B_m3.*m2.*m3.*t2.*t15.*t33.*2.0-l_AB.*l_B_m3.*m3.*m_motor.*t2.*t15.*t33.*2.0));
