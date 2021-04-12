function Mass_Rot_Op_Sp_inv = LLw_arm_op_inv(in1,in2)
%LLW_ARM_OP_INV
%    MASS_ROT_OP_SP_INV = LLW_ARM_OP_INV(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    15-Mar-2021 22:38:13

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
t7 = N.^2;
t8 = N.^3;
t10 = l_AB.^2;
t11 = l_A_m2.^2;
t12 = l_B_m3.^2;
t13 = l_OA.^2;
t14 = l_O_m1.^2;
t15 = m2.^2;
t16 = m3.^2;
t17 = m3.^3;
t18 = m_motor.^2;
t19 = I1.*I3;
t20 = I2.*I3;
t21 = I3.*I_motor;
t22 = I3.*Ir;
t9 = t7.^2;
t23 = t2.^2;
t24 = t3.^2;
t25 = I2.*t19;
t26 = I_motor.*t19;
t27 = I_motor.*t20;
t28 = Ir.*t19;
t29 = Ir.*t20;
t30 = N.*t22;
t31 = cos(t4);
t32 = Ir.*t21.*2.0;
t34 = I3.*t5;
t35 = I3.*t6;
t41 = Ir.*N.*t21.*4.0;
t42 = I1.*Ir.*t7;
t43 = I2.*Ir.*t7;
t44 = t7.*t22;
t45 = I_motor.*Ir.*t7;
t46 = I3.*m3.*t10;
t47 = I3.*m2.*t11;
t48 = I1.*m3.*t12;
t49 = I2.*m3.*t12;
t50 = I3.*m2.*t13;
t51 = I3.*m3.*t13;
t52 = I3.*m1.*t14;
t53 = I3.*m_motor.*t10;
t54 = I_motor.*m3.*t12;
t55 = Ir.*m3.*t12;
t56 = m3.*t10.*t19;
t57 = m2.*t11.*t19;
t59 = m2.*t13.*t20;
t60 = m3.*t13.*t20;
t61 = m1.*t14.*t20;
t62 = m_motor.*t10.*t19;
t63 = m3.*t10.*t21;
t64 = m2.*t11.*t21;
t67 = m2.*t13.*t21;
t68 = m3.*t13.*t21;
t69 = m1.*t14.*t21;
t70 = m3.*t10.*t22;
t71 = m2.*t11.*t22;
t74 = m2.*t13.*t22;
t75 = m3.*t13.*t22;
t76 = m1.*t14.*t22;
t77 = m_motor.*t10.*t21;
t78 = m_motor.*t10.*t22;
t85 = t6.*t7;
t86 = t6.*t8;
t87 = I3.*l_AB.*l_OA.*m3.*t2;
t88 = I3.*l_A_m2.*l_OA.*m2.*t2;
t89 = I3.*l_AB.*l_OA.*m_motor.*t2;
t90 = I3.*m_motor.*t13.*2.0;
t95 = Ir.*N.*l_AB.*l_B_m3.*m3.*t3;
t96 = Ir.*N.*l_AB.*l_OA.*m3.*t2;
t97 = Ir.*N.*l_A_m2.*l_OA.*m2.*t2;
t98 = Ir.*N.*l_AB.*l_OA.*m_motor.*t2;
t99 = m_motor.*t13.*t20.*2.0;
t100 = m_motor.*t13.*t21.*2.0;
t102 = m_motor.*t13.*t22.*2.0;
t104 = m3.*t5.*t12;
t105 = m3.*t6.*t12;
t108 = Ir.*t7.*t21.*3.0;
t114 = Ir.*t5.*t7;
t128 = Ir.*m3.*t7.*t10;
t129 = Ir.*m2.*t7.*t11;
t131 = Ir.*m2.*t7.*t13;
t132 = Ir.*m3.*t7.*t13;
t133 = Ir.*m1.*t7.*t14;
t134 = Ir.*m_motor.*t7.*t10;
t138 = m2.*m3.*t11.*t12;
t139 = m2.*m3.*t12.*t13;
t140 = m1.*m3.*t12.*t14;
t141 = m3.*m_motor.*t10.*t12;
t182 = Ir.*m_motor.*t7.*t13.*2.0;
t183 = l_A_m2.*l_OA.*m2.*m3.*t2.*t12;
t184 = l_AB.*l_OA.*m3.*m_motor.*t2.*t12;
t185 = t10.*t12.*t16;
t186 = t12.*t13.*t16;
t187 = m3.*m_motor.*t12.*t13.*2.0;
t209 = I3.*t10.*t13.*t16;
t210 = I3.*t11.*t13.*t15;
t226 = Ir.*l_AB.*l_B_m3.*m3.*t3.*t7.*2.0;
t227 = N.*l_AB.*l_B_m3.*m3.*t3.*t6.*2.0;
t228 = Ir.*l_AB.*l_OA.*m3.*t2.*t7.*2.0;
t229 = Ir.*l_A_m2.*l_OA.*m2.*t2.*t7.*2.0;
t230 = Ir.*l_AB.*l_OA.*m_motor.*t2.*t7.*2.0;
t235 = l_AB.*l_A_m2.*l_B_m3.*l_OA.*m2.*m3.*t2.*t3;
t236 = l_AB.*l_OA.*t2.*t12.*t16;
t256 = I3.*t10.*t13.*t18.*2.0;
t266 = l_B_m3.*l_OA.*m3.*m_motor.*t2.*t3.*t10;
t296 = l_B_m3.*l_OA.*t2.*t3.*t10.*t16;
t297 = t10.*t12.*t13.*t17;
t321 = Ir.*N.*l_AB.*l_B_m3.*t3.*t13.*t16.*2.0;
t334 = Ir.*t7.*t10.*t13.*t16;
t335 = Ir.*t7.*t11.*t13.*t15;
t349 = m3.*t11.*t12.*t13.*t15;
t372 = Ir.*t7.*t10.*t13.*t18.*2.0;
t373 = m3.*t10.*t12.*t13.*t18.*2.0;
t377 = Ir.*l_AB.*l_B_m3.*t3.*t7.*t13.*t16.*2.0;
t33 = t30.*2.0;
t36 = t31.^2;
t37 = -t30;
t39 = N.*t28.*2.0;
t40 = N.*t29.*2.0;
t58 = I2.*t48;
t65 = I_motor.*t48;
t66 = I_motor.*t49;
t72 = Ir.*t48;
t73 = Ir.*t49;
t79 = N.*t55;
t82 = -t41;
t83 = N.*t35.*4.0;
t84 = t45.*2.0;
t91 = I2.*t42;
t92 = t7.*t29;
t93 = I_motor.*t42;
t94 = I_motor.*t43;
t101 = Ir.*t54.*2.0;
t107 = t7.*t28.*2.0;
t109 = I2.*t85;
t110 = I1.*t6.*t9;
t111 = I2.*t6.*t9;
t112 = t9.*t35;
t113 = I_motor.*t85;
t115 = t95.*2.0;
t126 = Ir.*N.*t54.*4.0;
t127 = m_motor.*t13.*t30.*4.0;
t130 = t7.*t55;
t135 = I2.*l_B_m3.*l_OA.*m3.*t31;
t136 = I_motor.*l_B_m3.*l_OA.*m3.*t31;
t137 = Ir.*l_B_m3.*l_OA.*m3.*t31;
t142 = I2.*t86.*2.0;
t143 = t8.*t35.*4.0;
t144 = t7.*t35.*6.0;
t145 = I_motor.*t86.*2.0;
t146 = I_motor.*t6.*t9.*2.0;
t148 = m2.*t11.*t48;
t149 = m2.*t13.*t46;
t150 = m1.*t14.*t46;
t151 = m3.*t13.*t47;
t152 = m1.*t14.*t47;
t153 = m2.*t13.*t49;
t154 = m1.*t14.*t49;
t155 = m_motor.*t10.*t48;
t156 = m2.*t11.*t54;
t157 = m_motor.*t10.*t50;
t158 = m_motor.*t10.*t52;
t159 = m2.*t13.*t54;
t160 = m1.*t14.*t54;
t161 = m2.*t11.*t55;
t162 = m2.*t13.*t55;
t163 = m1.*t14.*t55;
t164 = m_motor.*t10.*t54;
t165 = m_motor.*t10.*t55;
t171 = m3.*t10.*t30.*-2.0;
t172 = m2.*t11.*t30.*-2.0;
t175 = m2.*t13.*t30.*-2.0;
t176 = m3.*t13.*t30.*-2.0;
t177 = m1.*t14.*t30.*-2.0;
t178 = m_motor.*t10.*t30.*-2.0;
t181 = N.*t105.*4.0;
t188 = m3.*t10.*t42;
t189 = m3.*t10.*t44;
t190 = m2.*t11.*t42;
t191 = m2.*t11.*t44;
t192 = m3.*t12.*t43;
t193 = m2.*t13.*t43;
t194 = m3.*t13.*t43;
t195 = m1.*t14.*t43;
t196 = m_motor.*t10.*t42;
t197 = m_motor.*t10.*t44;
t198 = m3.*t10.*t45;
t199 = m2.*t11.*t45;
t200 = m2.*t13.*t45;
t201 = m3.*t13.*t45;
t202 = m1.*t14.*t45;
t203 = m_motor.*t10.*t45;
t208 = I1.*t185;
t211 = I2.*t186;
t212 = I_motor.*t185;
t213 = I_motor.*t186;
t214 = Ir.*t185;
t215 = Ir.*t186;
t216 = m_motor.*t13.*t46.*3.0;
t217 = m_motor.*t13.*t47.*2.0;
t218 = m_motor.*t13.*t49.*2.0;
t219 = m_motor.*t13.*t54.*2.0;
t220 = m_motor.*t13.*t55.*2.0;
t221 = I1.*t95.*-2.0;
t222 = l_AB.*l_OA.*m3.*t2.*t30.*-2.0;
t223 = l_A_m2.*l_OA.*m2.*t2.*t30.*-2.0;
t224 = I_motor.*t95.*-2.0;
t225 = l_AB.*l_OA.*m_motor.*t2.*t30.*-2.0;
t237 = m3.*t12.*t42.*2.0;
t238 = m2.*t13.*t44.*2.0;
t239 = m3.*t13.*t44.*2.0;
t240 = m1.*t14.*t44.*2.0;
t241 = m3.*t12.*t45.*3.0;
t242 = m_motor.*t13.*t43.*2.0;
t243 = m_motor.*t13.*t44.*4.0;
t246 = m3.*t10.*t85;
t247 = m3.*t6.*t9.*t10;
t248 = m2.*t11.*t85;
t249 = m2.*t6.*t9.*t11;
t250 = t9.*t105;
t251 = m2.*t6.*t9.*t13;
t252 = m3.*t6.*t9.*t13;
t253 = m1.*t6.*t9.*t14;
t254 = m_motor.*t10.*t85;
t255 = m_motor.*t6.*t9.*t10;
t263 = -t227;
t264 = l_B_m3.*l_OA.*m2.*m3.*t11.*t31;
t265 = l_B_m3.*l_OA.*m3.*m_motor.*t10.*t31;
t267 = l_AB.*l_B_m3.*m3.*t3.*t42.*2.0;
t268 = l_AB.*l_OA.*m3.*t2.*t44.*2.0;
t269 = l_A_m2.*l_OA.*m2.*t2.*t44.*2.0;
t271 = l_AB.*l_OA.*m_motor.*t2.*t44.*2.0;
t273 = m3.*t10.*t86.*2.0;
t274 = m2.*t11.*t86.*2.0;
t275 = m3.*t12.*t86.*4.0;
t276 = m3.*t12.*t85.*6.0;
t277 = m_motor.*t10.*t86.*2.0;
t278 = m_motor.*t6.*t9.*t13.*2.0;
t279 = m1.*t14.*t138;
t280 = m_motor.*t10.*t139;
t281 = m_motor.*t10.*t140;
t287 = m_motor.*t13.*t95.*4.0;
t295 = l_B_m3.*l_OA.*t10.*t16.*t31;
t299 = m2.*t13.*t128;
t300 = m1.*t14.*t128;
t301 = m3.*t13.*t129;
t302 = m1.*t14.*t129;
t304 = m_motor.*t10.*t131;
t305 = m_motor.*t10.*t133;
t306 = l_AB.*l_B_m3.*m3.*t3.*t6.*t9.*2.0;
t307 = l_AB.*l_B_m3.*m3.*t3.*t85.*6.0;
t308 = l_AB.*l_B_m3.*m3.*t3.*t86.*6.0;
t309 = l_AB.*l_OA.*m3.*t2.*t86.*2.0;
t310 = l_AB.*l_OA.*m3.*t2.*t6.*t9.*2.0;
t311 = l_A_m2.*l_OA.*m2.*t2.*t86.*2.0;
t312 = l_A_m2.*l_OA.*m2.*t2.*t6.*t9.*2.0;
t313 = l_AB.*l_OA.*m_motor.*t2.*t86.*2.0;
t314 = l_AB.*l_OA.*m_motor.*t2.*t6.*t9.*2.0;
t319 = m_motor.*t13.*t138.*2.0;
t320 = Ir.*N.*t236.*2.0;
t322 = m2.*t13.*t95.*-2.0;
t323 = m1.*t14.*t95.*-2.0;
t329 = l_AB.*l_A_m2.*m3.*t23.*t50.*2.0;
t330 = l_AB.*l_A_m2.*m_motor.*t23.*t50.*2.0;
t338 = m_motor.*t13.*t128.*3.0;
t339 = m_motor.*t13.*t129.*2.0;
t345 = t24.*t185;
t346 = m2.*t13.*t185;
t347 = m1.*t14.*t185;
t348 = m2.*t11.*t186;
t351 = -t321;
t353 = t23.*t209;
t354 = t23.*t210;
t356 = I3.*t10.*t13.*t18.*t23;
t358 = m_motor.*t13.*t23.*t46.*2.0;
t359 = l_AB.*l_OA.*t3.*t12.*t16.*t31;
t363 = l_B_m3.*l_OA.*m3.*t31.*t85.*2.0;
t364 = l_B_m3.*l_OA.*m3.*t6.*t9.*t31.*2.0;
t365 = l_B_m3.*l_OA.*m3.*t31.*t86.*4.0;
t366 = l_AB.*l_B_m3.*m3.*t3.*t131.*2.0;
t367 = l_AB.*l_B_m3.*m3.*t3.*t133.*2.0;
t370 = l_AB.*l_B_m3.*m_motor.*t3.*t132.*4.0;
t374 = m_motor.*t13.*t185.*3.0;
t376 = Ir.*t7.*t236.*2.0;
t381 = l_AB.*l_A_m2.*m3.*t23.*t131.*2.0;
t382 = l_AB.*l_A_m2.*m_motor.*t23.*t131.*2.0;
t392 = l_AB.*l_A_m2.*m_motor.*t23.*t139.*2.0;
t395 = t23.*t297;
t396 = t24.*t297;
t397 = Ir.*N.*l_A_m2.*l_B_m3.*m2.*m3.*t2.*t13.*t31.*2.0;
t398 = Ir.*N.*l_AB.*l_B_m3.*m3.*m_motor.*t2.*t13.*t31.*2.0;
t403 = t23.*t349;
t404 = m3.*t10.*t12.*t13.*t18.*t23;
t405 = l_AB.*l_A_m2.*m2.*t23.*t186.*2.0;
t408 = Ir.*N.*l_AB.*l_B_m3.*t2.*t13.*t16.*t31.*2.0;
t411 = t23.*t334;
t412 = t23.*t335;
t413 = Ir.*t7.*t10.*t13.*t18.*t23;
t414 = m_motor.*t13.*t23.*t128.*2.0;
t418 = m_motor.*t13.*t23.*t185.*2.0;
t421 = l_A_m2.*l_B_m3.*m3.*t2.*t31.*t131.*2.0;
t422 = l_AB.*l_B_m3.*m_motor.*t2.*t31.*t132.*2.0;
t436 = Ir.*l_AB.*l_B_m3.*t2.*t7.*t13.*t16.*t31.*2.0;
t450 = l_AB.*l_A_m2.*m2.*t2.*t3.*t31.*t186.*2.0;
t452 = t2.*t3.*t31.*t297.*2.0;
t453 = m_motor.*t2.*t3.*t13.*t31.*t185.*2.0;
t38 = -t33;
t80 = -t39;
t81 = -t40;
t103 = t79.*2.0;
t106 = -t83;
t116 = -t79;
t118 = m3.*t10.*t33;
t119 = m2.*t11.*t33;
t120 = N.*t72.*2.0;
t121 = N.*t73.*2.0;
t122 = m2.*t13.*t33;
t123 = m3.*t13.*t33;
t124 = m1.*t14.*t33;
t125 = m_motor.*t10.*t33;
t147 = -t115;
t166 = I1.*t115;
t167 = l_AB.*l_OA.*m3.*t2.*t33;
t168 = l_A_m2.*l_OA.*m2.*t2.*t33;
t169 = I_motor.*t115;
t170 = l_AB.*l_OA.*m_motor.*t2.*t33;
t179 = -t126;
t180 = -t127;
t204 = N.*t137;
t205 = -t142;
t206 = -t143;
t207 = -t145;
t231 = -t181;
t232 = -t135;
t233 = -t136;
t234 = -t137;
t244 = m_motor.*t13.*t84;
t257 = t7.*t137;
t262 = m_motor.*t13.*t79.*4.0;
t270 = l_AB.*l_B_m3.*m3.*t3.*t84;
t283 = m2.*t13.*t115;
t284 = m1.*t14.*t115;
t288 = N.*t214.*2.0;
t289 = N.*t215.*2.0;
t290 = m2.*t11.*t79.*-2.0;
t291 = m2.*t13.*t79.*-2.0;
t292 = m1.*t14.*t79.*-2.0;
t293 = m_motor.*t10.*t79.*-2.0;
t298 = m2.*t11.*t130;
t303 = m_motor.*t10.*t130;
t315 = -t273;
t316 = -t274;
t317 = -t275;
t318 = -t277;
t324 = l_A_m2.*l_OA.*m2.*t2.*t79.*-2.0;
t325 = l_AB.*l_OA.*m_motor.*t2.*t79.*-2.0;
t326 = -t287;
t331 = -t264;
t332 = -t265;
t333 = t7.*t214;
t336 = m2.*t13.*t130.*2.0;
t337 = m1.*t14.*t130.*2.0;
t340 = m_motor.*t13.*t130.*4.0;
t341 = -t308;
t342 = -t309;
t343 = -t311;
t344 = -t313;
t350 = -t320;
t352 = t24.*t208;
t355 = t24.*t212;
t357 = t24.*t214;
t360 = -t329;
t361 = -t330;
t362 = -t295;
t368 = l_A_m2.*l_OA.*m2.*t2.*t130.*2.0;
t369 = l_AB.*l_OA.*m_motor.*t2.*t130.*2.0;
t371 = t7.*t215.*2.0;
t375 = -t358;
t378 = -t365;
t379 = t36.*t186;
t380 = -t345;
t383 = t36.*t211;
t384 = t36.*t213;
t385 = t36.*t215;
t387 = -t353;
t388 = -t354;
t390 = -t356;
t393 = -t359;
t399 = -t381;
t400 = -t382;
t401 = m2.*t13.*t345;
t402 = m1.*t14.*t345;
t406 = -t392;
t409 = Ir.*N.*t359.*2.0;
t419 = m_motor.*t13.*t345.*2.0;
t420 = -t405;
t424 = t36.*t297;
t425 = -t395;
t426 = -t396;
t427 = -t414;
t428 = t36.*t348;
t429 = m_motor.*t13.*t36.*t185;
t432 = -t403;
t433 = -t404;
t434 = -t418;
t437 = Ir.*t7.*t359.*2.0;
t438 = -t421;
t439 = -t422;
t441 = -t411;
t442 = -t412;
t443 = -t413;
t444 = -t436;
t117 = -t103;
t173 = -t120;
t174 = -t121;
t245 = t204.*2.0;
t258 = m2.*t11.*t103;
t259 = m2.*t13.*t103;
t260 = m1.*t14.*t103;
t261 = m_motor.*t10.*t103;
t282 = t257.*2.0;
t285 = l_A_m2.*l_OA.*m2.*t2.*t103;
t286 = l_AB.*l_OA.*m_motor.*t2.*t103;
t294 = -t262;
t327 = -t288;
t328 = -t289;
t386 = -t352;
t389 = -t355;
t391 = -t357;
t394 = t24.*t288;
t407 = -t379;
t410 = t24.*t333;
t415 = -t383;
t416 = -t384;
t417 = -t385;
t423 = t36.*t289;
t430 = m2.*t13.*t380;
t431 = m1.*t14.*t380;
t435 = -t419;
t445 = -t437;
t446 = -t424;
t447 = t36.*t371;
t448 = -t428;
t449 = -t429;
t451 = t7.*t385.*-2.0;
t272 = -t245;
t440 = -t410;
t454 = t25+t26+t27+t28+t29+t32+t34+t35+t56+t57+t58+t59+t60+t61+t62+t63+t64+t65+t66+t67+t68+t69+t70+t71+t72+t73+t74+t75+t76+t77+t78+t80+t81+t82+t91+t92+t93+t94+t99+t100+t101+t102+t104+t105+t106+t107+t108+t109+t110+t111+t112+t113+t114+t144+t146+t148+t149+t150+t151+t152+t153+t154+t155+t156+t157+t158+t159+t160+t161+t162+t163+t164+t165+t171+t172+t173+t174+t175+t176+t177+t178+t179+t180+t188+t189+t190+t191+t192+t193+t194+t195+t196+t197+t198+t199+t200+t201+t202+t203+t205+t206+t207+t208+t209+t210+t211+t212+t213+t214+t215+t216+t217+t218+t219+t220+t221+t222+t223+t224+t225+t231+t237+t238+t239+t240+t241+t242+t243+t244+t246+t247+t248+t249+t250+t251+t252+t253+t254+t255+t256+t263+t267+t268+t269+t270+t271+t276+t278+t279+t280+t281+t290+t291+t292+t293+t294+t297+t298+t299+t300+t301+t302+t303+t304+t305+t306+t307+t310+t312+t314+t315+t316+t317+t318+t319+t322+t323+t324+t325+t326+t327+t328+t333+t334+t335+t336+t337+t338+t339+t340+t341+t342+t343+t344+t346+t347+t348+t349+t350+t351+t360+t361+t363+t364+t366+t367+t368+t369+t370+t371+t372+t373+t374+t375+t376+t377+t378+t386+t387+t388+t389+t390+t391+t394+t397+t398+t399+t400+t406+t408+t409+t415+t416+t417+t420+t423+t425+t426+t427+t430+t431+t432+t433+t434+t435+t438+t439+t440+t441+t442+t443+t444+t445+t446+t448+t449+t450+t451+t452+t453;
t455 = 1.0./t454;
Mass_Rot_Op_Sp_inv = t455.*(t20+t21+t22+t38+t43+t44.*2.0+t45+t46+t47+t49+t53+t54+t55+t117+t128+t129+t130.*2.0+t134+t138+t141+t147+t185+t226+t380+t6.*t9)+t455.*(t19+t20+t21.*2.0+t22.*2.0+t38+t42+t43+t44+t46+t47+t48+t49+t50+t51+t52+t53+t54.*2.0+t55.*2.0+t84+t85+t87.*2.0+t88.*2.0+t89.*2.0+t90+t117+t128+t129+t130+t131+t132+t133+t134+t138+t139+t140+t141+t147+t182+t183.*2.0+t184.*2.0+t185+t186+t187+t226+t228+t229+t230+t236.*2.0+t272+t282-t359.*2.0+t380+t407)-t455.*(t20+t21+t22+t37+t43+t44+t45+t46+t47+t49+t53+t54+t55+t86+t87+t88+t89+t116+t128+t129+t130+t134+t138+t141+t147+t183+t184+t185-t204+t226+t236+t257+t380+t393+Ir.*l_AB.*l_OA.*m3.*t2.*t7+Ir.*l_A_m2.*l_OA.*m2.*t2.*t7+Ir.*l_AB.*l_OA.*m_motor.*t2.*t7).*2.0+t455.*(t30-t44+t79+t85-t86+t87+t88+t89+t95+t96+t97+t98-t130+t183+t184+t204+t232+t233+t234+t235+t236+t266+t296+t331+t332+t362+t393+t7.*t234-Ir.*l_AB.*l_B_m3.*m3.*t3.*t7).*2.0-t455.*(t19+t21+t22+t37+t48+t50+t51+t52+t54+t55-t85+t87+t88+t89+t90-t95+t96+t97+t98+t116+t139+t140+t183+t184+t186+t187+t232+t233+t234+t235+t236+t266+t296+t331+t332+t362+t393+t407+N.*t6+I1.*Ir.*N+I_motor.*Ir.*N+Ir.*N.*m1.*t14+Ir.*N.*m2.*t13+Ir.*N.*m3.*t13+Ir.*N.*m_motor.*t13.*2.0+I1.*l_AB.*l_B_m3.*m3.*t3+I_motor.*l_AB.*l_B_m3.*m3.*t3+Ir.*l_AB.*l_B_m3.*m3.*t3+l_AB.*l_B_m3.*t3.*t13.*t16-l_AB.*l_B_m3.*t2.*t13.*t16.*t31+l_AB.*l_B_m3.*m1.*m3.*t3.*t14+l_AB.*l_B_m3.*m2.*m3.*t3.*t13+l_AB.*l_B_m3.*m3.*m_motor.*t3.*t13.*2.0-l_A_m2.*l_B_m3.*m2.*m3.*t2.*t13.*t31-l_AB.*l_B_m3.*m3.*m_motor.*t2.*t13.*t31).*2.0+t455.*(t5+t6+t19+t21+t22+t38+t42+t43+t44+t48+t50+t51+t52+t54+t55+t84+t85+t90-t95.*4.0-t96.*2.0-t97.*2.0-t98.*2.0+t117+t128+t129+t130+t131+t132+t133+t134+t139+t140+t182+t186+t187+t226+t228+t229+t230+t272+t282+t407-N.*t6.*2.0+I1.*I2+I1.*I_motor+I2.*I_motor+I1.*Ir+I2.*Ir+I_motor.*Ir.*2.0-I2.*Ir.*N.*2.0-I_motor.*Ir.*N.*2.0+I1.*m2.*t11+I1.*m3.*t10+I2.*m1.*t14+I2.*m2.*t13+I2.*m3.*t13+I1.*m_motor.*t10+I_motor.*m2.*t11+I_motor.*m3.*t10+I2.*m_motor.*t13.*2.0+I_motor.*m1.*t14+I_motor.*m2.*t13+I_motor.*m3.*t13+Ir.*m2.*t11+Ir.*m3.*t10+Ir.*m1.*t14+Ir.*m2.*t13+Ir.*m3.*t13+I_motor.*m_motor.*t10+I_motor.*m_motor.*t13.*2.0+Ir.*m_motor.*t10+Ir.*m_motor.*t13.*2.0+t10.*t13.*t16+t11.*t13.*t15+t10.*t13.*t18.*2.0+m1.*m2.*t11.*t14+m1.*m3.*t10.*t14+m2.*m3.*t10.*t13+m2.*m3.*t11.*t13+m1.*m_motor.*t10.*t14+m2.*m_motor.*t10.*t13+m2.*m_motor.*t11.*t13.*2.0+m3.*m_motor.*t10.*t13.*3.0-t10.*t13.*t16.*t23-t11.*t13.*t15.*t23-t10.*t13.*t18.*t23-Ir.*N.*m2.*t11.*2.0-Ir.*N.*m3.*t10.*2.0-Ir.*N.*m_motor.*t10.*2.0+I1.*l_AB.*l_B_m3.*m3.*t3.*2.0+I_motor.*l_AB.*l_B_m3.*m3.*t3.*2.0+Ir.*l_AB.*l_B_m3.*m3.*t3.*2.0+l_AB.*l_B_m3.*t3.*t13.*t16.*2.0-m3.*m_motor.*t10.*t13.*t23.*2.0-l_AB.*l_B_m3.*t2.*t13.*t16.*t31.*2.0-l_AB.*l_A_m2.*m2.*m3.*t13.*t23.*2.0+l_AB.*l_B_m3.*m1.*m3.*t3.*t14.*2.0+l_AB.*l_B_m3.*m2.*m3.*t3.*t13.*2.0-l_AB.*l_A_m2.*m2.*m_motor.*t13.*t23.*2.0+l_AB.*l_B_m3.*m3.*m_motor.*t3.*t13.*4.0-l_A_m2.*l_B_m3.*m2.*m3.*t2.*t13.*t31.*2.0-l_AB.*l_B_m3.*m3.*m_motor.*t2.*t13.*t31.*2.0);
