function Mass_Joint_Sp_inv = A_inv_arm_rotors(in1,in2)
%A_INV_ARM_ROTORS
%    MASS_JOINT_SP_INV = A_INV_ARM_ROTORS(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    30-May-2021 10:56:50

I1 = in2(4,:);
I2 = in2(5,:);
I3 = in2(6,:);
Ir1 = in2(7,:);
Ir2 = in2(8,:);
Ir3 = in2(9,:);
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
t4 = th1+th2;
t5 = l_AB.^2;
t6 = l_A_m2.^2;
t7 = l_B_m3.^2;
t8 = l_OA.^2;
t9 = l_O_m1.^2;
t10 = m2.^2;
t11 = m3.^2;
t12 = I1.*I3;
t13 = I2.*I3;
t14 = I1.*Ir3;
t15 = I3.*Ir1;
t16 = I2.*Ir3;
t17 = I3.*Ir2;
t18 = Ir1.*Ir3;
t19 = Ir2.*Ir3;
t20 = t2.^2;
t21 = I2.*t12;
t22 = I2.*t14;
t23 = Ir2.*t12;
t24 = Ir1.*t13;
t25 = Ir2.*t14;
t26 = Ir1.*t16;
t27 = Ir2.*t15;
t28 = Ir2.*t18;
t29 = t3.^2;
t30 = cos(t4);
t31 = sin(t4);
t32 = t4+th3;
t33 = cos(t32);
t34 = sin(t32);
t35 = t30.^2;
t36 = t31.^2;
t39 = I3.*m2.*t8.*t20;
t40 = I3.*m3.*t8.*t20;
t41 = I3.*m1.*t9.*t20;
t42 = Ir3.*m2.*t8.*t20;
t43 = Ir3.*m3.*t8.*t20;
t44 = Ir3.*m1.*t9.*t20;
t45 = I3.*m2.*t8.*t29;
t46 = I3.*m3.*t8.*t29;
t47 = I3.*m1.*t9.*t29;
t48 = Ir3.*m2.*t8.*t29;
t49 = Ir3.*m3.*t8.*t29;
t50 = Ir3.*m1.*t9.*t29;
t51 = m2.*t8.*t13.*t20;
t52 = m3.*t8.*t13.*t20;
t53 = m1.*t9.*t13.*t20;
t54 = m2.*t8.*t16.*t20;
t55 = m2.*t8.*t17.*t20;
t56 = m3.*t8.*t16.*t20;
t57 = m3.*t8.*t17.*t20;
t58 = m1.*t9.*t16.*t20;
t59 = m1.*t9.*t17.*t20;
t60 = m2.*t8.*t19.*t20;
t61 = m3.*t8.*t19.*t20;
t62 = m1.*t9.*t19.*t20;
t63 = m2.*t8.*t13.*t29;
t64 = m3.*t8.*t13.*t29;
t65 = m1.*t9.*t13.*t29;
t66 = I3.*l_AB.*l_OA.*m3.*t2.*t30;
t67 = I3.*l_A_m2.*l_OA.*m2.*t2.*t30;
t68 = m2.*t8.*t16.*t29;
t69 = m2.*t8.*t17.*t29;
t70 = m3.*t8.*t16.*t29;
t71 = m3.*t8.*t17.*t29;
t72 = m1.*t9.*t16.*t29;
t73 = m1.*t9.*t17.*t29;
t74 = Ir3.*l_AB.*l_OA.*m3.*t2.*t30;
t75 = Ir3.*l_A_m2.*l_OA.*m2.*t2.*t30;
t76 = m2.*t8.*t19.*t29;
t77 = m3.*t8.*t19.*t29;
t78 = m1.*t9.*t19.*t29;
t79 = I3.*l_AB.*l_OA.*m3.*t3.*t31;
t80 = I3.*l_A_m2.*l_OA.*m2.*t3.*t31;
t81 = Ir3.*l_AB.*l_OA.*m3.*t3.*t31;
t82 = Ir3.*l_A_m2.*l_OA.*m2.*t3.*t31;
t231 = I3.*l_AB.*l_A_m2.*m2.*m3.*t2.*t3.*t8.*t30.*t31.*4.0;
t232 = Ir3.*l_AB.*l_A_m2.*m2.*m3.*t2.*t3.*t8.*t30.*t31.*4.0;
t271 = I3.*t2.*t3.*t5.*t8.*t11.*t30.*t31.*2.0;
t272 = I3.*t2.*t3.*t6.*t8.*t10.*t30.*t31.*2.0;
t273 = Ir3.*t2.*t3.*t5.*t8.*t11.*t30.*t31.*2.0;
t274 = Ir3.*t2.*t3.*t6.*t8.*t10.*t30.*t31.*2.0;
t37 = t33.^2;
t38 = t34.^2;
t83 = I2.*l_B_m3.*l_OA.*m3.*t3.*t34;
t84 = Ir2.*l_B_m3.*l_OA.*m3.*t3.*t34;
t85 = I2.*l_B_m3.*l_OA.*m3.*t2.*t33;
t86 = Ir2.*l_B_m3.*l_OA.*m3.*t2.*t33;
t87 = I3.*m3.*t5.*t35;
t88 = I3.*m2.*t6.*t35;
t89 = Ir3.*m3.*t5.*t35;
t90 = Ir3.*m2.*t6.*t35;
t91 = I3.*m3.*t5.*t36;
t92 = I3.*m2.*t6.*t36;
t93 = Ir3.*m3.*t5.*t36;
t94 = Ir3.*m2.*t6.*t36;
t95 = m3.*t5.*t12.*t35;
t96 = m2.*t6.*t12.*t35;
t97 = m3.*t5.*t14.*t35;
t98 = m3.*t5.*t15.*t35;
t99 = m2.*t6.*t14.*t35;
t100 = m2.*t6.*t15.*t35;
t101 = m3.*t5.*t18.*t35;
t102 = m2.*t6.*t18.*t35;
t103 = m3.*t5.*t12.*t36;
t104 = m2.*t6.*t12.*t36;
t105 = m3.*t5.*t14.*t36;
t106 = m3.*t5.*t15.*t36;
t107 = m2.*t6.*t14.*t36;
t108 = m2.*t6.*t15.*t36;
t109 = m3.*t5.*t18.*t36;
t110 = m2.*t6.*t18.*t36;
t126 = I1.*l_AB.*l_B_m3.*m3.*t30.*t33;
t129 = Ir1.*l_AB.*l_B_m3.*m3.*t30.*t33;
t131 = I1.*l_AB.*l_B_m3.*m3.*t31.*t34;
t132 = Ir1.*l_AB.*l_B_m3.*m3.*t31.*t34;
t135 = l_AB.*l_A_m2.*l_B_m3.*l_OA.*m2.*m3.*t3.*t34.*t36;
t136 = m3.*t5.*t35.*t39;
t137 = m3.*t5.*t35.*t41;
t138 = m3.*t6.*t35.*t39;
t139 = m2.*t6.*t35.*t41;
t140 = m3.*t5.*t35.*t42;
t141 = m3.*t5.*t35.*t44;
t142 = m3.*t6.*t35.*t42;
t143 = m2.*t6.*t35.*t44;
t144 = l_AB.*l_A_m2.*m3.*t35.*t39.*2.0;
t145 = l_AB.*l_A_m2.*m3.*t35.*t42.*2.0;
t146 = m3.*t5.*t35.*t45;
t147 = m3.*t5.*t36.*t39;
t148 = m3.*t5.*t35.*t47;
t149 = m3.*t5.*t36.*t41;
t150 = m3.*t6.*t35.*t45;
t151 = m3.*t6.*t36.*t39;
t152 = m2.*t6.*t35.*t47;
t153 = m2.*t6.*t36.*t41;
t154 = m3.*t5.*t35.*t48;
t155 = m3.*t5.*t36.*t42;
t156 = m3.*t5.*t35.*t50;
t157 = m3.*t5.*t36.*t44;
t158 = m3.*t6.*t35.*t48;
t159 = m3.*t6.*t36.*t42;
t160 = m2.*t6.*t35.*t50;
t161 = m2.*t6.*t36.*t44;
t162 = m3.*t5.*t36.*t45;
t163 = m3.*t5.*t36.*t47;
t164 = m3.*t6.*t36.*t45;
t165 = m2.*t6.*t36.*t47;
t166 = m3.*t5.*t36.*t48;
t167 = m3.*t5.*t36.*t50;
t168 = m3.*t6.*t36.*t48;
t169 = m2.*t6.*t36.*t50;
t172 = l_AB.*l_A_m2.*m3.*t36.*t45.*2.0;
t173 = l_AB.*l_A_m2.*m3.*t36.*t48.*2.0;
t180 = l_AB.*l_A_m2.*l_B_m3.*l_OA.*m2.*m3.*t2.*t33.*t35;
t181 = l_AB.*l_B_m3.*m2.*m3.*t8.*t20.*t31.*t34;
t182 = l_AB.*l_B_m3.*m1.*m3.*t9.*t20.*t31.*t34;
t184 = l_B_m3.*l_OA.*m2.*m3.*t3.*t6.*t34.*t35;
t185 = l_AB.*l_B_m3.*m2.*m3.*t8.*t29.*t31.*t34;
t186 = l_AB.*l_B_m3.*m1.*m3.*t9.*t29.*t31.*t34;
t187 = l_A_m2.*l_B_m3.*m2.*m3.*t8.*t29.*t31.*t34;
t189 = l_B_m3.*l_OA.*m2.*m3.*t3.*t6.*t34.*t36;
t190 = I3.*t5.*t8.*t11.*t29.*t35;
t191 = I3.*t5.*t8.*t11.*t20.*t36;
t192 = I3.*t6.*t8.*t10.*t29.*t35;
t193 = I3.*t6.*t8.*t10.*t20.*t36;
t194 = Ir3.*t5.*t8.*t11.*t29.*t35;
t195 = Ir3.*t5.*t8.*t11.*t20.*t36;
t196 = Ir3.*t6.*t8.*t10.*t29.*t35;
t197 = Ir3.*t6.*t8.*t10.*t20.*t36;
t200 = l_A_m2.*l_B_m3.*m2.*m3.*t2.*t3.*t8.*t31.*t33;
t201 = l_A_m2.*l_B_m3.*m2.*m3.*t2.*t3.*t8.*t30.*t34;
t214 = l_AB.*l_B_m3.*m2.*m3.*t8.*t20.*t30.*t33;
t215 = l_AB.*l_B_m3.*m1.*m3.*t9.*t20.*t30.*t33;
t216 = l_A_m2.*l_B_m3.*m2.*m3.*t8.*t20.*t30.*t33;
t218 = l_B_m3.*l_OA.*m2.*m3.*t2.*t6.*t33.*t35;
t225 = l_AB.*l_B_m3.*m2.*m3.*t8.*t29.*t30.*t33;
t226 = l_AB.*l_B_m3.*m1.*m3.*t9.*t29.*t30.*t33;
t228 = l_B_m3.*l_OA.*m2.*m3.*t2.*t6.*t33.*t36;
t233 = l_AB.*l_A_m2.*l_B_m3.*l_OA.*m2.*m3.*t3.*t30.*t31.*t33;
t234 = l_AB.*l_A_m2.*l_B_m3.*l_OA.*m2.*m3.*t2.*t30.*t31.*t34;
t238 = l_AB.*l_B_m3.*t2.*t3.*t8.*t11.*t31.*t33;
t239 = l_AB.*l_B_m3.*t2.*t3.*t8.*t11.*t30.*t34;
t245 = l_AB.*l_B_m3.*t8.*t11.*t29.*t30.*t33;
t247 = l_B_m3.*l_OA.*t2.*t5.*t11.*t33.*t36;
t248 = l_AB.*l_B_m3.*t8.*t11.*t20.*t31.*t34;
t250 = l_B_m3.*l_OA.*t3.*t5.*t11.*t34.*t35;
t254 = -t231;
t255 = -t232;
t275 = l_B_m3.*l_OA.*t3.*t5.*t11.*t30.*t31.*t33;
t276 = l_B_m3.*l_OA.*t2.*t5.*t11.*t30.*t31.*t34;
t285 = -t271;
t286 = -t272;
t287 = -t273;
t288 = -t274;
t289 = l_AB.*l_OA.*t3.*t7.*t11.*t30.*t33.*t34;
t290 = l_AB.*l_OA.*t2.*t7.*t11.*t31.*t33.*t34;
t291 = t2.*t3.*t7.*t8.*t11.*t33.*t34.*2.0;
t295 = I2.*t2.*t3.*t7.*t8.*t11.*t33.*t34.*-2.0;
t296 = Ir2.*t2.*t3.*t7.*t8.*t11.*t33.*t34.*-2.0;
t299 = t5.*t7.*t11.*t30.*t31.*t33.*t34.*2.0;
t311 = I1.*t5.*t7.*t11.*t30.*t31.*t33.*t34.*-2.0;
t312 = Ir1.*t5.*t7.*t11.*t30.*t31.*t33.*t34.*-2.0;
t345 = l_AB.*l_A_m2.*m2.*t7.*t8.*t11.*t29.*t30.*t31.*t33.*t34.*2.0;
t346 = m2.*t2.*t3.*t6.*t7.*t8.*t11.*t33.*t34.*t35.*-2.0;
t347 = m2.*t2.*t3.*t6.*t7.*t8.*t11.*t33.*t34.*t36.*-2.0;
t348 = l_AB.*l_A_m2.*m2.*t7.*t8.*t11.*t20.*t30.*t31.*t33.*t34.*2.0;
t353 = m2.*t5.*t7.*t8.*t11.*t20.*t30.*t31.*t33.*t34.*-2.0;
t354 = m1.*t5.*t7.*t9.*t11.*t20.*t30.*t31.*t33.*t34.*-2.0;
t355 = m2.*t5.*t7.*t8.*t11.*t29.*t30.*t31.*t33.*t34.*-2.0;
t356 = m1.*t5.*t7.*t9.*t11.*t29.*t30.*t31.*t33.*t34.*-2.0;
t111 = I1.*m3.*t7.*t37;
t112 = I2.*m3.*t7.*t37;
t113 = Ir1.*m3.*t7.*t37;
t114 = Ir2.*m3.*t7.*t37;
t115 = I1.*m3.*t7.*t38;
t116 = I2.*m3.*t7.*t38;
t117 = Ir1.*m3.*t7.*t38;
t118 = Ir2.*m3.*t7.*t38;
t119 = -t83;
t120 = -t84;
t133 = -t85;
t134 = -t86;
t170 = m2.*m3.*t7.*t8.*t20.*t37;
t171 = m1.*m3.*t7.*t9.*t20.*t37;
t174 = m2.*m3.*t7.*t8.*t29.*t37;
t175 = m2.*m3.*t7.*t8.*t20.*t38;
t176 = m1.*m3.*t7.*t9.*t29.*t37;
t177 = m1.*m3.*t7.*t9.*t20.*t38;
t178 = m2.*m3.*t7.*t8.*t29.*t38;
t179 = m1.*m3.*t7.*t9.*t29.*t38;
t183 = l_A_m2.*l_OA.*m2.*m3.*t3.*t7.*t31.*t37;
t188 = l_A_m2.*l_OA.*m2.*m3.*t3.*t7.*t31.*t38;
t198 = -t144;
t199 = -t145;
t202 = -t172;
t203 = -t173;
t204 = t7.*t8.*t11.*t29.*t37;
t205 = t7.*t8.*t11.*t20.*t38;
t217 = l_A_m2.*l_OA.*m2.*m3.*t2.*t7.*t30.*t37;
t227 = l_A_m2.*l_OA.*m2.*m3.*t2.*t7.*t30.*t38;
t235 = m2.*m3.*t6.*t7.*t35.*t37;
t236 = m2.*m3.*t6.*t7.*t36.*t37;
t237 = m2.*m3.*t6.*t7.*t35.*t38;
t240 = m2.*m3.*t6.*t7.*t36.*t38;
t246 = l_AB.*l_OA.*t2.*t7.*t11.*t30.*t38;
t249 = l_AB.*l_OA.*t3.*t7.*t11.*t31.*t37;
t251 = -t184;
t252 = -t187;
t253 = -t189;
t256 = t5.*t7.*t11.*t36.*t37;
t257 = t5.*t7.*t11.*t35.*t38;
t258 = -t200;
t259 = -t201;
t268 = -t216;
t269 = -t218;
t270 = -t228;
t277 = -t238;
t278 = -t239;
t283 = -t247;
t284 = -t250;
t292 = I2.*t291;
t293 = Ir2.*t291;
t294 = -t291;
t297 = -t289;
t298 = -t290;
t300 = -t299;
t301 = I1.*t299;
t302 = Ir1.*t299;
t319 = m3.*t6.*t7.*t8.*t10.*t29.*t35.*t37;
t320 = m3.*t6.*t7.*t8.*t10.*t20.*t36.*t37;
t328 = m3.*t6.*t7.*t8.*t10.*t29.*t35.*t38;
t329 = m3.*t6.*t7.*t8.*t10.*t20.*t36.*t38;
t333 = l_AB.*l_A_m2.*m2.*t2.*t3.*t7.*t8.*t11.*t30.*t31.*t37.*2.0;
t334 = l_AB.*l_A_m2.*m2.*t2.*t3.*t7.*t8.*t11.*t30.*t31.*t38.*2.0;
t335 = l_AB.*l_A_m2.*m2.*t35.*t291;
t336 = l_AB.*l_A_m2.*m2.*t36.*t291;
t337 = m3.*t2.*t3.*t6.*t7.*t8.*t10.*t30.*t31.*t37.*2.0;
t339 = m3.*t2.*t3.*t6.*t7.*t8.*t10.*t30.*t31.*t38.*2.0;
t341 = m2.*t6.*t35.*t291;
t342 = m2.*t6.*t36.*t291;
t349 = m2.*t8.*t20.*t299;
t350 = m1.*t9.*t20.*t299;
t351 = m2.*t8.*t29.*t299;
t352 = m1.*t9.*t29.*t299;
t121 = I2.*t111;
t122 = Ir2.*t111;
t123 = Ir1.*t112;
t124 = Ir2.*t113;
t125 = I2.*t115;
t127 = Ir2.*t115;
t128 = Ir1.*t116;
t130 = Ir2.*t117;
t206 = m2.*t8.*t20.*t112;
t207 = m1.*t9.*t20.*t112;
t208 = m2.*t8.*t20.*t114;
t209 = m1.*t9.*t20.*t114;
t210 = m2.*t8.*t29.*t112;
t211 = m2.*t8.*t20.*t116;
t212 = m1.*t9.*t29.*t112;
t213 = m1.*t9.*t20.*t116;
t219 = m2.*t8.*t29.*t114;
t220 = m2.*t8.*t20.*t118;
t221 = m1.*t9.*t29.*t114;
t222 = m1.*t9.*t20.*t118;
t223 = m2.*t8.*t29.*t116;
t224 = m1.*t9.*t29.*t116;
t229 = m2.*t8.*t29.*t118;
t230 = m1.*t9.*t29.*t118;
t241 = I2.*t204;
t242 = I2.*t205;
t243 = Ir2.*t204;
t244 = Ir2.*t205;
t260 = m2.*t6.*t35.*t111;
t261 = m2.*t6.*t35.*t113;
t262 = m2.*t6.*t36.*t111;
t263 = m2.*t6.*t35.*t115;
t264 = m2.*t6.*t36.*t113;
t265 = m2.*t6.*t35.*t117;
t266 = m2.*t6.*t36.*t115;
t267 = m2.*t6.*t36.*t117;
t279 = I1.*t256;
t280 = I1.*t257;
t281 = Ir1.*t256;
t282 = Ir1.*t257;
t303 = m2.*t6.*t35.*t171;
t304 = m2.*t6.*t35.*t176;
t305 = m2.*t6.*t36.*t171;
t306 = m2.*t6.*t35.*t177;
t307 = m2.*t6.*t36.*t176;
t308 = m2.*t6.*t35.*t179;
t309 = m2.*t6.*t36.*t177;
t310 = m2.*t6.*t36.*t179;
t313 = m2.*t8.*t20.*t256;
t314 = m2.*t5.*t35.*t205;
t315 = m1.*t9.*t20.*t256;
t316 = m1.*t9.*t20.*t257;
t317 = m2.*t6.*t35.*t204;
t318 = m2.*t6.*t35.*t205;
t321 = l_AB.*l_A_m2.*m2.*t35.*t205.*2.0;
t322 = m2.*t5.*t36.*t204;
t323 = m2.*t8.*t29.*t257;
t324 = m1.*t9.*t29.*t256;
t325 = m1.*t9.*t29.*t257;
t326 = m2.*t6.*t36.*t204;
t327 = m2.*t6.*t36.*t205;
t330 = l_AB.*l_A_m2.*m2.*t36.*t204.*2.0;
t338 = -t333;
t340 = -t334;
t343 = -t337;
t344 = -t339;
t357 = t13+t16+t17+t19+t66+t67+t74+t75+t79+t80+t81+t82+t87+t88+t89+t90+t91+t92+t93+t94+t112+t114+t116+t118+t183+t188+t217+t227+t235+t236+t237+t240+t246+t249+t256+t257+t297+t298+t300;
t358 = t66+t67+t74+t75+t79+t80+t81+t82+t119+t120+t133+t134+t135+t180+t183+t188+t217+t227+t233+t234+t246+t249+t251+t253+t269+t270+t275+t276+t283+t284+t297+t298;
t331 = -t321;
t332 = -t330;
t359 = t12+t14+t15+t18+t39+t40+t41+t42+t43+t44+t45+t46+t47+t48+t49+t50+t111+t113+t115+t117+t126+t129+t131+t132+t170+t171+t174+t175+t176+t177+t178+t179+t181+t182+t185+t186+t204+t205+t214+t215+t225+t226+t245+t248+t252+t258+t259+t268+t277+t278+t294+t358;
t360 = t21+t22+t23+t24+t25+t26+t27+t28+t51+t52+t53+t54+t55+t56+t57+t58+t59+t60+t61+t62+t63+t64+t65+t68+t69+t70+t71+t72+t73+t76+t77+t78+t95+t96+t97+t98+t99+t100+t101+t102+t103+t104+t105+t106+t107+t108+t109+t110+t121+t122+t123+t124+t125+t127+t128+t130+t136+t137+t138+t139+t140+t141+t142+t143+t146+t147+t148+t149+t150+t151+t152+t153+t154+t155+t156+t157+t158+t159+t160+t161+t162+t163+t164+t165+t166+t167+t168+t169+t190+t191+t192+t193+t194+t195+t196+t197+t198+t199+t202+t203+t206+t207+t208+t209+t210+t211+t212+t213+t219+t220+t221+t222+t223+t224+t229+t230+t241+t242+t243+t244+t254+t255+t260+t261+t262+t263+t264+t265+t266+t267+t279+t280+t281+t282+t285+t286+t287+t288+t295+t296+t303+t304+t305+t306+t307+t308+t309+t310+t311+t312+t313+t314+t315+t316+t317+t318+t319+t320+t322+t323+t324+t325+t326+t327+t328+t329+t331+t332+t335+t336+t338+t340+t343+t344+t345+t346+t347+t348+t353+t354+t355+t356;
t361 = 1.0./t360;
t362 = t357.*t361;
t364 = t358.*t361;
t365 = t359.*t361;
t363 = -t362;
t366 = -t365;
Mass_Joint_Sp_inv = reshape([t361.*(t13+t16+t17+t19+t87+t88+t89+t90+t91+t92+t93+t94+t112+t114+t116+t118+t235+t236+t237+t240+t256+t257+t300),t363,t364,t363,t361.*(t12+t13+t14+t15+t16+t17+t18+t19+t39+t40+t41+t42+t43+t44+t45+t46+t47+t48+t49+t50+t66.*2.0+t67.*2.0+t74.*2.0+t75.*2.0+t79.*2.0+t80.*2.0+t81.*2.0+t82.*2.0+t87+t88+t89+t90+t91+t92+t93+t94+t111+t112+t113+t114+t115+t116+t117+t118+t170+t171+t174+t175+t176+t177+t178+t179+t183.*2.0+t188.*2.0+t204+t205+t217.*2.0+t227.*2.0+t235+t236+t237+t240+t246.*2.0+t249.*2.0+t256+t257-t289.*2.0-t290.*2.0+t294+t300),t366,t364,t366,t361.*(t12+t14+t15+t18+t39+t40+t41+t42+t43+t44+t45+t46+t47+t48+t49+t50+t111+t113+t115+t117+t126.*2.0+t129.*2.0+t131.*2.0+t132.*2.0+t170+t171+t174+t175+t176+t177+t178+t179+t181.*2.0+t182.*2.0+t185.*2.0+t186.*2.0-t187.*2.0-t200.*2.0-t201.*2.0+t204+t205+t214.*2.0+t215.*2.0-t216.*2.0+t225.*2.0+t226.*2.0-t238.*2.0-t239.*2.0+t245.*2.0+t248.*2.0+t294+I1.*I2+I1.*Ir2+I2.*Ir1+Ir1.*Ir2+I2.*m1.*t9.*t20+I2.*m2.*t8.*t20+I2.*m3.*t8.*t20+I2.*m1.*t9.*t29+I2.*m2.*t8.*t29+I2.*m3.*t8.*t29+I1.*m2.*t6.*t35+I1.*m3.*t5.*t35+I1.*m2.*t6.*t36+I1.*m3.*t5.*t36+Ir2.*m1.*t9.*t20+Ir2.*m2.*t8.*t20+Ir2.*m3.*t8.*t20+Ir2.*m1.*t9.*t29+Ir2.*m2.*t8.*t29+Ir2.*m3.*t8.*t29+Ir1.*m2.*t6.*t35+Ir1.*m3.*t5.*t35+Ir1.*m2.*t6.*t36+Ir1.*m3.*t5.*t36+t5.*t8.*t11.*t20.*t36+t6.*t8.*t10.*t20.*t36+t5.*t8.*t11.*t29.*t35+t6.*t8.*t10.*t29.*t35+m1.*m2.*t6.*t9.*t20.*t35+m1.*m3.*t5.*t9.*t20.*t35+m2.*m3.*t5.*t8.*t20.*t35+m1.*m2.*t6.*t9.*t20.*t36+m1.*m3.*t5.*t9.*t20.*t36+m2.*m3.*t5.*t8.*t20.*t36+m2.*m3.*t6.*t8.*t20.*t35+m2.*m3.*t6.*t8.*t20.*t36+m1.*m2.*t6.*t9.*t29.*t35+m1.*m3.*t5.*t9.*t29.*t35+m2.*m3.*t5.*t8.*t29.*t35+m1.*m2.*t6.*t9.*t29.*t36+m1.*m3.*t5.*t9.*t29.*t36+m2.*m3.*t5.*t8.*t29.*t36+m2.*m3.*t6.*t8.*t29.*t35+m2.*m3.*t6.*t8.*t29.*t36-l_AB.*l_A_m2.*m2.*m3.*t8.*t20.*t35.*2.0-l_AB.*l_A_m2.*m2.*m3.*t8.*t29.*t36.*2.0-t2.*t3.*t5.*t8.*t11.*t30.*t31.*2.0-t2.*t3.*t6.*t8.*t10.*t30.*t31.*2.0-l_AB.*l_A_m2.*m2.*m3.*t2.*t3.*t8.*t30.*t31.*4.0)],[3,3]);
