/*********************************************************************
 *
 *    Step fiber waveguide
 *
 *********************************************************************/
DefineConstant[scale = 4053667.9401158620603];
DefineConstant[factor = 1];
DefineConstant[r1 = 0.055];
DefineConstant[r2 = 0.005];
DefineConstant[r3 = 0.050];

r1*=scale;
r2*=scale;
r3*=scale;

dist1=r1;
dist2=r2+dist1;
dist3=r3+dist2;
pi_val = 3.14159265358979323846264338327950288;
bound = dist1 + 5 * 2 * pi_val;//dist1+2lambda
pml_d =  2 * 3 * 2 * pi_val;
nlayers_pml = 3 * factor;
ndiv_circle_arc = 20 * factor;
ndiv_circle_radii = 6 * factor;
ndiv_bound_diag = 4 * factor;

minR = ((r1 < r2) ? r1 : r2 ) < r3 ? ((r1 < r2) ? r1 : r2 ) : r3;
lcInner = minR/2.;
p_0 = newp; Point(p_0) = { 0 , 0 , 0 , lcInner};//center

//inner circle
p_c1_1 = newp; Point(p_c1_1) = { dist1 * Cos(pi_val/4) , dist1 * Sin(pi_val/4) , 0 , lcInner};
p_c1_2 = newp; Point(p_c1_2) = {-dist1 * Cos(pi_val/4) , dist1 * Sin(pi_val/4) , 0 , lcInner};
p_c1_3 = newp; Point(p_c1_3) = {-dist1 * Cos(pi_val/4) ,-dist1 * Sin(pi_val/4) , 0 , lcInner};
p_c1_4 = newp; Point(p_c1_4) = { dist1 * Cos(pi_val/4) ,-dist1 * Sin(pi_val/4) , 0 , lcInner};

//inner square
prop=0.7071067812;
p_s1_1 = newp; Point(p_s1_1) = { dist1 * prop / Sqrt(2) , dist1 * prop / Sqrt(2) , 0 , lcInner};
p_s1_2 = newp; Point(p_s1_2) = {-dist1 * prop / Sqrt(2) , dist1 * prop / Sqrt(2) , 0 , lcInner};
p_s1_3 = newp; Point(p_s1_3) = {-dist1 * prop / Sqrt(2) ,-dist1 * prop / Sqrt(2) , 0 , lcInner};
p_s1_4 = newp; Point(p_s1_4) = { dist1 * prop / Sqrt(2) ,-dist1 * prop / Sqrt(2) , 0 , lcInner};

//inner square sides
t_s_1 = newl; Line(t_s_1) = {p_s1_1,p_s1_2};
t_s_2 = newl; Line(t_s_2) = {p_s1_2,p_s1_3};
t_s_3 = newl; Line(t_s_3) = {p_s1_3,p_s1_4};
t_s_4 = newl; Line(t_s_4) = {p_s1_4,p_s1_1};

//radii
t_c1_1 = newl; Line(t_c1_1) = {p_s1_1,p_c1_1};
t_c1_2 = newl; Line(t_c1_2) = {p_s1_2,p_c1_2};
t_c1_3 = newl; Line(t_c1_3) = {p_s1_3,p_c1_3};
t_c1_4 = newl; Line(t_c1_4) = {p_s1_4,p_c1_4};


//arcs
a_c1_1 = newl; Circle(a_c1_1) = {p_c1_1 , p_0 , p_c1_2};
a_c1_2 = newl; Circle(a_c1_2) = {p_c1_2 , p_0 , p_c1_3};
a_c1_3 = newl; Circle(a_c1_3) = {p_c1_3 , p_0 , p_c1_4};
a_c1_4 = newl; Circle(a_c1_4) = {p_c1_4 , p_0 , p_c1_1};

innerBorder[] = {a_c1_1,a_c1_2,a_c1_3,a_c1_4};
innerRadii[] = {t_c1_1,t_c1_2,t_c1_3,t_c1_4};

//surfaces
c_s_1 = newll; Line Loop(c_s_1) = {t_s_1,t_s_2,t_s_3,t_s_4};
s_s_1 = news; Surface(s_s_1) = {c_s_1};

c_c1_1 = newll; Line Loop(c_c1_1) = {t_c1_1,a_c1_1,-t_c1_2, -t_s_1};
s_c1_1 = news; Surface (s_c1_1) = {c_c1_1};
c_c1_2 = newll; Line Loop(c_c1_2) = {t_c1_2,a_c1_2,-t_c1_3, - t_s_2};
s_c1_2 = news; Surface (s_c1_2) = {c_c1_2};
c_c1_3 = newll; Line Loop(c_c1_3) = {t_c1_3,a_c1_3,-t_c1_4, - t_s_3};
s_c1_3 = news; Surface (s_c1_3) = {c_c1_3};
c_c1_4 = newll; Line Loop(c_c1_4) = {t_c1_4,a_c1_4,-t_c1_1, - t_s_4};
s_c1_4 = news; Surface (s_c1_4) = {c_c1_4};



///circle 2
//inner circle
p_c2_1 = newp; Point(p_c2_1) = { dist2 * Cos(pi_val/4) , dist2 * Sin(pi_val/4) , 0 , lcInner};
p_c2_2 = newp; Point(p_c2_2) = {-dist2 * Cos(pi_val/4) , dist2 * Sin(pi_val/4) , 0 , lcInner};
p_c2_3 = newp; Point(p_c2_3) = {-dist2 * Cos(pi_val/4) ,-dist2 * Sin(pi_val/4) , 0 , lcInner};
p_c2_4 = newp; Point(p_c2_4) = { dist2 * Cos(pi_val/4) ,-dist2 * Sin(pi_val/4) , 0 , lcInner};

//radii
t_c2_1 = newl; Line(t_c2_1) = {p_c1_1,p_c2_1};
t_c2_2 = newl; Line(t_c2_2) = {p_c1_2,p_c2_2};
t_c2_3 = newl; Line(t_c2_3) = {p_c1_3,p_c2_3};
t_c2_4 = newl; Line(t_c2_4) = {p_c1_4,p_c2_4};


//arcs
a_c2_1 = newl; Circle(a_c2_1) = {p_c2_1 , p_0 , p_c2_2};
a_c2_2 = newl; Circle(a_c2_2) = {p_c2_2 , p_0 , p_c2_3};
a_c2_3 = newl; Circle(a_c2_3) = {p_c2_3 , p_0 , p_c2_4};
a_c2_4 = newl; Circle(a_c2_4) = {p_c2_4 , p_0 , p_c2_1};


c_c2_1 = newll; Line Loop(c_c2_1) = {t_c2_1,a_c2_1,-t_c2_2, -a_c1_1};
s_c2_1 = news; Surface (s_c2_1) = {c_c2_1};
c_c2_2 = newll; Line Loop(c_c2_2) = {t_c2_2,a_c2_2,-t_c2_3, - a_c1_2};
s_c2_2 = news; Surface (s_c2_2) = {c_c2_2};
c_c2_3 = newll; Line Loop(c_c2_3) = {t_c2_3,a_c2_3,-t_c2_4, - a_c1_3};
s_c2_3 = news; Surface (s_c2_3) = {c_c2_3};
c_c2_4 = newll; Line Loop(c_c2_4) = {t_c2_4,a_c2_4,-t_c2_1, - a_c1_4};
s_c2_4 = news; Surface (s_c2_4) = {c_c2_4};


///circle 3
//inner circle
p_c3_1 = newp; Point(p_c3_1) = { dist3 * Cos(pi_val/4) , dist3 * Sin(pi_val/4) , 0 , lcInner};
p_c3_2 = newp; Point(p_c3_2) = {-dist3 * Cos(pi_val/4) , dist3 * Sin(pi_val/4) , 0 , lcInner};
p_c3_3 = newp; Point(p_c3_3) = {-dist3 * Cos(pi_val/4) ,-dist3 * Sin(pi_val/4) , 0 , lcInner};
p_c3_4 = newp; Point(p_c3_4) = { dist3 * Cos(pi_val/4) ,-dist3 * Sin(pi_val/4) , 0 , lcInner};

//radii
t_c3_1 = newl; Line(t_c3_1) = {p_c2_1,p_c3_1};
t_c3_2 = newl; Line(t_c3_2) = {p_c2_2,p_c3_2};
t_c3_3 = newl; Line(t_c3_3) = {p_c2_3,p_c3_3};
t_c3_4 = newl; Line(t_c3_4) = {p_c2_4,p_c3_4};


//arcs
a_c3_1 = newl; Circle(a_c3_1) = {p_c3_1 , p_0 , p_c3_2};
a_c3_2 = newl; Circle(a_c3_2) = {p_c3_2 , p_0 , p_c3_3};
a_c3_3 = newl; Circle(a_c3_3) = {p_c3_3 , p_0 , p_c3_4};
a_c3_4 = newl; Circle(a_c3_4) = {p_c3_4 , p_0 , p_c3_1};


c_c3_1 = newll; Line Loop(c_c3_1) = {t_c3_1,a_c3_1,-t_c3_2, -a_c2_1};
s_c3_1 = news; Surface (s_c3_1) = {c_c3_1};
c_c3_2 = newll; Line Loop(c_c3_2) = {t_c3_2,a_c3_2,-t_c3_3, - a_c2_2};
s_c3_2 = news; Surface (s_c3_2) = {c_c3_2};
c_c3_3 = newll; Line Loop(c_c3_3) = {t_c3_3,a_c3_3,-t_c3_4, - a_c2_3};
s_c3_3 = news; Surface (s_c3_3) = {c_c3_3};
c_c3_4 = newll; Line Loop(c_c3_4) = {t_c3_4,a_c3_4,-t_c3_1, - a_c2_4};
s_c3_4 = news; Surface (s_c3_4) = {c_c3_4};

Transfinite Surface{s_c1_1,s_c1_2,s_c1_3,s_c1_4};//inner well
Transfinite Surface{s_s_1};//inner square
Transfinite Line{a_c1_1,a_c1_2,a_c1_3,a_c1_4} = ndiv_circle_arc;//inner well
Transfinite Line{t_s_1,t_s_2,t_s_3,t_s_4} = ndiv_circle_arc;//inner square
//Transfinite Line{t_c1_1,t_c1_2,t_c1_3,t_c1_4} = ndiv_circle_radii;//radii inner well

Transfinite Surface{s_c2_1,s_c2_2,s_c2_3,s_c2_4};//casing
Transfinite Line{a_c2_1,a_c2_2,a_c2_3,a_c2_4} = ndiv_circle_arc;//casing
//Transfinite Line{t_c2_1,t_c2_2,t_c2_3,t_c2_4} = ndiv_circle_radii-1;//radii casing well

Transfinite Surface{s_c3_1,s_c3_2,s_c3_3,s_c3_4};//casing
Transfinite Line{a_c3_1,a_c3_2,a_c3_3,a_c3_4} = ndiv_circle_arc;//casing
//Transfinite Line{t_c3_1,t_c3_2,t_c3_3,t_c3_4} = ndiv_circle_radii-1;//radii casing well

//boundary
//p_b_1 = 101; Point(p_b_1) = { bound , bound , 0 , lcBound};
//p_b_2 = 102; Point(p_b_2) = {-bound , bound , 0 , lcBound};
//p_b_3 = 103; Point(p_b_3) = {-bound ,-bound , 0 , lcBound};
//p_b_4  = 100; Point(p_b_4)  = { bound ,-bound , 0 , lcBound};


////boundary region delimiters
//t_b_1 = newl; Line(t_b_1) = {p_c1_1,p_b_1};
//t_b_2 = newl; Line(t_b_2) = {p_c1_2,p_b_2};
//t_b_3 = newl; Line(t_b_3) = {p_c1_3,p_b_3};
//t_b_4 = newl; Line(t_b_4) = {p_c1_4,p_b_4};

//a_b_1 = newl; Line(a_b_1) = {p_b_1 , p_b_2};
//a_b_2 = newl; Line(a_b_2) = {p_b_2 , p_b_3};
//a_b_3 = newl; Line(a_b_3) = {p_b_3 , p_b_4 };
//a_b_4  = newl; Line(a_b_4)  = {p_b_4  , p_b_1};

//boundaryOrig[] = {a_b_1,a_b_2,a_b_3,a_b_4};

//Printf ("now?");
//c_b_1 = newll; Line Loop(c_b_1) = {t_b_1,a_b_1 , -t_b_2 , -a_c1_1};
//c_b_2 = newll; Line Loop(c_b_2) = {t_b_2,a_b_2 , -t_b_3 , -a_c1_2};
//c_b_3 = newll; Line Loop(c_b_3) = {t_b_3,a_b_3 , -t_b_4 , -a_c1_3};
//c_b_4 = newll; Line Loop(c_b_4) = {t_b_4,a_b_4 , -t_b_1 , -a_c1_4};
//Printf ("yes!");
//s_b_1 = news; Surface (s_b_1) = {c_b_1};
//s_b_2 = news; Surface (s_b_2) = {c_b_2};
//s_b_3 = news; Surface (s_b_3) = {c_b_3};
//s_b_4 = news; Surface (s_b_4) = {c_b_4};


////+x
//pml_xP_1  = newp; Point(pml_xP_1)  = { bound + pml_d,-bound , 0 , pml_d/nlayers_pml};
//pml_xP_2  = newp; Point(pml_xP_2)  = { bound + pml_d,+bound , 0 , pml_d/nlayers_pml};
////+y
//pml_yP_1  = newp; Point(pml_yP_1)  = { bound,bound + pml_d , 0 , pml_d/nlayers_pml};
//pml_yP_2  = newp; Point(pml_yP_2)  = {-bound,bound + pml_d , 0 , pml_d/nlayers_pml};
////-x
//pml_xM_1  = newp; Point(pml_xM_1)  = {-bound - pml_d,+bound , 0 , pml_d/nlayers_pml};
//pml_xM_2  = newp; Point(pml_xM_2)  = {-bound - pml_d,-bound , 0 , pml_d/nlayers_pml};
////-y
//pml_yM_1  = newp; Point(pml_yM_1)  = {-bound, -bound - pml_d , 0 , pml_d/nlayers_pml};
//pml_yM_2  = newp; Point(pml_yM_2)  = { bound, -bound - pml_d , 0 , pml_d/nlayers_pml};
////+x-y
//pml_co_1  = newp; Point(pml_co_1)  = { bound + pml_d, -bound - pml_d, 0 , pml_d/nlayers_pml};
////+x+y
//pml_co_2  = newp; Point(pml_co_2)  = { bound + pml_d,  bound + pml_d, 0 , pml_d/nlayers_pml};
////-x+y
//pml_co_3  = newp; Point(pml_co_3)  = {-bound - pml_d,  bound + pml_d, 0 , pml_d/nlayers_pml};
////-x-y
//pml_co_4  = newp; Point(pml_co_4)  = {-bound - pml_d, -bound - pml_d, 0 , pml_d/nlayers_pml};


////+x
//a_xp_1  = newl; Line(a_xp_1)  = {p_b_1,pml_xP_2};
//a_xp_2  = newl; Line(a_xp_2)  = {pml_xP_2,pml_xP_1};
//a_xp_3  = newl; Line(a_xp_3)  = {pml_xP_1,p_b_4};
//a_xp    = newll; Line Loop (a_xp) = {-a_b_4,-a_xp_3,-a_xp_2,-a_xp_1};
//s_xp    = news; Surface (s_xp) = {a_xp};
////+y
//a_yp_1  = newl; Line(a_yp_1)  = {p_b_2,pml_yP_2};
//a_yp_2  = newl; Line(a_yp_2)  = {pml_yP_2,pml_yP_1};
//a_yp_3  = newl; Line(a_yp_3)  = {pml_yP_1,p_b_1};
//a_yp    = newll; Line Loop (a_yp) = {a_b_1,a_yp_1,a_yp_2,a_yp_3};
//s_yp    = news; Surface (s_yp) = {-a_yp};
////-x
//a_xm_1  = newl; Line(a_xm_1)  = {p_b_3,pml_xM_2};
//a_xm_2  = newl; Line(a_xm_2)  = {pml_xM_2,pml_xM_1};
//a_xm_3  = newl; Line(a_xm_3)  = {pml_xM_1,p_b_2};
//a_xm    = newll; Line Loop (a_xm) = {a_b_2,a_xm_1,a_xm_2,a_xm_3};
//s_xm    = news; Surface (s_xm) = {-a_xm};
////-y
//a_ym_1  = newl; Line(a_ym_1)  = {p_b_4,pml_yM_2};
//a_ym_2  = newl; Line(a_ym_2)  = {pml_yM_2,pml_yM_1};
//a_ym_3  = newl; Line(a_ym_3)  = {pml_yM_1,p_b_3};
//a_ym    = newll; Line Loop (a_ym) = {a_b_3,a_ym_1,a_ym_2,a_ym_3};
//s_ym    = news; Surface (s_ym) = {-a_ym};
////+x-y
//a_xpym_1  = newl; Line(a_xpym_1)  = {pml_yM_2,pml_co_1};
//a_xpym_2  = newl; Line(a_xpym_2)  = {pml_co_1,pml_xP_1};
//a_xpym    = news; Line Loop (a_xpym) = {a_ym_1,a_xpym_1,a_xpym_2,a_xp_3};
//s_xpym    = news; Ruled Surface (s_xpym) = { a_xpym };
////+x+y
//a_xpyp_1  = newl; Line(a_xpyp_1)  = {pml_xP_2,pml_co_2};
//a_xpyp_2  = newl; Line(a_xpyp_2)  = {pml_co_2,pml_yP_1};
//a_xpyp    = news; Line Loop (a_xpyp) = {a_xp_1,a_xpyp_1,a_xpyp_2,a_yp_3};
//s_xpyp    = news; Ruled Surface (s_xpyp) = { a_xpyp };
////-x+y
//a_xmyp_1  = newl; Line(a_xmyp_1)  = {pml_yP_2,pml_co_3};
//a_xmyp_2  = newl; Line(a_xmyp_2)  = {pml_co_3,pml_xM_1};
//a_xmyp    = news; Line Loop (a_xmyp) = {a_yp_1,a_xmyp_1,a_xmyp_2,a_xm_3};
//s_xmyp    = news; Ruled Surface (s_xmyp) = { a_xmyp };
////-x-y
//a_xmym_1  = newl; Line(a_xmym_1)  = {pml_xM_2,pml_co_4};
//a_xmym_2  = newl; Line(a_xmym_2)  = {pml_co_4,pml_yM_1};
//a_xmym    = news; Line Loop (a_xmym) = {a_xm_1,a_xmym_1,a_xmym_2,a_ym_3};
//s_xmym    = news; Ruled Surface (s_xmym) = { a_xmym };

////totalboundary
//boundary[] = {-a_xp_2,a_xpyp_1,a_xpyp_2,-a_yp_2,a_xmyp_1,a_xmyp_2,-a_xm_2,a_xmym_1,a_xmym_2,-a_ym_2,a_xpym_1,a_xpym_2};



////Field[1] = Attractor;
////Field[1].NodesList = {p_0};

////Field[2] = Threshold;
////Field[2].IField = 1;
////Field[2].LcMin = lcInner;
////Field[2].LcMax = lcBound;
////Field[2].DistMin = dist1;
////Field[2].DistMax = bound;
////Field[2].StopAtDistMax = 1;
////Background Field = 2;
//Transfinite Surface{s_xp,s_yp,s_xm,s_ym,s_xpym,s_xpyp,s_xmyp,s_xmym};
//Transfinite Surface{s_b_1,s_b_2,s_b_3,s_b_4};

//Transfinite Line{a_b_4,a_b_1,a_b_2,a_b_3,a_xp_2,a_yp_2,a_xm_2,a_ym_2} = ndiv_circle_arc;
//Transfinite Line{t_b_1,t_b_2,t_b_3,t_b_4} = ndiv_bound_diag;
//Transfinite Line{a_xp_1,a_xp_3,a_yp_1,a_yp_3,a_xm_1,a_xm_3,a_ym_1,a_ym_3} = nlayers_pml + 1;

////////physical entities
//Physical Surface(1) = {s_s_1,s_c_1,s_c_2,s_c_3,s_c_4};//core domain
//Physical Surface(2) = {s_b_1,s_b_2,s_b_3,s_b_4};//air domain
//Physical Surface(10) = {s_xp};//xp
//Physical Surface(11) = {s_yp};//yp
//Physical Surface(12) = {s_xm};//xm
//Physical Surface(13) = {s_ym};//ym
//Physical Surface(14) = {s_xpym};//xpym
//Physical Surface(15) = {s_xpyp};//xpyp
//Physical Surface(16) = {s_xmyp};//xmyp
//Physical Surface(17) = {s_xmym};//xmym
//Physical Line(18) = boundary[] ;//dirichlet boundary condition

Coherence;