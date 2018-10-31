/*********************************************************************
 *
 *    Step fiber waveguide
 *
 *********************************************************************/
DefineConstant[scale = 4053667.9401158620603];
DefineConstant[factor = 1];
DefineConstant[el_size = 0.008333333333];
DefineConstant[r1 = 0.055];
DefineConstant[r2 = 0.005];
DefineConstant[r3 = 0.050];

r1*=scale;
r2*=scale;
r3*=scale;
el_size*=scale;

dist1=r1;
dist2=r2+dist1;
dist3=r3+dist2;
pi_val = 3.14159265358979323846264338327950288;

//ndiv_circle_arc = 20 * factor;
//ndiv_circle_radii = 6 * factor;
//ndiv_bound_diag = 4 * factor;

ndiv_circle_arc = Ceil(pi_val/2 * dist3 / el_size)+1;

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
//Transfinite Line{t_c2_1,t_c2_2,t_c2_3,t_c2_4} = 2;//ndiv_circle_radii-1;//radii casing well

Transfinite Surface{s_c3_1,s_c3_2,s_c3_3,s_c3_4};//casing
Transfinite Line{a_c3_1,a_c3_2,a_c3_3,a_c3_4} = ndiv_circle_arc;//casing
//Transfinite Line{t_c3_1,t_c3_2,t_c3_3,t_c3_4} = 3;//radii casing well

//boundary


allsurfaces[] = Surface '*'; // For recovering all volumes just after the box for the air

boundary[] = CombinedBoundary{ Surface{allsurfaces[]}; } ; // Boundary of all the volumes for your surface loop

Physical Surface("water_inner",1) = {s_s_1,s_c1_1,s_c1_2,s_c1_3,s_c1_4};
Physical Surface("casing",2) 	  = {s_c2_1,s_c2_2,s_c2_3,s_c2_4} ;//pml domain
Physical Surface("water_outer",3) = {s_c3_1,s_c3_2,s_c3_3,s_c3_4} ;//pml domain
Physical Surface("boundary", 4)   = boundary[];
Coherence;