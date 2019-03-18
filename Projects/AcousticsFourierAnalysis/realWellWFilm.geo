DefineConstant[scale = 1];
DefineConstant[factor = 1];
DefineConstant[r1 = 0.09525];
DefineConstant[r2 = 0.01750];
DefineConstant[r3 = 0.03725];
DefineConstant[r4 = 0.01750];
DefineConstant[r5 = 0.04950];
DefineConstant[r6 = 0.00050];
DefineConstant[r7 = 0.08500];

DefineConstant[el_size_1 = 0.0125];
DefineConstant[el_size_2 = 0.005];
DefineConstant[el_size_3 = 0.0175];
DefineConstant[el_size_4 = 0.025];

r1*=scale;
r2*=scale;
r3*=scale;
r4*=scale;
r5*=scale;
r6*=scale;
r7*=scale;


el_size_1*=scale;
el_size_2*=scale;
el_size_3*=scale;
el_size_4*=scale;

dist1=r1;
dist2=r2+dist1;
dist3=r3+dist2;
dist4=r4+dist3;
dist5=r5+dist4;
dist6=r6+dist5;
dist7=r7+dist6;
pi_val = 3.14159265358979323846264338327950288;


p_0 = newp; Point(p_0) = { 0 , 0 , 0 };//center

////////////////////////////INNER CIRCLE
//inner circle
p_c1_1 = newp; Point(p_c1_1) = { dist1 * Cos(pi_val/4) , dist1 * Sin(pi_val/4) , 0 };
p_c1_2 = newp; Point(p_c1_2) = {-dist1 * Cos(pi_val/4) , dist1 * Sin(pi_val/4) , 0 };
p_c1_3 = newp; Point(p_c1_3) = {-dist1 * Cos(pi_val/4) ,-dist1 * Sin(pi_val/4) , 0 };
p_c1_4 = newp; Point(p_c1_4) = { dist1 * Cos(pi_val/4) ,-dist1 * Sin(pi_val/4) , 0 };

//arcs
a_c1_1 = newl; Circle(a_c1_1) = {p_c1_1 , p_0 , p_c1_2};
a_c1_2 = newl; Circle(a_c1_2) = {p_c1_2 , p_0 , p_c1_3};
a_c1_3 = newl; Circle(a_c1_3) = {p_c1_3 , p_0 , p_c1_4};
a_c1_4 = newl; Circle(a_c1_4) = {p_c1_4 , p_0 , p_c1_1};

//surfaces
edges_list_1[] = {a_c1_1,a_c1_2,a_c1_3,a_c1_4};
ll_c1_1 = newll; Line Loop(ll_c1_1) = edges_list_1[];
s_c1_1 = news; Plane Surface(s_c1_1) = {ll_c1_1};


///circle 2
//inner circle
p_c2_1 = newp; Point(p_c2_1) = { dist2 * Cos(pi_val/4) , dist2 * Sin(pi_val/4) , 0 };
p_c2_2 = newp; Point(p_c2_2) = {-dist2 * Cos(pi_val/4) , dist2 * Sin(pi_val/4) , 0 };
p_c2_3 = newp; Point(p_c2_3) = {-dist2 * Cos(pi_val/4) ,-dist2 * Sin(pi_val/4) , 0 };
p_c2_4 = newp; Point(p_c2_4) = { dist2 * Cos(pi_val/4) ,-dist2 * Sin(pi_val/4) , 0 };

//arcs
a_c2_1 = newl; Circle(a_c2_1) = {p_c2_1 , p_0 , p_c2_2};
a_c2_2 = newl; Circle(a_c2_2) = {p_c2_2 , p_0 , p_c2_3};
a_c2_3 = newl; Circle(a_c2_3) = {p_c2_3 , p_0 , p_c2_4};
a_c2_4 = newl; Circle(a_c2_4) = {p_c2_4 , p_0 , p_c2_1};

//surfaces
edges_list_2[] = {a_c2_1,a_c2_2,a_c2_3,a_c2_4};
ll_c2_1 = newll; Line Loop(ll_c2_1) = edges_list_2[];
s_c2_1 = news; Plane Surface(s_c2_1) = {ll_c2_1,ll_c1_1};


///circle 3
//inner circle
p_c3_1 = newp; Point(p_c3_1) = { dist3 * Cos(pi_val/4) , dist3 * Sin(pi_val/4) , 0 };
p_c3_2 = newp; Point(p_c3_2) = {-dist3 * Cos(pi_val/4) , dist3 * Sin(pi_val/4) , 0 };
p_c3_3 = newp; Point(p_c3_3) = {-dist3 * Cos(pi_val/4) ,-dist3 * Sin(pi_val/4) , 0 };
p_c3_4 = newp; Point(p_c3_4) = { dist3 * Cos(pi_val/4) ,-dist3 * Sin(pi_val/4) , 0 };

//arcs
a_c3_1 = newl; Circle(a_c3_1) = {p_c3_1 , p_0 , p_c3_2};
a_c3_2 = newl; Circle(a_c3_2) = {p_c3_2 , p_0 , p_c3_3};
a_c3_3 = newl; Circle(a_c3_3) = {p_c3_3 , p_0 , p_c3_4};
a_c3_4 = newl; Circle(a_c3_4) = {p_c3_4 , p_0 , p_c3_1};

//surfaces
edges_list_3[] = {a_c3_1,a_c3_2,a_c3_3,a_c3_4};
ll_c3_1 = newll; Line Loop(ll_c3_1) = edges_list_3[];
s_c3_1 = news; Plane Surface(s_c3_1) = {ll_c3_1,ll_c2_1};

///circle 4
//inner circle
p_c4_1 = newp; Point(p_c4_1) = { dist4 * Cos(pi_val/4) , dist4 * Sin(pi_val/4) , 0 };
p_c4_2 = newp; Point(p_c4_2) = {-dist4 * Cos(pi_val/4) , dist4 * Sin(pi_val/4) , 0 };
p_c4_3 = newp; Point(p_c4_3) = {-dist4 * Cos(pi_val/4) ,-dist4 * Sin(pi_val/4) , 0 };
p_c4_4 = newp; Point(p_c4_4) = { dist4 * Cos(pi_val/4) ,-dist4 * Sin(pi_val/4) , 0 };

//arcs
a_c4_1 = newl; Circle(a_c4_1) = {p_c4_1 , p_0 , p_c4_2};
a_c4_2 = newl; Circle(a_c4_2) = {p_c4_2 , p_0 , p_c4_3};
a_c4_3 = newl; Circle(a_c4_3) = {p_c4_3 , p_0 , p_c4_4};
a_c4_4 = newl; Circle(a_c4_4) = {p_c4_4 , p_0 , p_c4_1};

//surfaces
edges_list_4[] = {a_c4_1,a_c4_2,a_c4_3,a_c4_4};
ll_c4_1 = newll; Line Loop(ll_c4_1) = edges_list_4[];
s_c4_1 = news; Plane Surface(s_c4_1) = {ll_c4_1,ll_c3_1};

///circle 5
//inner circle
p_c5_1 = newp; Point(p_c5_1) = { dist5 * Cos(pi_val/4) , dist5 * Sin(pi_val/4) , 0 };
p_c5_2 = newp; Point(p_c5_2) = {-dist5 * Cos(pi_val/4) , dist5 * Sin(pi_val/4) , 0 };
p_c5_3 = newp; Point(p_c5_3) = {-dist5 * Cos(pi_val/4) ,-dist5 * Sin(pi_val/4) , 0 };
p_c5_4 = newp; Point(p_c5_4) = { dist5 * Cos(pi_val/4) ,-dist5 * Sin(pi_val/4) , 0 };

//arcs
a_c5_1 = newl; Circle(a_c5_1) = {p_c5_1 , p_0 , p_c5_2};
a_c5_2 = newl; Circle(a_c5_2) = {p_c5_2 , p_0 , p_c5_3};
a_c5_3 = newl; Circle(a_c5_3) = {p_c5_3 , p_0 , p_c5_4};
a_c5_4 = newl; Circle(a_c5_4) = {p_c5_4 , p_0 , p_c5_1};

//surfaces
edges_list_5[] = {a_c5_1,a_c5_2,a_c5_3,a_c5_4};
ll_c5_1 = newll; Line Loop(ll_c5_1) = edges_list_5[];
s_c5_1 = news; Plane Surface(s_c5_1) = {ll_c5_1,ll_c4_1};

///circle 6
//inner circle
p_c6_1 = newp; Point(p_c6_1) = { dist6 * Cos(pi_val/4) , dist6 * Sin(pi_val/4) , 0 };
p_c6_2 = newp; Point(p_c6_2) = {-dist6 * Cos(pi_val/4) , dist6 * Sin(pi_val/4) , 0 };
p_c6_3 = newp; Point(p_c6_3) = {-dist6 * Cos(pi_val/4) ,-dist6 * Sin(pi_val/4) , 0 };
p_c6_4 = newp; Point(p_c6_4) = { dist6 * Cos(pi_val/4) ,-dist6 * Sin(pi_val/4) , 0 };

//arcs
a_c6_1 = newl; Circle(a_c6_1) = {p_c6_1 , p_0 , p_c6_2};
a_c6_2 = newl; Circle(a_c6_2) = {p_c6_2 , p_0 , p_c6_3};
a_c6_3 = newl; Circle(a_c6_3) = {p_c6_3 , p_0 , p_c6_4};
a_c6_4 = newl; Circle(a_c6_4) = {p_c6_4 , p_0 , p_c6_1};

//surfaces
edges_list_6[] = {a_c6_1,a_c6_2,a_c6_3,a_c6_4};
ll_c6_1 = newll; Line Loop(ll_c6_1) = edges_list_6[];
s_c6_1 = news; Plane Surface(s_c6_1) = {ll_c6_1,ll_c5_1};

///circle 7
//inner circle
p_c7_1 = newp; Point(p_c7_1) = { dist7 * Cos(pi_val/4) , dist7 * Sin(pi_val/4) , 0 };
p_c7_2 = newp; Point(p_c7_2) = {-dist7 * Cos(pi_val/4) , dist7 * Sin(pi_val/4) , 0 };
p_c7_3 = newp; Point(p_c7_3) = {-dist7 * Cos(pi_val/4) ,-dist7 * Sin(pi_val/4) , 0 };
p_c7_4 = newp; Point(p_c7_4) = { dist7 * Cos(pi_val/4) ,-dist7 * Sin(pi_val/4) , 0 };

//arcs
a_c7_1 = newl; Circle(a_c7_1) = {p_c7_1 , p_0 , p_c7_2};
a_c7_2 = newl; Circle(a_c7_2) = {p_c7_2 , p_0 , p_c7_3};
a_c7_3 = newl; Circle(a_c7_3) = {p_c7_3 , p_0 , p_c7_4};
a_c7_4 = newl; Circle(a_c7_4) = {p_c7_4 , p_0 , p_c7_1};

//surfaces
edges_list_7[] = {a_c7_1,a_c7_2,a_c7_3,a_c7_4};
ll_c7_1 = newll; Line Loop(ll_c7_1) = edges_list_7[];
s_c7_1 = news; Plane Surface(s_c7_1) = {ll_c7_1,ll_c6_1};


Coherence;

allsurfaces[] = Surface '*';

boundary[] = CombinedBoundary{ Surface{allsurfaces[]}; } ;

Physical Surface("water",1) = {s_c1_1, s_c3_1,s_c6_1};
Physical Surface("steel",2) 	  = {s_c2_1, s_c4_1} ;
Physical Surface("cement",3) = {s_c5_1} ;
Physical Surface("rock",4) = {s_c7_1} ;
Physical Line("boundary", 5)   = boundary[];

waterMaterials[] = Physical Surface(1);
waterMaterialsBound[] = CombinedBoundary{ Surface{waterMaterials[]}; };

steelMaterials[] = Physical Surface(2);
steelMaterialsBound[] = CombinedBoundary{ Surface{steelMaterials[]}; };

cementMaterials[] = Physical Surface(3);
cementMaterialsBound[] = CombinedBoundary{ Surface{cementMaterials[]}; };

rockMaterials[] = Physical Surface(4);
rockMaterialsBound[] = CombinedBoundary{ Surface{rockMaterials[]}; };



Field[1] = MathEval;
Field[1].F = Sprintf("%g", el_size_1);

Field[2] = Restrict;
Field[2].IField = 1;
Field[2].FacesList = {Abs(waterMaterials[])};
Field[2].EdgesList = {Abs(waterMaterialsBound[])};

Field[3] = MathEval;
Field[3].F = Sprintf("%g", el_size_2);

Field[4] = Restrict;
Field[4].IField = 3;
Field[4].FacesList = {Abs(steelMaterials[])};
Field[4].EdgesList = {Abs(steelMaterialsBound[])};

Field[5] = MathEval;
Field[5].F = Sprintf("%g", el_size_3);

Field[6] = Restrict;
Field[6].IField = 5;
Field[6].FacesList = {Abs(cementMaterials[])};
Field[6].EdgesList = {Abs(cementMaterialsBound[])};

Field[7] = MathEval;
Field[7].F = Sprintf("%g", el_size_4);

Field[8] = Restrict;
Field[8].IField = 7;
Field[8].FacesList = {Abs(rockMaterials[])};
Field[8].EdgesList = {Abs(rockMaterialsBound[])};

Field[10] = Min;
Field[10].FieldsList = {2,4,6,8};

Background Field = 10;

Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.ElementOrder=2;
