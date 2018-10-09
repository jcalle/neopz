/*********************************************************************
 *
 *    homogenenous waveguide(for testing purposes)
 *
 *********************************************************************/
SetFactory("OpenCASCADE");

DefineConstant[length = 60];
DefineConstant[height = 20];
DefineConstant[el_size = 4.];
DefineConstant[pml_length = 20];



ndiv_x=Ceil(length/el_size) + 1;
ndiv_y=Ceil(height/el_size) + 1;
ndiv_pml=Ceil(pml_length/el_size);
lc = 1e6;

//inner domain
p_1 = newp; Point(p_1) = { length, 0, 0, lc} ;
p_2 = newp; Point(p_2) = { length, height, 0, lc} ;
p_3 = newp; Point(p_3) = { 0, height, 0, lc} ;
p_4 = newp; Point(p_4) = { 0, 0, 0, lc} ;

l_1 = newl; Line(l_1) = {p_1,p_2} ;
l_2 = newl; Line(l_2) = {p_2,p_3} ;
l_3 = newl; Line(l_3) = {p_3,p_4} ;
l_4 = newl; Line(l_4) = {p_4,p_1} ;

//pml domain
p_5 = newp; Point(p_5) = { length + pml_length, 0, 0, lc} ; //right
p_6 = newp; Point(p_6) = { length + pml_length, height, 0, lc} ; //right

p_7 = newp; Point(p_7) = { length + pml_length, height + pml_length, 0, lc} ; //upright

p_8= newp; Point(p_8) = { length, height + pml_length, 0, lc} ; //up
p_9 = newp; Point(p_9) = { 0     , height + pml_length, 0, lc} ; //up

p_10 = newp; Point(p_10) = { -pml_length, height + pml_length, 0, lc} ; //upleft


p_11 = newp; Point(p_11) = { -pml_length, height, 0, lc} ; //left
p_12 = newp; Point(p_12) = { -pml_length, 0, 0, lc} ; //left

p_13 = newp; Point(p_13) = { 0 - pml_length, 0 - pml_length, 0, lc} ; //downleft

p_14 = newp; Point(p_14) = { 0, 0 - pml_length, 0, lc} ; //down
p_15 = newp; Point(p_15) = { length, 0 - pml_length, 0, lc} ; //down

p_16 = newp; Point(p_16) = { length + pml_length, 0 - pml_length, 0, lc} ; //downright

l_5 = newl; Line(l_5) = {p_1,p_5} ;
l_6 = newl; Line(l_6) = {p_5,p_6} ;
l_7 = newl; Line(l_7) = {p_6,p_2} ;

l_8 = newl; Line(l_8) = {p_6,p_7} ;
l_9 = newl; Line(l_9) = {p_7,p_8} ;
l_10 = newl; Line(l_10) = {p_8,p_2} ;

l_11 = newl; Line(l_11) = {p_8,p_9} ;
l_12 = newl; Line(l_12) = {p_9,p_3} ;

l_13 = newl; Line(l_13) = {p_9,p_10} ;
l_14 = newl; Line(l_14) = {p_10,p_11} ;
l_15 = newl; Line(l_15) = {p_11,p_3} ;


l_16 = newl; Line(l_16) = {p_11,p_12} ;
l_17 = newl; Line(l_17) = {p_12,p_4} ;

l_18 = newl; Line(l_18) = {p_12,p_13} ;
l_19 = newl; Line(l_19) = {p_13,p_14} ;
l_20 = newl; Line(l_20) = {p_14,p_4} ;

l_21 = newl; Line(l_21) = {p_14,p_15} ;
l_22 = newl; Line(l_22) = {p_15,p_1} ;

l_23 = newl; Line(l_23) = {p_15,p_16} ;
l_24 = newl; Line(l_24) = {p_16,p_5} ;


ll_1 = newll; Line Loop(ll_1) = {l_1,l_2,l_3,l_4} ;

ll_2 = newll; Line Loop(ll_2) = {l_6,l_7,-l_1,l_5} ;
ll_3 = newll; Line Loop(ll_3) = {l_8,l_9,l_10,l_7} ;
ll_4 = newll; Line Loop(ll_4) = {l_11,l_12,-l_2,-l_10} ;
ll_5 = newll; Line Loop(ll_5) = {l_13,l_14,l_15,-l_12} ;
ll_6 = newll; Line Loop(ll_6) = {l_16,l_17,-l_3,-l_15} ;
ll_7 = newll; Line Loop(ll_7) = {l_18,l_19,l_20,-l_17} ;
ll_8 = newll; Line Loop(ll_8) = {l_21,l_22,-l_4,-l_20} ;
ll_9 = newll; Line Loop(ll_9) = {l_23,l_24,-l_5,-l_22} ;

s_1 = news; Plane Surface(s_1) = {ll_1} ;

s_2 = news; Plane Surface(s_2) = {ll_2} ;
s_3 = news; Plane Surface(s_3) = {ll_3} ;
s_4 = news; Plane Surface(s_4) = {ll_4} ;
s_5 = news; Plane Surface(s_5) = {ll_5} ;
s_6 = news; Plane Surface(s_6) = {ll_6} ;
s_7 = news; Plane Surface(s_7) = {ll_7} ;
s_8 = news; Plane Surface(s_8) = {ll_8} ;
s_9 = news; Plane Surface(s_9) = {ll_9} ;


bound[] = {l_6,l_8,l_9,l_11,l_13,l_14,l_16,l_18,l_19,l_21,l_23,l_24};
Transfinite Line {l_2,l_4,l_11,l_21} = ndiv_x; 
Transfinite Line {l_1,l_3,l_6,l_16} = ndiv_y;
Transfinite Line {l_5,l_7,l_8,l_9,l_10,l_12,l_13,l_14,l_15,l_17,l_18,l_19,l_20,l_22,l_23,l_24} = ndiv_pml; 

Transfinite Surface{s_1};
Transfinite Surface{s_2};
Transfinite Surface{s_3};
Transfinite Surface{s_4};
Transfinite Surface{s_5};
Transfinite Surface{s_6};
Transfinite Surface{s_7};
Transfinite Surface{s_8};
Transfinite Surface{s_9};

////physical entities
Physical Surface("water",1) = {s_1} ;//regular domain
Physical Surface("pmlXP",2) = {s_2} ;//pml domain
Physical Surface("pmlXPYP",3) = {s_3} ;//pml domain
Physical Surface("pmlYP",4) = {s_4} ;//pml domain
Physical Surface("pmlXMYP",5) = {s_5} ;//pml domain
Physical Surface("pmlXM",6) = {s_6} ;//pml domain
Physical Surface("pmlXMYM",7) = {s_7} ;//pml domain
Physical Surface("pmlYM",8) = {s_8} ;//pml domain
Physical Surface("pmlXPYM",9) = {s_9} ;//pml domain
Physical Line("bound",10) = bound[] ;//dirichlet boundary condition

