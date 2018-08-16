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
p_5 = newp; Point(p_5) = { length + pml_length, 0, 0, lc} ;
p_6 = newp; Point(p_6) = { length + pml_length, height, 0, lc} ;
p_7 = newp; Point(p_7) = { -pml_length, height, 0, lc} ;
p_8 = newp; Point(p_8) = { -pml_length, 0, 0, lc} ;

l_5 = newl; Line(l_5) = {p_1,p_5} ;
l_6 = newl; Line(l_6) = {p_5,p_6} ;
l_7 = newl; Line(l_7) = {p_6,p_2} ;

l_8 = newl; Line(l_8) = {p_3,p_7} ;
l_9 = newl; Line(l_9) = {p_7,p_8} ;
l_10 = newl; Line(l_10) = {p_8,p_4} ;

ll_1 = newll; Line Loop(ll_1) = {l_1,l_2,l_3,l_4} ;
ll_2 = newll; Line Loop(ll_2) = {l_6,l_7,l_1,l_5} ;
ll_3 = newll; Line Loop(ll_3) = {l_3,l_8,l_9,l_10} ;

s_1 = news; Plane Surface(s_1) = {ll_1} ;
s_2 = news; Plane Surface(s_2) = {ll_2} ;
s_3 = news; Plane Surface(s_3) = {ll_3} ;

bound[] = {l_7,l_2,l_8,l_9,l_10,l_4,l_5,l_6};
Transfinite Line {l_2,l_4} = ndiv_x; 
Transfinite Line {l_1,l_3,l_6,l_9} = ndiv_y;
Transfinite Line {l_5,l_7,l_8,l_10} = ndiv_pml; 

Transfinite Surface{s_1};
Transfinite Surface{s_2};
Transfinite Surface{s_3};

////physical entities
Physical Surface("water",1) = {s_1} ;//regular domain
Physical Surface("pmlLeft",2) = {s_2} ;//pml domain
Physical Surface("pmlRight",3) = {s_3} ;//pml domain
Physical Line("bound",4) = bound[] ;//dirichlet boundary condition

