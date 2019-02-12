/*********************************************************************
 *
 *    homogenenous waveguide(for testing purposes)
 *
 *********************************************************************/
SetFactory("OpenCASCADE");

DefineConstant[length = 10];
DefineConstant[height = 20];
DefineConstant[el_size_1 = 4.];



ndiv_x=Ceil(length/el_size_1) + 1;
ndiv_y=Ceil(height/el_size_1) + 1;
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

ll_1 = newll; Line Loop(ll_1) = {l_1,l_2,l_3,l_4} ;

s_1 = news; Plane Surface(s_1) = {ll_1} ;

Transfinite Line {l_2,l_4} = ndiv_x;


bound[] = {l_1,l_2,l_3,l_4};
Transfinite Line {l_1,l_3} = ndiv_y;
Transfinite Surface{s_1};


////physical entities
Physical Surface("water",1) = {s_1} ;//regular domain
Physical Line("bound1",2) = {l_3};//neumann boundary condition
Physical Line("bound2",3) = {l_1,l_2,l_4} ;//dirichlet boundary condition

