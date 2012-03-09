#include "matrix.h"

int main()
{
     Matrix a;
     Matrix G(5,5),Gdiv(5,5);
     G.iran();
    // G.output();
     Gdiv = G / 2;
     //Gdiv.output();
     Matrix fx;
     Matrix fy;
     a.load("excuse.dat"); 
     //a.loadImage("/home/pradeep/class/poems/leisure/project/dataset/iso_sam/ru/1000.png"); 
     int r,c;
     a.size(r,c);
     Matrix phi((double)1,r,c),u(r,c),u1(r,c),u2(r,c),u3(r,c),spf(r,c);
phi.mask();
phi.ineg();
u = phi;
//u.dump();
//a.dump();
     //return 0;
double del=1,mu=25;
//u.output();
     //a.iran();

     //a.output();
     //fx.output();
     short iter=0,i;
     double  c1,c2,me;
FOR(i,iter){
	fx = u.gradX();
	fy = u.gradY();
	u1 = u < 0;
	u2 = a.elemul(u1);
	c1 = u2.isum() / u1.isum();
	u3 = u >= 0; 
	u2 = a.elemul(u3);
	c2 = u2.isum() / u3.isum();
	me = (c1 + c2 )/2;
	spf = a.iminus(me);
	spf = spf / spf.imax();
	u= spf.elemul((fx.ipow(2) + fy.ipow(2)).isqrt()).imul(mu).imul(del) + u;
	u = (u >= 0) - (u < 0);
	u.conv2d();
}
//u.output();
//u.dump();
//a.conv2d();
a.Contour();
//a.dump();

//u.dim();
}