#ifndef MY_MATRIX_H
#define MY_MATRIX_H
#include <iomanip>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <time.h>
#include <math.h>
#include <algorithm>
#include <string>
#include <Imlib2.h> 
#include "bitmaplib.h"
using namespace std;
const double tol = .00001;
#define FOR(i,n) for(i = 0;i < n; i++)
#define FOR3(i,x,n) for(i = x;i < n; i++)
#define MAX(i,n) max(i,n)
#define MIN(i,n) max(i,n)
class Matrix
{
	private:

		double **m;

		short r;
		
		short c;
		short i;
		short j;
	long long nump;
	int NX;
int NY;
		BITMAP4 *image1;
	public:
		        Matrix (  ) { r = 0 ; c = 0 ; m = NULL;
	NX = 34;
NY = 32;
nump = 0;
	}
		void loadImage( const string & file ){
			ifstream fin;
		Imlib_Image *image;

image =(Imlib_Image *)imlib_load_image(file.c_str());
imlib_context_set_image(image);
				imlib_context_set_anti_alias(1);
c=imlib_image_get_width(),r=imlib_image_get_height();
DATA32 *data;
data = imlib_image_get_data_for_reading_only();
	
    m=new double * [ r ];
	
    for ( i=0; i < r; i++ )
		m[ i ] = new double [ c ];
int ir,ig,ib;
FOR(i,r){
FOR(j,c){
 ir = ((*data >> 16) & 0xff)* 0.299;
 ig = ((*data >> 8) & 0xff)* 0.587;
 ib = (*data & 0xff) *0.114;
m[i][j]=((ir + ig + ib) & 0xff);
data++;
}
}
		}
		void mask(){
			FOR3(i,10,r-10){
FOR3(j,10,c-10){
            m[i][j] = -1;
}
}
		}
	void dump(){
	string filename="edump.dat";
	ofstream out(filename.c_str (  ), ios::out );
    //out << r<<" ";
	//out << c<<endl;
	 //double max = imax();
	FOR(i,r){
FOR(j,c){
            //out <<setprecision(1)<< (double)(m[i][j] / max)<<" ";
            out << m[i][j]<<" ";
}
out<<endl;
}	
	}
void size(int& w,int& u)	
{
w = r;
u = c;
}
void dim()	
{
cout<<r<<" "<<c<<endl;
}
Matrix elemul(Matrix  mat)
{
	    Matrix tmp(r,c);
	FOR(i,r){
FOR(j,c){
            tmp.m[i][j] = m[i][j]*mat.m[i][j];
}
}
        return tmp;
}

Matrix imul(double x)
{
	    Matrix tmp(r,c);
	FOR(i,r){
FOR(j,c){
            tmp.m[i][j] = m[i][j]*x;
}
}
        return tmp;
}

Matrix operator+(Matrix  mat)
{
	    Matrix tmp(r,c);
	FOR(i,r){
FOR(j,c){
            tmp.m[i][j] = m[i][j] + mat.m[i][j];
}
}
        return tmp;
}
Matrix operator-(Matrix  mat)
{
	    Matrix tmp(r,c);
	FOR(i,r){
FOR(j,c){
            tmp.m[i][j] = m[i][j] - mat.m[i][j];
}
}
        return tmp;
}
//~ void conv2d(Matrix x)
	//~ {
		//~ FOR(i,r){
//~ FOR(j,c){
	//~ }
//~ }
//~ }
Matrix gradX()
	{
	    Matrix fx(r,c);
	short jN,iE;
	double Cval,Eval,Nval;
	FOR(j,r){// foreach row
jN = j+1;
if(jN>r-1) jN=r-1;
FOR(i,c){
iE = i+1;
if(iE>c-1) iE=c-1;
 Cval = m[j][i];

 Eval = m[j][iE];
 Nval = m[jN][i];
fx.m[j][i] = (Eval-Cval); //gradient computation with forward difference
	}
}
return fx;
}
Matrix gradY()
	{
		    Matrix fy(r,c);
	short jN,iE;
	double Cval,Eval,Nval;
	FOR(j,r){// foreach row
jN = j+1;
if(jN>r-1) jN=r-1;
FOR(i,c){
iE = i+1;
if(iE>c-1) iE=c-1;
 Cval = m[j][i];
 Eval = m[j][iE];
 Nval = m[jN][i];
fy.m[j][i] = (Nval-Cval); 
	}
}
return fy;
}
double imax(){
double max = m[0][0];
FOR(i,r){
FOR(j,c){
if(m[i][j] > max){
max = m[i][j];
}
}
}
return max;
}

double imin(){
double min = m[0][0];
FOR(i,r){
FOR(j,c){
if(m[i][j] < min){
min = m[i][j];
}
}
}
return min;
}

void iran(){
srand(time(NULL));
FOR(i,r){
FOR(j,c){
m[i][j] = rand() * 1e-10;
}
}
}
Matrix ipow(double x){
double max = m[0][0];
FOR(i,r){
FOR(j,c){
m[i][j] = pow(m[i][j],x);
}
}
return *this;
}

void ineg(){
double max = m[0][0];
FOR(i,r){
FOR(j,c){
m[i][j] = -m[i][j];
}
}
}
Matrix iminus(double x){
FOR(i,r){
FOR(j,c){
m[i][j] -= x;
}
}
return *this;
}
Matrix isqrt(){
double max = m[0][0];
FOR(i,r){
FOR(j,c){
m[i][j] = sqrt((m[i][j]));
}
}
return *this;
}
double isum(){
double  sum;
FOR(i,r){
FOR(j,c){
sum += m[i][j] ;
}
}
return sum;
}
Matrix operator/ (double x)
{
Matrix tmp(r,c);
FOR(i,r){
FOR(j,c){
tmp.m[i][j] = m[i][j] / x ;
}
}
return tmp;
}

Matrix operator>=(int w)
{
	Matrix tmp(r,c);
	for (int i = 0; i < r; i++)
	{
		for (int j = 0; j < c; j++)
		{
			if (m[i][j] >= 0)
				tmp.m[i][j] = 1;
			else tmp.m[i][j] = 0;
		}
	}
	return tmp;
}

Matrix operator<(int w)
{
	Matrix tmp(r,c);
	for (int i = 0; i < r; i++)
	{
		for (int j = 0; j < c; j++)
		{
			if (m[i][j] < 0)
				tmp.m[i][j] = 1;
			else tmp.m[i][j] = 0;
		}
	}
	return tmp;
}

void conv2d(){
double sum;
	Matrix tmp(r,c);
int M = r / 2, N = c /2;
FOR(i,M){
FOR(j,N){
sum += m[i][j] ;
m[i][j] = sum / (r*c);
}
}
}

		Matrix ( int, int) ;
		Matrix ( double, int, int ) ;
		// file constructor and function
		Matrix ( const std::string & ) ;
		void load ( const std::string & ) ;
		// random Matrix
		Matrix ( int, int, char ) ;
		// copy constructors
		Matrix ( const Matrix & ) ;
		// agumentated matrix
		Matrix ( const Matrix &, const Matrix & );
	    // deconstructor
        ~Matrix();
		//  equal operator
		Matrix & operator= ( const Matrix & );
		// operators
		Matrix & operator+= ( Matrix & ) ;
		Matrix & operator-= ( Matrix & ) ;
		Matrix & operator*= ( Matrix & ) ;
		//Matrix operator - ( Matrix  n);
       // Matrix& operator + ( Matrix & );
        Matrix& operator * ( Matrix & );
        // friend operators
        friend ostream& operator<<(ostream & out,const Matrix&);
		// bool operators
		bool operator== ( Matrix & ) const;
		//bool operator> (Matrix & mat);
		Matrix operator> (int x);
		// manipulation functions
        Matrix vecmat (Matrix, int);
		Matrix matvec (Matrix, int);
		Matrix matmat (Matrix);
		Matrix subvec (Matrix, int);
		Matrix diagtimes(Matrix, int);
		Matrix nummat(double num);
  	    void print();
  	    double DistanceFormula();
  	    double inner(Matrix &, int);
  	    double minus(int, int);
  	    void isRight(Matrix x,Matrix b);
		// output
		void output (  ) const ;
		// get/put functions
		int getR (  ) const { return r ; }
		int getC (  ) const { return c ; }
		double getOne ( int row, int col ) const { return m [ row ] [ col ] ; }
		void putOne ( double d, int row, int col ) { m [ row ] [ col ] = d ; }
		void putColumn ( const Matrix, int j  ) ;
		void putMatrix ( const Matrix ) ;
        // Inverse matrix
        Matrix Inverse();
        // transposed matrix
        Matrix Transposed();
        // row-echolon form
        void REF ( Matrix & );
        // reduced row-echelon form
        void RREF ( Matrix & );
        // matrix utility
        void lu (Matrix & );
        Matrix BackSolve (Matrix & );
        Matrix ForwardSolve (Matrix &);
        Matrix BackwardSolve ( Matrix & );
        Matrix FrontSolve (Matrix & );
        void jacobi ( Matrix & );
        void gs ( Matrix & );
        void sor ( Matrix & );
        void pm ( Matrix &,double );
        void sur ( Matrix & );
        double norm (  );
        void ods ( Matrix & );
        void rotatem (  );
	/*
   Derivation from the fortran version of CONREC by Paul Bourke
   d               ! matrix of data to contour
   ilb,iub,jlb,jub ! index bounds of data matrix
   x               ! data matrix column coordinates
   y               ! data matrix row coordinates
   nc              ! number of contour levels
   z               ! contour levels in increasing order
*/
#define SCALE 5
#define NCONTOUR 1
void Contour(){

double contours[NCONTOUR];
double **data;
double zmin=1e32,zmax=-1e32;
   int i,j,ii,jj;

data = new double*[SCALE*NX];
FOR(i,SCALE*NX){
data[i] = new double[SCALE*NY];
}
FOR(i,r){
FOR(j,c){
FOR(ii,SCALE){
FOR(jj,SCALE){
            data[SCALE*(int)i+ii][SCALE*(int)j+jj] = m[i][j];
}
}
}
}
int sum,n;
   for (i=0;i<SCALE*NX;i++) {
      for (j=0;j<NY*SCALE;j++) {
         n = 0;
         sum = 0;
         for (ii=-4;ii<=4;ii++) {
            for (jj=-4;jj<=4;jj++) {
               if (i + ii < 0 || i + ii >= SCALE*NX)
                  continue;
               if (j + jj < 0 || j + jj >= SCALE*NY)
                  continue;
               sum += data[i+ii][j+jj];
               n++;
            }
         }
         if (n <= 0) {
            fprintf(stderr,"No cells averaged, this shouldn't happen!\n");
            exit(-1);
         }
         data[i][j] = sum / n;
      }
   }

double *xaxis = new double[SCALE*NX];
double *yaxis = new double[SCALE*NY];

	image1 = Create_Bitmap(NX*SCALE,NY*SCALE);

FOR(i,SCALE*NX){
      xaxis[i] = i;
}
FOR(i,SCALE*NY){
      yaxis[i] = i;
}
zmax = imax();
zmin = imin();
   contours[0] = (zmax + zmin) / 2;
   //~ contours[2] = zmin + (zmax - zmin) / 8;
//~ //   contours[0] = 92;
   //~ contours[1] = zmin + (zmax - zmin) / 4;

   //~ contours[3] = zmax - (zmax - zmin) / 4;
   //~ contours[4] = zmax - (zmax - zmin) / 8;
//~ cout<<zmax<<" "<<zmin<<endl;
   CONREC(data,0,SCALE*NX-1,0,SCALE*NY-1,
      xaxis,yaxis,NCONTOUR,contours);
FILE* fptr;
   if ((fptr = fopen("image.tga","w")) == NULL) {
      fprintf(stderr,"Failed to open output image\n");
      exit(-1);
   }
   Write_Bitmap(fptr,image1,SCALE*NX,SCALE*NY,12);
   fclose(fptr);
cout<<nump<<endl;
	//exit(0);
}

 void ConrecLine(double x1,double y1,double x2,double y2,double z){
    BITMAP4 black = {255,255,255,255};
 //cout<<x1<<" "<<x2<<" "<<y1<<" "<<y2<<endl;
  Draw_Line(image1,SCALE*NX,SCALE*NY,(int)x1,(int)y1,(int)x2,(int)y2,black);
nump++;
 }

void CONREC(double **d,int ilb,int iub,int jlb,int jub,
   double *x,double *y,int nc,double *z)
{
#define xsect(p1,p2) (h[p2]*xh[p1]-h[p1]*xh[p2])/(h[p2]-h[p1])
#define ysect(p1,p2) (h[p2]*yh[p1]-h[p1]*yh[p2])/(h[p2]-h[p1])
//cout<<"blah";
   int m1,m2,m3,case_value;
   double dmin,dmax,x1,x2,y1,y2;
   int i,j,k,m;
   double h[5];
   int sh[5];
   double xh[5],yh[5];
   int im[4] = {0,1,1,0},jm[4]={0,0,1,1};
   int castab[3][3][3] = {
     { {0,0,8},{0,2,5},{7,6,9} },
     { {0,3,4},{1,3,1},{4,3,0} },
     { {9,6,7},{5,2,0},{8,0,0} }
   };
   double temp1,temp2;

   for (j=(jub-1);j>=jlb;j--) {
      for (i=ilb;i<=iub-1;i++) {
         temp1 = MIN(d[i][j],d[i][j+1]);
         temp2 = MIN(d[i+1][j],d[i+1][j+1]);
         dmin  = MIN(temp1,temp2);
         temp1 = MAX(d[i][j],d[i][j+1]);
         temp2 = MAX(d[i+1][j],d[i+1][j+1]);
         dmax  = MAX(temp1,temp2);
         //if (dmax < z[0])
         if (!(dmax < z[0] || dmin > z[nc-1]))
		 {//cout<<"blah";
		 //cout<< dmax<<endl;
            continue;
		 }
	 //cout<<"blah";
         for (k=0;k<nc;k++) {
            if (!(z[k] < dmin || z[k] > dmax)){
	    	// cout<<"blah";
               continue;
	    }
            for (m=4;m>=0;m--) {
               if (m > 0) {
                  h[m]  = d[i+im[m-1]][j+jm[m-1]]-z[k];
                  xh[m] = x[i+im[m-1]];
                  yh[m] = y[j+jm[m-1]];
               } else {
                  h[0]  = 0.25 * (h[1]+h[2]+h[3]+h[4]);
                  xh[0] = 0.50 * (x[i]+x[i+1]);
                  yh[0] = 0.50 * (y[j]+y[j+1]);
               }
               if (h[m] > 0.0)
                  sh[m] = 1;
               else if (h[m] < 0.0)
                  sh[m] = -1;
               else
                  sh[m] = 0;
            }

            /*
               Note: at this stage the relative heights of the corners and the
               centre are in the h array, and the corresponding coordinates are
               in the xh and yh arrays. The centre of the box is indexed by 0
               and the 4 corners by 1 to 4 as shown below.
               Each triangle is then indexed by the parameter m, and the 3
               vertices of each triangle are indexed by parameters m1,m2,and m3.
               It is assumed that the centre of the box is always vertex 2
               though this isimportant only when all 3 vertices lie exactly on
               the same contour level, in which case only the side of the box
               is drawn.
                  vertex 4 +-------------------+ vertex 3
                           | \               / |
                           |   \    m-3    /   |
                           |     \       /     |
                           |       \   /       |
                           |  m=2    X   m=2   |       the centre is vertex 0
                           |       /   \       |
                           |     /       \     |
                           |   /    m=1    \   |
                           | /               \ |
                  vertex 1 +-------------------+ vertex 2
            */
            /* Scan each triangle in the box */
	  //  cout<<"blah";
            for (m=1;m<=4;m++) {
               m1 = m;
               m2 = 0;
               if (m != 4)
                  m3 = m + 1;
               else
                  m3 = 1;
               if ((case_value = castab[sh[m1]+1][sh[m2]+1][sh[m3]+1]) == 0)
		       {
		   //    cout<<"blah";
                  continue;
		       }
               switch (case_value) {
               case 1: /* Line between vertices 1 and 2 */
                   x1 = xh[m1];
                   y1 = yh[m1];
                   x2 = xh[m2];
                   y2 = yh[m2];
                   break;
               case 2: /* Line between vertices 2 and 3 */
                   x1 = xh[m2];
                   y1 = yh[m2];
                   x2 = xh[m3];
                   y2 = yh[m3];
                   break;
               case 3: /* Line between vertices 3 and 1 */
                   x1 = xh[m3];
                   y1 = yh[m3];
                   x2 = xh[m1];
                   y2 = yh[m1];
                   break;
               case 4: /* Line between vertex 1 and side 2-3 */
                   x1 = xh[m1];
                   y1 = yh[m1];
                   x2 = xsect(m2,m3);
                   y2 = ysect(m2,m3);
                   break;
               case 5: /* Line between vertex 2 and side 3-1 */
                   x1 = xh[m2];
                   y1 = yh[m2];
                   x2 = xsect(m3,m1);
                   y2 = ysect(m3,m1);
                   break;
               case 6: /* Line between vertex 3 and side 1-2 */
                   x1 = xh[m1];
                   y1 = yh[m1];
                   x2 = xsect(m1,m2);
                   y2 = ysect(m1,m2);
                   break;
               case 7: /* Line between sides 1-2 and 2-3 */
                   x1 = xsect(m1,m2);
                   y1 = ysect(m1,m2);
                   x2 = xsect(m2,m3);
                   y2 = ysect(m2,m3);
                   break;
               case 8: /* Line between sides 2-3 and 3-1 */
                   x1 = xsect(m2,m3);
                   y1 = ysect(m2,m3);
                   x2 = xsect(m3,m1);
                   y2 = ysect(m3,m1);
                   break;
               case 9: /* Line between sides 3-1 and 1-2 */
                   x1 = xsect(m3,m1);
                   y1 = ysect(m3,m1);
                   x2 = xsect(m1,m2);
                   y2 = ysect(m1,m2);
                   break;
               default:
                   break;
               }
               /* Finally draw the line */
               ConrecLine(x1,y1,x2,y2,z[k]);
            } /* m */
         } /* k - contour */
      } /* i */
   } /* j */
}
};


#endif

