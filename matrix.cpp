

#include "matrix.h"

// identity constructor

Matrix::Matrix ( int ro, int co )
{
  	r = ro;
	c = co;
	m = new double * [ ro ];

	for ( int i = 0; i != ro; i++ )
		m [ i ] = new double [ co ];
	for ( int i = 0; i != ro; i++ )
	{
		for(int j = 0; j != co; j++)
		{
			if ( i == j )
				m [ i ] [ j ] = 1.0;
			else
				m [ i ] [ j ] = 0.0;
		}
	}
}

// constructors that makes a matrix of nums

Matrix::Matrix ( double num, int ro, int co )
{
	int i;
	r = ro;
	c = co;
	m = new double * [ ro ];

	for ( i = 0; i != ro; i++ )
		m [ i ] = new double [ co ];
	
    for ( i = 0; i != ro; i++ )
	{
		for(int j = 0; j != co; j++)
		{
				m [ i ] [ j ] = num;
		}
	}
}

// file constructor and function

Matrix::Matrix ( const string& filename )
{
	ifstream input;
	input.open ( filename.c_str (  ), ios::in );
    input >> r >> c;
	m = new double * [ r ];
	for ( int i = 0; i < r; i++ )
	{
		m [ i ] = new double [ c ];
    }
	for ( int i = 0; i < r; i++ )
	{
		for ( int j = 0; j < c; j++ )
		{
			input >> m [ i ] [ j ];
		}
	}
	  
   input.close();
}

void Matrix::load ( const string & fname )
{

	ifstream fin;
	fin.open ( fname.c_str(  ), ios::in );
	fin >> r >> c;
	
    m=new double * [ r ];
	
    for ( int i=0; i < r; i++ )
		m[ i ] = new double [ c ];

	for ( int i = 0; i < r; i++ )
	{	for ( int j = 0; j < c; j++)
		{	fin >> m [ i ] [ j ];}
	}
    fin . close (  );
}

//random and diagnol constructors

Matrix::Matrix ( int ro, int co, char o )
{
	r = ro;
	c = co;
	double sum;
	int Rvalue, coin;
	
	m = new double *[r];
	for(int i=0; i<r; i++)
	{
		m[i]=new double [c];
	}
	if(o=='r')
	{

		for(int i=0; i<r; i++)
		{
			for(int j=0; j<c; j++)
			{
				Rvalue = rand()%100000;
				coin=rand()%2;
				if(coin==0)
				{
					Rvalue=-1*Rvalue;
				}	
				m[i][j]=(double)Rvalue/100000;
			}
		}
	}
	else if(o=='d')
	{
		for ( int i=0; i < r; i++ )
		{
			for(int j=0; j<c; j++)
            {
				Rvalue = rand();//100000
				coin=rand()%2;
				if(coin==0)
				{
					Rvalue=-1*Rvalue;
				}	
				m[i][j]=(double)Rvalue/1;
            }
        }
		for ( int i=0; i < r; i++ )
        {
			sum=0.0;
			for(int j=0; j<c; j++)
			{
				if(i==j)
				{
					sum=sum+fabs(m[i][j]);
				}
          }
			m[i][i]=sum+1.0;   
        }
	}
}

// copy constructor

Matrix::Matrix ( const Matrix& mat )
{
	r = mat.r;
	c = mat.c;
    m = mat.m;
}

// deconstructor

Matrix::~Matrix()
{
}

// equals operator

Matrix& Matrix::operator=( const Matrix& mat )
{
    if(this != &mat)
    {
	        r = mat.r;
	        c = mat.c;
            m = mat.m;
    }
	return *this;
}

// operators // changes a b, use this +,-,*= 

Matrix& Matrix::operator+= ( Matrix & adder )
{
    if ( r != adder.r || c != adder.c )
    {
		cout << "matricies are not equivilant size." << endl;
		output();
        cout << endl;
		adder.output();
		return *this;
    }
	  
      for ( int i = 0; i != r; i++ )
      {
          for ( int j = 0; j != c; j++ )
          {
              m [ i ] [ j ] += adder.m [ i ] [ j ];
          }
      }
	  return *this;
}
 
Matrix& Matrix::operator*= ( Matrix& multiplier )
{
      if ( c != multiplier.r )
      {
           cout << "columns in first Matrix does not equal that in second Matrix" << endl ;
		   output (  ) ;
		   cout << endl ;
		   return *this ;
	  }
	  else
	  {
		  double temp = 0 ;
		  for ( int i = 0; i < r; i++ )
		  {
				  for ( int j = 0; j < multiplier.c; j++ )
				  {
						  for ( int k = 0 ; k < ( c+multiplier.r )/2 ; k++ )
						  {
								  temp += m [ i ] [ k ] * multiplier.m [ k ] [ j ] ;
						  }
						  m [ i ] [ j ] = temp ;
						  temp = 0 ;
				  }
		  }
	  }
      return *this ;
}
 
Matrix& Matrix::operator-= ( Matrix & adder )
{
      if ( r != adder.r || c != adder.c )
      {
           cout << "matricies are not equivilant size." << endl;
           output();
		   cout << endl;
		   adder.output();
		   return *this;
      }
	  
      for ( int i = 0; i != r; i++ )
      {
          for ( int j = 0; j != c; j++ )
          {
              m [ i ] [ j ] -= adder.m [ i ] [ j ];
          }
      }
	  return *this;
}

//~ Matrix& Matrix::operator+ ( Matrix & adder )
//~ {
    //~ if ( r != adder.r || c != adder.c )
    //~ {
		//~ cout << "matricies are not equivilant size." << endl;
		//~ output();
        //~ cout << endl;
		//~ adder.output();
		//~ return *this;
    //~ }
	  
    //~ for ( int i = 0; i != r; i++ )
    //~ {
          //~ for ( int j = 0; j != c; j++ )
          //~ {
              //~ m [ i ] [ j ] = m [ i ] [ j ] + adder.m [ i ] [ j ];
          //~ }
    //~ }
	//~ return *this;
//~ }
 
Matrix& Matrix::operator* ( Matrix& multiplier )
{
        double**mat;
        mat = new double*[multiplier.r];
        for(int i = 0; i < r; i++)
                mat[i] = new double[c];
      if ( c != multiplier.r )
      {
           cout << "columns in first Matrix does not equal that in second Matrix" << endl ;
		   output (  ) ;
		   cout << endl ;
		   multiplier.output();
		   return *this ;
	  }
	  else
	  {
		  double temp = 0 ;
		  for ( int i = 0; i < r; i++ )
		  {
				  for ( int j = 0; j < multiplier.c; j++ )
				  {
						  for ( int k = 0 ; k < r ; k++ )
						  {
								  temp += m [ i ] [ k ] * multiplier.m [ k ] [ j ] ;
						  }
						  mat [ i ] [ j ] = temp ;
						  temp = 0 ;
				  }
		  }
	  }
	  m = mat;
	  c = multiplier.c;
      return *this ;
}

//~ Matrix Matrix::operator- ( Matrix  adder )
//~ {
      //~ if ( r != adder.r || c != adder.c )
      //~ {
           //~ cout << "matricies are not equivilant size." << endl;
           //~ output();
		   //~ cout << endl;
		   //~ adder.output();
		   //~ return *this;
      //~ }
	  
      //~ for ( int i = 0; i != r; i++ )
      //~ {
          //~ for ( int j = 0; j != c; j++ )
          //~ {
              //~ m [ i ] [ j ] = m [i][j] - adder.m [ i ] [ j ];
          //~ }
      //~ }
	  //~ return *this;
//~ }

// friend functions

ostream& operator<<(ostream & out,const Matrix & mat)
{
         mat.output();
         return out;
}

// bool functions

bool Matrix::operator== ( Matrix & mat ) const
{
	if ( r != mat.r && c != mat.c )
		return false;
	for ( int i = 0; i < r; i++ )
	{
		for ( int j = 0; j < c; j++)
		{
			if ( fabs(m [ i ] [ j ] - mat.m [ i ] [ j ]) > .0001)
				return false;
		}
	}
	return true;
}

//~ bool Matrix::operator>(Matrix & mat)
//~ {
	//~ for (int i = 0; i < mat.getR (); i++)
	//~ {
		//~ for (int j = 0; j < mat.getC (); j++)
		//~ {
			//~ if (m[i][j]>mat.m[i][j])
				//~ return false;
		//~ }
	//~ }
	//~ return false;
//~ }

void Matrix::print(){
for (int i = 0; i < r; i++)
	{
		for (int j = 0; j < c; j++)
		{
			cout<<m[i][j] <<" ";
		}
	cout<<endl;
	}
}
Matrix Matrix::operator>(int w)
{
	Matrix tmp(r,c);
	for (int i = 0; i < r; i++)
	{
		for (int j = 0; j < c; j++)
		{
			if (m[i][j] > 0)
				tmp.m[i][j] = 1;
			else tmp.m[i][j] = 0;
		}
	}
	return tmp;
}

// manipulation functions

Matrix Matrix::diagtimes(Matrix x,int n)
{
       Matrix sol(n,1);
       
       for(int i=0; i<n; i++)
       {
               sol.m[i][0]=m[i][0]*x.m[i][0];
       }
       
       return sol;

}


Matrix Matrix::vecmat(Matrix a, int n)
{
       cout << c;
       Matrix sol(r, 1);
double sum;
    for(int i=0; i<r; i++)
    {
        sum=0.0;
        for(int j=0; j<c; j++)
        {
            sum = sum+m[j][0]*a.m[i][j];
        }
        sol.m[i][0]=sum;
    }
    return sol;
}
Matrix Matrix::matvec(Matrix x, int n)
{
    double sum;
    
    Matrix y (n,1);
    
    for(int i=0; i<n; i++)
    {
        sum=0.0;
        for(int c=0; c<n; c++)
        {
            sum = sum+m[i][c]*x.m[c][0];
        }
        y.m[i][0]=sum;
    }
    
    return y;
}

Matrix Matrix::nummat(double num)
{
       double sum;
       
       Matrix y (r,1);
       
       for(int i = 0; i<r; i++)
       {
        y.m[i][0] = num*m[i][0];
       }
       
       return y;
}

double Matrix::inner(Matrix & mat, int n)
{
       double dot=0;
        for (int i = 0; i < n; i++)
            dot = dot + m[i][0]*mat.m[i][0];
        return dot;
}

Matrix Matrix::matmat(Matrix x)
{
    double sum;
    
    Matrix y (r, x.c);
    
    for ( int i = 0; i != r; i++ )
	  {
		  for ( int j = 0; j != r; j++ )
		  {
				  for ( int k = 0; k != c; k++ )
						  sum = sum + m [ i ] [ k ] * x.m [ k ] [ j ];
				  y.m [ i ] [ j ] = sum;
				  sum = 0;
		  }
	  }
    
    return y; 
}

double Matrix::DistanceFormula()
{
       double sum=0.0;
       
       for ( int i=1; i < r; i++ )
       {
               sum = sum+ pow(m[i][0],2.0);
       }
       
       sum=sqrt(sum);
       
       return sum;
}

Matrix Matrix::subvec(Matrix v, int n)
{
       Matrix sol(n,1);
       
       for(int i=0; i<n; i++)
       {
               sol.m[i][0]=m[i][0]-v.m[i][0];
       }
       
       return sol;
}

double Matrix::minus(int i, int j)
{
       return - m[i][j];
}

//get/put functions

void Matrix::putColumn(const Matrix mat, int j)
{
	for( int i = 0; i < r; i++)
	{
		m [ i ] [ j ]=mat.m[i][j];
	}
}

void Matrix::putMatrix(const Matrix mat)
{
     Matrix matr(mat.r,mat.c);
     
     for(int i = 0; i < mat.r; i++)
             for(int j = 0; j < mat.c; j++)
                     matr.m [i][j]=mat.m[i][j];
     m = matr.m;
     
     c++;
}

// output functions

void Matrix::output() const
{
	int i,j;
    //ofstream file;

    //file.open("",fstream::out);
	for ( i = 0; i < r; i++ )
	{
		for ( j = 0; j < c; j++ )
	    {
			cout << fixed << setw(10) << setprecision(2) << m [ i ] [ j ];
		}
		cout << endl;
    }
    //file.close();
}

// Inverse

Matrix Matrix::Inverse (  )
{
     Matrix mat(r,1);
     
     for(int i = 0; i < r; i++)
     {
             for(int j = 0; j<c;j++)
             {
                     if (i == j)
                        mat.m[j][0] = 1/m[i][j];
             }
     }
     
     return mat;
}

// Transposed

Matrix Matrix::Transposed (  )
{
     Matrix mat(c,r);
     
     for ( int i = 0; i < c; i++)
     {
         for ( int j  = 0; j < r; j++ )
         {
               mat.m [ i ] [ j ] = m [ j ] [ i ];
         }  
     }
     
     return mat;
}

// ref

void Matrix::REF(Matrix & b)
{
     Matrix x;
     x = FrontSolve(b);
     
     cout << "ref:"<< endl;
     
     cout << "without backsolve x:" << endl;
     
     cout << x;
     
     x = BackwardSolve(x);
     
     cout << "with backsolve x:"<<endl;
     
     x.output();
     
     isRight(x,b);
}

// rref

void Matrix::RREF(Matrix & b)
{
    Matrix x;
    
    x = FrontSolve(b);
    
    x = BackSolve(x);
    
    cout << " x:" << endl;
    
    x.output();
    
    isRight(x,b);
}


void Matrix::lu ( Matrix & b )
{
    double sum;
    cout << *this<<endl;
    cout << b<<endl;
    
    Matrix l(r,c);
    
    Matrix u(r,c);
    
    for ( int k = 0 ; k < r ; k++ )
	{
		l.m [ k ] [ 0 ] = m [ k ] [ 0 ];
	}

	for ( int l = 1 ; l < c ; l++ )
	{
        if (m[0][0]!= 0)
		   u.m[0][l] = m[0][l]/m[0][0];
	}
    double tempd;
    for ( int i = 1; i < c; i++) // lu facterization
	{
		for(int j = i; j < r; j++)//vary over the rows
		{
                sum = 0.0;
	            for( int k = 0; k < j ; k++)
	            {
			             sum +=l.m[j][k] * u.m[k][i];
                }
                tempd=m[j][i];
                l.m[j][i]=tempd-sum;
		}
		for(int j = i+1; j < c; j++)//generate various cols
		{
                sum = 0.0;
	            for( int k = 0; k < j-1; k++)
                {
	                 sum += l.m[i][k] * u.m[k][j];
                } 
                tempd=l.m[i][i];

                u.m[i][j] = (1.0/tempd)*(m[i][j]-sum);
		}
    }
    cout << l<<endl;
    cout << u<< endl;
    Matrix y = l.ForwardSolve(b);
    Matrix x = u.BackwardSolve(y); 
    cout << "x" << endl;
    cout << x;
    
    isRight(x,b);
}

Matrix Matrix::BackwardSolve(Matrix & b)
{
       Matrix y(b.r,b.c);
       double sum;
       for (int k = c-1; k >= 0; k--)
       {
                sum = 0.0;
                for (int i = k+1; i < c; i++)
                {
                     sum += (m[k][i]*y.m[i][0]);
                } 
                y.m[k][0] = ((b.m[k][0]-sum)/m[k][k]);
       }
       return y;
}
void Matrix::isRight(Matrix x, Matrix b)
{
    if(matvec(x,x.r)==b)
         cout << "Right!!!!!!!!" << endl;
    else
         cout << "Wrong!!!!!!!!" << endl;
}

Matrix Matrix::ForwardSolve(Matrix & b)
{
       Matrix y(b.r,b.c);
       double sum;
       for (int k = 0; k < c; k++)
       {
                sum = 0.0;
                for (int i = 0; i <= (k-1); i++)
                {
                     sum -= (m[k][i]*y.m[i][0]);
                } 
                y.m[k][0] = ((b.m[k][0]+sum)/m[k][k]);
                
       } 
       return y;
}

Matrix Matrix::FrontSolve(Matrix & b)
{
    double temp = 0;
    Matrix y;
    y = b;
    for(int k = 0; k < r-1; k++)
	{
		temp = m [ k ] [ k ];
        
        for ( int j = 0; j < c; j++)
		{
            if ( m [k][j] != 0 )
               m [ k ] [ j ] /= temp;
        }
        if ( y.m [ k ] [ 0 ] != 0 )
             y.m [ k ] [ 0 ] /= temp;
        
        for(int j = k+1; j < r; j++)
		{
				temp = m [ j ] [ k ];

				for(int i = 0; i < c; i++)
		        {
                        m [ k ] [ i ] *= temp;
                }
                y.m [ k ] [ 0 ] *= temp;
				
                for(int i = 0; i < c; i++)
		        {
                        m [ j ] [ i ] -= m [ k ] [ i ];
                }
                y.m [ j ] [ 0 ] -= y.m [ k ] [ 0 ];
                
                for ( int i = 0; i < c; i++ )
		        {
                        if ( m [ k ] [ i ] != 0 )
                           m [ k ] [ i ] /= temp;
                }
                if ( y.m [ k ] [ 0 ] != 0 )
                     y.m [ k ] [ 0 ] /= temp;
		}
    }
    
    temp = m [ r - 1 ] [ c - 1 ];

	for ( int i = 0; i < c; i++ )
    {
        if ( m [ r - 1 ] [ i ] != 0 )
             m [ r - 1 ] [ i ] /= temp;
    }
    if ( y.m [ r - 1 ] [ 0 ] != 0 )
         y.m [ r - 1 ] [ 0 ] /= temp;
    
    return y;
}

Matrix Matrix::BackSolve(Matrix & b)
{
    double temp = 0.0;
    Matrix y;
    y = b;
	for(int k = r - 1; k > 0; k--)
	{
		for(int j = k-1; j >= 0; j--)
		{
			temp = m [ j ] [ k ];

				for(int i = 0; i < c; i++)
		        {
                        m [ k ] [ i ] *= temp;
                }
                y.m [ k ] [ 0 ] *= temp;
				
                for(int i = 0; i < c; i++)
		        {
                        m [ j ] [ i ] -= m [ k ] [ i ];
                }
                y.m [ j ] [ 0 ] -= b.m [ k ] [ 0 ];
                
                for ( int i = 0; i < c; i++ )
		        {
                    if ( m [ k ] [ i ] != 0 )

                         m [ k ] [ i ] /= temp;
                }

                if ( y.m [ k ] [ 0 ] != 0 )

                     y.m [ k ] [ 0 ] /= temp;
		}
    }
    
    return y;
}

void Matrix::jacobi(Matrix & b)
{
       Matrix LU ( r, c ) ;
       double num = 0.0;
       for ( int i = 0 ; i < r ; i++ )
       {
               for ( int j = 0 ; j < c ; j++ )
               {
                       if ( i != j )
                       {
                           LU.m[i][j]=m[i][j];
                       }
                       else
                       {
                           LU.m[i][j]=0;
                       }
               }
       }

       cout << "what number do you want to try?hint try 0" << endl ;
       cin >> num ;
       
       Matrix xold ( num , b . r , b . c ) ;
       Matrix xnew ( num , b . r , b . c ) ;
       Matrix D ;
       Matrix temp (b.r,b.c);
       
       double normal;

       int k = 0;     
  
       D = Inverse();
	   
       do
       {
           xold = xnew;
           
           temp = LU.matvec ( xold , b . r );
           temp = b.subvec ( temp , b . r );
              
           xnew = D.diagtimes ( temp , b . r );

           temp = matvec ( xnew , b . r ) ;
           temp = b.subvec ( temp , b . r );
           
           normal = temp.norm();
           
           k++;
       }
       while ( normal > tol && k < 100 ) ;
       
       if ( k == 100 )
          cout << "Diverges" << endl;
       cout << "iterations: "<<k<<endl;
       cout << xnew;
       isRight(xnew, b);
}

void Matrix::gs(Matrix & b)
{
       Matrix LU ( r, c ) ;
       
       double num = 0.0;
       
       for ( int i = 0 ; i < r ; i++ )
       {
               for ( int j = 0 ; j < c ; j++ )
               {
                       if ( i != j )
                       {
                           LU.m[i][j]=m[i][j];
                       }
                       else
                       {
                           LU.m[i][j]=0;
                       }
               }
       }
       
       cout << "what number do you want to try?hint try 0" << endl ;
       cin >> num ;
       
       Matrix xold ( num , b . r , 1 ) ;
       Matrix xnew ( num , b . r , 1 ) ;
       Matrix D ;
       Matrix temp ( b . r , 1 );
       
       double normal, sum = 0;
       int k = 0;
       D = Inverse();
	   
       do
       {
           xold = xnew;

           for ( int i = 0; i < b.r; i++)
           {
                for( int j = 0; j < b.r; j++)   
                    sum = sum +LU.m[i][j] * xold .m[j][0];
                    sum = b.m[i][0] - sum;
                    xnew.m[i][0] = D.m[i][0] * sum;
                    sum = 0;
           }
           
           temp = matvec ( xnew , b . r ) ;
           temp = b.subvec ( temp , b . r );
           
           normal = temp.norm();
           
           k++;
       }
       while ( normal > tol && k < 100) ;
       
       if ( k == 100)
            cout << "Diverges" << endl;
       cout << "iterations: "<<k<<endl;
       cout << xnew;
       isRight(xnew,b);
}

void Matrix::sor(Matrix & b)
{
       Matrix LU ( r, c ) ;
       
       double num = 0.0;
       
       for ( int i = 0 ; i < r ; i++ )
       {
               for ( int j = 0 ; j < c ; j++ )
               {
                       if ( i != j )
                       {
                           LU.m[i][j]=m[i][j];
                       }
                       else
                       {
                           LU.m[i][j]=0;
                       }
               }
       }

       cout << "what number do you want to try?hint try 0" << endl ;
       cin >> num ;
       
       Matrix xold  ( num , b . r , 1 ) ;
       Matrix xnew  ( num , b . r , 1 ) ;
       Matrix D ;
       Matrix temp  ( num , b . r , 1 ) ;
    
       double normal, sum = 0, deltaB;
       int k = 0;
       D = Inverse();
	   
       do
       {
           xold = xnew;
           
           for ( int i = 0; i < b.r; i++)
           {    for( int j = 0; j < b.r; j++)   
                    sum = sum +LU.m[i][j] * xold .m[j][0];
                    sum = b.m[i][0] - sum;
                    xnew.m[i][0] = D.m[i][0] * sum;
                    sum = 0;
           }
           
           for(int i = 0; i < b.r; i++)
           {
                   deltaB = xnew.m[i][0] - xold.m[i][0];
                   deltaB = 1.2*deltaB;
                   xnew.m[i][0]=xold.m[i][0]+deltaB;
           }           
           
           temp = matvec ( xnew , b . r ) ;
           temp = b.subvec ( temp , b . r );
           
           normal = temp.norm();
           k++;
       }
       while ( normal > tol && k < 100 ) ;
       
       if ( k == 100 )
          cout << "Diverges" << endl;
       cout << "iterations: "<<k<<endl;
       cout << xnew;
       isRight(xnew, b);
}

void Matrix::sur(Matrix&b)
{
     Matrix LU ( r, c ) ;
       
       double num = 0.0;
       
       for ( int i = 0 ; i < r ; i++ )
       {
               for ( int j = 0 ; j < c ; j++ )
               {
                       if ( i != j )
                       {
                           LU.m[i][j]=m[i][j];
                       }
                       else
                       {
                           LU.m[i][j]=0;
                       }
               }
       }

       cout << "what number do you want to try?hint try 0" << endl ;
       cin >> num ;
       
       Matrix xold  ( num , b . r , 1 ) ;
       Matrix xnew  ( num , b . r , 1 ) ;
       Matrix D ;
       Matrix temp  ( num , b . r , 1 ) ;
    
       double normal, sum = 0, deltaB;
       int k = 0;
       D = Inverse();
	   
       do
       {
           xold = xnew;
           
           for ( int i = 0; i < b.r; i++)
           {    for( int j = 0; j < b.r; j++)   
                    sum = sum +LU.m[i][j] * xold .m[j][0];
                sum = b.m[i][0] - sum;
                xnew.m[i][0] = D.m[i][0] * sum;
                sum = 0;
           }
           
           for(int i = 0; i < b.r; i++)
           {
                   deltaB = xnew.m[i][0] - xold.m[i][0];
                   deltaB = .8*deltaB;
                   xnew.m[i][0]=xold.m[i][0]+deltaB;
           }           
           
           temp = matvec ( xnew , b . r ) ;
           temp = b.subvec ( temp , b . r );
           
           normal = temp.norm();
           k++;
       }
       while ( normal > tol && k < 100 ) ;
       
       if ( k == 100 )
          cout << "Diverges" << endl;
       cout << "iterations: "<<k<<endl;
       cout << xnew;
       isRight(xnew,b);
}

double Matrix::norm(  )
{
       double sum=DistanceFormula();
       
       return sum;
}

void Matrix::pm( Matrix & b, double m)
{
     Matrix	y(0.0,b.r,1); 
     Matrix r(0.0,b.r,1); 
     double s, d, t, k = 0;
     
     do 
     {
         s = b.norm ( );
         for (int i = 0; i < b.r; i++)
            y.m[i][0] = b.m[i][0]/s;
         matvec (b,b.r);
         d = b.inner(y,b.r)/y.inner(y,b.r);
         for (int i = 0; i < b.r; i++)
            r.m[i][0] = (d)*y.m[i][0] - b.m[i][0];
         t = r.norm ();
         k++;   
         cout << "dominant eigenvalue:" << d <<endl;
         cout << "dominant eigenvector:" << endl;
         for (int i = 0; i < b.r; i++)
                     cout << r.m[i][0] << endl;
     }
     while ((t > tol) && (k < m));
     
}

void Matrix::ods(Matrix & b)
{
     Matrix A;
     Matrix B;
     Matrix x;
     A = Transposed();
     B = A.matmat(*this);
     x = A.matvec(b,c);
     B.jacobi(x);
}

void Matrix::rotatem()
{
     int n = r, p, q, k = 1 ;
     Matrix R ( n , n ) ;
     Matrix P ( n , n ) ;
     Matrix A=*this;
     Matrix temp;
     double* eigenvalue = new double [ n ] ;
     double angle , magnitude , sum ;
     ofstream file;
     
     file.open("q4 solutions.txt", ios::out);
     
     for ( int i = 0; i < n-1; i++ )
     {
         for ( int j = i; j < n-1; j++ )
         {
                 if ( m[i+1][i+1] > m [ i ] [ j ] )
                 {   
                         p = i+1;
                         q = j+1;
                 }
         }
     }
     while ( (abs(A.m[p][q]) > 1) )
     {
           cout << A.m[p][q]<<endl;
           
           sum = (2*A.m[p][q]/(A.m[q][p]));
           angle = (1/2)*std::atan(sum);
           
           R.m [ p ] [ p ] = cos(angle);
           R.m [ q ] [ q ] = cos(angle);
           R.m [ p ] [ q ] = -sin(angle);
           R.m [ q ] [ p ] = -sin(angle);
           P = P.matmat(R);
           temp = A.matmat(R);
           A = R . Transposed (  ).matmat( temp );
           
           for (int i = 0; i < n-1; i++)
           {
               for ( int j = i; j < n-1; j++ )
               {
                     if ( A.m[i+1][j+1] > A.m [ i ] [ j ] )
                     {
                         p = i+1;
                         q = j+1;
                     }
               }
           }
           k++; cout<< k;
     }
     file << "eigenvalues:  ";
     for(int i = 0; i < n; i++)
     {
           file<< A.m[i][i] << ", ";
     }
     file << endl;
     file << "eigenvector:  " << endl;
     for(int j = 0; j < n; j++)
     {
             file <<j+1 <<". " << endl;
             for(int i = 0; i < n; i++)
             {
                     file << "["<<R.m[i][j] <<"]"<< endl;
             }
             file << endl;
     }
     file.close();
     cout << "eofrm" << endl;
}

