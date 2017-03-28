#include <stdio.h>
#include <math.h>


double L(double x, double a, double b, double c, double alp) {
	double ex = x*x*x - (a*b*c*x*x) + ((a*b+b*c+c*a)*x) - (a+b+c);
	double res = (double)ex/alp;
  	printf("x %lf , a %lf , b %lf, c %lf, alp %lf , L %lf \n", x , a, b , c , alp , res);
	return res;
}

double del_L(double x, double a, double b, double c, double alp){
	double ex = 3*x*x - (2*a*b*c*x) + (a*b+b*c+c*a);
	double res = (double)ex/alp;
	printf("DEL   x %lf , a %lf , b %lf, c %lf, alp %lf , L %lf \n", x , a, b , c , alp , res);
	return res;
}

double sub_poly(double y, double y0, double y1, double y2, double y3, double f[4][4], int a, int b ) {
	//printf("Sub_poly_open\n");
	double y_alp0=(y0-y1)*(y0-y2)*(y0-y3);
	double y_alp1=(y1-y0)*(y1-y2)*(y1-y3);
	double y_alp2=(y2-y0)*(y2-y1)*(y2-y3);
	double y_alp3=(y3-y0)*(y3-y1)*(y3-y2);

	double ex0 = L(y, y1, y2, y3, y_alp0)*f[a][b];
	double ex1 = L(y, y0, y2, y3, y_alp1)*f[a][b+1];
	double ex2 = L(y, y0, y1, y3, y_alp2)*f[a][b+2];
	double ex3 = L(y, y0, y1, y2, y_alp3)*f[a][b+3];
	//printf("Sub_poly_close\n");
	return ex0+ex1+ex2+ex3;
}

double poly(double x, double y, double x0, double x1, double x2, double x3, double y0, double y1, double y2, double y3, double f[4][4]){
	//printf("poly_open\n");
	double x_alp0=(x0-x1)*(x0-x2)*(x0-x3);
	double x_alp1=(x1-x0)*(x1-x2)*(x1-x3);
	double x_alp2=(x2-x0)*(x2-x1)*(x2-x3);
	double x_alp3=(x3-x0)*(x3-x1)*(x3-x2);

	double P_x_y_0 = L(x, x1, x2, x3, x_alp0) * sub_poly(y, y0, y1, y2, y3, f, 0, 0);
	double P_x_y_1 = L(x, x0, x2, x3, x_alp1) * sub_poly(y, y0, y1, y2, y3, f, 1, 0);
	double P_x_y_2 = L(x, x0, x1, x3, x_alp2) * sub_poly(y, y0, y1, y2, y3, f, 2, 0);
	double P_x_y_3 = L(x, x0, x1, x2, x_alp3) * sub_poly(y, y0, y1, y2, y3, f, 3, 0);
	double P_x_y = P_x_y_0 + P_x_y_1 + P_x_y_2 + P_x_y_3;
	//printf("poly_close\n");
	return P_x_y;
}

double sub_del_poly_y(double y, double y0, double y1, double y2, double y3, double f[4][4], int a, int b) {
	//printf("sub_del_poly_y_open\n");
	double y_alp0=(y0-y1)*(y0-y2)*(y0-y3);
	double y_alp1=(y1-y0)*(y1-y2)*(y1-y3);
	double y_alp2=(y2-y0)*(y2-y1)*(y2-y3);
	double y_alp3=(y3-y0)*(y3-y1)*(y3-y2);

	double ex0 = del_L(y, y1, y2, y3, y_alp0)*f[a][b];
	double ex1 = del_L(y, y0, y2, y3, y_alp1)*f[a][b+1];
	double ex2 = del_L(y, y0, y1, y3, y_alp2)*f[a][b+2];
	double ex3 = del_L(y, y0, y1, y2, y_alp3)*f[a][b+3];
	//printf("sub_del_poly_y_close\n");

	return ex0+ex1+ex2+ex3;
}

double del_poly_x(double x, double y, double x0, double x1, double x2, double x3, double y0, double y1, double y2, double y3 , double f[4][4]){
	//printf("del_poly_x_open\n");

	double x_alp0=(x0-x1)*(x0-x2)*(x0-x3);
	double x_alp1=(x1-x0)*(x1-x2)*(x1-x3);
	double x_alp2=(x2-x0)*(x2-x1)*(x2-x3);
	double x_alp3=(x3-x0)*(x3-x1)*(x3-x2);

	double P_x_y_0 = del_L(x, x1, x2, x3, x_alp0) * sub_poly(y, y0, y1, y2, y3, f, 0, 0);
	double P_x_y_1 = del_L(x, x0, x2, x3, x_alp1) * sub_poly(y, y0, y1, y2, y3, f, 1, 0);
	double P_x_y_2 = del_L(x, x0, x1, x3, x_alp2) * sub_poly(y, y0, y1, y2, y3, f, 2, 0);
	double P_x_y_3 = del_L(x, x0, x1, x2, x_alp3) * sub_poly(y, y0, y1, y2, y3, f, 3, 0);
	double P_x_y = P_x_y_0 + P_x_y_1 + P_x_y_2 + P_x_y_3;
	//printf("del_poly_x_close\n");
	return P_x_y;
}

double del_poly_y(double x, double y, double x0, double x1, double x2, double x3, double y0, double y1, double y2, double y3, double f[4][4]){
	//printf("del_poly_y_open\n");
	double x_alp0=(x0-x1)*(x0-x2)*(x0-x3);
	double x_alp1=(x1-x0)*(x1-x2)*(x1-x3);
	double x_alp2=(x2-x0)*(x2-x1)*(x2-x3);
	double x_alp3=(x3-x0)*(x3-x1)*(x3-x2);

	double P_x_y_0 = L(x, x1, x2, x3, x_alp0) * sub_del_poly_y(y, y0, y1, y2, y3, f, 0, 0);
	double P_x_y_1 = L(x, x0, x2, x3, x_alp1) * sub_del_poly_y(y, y0, y1, y2, y3, f, 1, 0);
	double P_x_y_2 = L(x, x0, x1, x3, x_alp2) * sub_del_poly_y(y, y0, y1, y2, y3, f, 2, 0);
	double P_x_y_3 = L(x, x0, x1, x2, x_alp3) * sub_del_poly_y(y, y0, y1, y2, y3, f, 3, 0);
	double P_x_y = P_x_y_0 + P_x_y_1 + P_x_y_2 + P_x_y_3;
	//printf("del_poly_y_close\n");
	return P_x_y;
}

void multiplyMatrices(double firstMatrix[2][2], double secondMatrix[2][1], double mult[2][1], int rowFirst, int columnFirst, int rowSecond, int columnSecond)
{
	int i, j, k;

	// Initializing elements of matrix mult to 0.
	for(i = 0; i < rowFirst; ++i)
	{
		for(j = 0; j < columnSecond; ++j)
		{
			mult[i][j] = 0;
		}
	}

	// Multiplying matrix firstMatrix and secondMatrix and storing in array mult.
	for(i = 0; i < rowFirst; ++i)
	{
		for(j = 0; j < columnSecond; ++j)
		{
			for(k=0; k<columnFirst; ++k)
			{
				mult[i][j] += firstMatrix[i][k] * secondMatrix[k][j];
			}
		}
	}

	//return mult;
}


void display(double mult[][2], int rowFirst, int columnSecond)
{
	int i, j;
	printf("\nOutput Matrix:\n");
	for(i = 0; i < rowFirst; ++i)
	{
		for(j = 0; j < columnSecond; ++j)
		{
			printf("%lf  ", mult[i][j]);
			if(j == columnSecond - 1)
				printf("\n\n");
		}
	}
}


void inv_2_2(double A[2][2], double inv[2][2]){
    double a = A[0][0];
    double b = A[0][1];
    double c = A[1][0];
    double d = A[1][1];

    double det = a*d - b*c ; 

    inv[0][0] = d/det;
    inv[0][1] = -b/det;
    inv[1][0] = -c/det;
    inv[1][1] = a/det;

}

void display_inv(double mult[][2], int rowFirst, int columnSecond)
{
    int i, j;
    printf("\nOutput Matrix:\n");
    for(i = 0; i < rowFirst; ++i)
    {
        for(j = 0; j < columnSecond; ++j)
        {
            printf("%lf  ", mult[i][j]);
            if(j == columnSecond - 1)
                printf("\n\n");
        }
    }
}


int newton_raphson(double result[2], double x0, double x1, double x2, double x3, double y0, double y1, double y2, double y3, double f[4][4] , double g[4][4]){

	int max_iterations = 100000 ; 
	double TOL = 0.001; 

	double X0[2][1];
	double X1[2][1];

	X0[0][0] = 1;
	X0[1][0] = 2;
	int i;
	for(i = 0 ; i<=max_iterations ; i++ ){

		double J[2][2];

		J[0][0] = del_poly_x(X0[0][0] , X0[1][0] , x0 , x1 , x2 , x3 , y0 , y1 , y2 , y3 , f );  // at ( U xy[0] , xy[1]) 
		J[0][1] = del_poly_y(X0[0][0] , X0[1][0] , x0 , x1 , x2 , x3 , y0 , y1 , y2 , y3 , f ); // at ( U xy[0] , xy[1])
		J[1][0] = del_poly_x(X0[0][0] , X0[1][0] , x0 , x1 , x2 , x3 , y0 , y1 , y2 , y3 , g );// at ( V xy[0] , xy[1])
		J[1][1] = del_poly_y(X0[0][0] , X0[1][0] , x0 , x1 , x2 , x3 , y0 , y1 , y2 , y3 , g );	// at ( V xy[0] , xy[1])

		double J_inv[2][2];
		inv_2_2(J, J_inv);

		double fun[2][1];

		fun[0][0] = poly(X0[0][0] , X0[1][0] , x0 , x1 , x2 , x3 , y0 , y1 , y2 , y3 , f );  // at U  ( xy[0][0] , xy[1][0])
		fun[1][0] = poly(X0[0][0] , X0[1][0] , x0 , x1 , x2 , x3 , y0 , y1 , y2 , y3 , g ); // at V  ( xy[0][0] , xy[1][0])

		//printf("Poly for F %lf \n", fun[0][0]);
		//printf("Poly for G %lf  \n", fun[1][0]);

		double res[2][1];

		multiplyMatrices(J_inv , fun , res , 2 , 2 , 2 , 2);
		printf("X1 : %f , %f       Res : %f ,  %f\n" , X1[0][0] , X1[1][0] , res[0][0] , res[1][0]);

		X1[0][0] = X0[0][0] - res[0][0];
		X1[1][0] = X0[1][0] - res[1][0];

		if(fabs(X1[0][0] - X0[0][0]) < TOL && fabs(X1[1][0] - X0[1][0]) < TOL  ){
			result[0] = X1[0][0];
			result[1] = X1[1][0];
			return 1; 
		}
		else{
			X0[0][0] = X1[0][0];
			X0[1][0] = X1[1][0];
		}


	}
	return 0;

}




int main(){
	

	double X[4] = { 1, 3 , 5 , 7 };
	double Y[4] = { 2 , 4 , 6, 8 };
	double U[4][4] = {
		{0.5 , 0.4 , 0.7, 1.25 },
		{0.3 , 0.45 , 0.5, 0.15 },
		{1.5 , 0.6 ,0.35 , 0.2 },
		{0.95 , 0.55 , 0.4 , 0.65 }
	};

	double V[4][4] = {
		{0.7 , 0.4 , 0.45, 0.25 },
		{0.15 , 1.0 , 0.9, 0.10 },
		{0.35 , 0.60 , 0.70 , 0.83 },
		{0.5 , 0.80 , 0.75, 0.95 }
	};

	double result[2];

	int response =  newton_raphson(result , X[0] , X[1] , X[2] , X[3] , Y[0] ,Y[1] , Y[2] , Y[3] , U , V) ; 
	if(response){
		printf("x = %lf  ", result[0]);	
		printf("y = %lf  ", result[1]);
	}

	else{
		printf("Max iterations reached");
	}
		
}