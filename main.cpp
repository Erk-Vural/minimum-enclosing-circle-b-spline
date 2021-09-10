// C program for find minimum-enclosing-circle
// and draw B-spline with given points using de boor's algorithm
// Author @Bugrahan Erk Vural
// 13/11/2020

#include <stdio.h>
#include <graphics.h>
#include <math.h>
#include <conio.h>

// Defining infinity
const double INF = 1e18;
// Coordinate plane size
int pln_sz = 500;

// Structure to represent a 2D point
struct Point {
	double X, Y;
};

// Structure to represent a 2D circle
struct Circle {
	Point C;
	double R;
};

// Check points both x,y values for biggest number
int find_biggest_pnt(struct Point P[ ], int counter){
    int biggest_pnt_X = P[0].X;
    int biggest_pnt_Y = P[0].Y;

    for(int i = 0; i < counter; i++){
        if(P[i].X > biggest_pnt_X){
            biggest_pnt_X = P[i].X;
        }
        if(P[i].Y > biggest_pnt_Y){
            biggest_pnt_Y = P[i].Y;
        }
    }
    if(biggest_pnt_X > biggest_pnt_Y){
        return biggest_pnt_X;
    }
    else if(biggest_pnt_Y > biggest_pnt_X){
        return biggest_pnt_Y;
    }else{
        return biggest_pnt_X;
    }
}

// find the fitting interval for displaying numbers
int set_interval(struct Point P[ ], int counter){
    int biggest = find_biggest_pnt(P, counter);

    if(biggest <= 25){
        return 1;
    }
    else if(biggest > 25 && biggest <= 100){
        return 4;
    }
    else if(biggest > 100 && biggest <= 250){
        return 10;
    }
    else{
        return 40;
    }
}

// Draw coordinate plane
void draw_coordinade_plane(int factor){

    int x = -25 * factor;
    int y = -25 * factor;
    char k[50];
    char s[50];

    line(0, pln_sz, pln_sz*2, pln_sz);// X-axis line
    line(pln_sz,0,pln_sz,pln_sz*2);// Y-axis line

    // X-axis numbers and lines
    for(int i=0; i<pln_sz*2; i=i+20){
        sprintf(s, "%d",x);
        if(i!=pln_sz){
            outtextxy(i-10, pln_sz+10, s);
        }
        line(i, pln_sz, i, pln_sz+5);
        x=x+factor;
    }
    // Y-axis numbers and lines
    for(int j=0; j<pln_sz*2; j=j+20){
        sprintf(k, "%d",-y);
        if(j!=pln_sz){
            outtextxy(pln_sz+10, j-10, k);
        }
        line(pln_sz, j, pln_sz+5, j);
        y=y+factor;
    }
}

// Show points on screen
void show_points(struct Point P[ ], int counter, int factor){

    while(counter >= 0){
        int point_x = (pln_sz + (P[counter].X * 20)/factor);
        int point_y = (pln_sz - (P[counter].Y * 20)/factor);
        char coordinate_text[20];
        sprintf(coordinate_text,"{%d, %d}",
            int(P[counter].X),int(P[counter].Y));

        circle(point_x, point_y,3);
        outtextxy(point_x + 5, point_y + 5,coordinate_text);
        counter--;
    }
}

// MEC Functions

// Function to return the euclidean distance
// between two points
double cal_dist(const Point& a, const Point& b)
{
	return sqrt(pow(a.X - b.X, 2) + pow(a.Y - b.Y, 2));
}

// Function to check whether a point lies inside
// or on the boundaries of the circle
bool is_pnt_inside(const Circle& c, const Point& p)
{
	return cal_dist(c.C, p) <= c.R;
}

// The following two functions are the functions used
// To find  circle when three points are given.

// Helper method to get a circle defined by 3 points
Point get_circle_center(double ax, double ay,
						double bx, double by)
{
	double A = ax * ax + ay * ay;
	double B = bx * bx + by * by;
	double C = ax * by - ay * bx;
	return { (by * A - ay * B) / (2 * C),
			(ax * B - bx * A) / (2 * C) };
}

// Function to return circle that intersects three points
Circle circle_from(const Point& A, const Point& B,
				const Point& C){
	Point I = get_circle_center(B.X - A.X, B.Y - A.Y,
								C.X - A.X, C.Y - A.Y);
	I.X += A.X;
	I.Y += A.Y;
	return { I, cal_dist(I, A) };
}

// Function to return circle that intersects 2 points
Circle circle_from(const Point& A, const Point& B){
	// Center is at middle point of A and B
	Point C = { (A.X + B.X) / 2.0, (A.Y + B.Y) / 2.0 };

	// Radius is the distance between AB
	return { C, cal_dist(A, B) / 2.0 };
}

// Function to check that a circle encloses the given points
bool is_valid_circle(const Circle& c, struct Point P[ ],int p_count){

	// Iterating through all the points to check
	// whether the points lie inside the circle or not
	while( p_count >= 0){
        if (!is_pnt_inside(c, P[p_count])){
            p_count --;
            return false;
        }
		p_count --;
	}
	return true;
}

// Function to return find the minimum enclosing
// circle from the given set of points
Circle minimum_enclosing_circle( struct Point P[ ],int p_count){
	if (p_count == 0)
		return { { 0, 0 }, 0 };
	if (p_count == 1)
		return { P[0], 0 };

	// Set initial MEC to have infinity radius
	Circle mec = { { 0, 0 }, INF };

	// Calculations to find the circle

	// Go over all pair of points
	for (int i = 0; i < p_count; i++) {
		for (int j = i + 1; j < p_count; j++) {

			// Get the smallest circle that
			// intersects P[i] and P[j]
			Circle tmp = circle_from(P[i], P[j]);

			// Update MEC if tmp encloses all points
			// and has a smaller radius
			if (tmp.R < mec.R && is_valid_circle(tmp, P, p_count)){
                mec = tmp;
			}
		}
	}

	// Go over all triples of points
	for (int i = 0; i < p_count; i++) {
		for (int j = i + 1; j < p_count; j++) {
			for (int k = j + 1; k < p_count; k++) {

				// Get the circle that intersects P[i], P[j], P[k]
				Circle tmp = circle_from(P[i], P[j], P[k]);

				// Update MEC if tmp encloses all points
				// and has smaller radius
				if (tmp.R < mec.R && is_valid_circle(tmp, P, p_count))
					mec = tmp;
			}
		}
	}
	return mec;
}

// Draw Circle & print values at console
void draw_circle(Circle mec,struct Point P[ ], int counter){
    int factor = set_interval(P, counter);

    printf("\n\nCenter: {%.1f, %.1f} Radius: %.1f\n", mec.C.X, mec.C.Y, mec.R);
	setcolor(RED);
	int X = pln_sz + (int(mec.C.X * 20)/factor);
	int Y = pln_sz - (int(mec.C.Y * 20)/factor);
	int R = (int(mec.R * 20)/factor);
	circle(X, Y, R);

	setcolor(YELLOW);
	char coordinate_text[20];
        sprintf(coordinate_text,"R{%.1f, %.1f}", mec.C.X, mec.C.Y);
        circle(X, Y,3);
        outtextxy(X + 5, Y + 5,coordinate_text);
}

// B-spline Functions

// Arrange points in ascending order
void arrange_pnt(Point P[ ], int counter){

    int i, j;
    struct Point tmp;

    for (i = 0; i < counter - 1; i++)
    {
        for (j = 0; j < (counter - 1-i); j++)
        {
            if ((P[j].X + P[j].Y) < (P[j+1].X + P[j+1].Y))
            {
                tmp = P[j];
                P[j] = P[j + 1];
                P[j + 1] = tmp;
            }
        }
    }
    /*//Code to test ascending is true by printing new point set
    printf("\n");
    while(counter >= 0){
        printf("{%d,%d} ",int(P[counter].X),int(P[counter].Y));
        counter--;
    }*/
}

//De boor's algorithm
Point deBoor (int k, int h, int i, double t, double *knots, Point cnt[ ]){
	Point P;
	P.X=0; P.Y=0;
	if (h==0){
		P=cnt[i];
	}
	else{
		P.X=( 1-(t-knots[i])/(knots[i+k+1-h]-knots[i]) )*deBoor(k,h-1,i-1,t,knots,cnt).X
		 + (t-knots[i])/(knots[i+k+1-h]-knots[i])*deBoor(k,h-1,i,t,knots, cnt).X;

		P.Y=( 1-(t-knots[i])/(knots[i+k+1-h]-knots[i]) )*deBoor(k,h-1,i-1,t,knots,cnt).Y
		 + (t-knots[i])/(knots[i+k+1-h]-knots[i])*deBoor(k,h-1,i,t,knots, cnt).Y;
	}
	return P;
}

// Function to find B-spline
void bspline(struct Point P[ ], int counter){
    int factor = set_interval(P, counter);
    arrange_pnt(P, counter);

    Point p;

    p.X = 0;
    p.Y = 0;

	double t;
    int degree = 2;
    int multi = 10*(20*counter*factor*degree);//multiplier for line drawing

    double knots [(counter+degree+1)];
    for(int i=0; i<counter+degree+1; i++) {
      knots[i] = i;
    }

    //determine the t values of the first & last point
	double t0 = knots[0];
	double tmax = knots[counter];

   	for (int index=0; index < multi; index++){
		t = t0 + index*(tmax-t0)/ multi;
		for(int j=degree;j<counter;j++){
			if (t>=knots[j] && t<=knots[j+1]){
				p = deBoor(degree,degree,j,t,knots, P);
				setcolor(BLUE);
                circle(pln_sz +(p.X * 20)/factor, pln_sz - (p.Y * 20)/factor, 1);
				break;
			}
		}
	}
}
int main()
{
    initwindow(pln_sz*2, pln_sz*2);

    int pnt_cntr = 0;
    struct Point points[32];

    //assign array to 0 to avoid runtime errors
    for(int j = 31; j >= 0;j--){
    points[j].X = 0;
    points[j].Y = 0;
    };

    //Read points from file
    FILE *read_points = fopen("coordinates2.txt", "r");

    if(read_points == NULL){
    printf("file is empty");
    exit(EXIT_FAILURE);
    }
    while(!feof(read_points)){
        fscanf(read_points,"%lf %lf",&points[pnt_cntr].X, &points[pnt_cntr].Y);
        printf("{%d,%d} ",int(points[pnt_cntr].X),int(points[pnt_cntr].Y));

        pnt_cntr++;
    }

    pnt_cntr--;
    fclose(read_points);

    // Draw coordinate plane and show points
    int factor = set_interval(points, pnt_cntr);

    draw_coordinade_plane(factor);

    show_points(points, pnt_cntr, factor);

    // Find & Draw MEC
	Circle mec = minimum_enclosing_circle(points, pnt_cntr);

	draw_circle(mec, points, pnt_cntr);

	// Find & Draw B-spline
    bspline(points, pnt_cntr);

    getch();
    closegraph();
    return 0;
}
