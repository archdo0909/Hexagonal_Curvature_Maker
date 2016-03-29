// air_bag_2014_10_21.cpp : コンソール アプリケーションのエントリ ポイントを定義します。
//

// finger_contact_model.cpp : コンソール アプリケーションのエントリ ポイントを定義します。
//

#include <Windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <GL/glut.h>
#include <math.h>
#include <string.h>
#include <vector>

using namespace std;
#define Hexa_h 78
#define Hexa_w 40
#define judge 0.0001
#define e_judge 0.000000001
int node_Num_m[2];
int node_Num_n[2];
bool up_flag = false;
bool down_flag = false;
bool close_flag = false;
bool open_flag = false;
bool close_flag_n = false;
bool open_flag_n = false;
double c[3];
int window_size_x = 1100;
int window_size_y = window_size_x;
GLfloat orange[] = { 255.0 / 256.0, 153.0 / 256.0, 0.0 / 256.0, 0.9 };
GLfloat blue[] = { 0.0 / 256.0, 65.0 / 256.0, 255.0 / 256.0, 0.4 };
GLfloat white[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat blue2[] = { 102.0 / 256.0, 204.0 / 256.0, 255.0 / 256.0, 0.9 };
GLfloat blue_node[] = { 0.5, 0.5, 1.0, 1.0 };
int w_view;
int h_view;
int first_count = 1;
double wall_z = 0.5;
double wall_n = 5.0;
int num_count = 0;
int con_count = 0;
int tri_count = 0;
double damp_k = 1000.0;
double damp_k_normal = 10;
double dv = 3.5;
double node_Radius = 0.08;
double View_from[3] = { 0.0, 120.0, 120.0 };
//double View_from[3] = { 0.0, 20.0, 70.0 };  //side
//double View_from[3] = { 0.0, 70.0, 20.0 };  //above
double View_to[3] = { 0.0, -15.0, 0.0 };
double View_from2[3] = { 0.0, 15.0, 0.01 };
double View_to2[3] = { 0.0, -10.0, 0.0 };
double View_from3[3] = { 0.0, 13.0, 0.01 };
double View_to3[3] = { 0.0, -10.0, 0.0 };
bool MouseFlagRight = false;
bool MouseFlagLeft = false;
bool MouseFlagMiddle = false;
bool View_point_flag = false;
GLUnurbsObj *theNurb; 
typedef struct{
	double x[3];
}position;
typedef struct{
	int number;
	int node_Num_w;
	int node_Num_h;
}face;
typedef struct{
	int t[3];
	int color;
	double normal[3];
}triangle_d;
typedef struct{
	bool torf;
	double len;
	//int color;
}edge_d;
typedef struct{
	int number;
	int edge_flag;			//ノード番号
	int none_flag;
	position pos;		//位置
}node;
typedef struct{
	int number;
	int edge_flag;			//ノード番号
	int none_flag;
	int N[6];
	int T[6];
	int index_n;
	int index_t;
	double m_normal[3];
	double n_normal[3];
	double cosa[6];
	double cosb[6];
	double cota[6];
	double cotb[6];
	double K;
	position pos;		//位置
	position del_pos;	//速度
	position acc;		//加速度
	double color_grad;
}node2;
typedef struct{
	bool flag;
	position ini_position;
}ini_flag;
typedef struct{
	int number;
	position pos;
	face face;
}point;
static node2 node_surface2[50000];
static node node_surface[1000][1000][3];
static edge_d edge[10000][10000];
static triangle_d triangle_data[50000];

void mult_matrix3x1(double *c, double *a, double *b){
	for (int i = 0; i < 3; i++){
		//for(int j=0;j<3;j++){
		c[i] = 0.0;
		for (int k = 0; k < 3; k++){
			c[i] += a[3 * i + k] * b[k];
		}
		//}
	}
}
void gaiseki_9_3(double *a, double *b){
	a[0] = 0.0;
	a[1] = -b[2];
	a[2] = b[1];
	a[3] = b[2];
	a[4] = 0.0;
	a[5] = -b[0];
	a[6] = -b[1];
	a[7] = b[0];
	a[8] = 0.0;
}
void sphere(double R, double precise, GLfloat sph_col[10]){

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, sph_col);
	GLUquadricObj *sphere;

	sphere = gluNewQuadric();
	gluQuadricDrawStyle(sphere, GLU_FILL);

	gluSphere(sphere, R, precise, precise);
}
void View_control(bool vector_flag){
	double View_distance;
	double temp[5];
	temp[2] = View_from[2] - View_to[2];
	temp[1] = View_from[0] - View_to[0];
	temp[0] = pow(temp[1], 2.0) + pow(temp[2], 2.0);
	View_distance = pow(temp[0], 0.5);
		//printf("%f\n", View_distance);
	temp[0] = View_from[2] - View_to[2];
	temp[3] = temp[0] / View_distance;
	temp[1] = View_from[0] - View_to[0];
	temp[4] = temp[1] / View_distance;
	temp[2] = atan2(temp[4], temp[3]);
	//temp[2] = acos(temp[1]);
	if (vector_flag) temp[2] = temp[2] + 0.01;
	else temp[2] = temp[2] - 0.01;
	temp[0] = View_distance * cos(temp[2]);
	temp[1] = View_distance * sin(temp[2]);
	View_from[2] = View_to[2] + temp[0];
	View_from[0] = View_to[0] + temp[1];
}
void View_control_up_down(bool vector_flag){
	double View_distance;
	double temp[5];
	temp[2] = View_from[1] - View_to[1];
	temp[1] = View_from[0] - View_to[0];
	temp[0] = pow(temp[1], 2.0) + pow(temp[2], 2.0);
	View_distance = pow(temp[0], 0.5);
	//	printf("%f\n", View_distance);
	temp[0] = View_from[1] - View_to[1];
	temp[3] = temp[0] / View_distance;
	temp[1] = View_from[0] - View_to[0];
	temp[4] = temp[1] / View_distance;
	temp[2] = atan2(temp[4], temp[3]);
	//temp[2] = acos(temp[1]);
	if (vector_flag) temp[2] = temp[2] + 0.01;
	else temp[2] = temp[2] - 0.01;
	temp[0] = View_distance * cos(temp[2]);
	temp[1] = View_distance * sin(temp[2]);
	View_from[1] = View_to[1] + temp[0];
	View_from[0] = View_to[0] + temp[1];
}
void View_control2(bool vector_flag){
	double View_distance;
	double temp[5];
	temp[2] = View_from2[2] - View_to2[2];
	temp[1] = View_from2[0] - View_to2[0];
	temp[0] = pow(temp[1], 2.0) + pow(temp[2], 2.0);
	View_distance = pow(temp[0], 0.5);
	temp[0] = View_from2[2] - View_to2[2];
	temp[3] = temp[0] / View_distance;
	temp[1] = View_from2[0] - View_to2[0];
	temp[4] = temp[1] / View_distance;
	temp[2] = atan2(temp[4], temp[3]);
	if (vector_flag) temp[2] = temp[2] + 0.01;
	else temp[2] = temp[2] - 0.01;
	temp[0] = View_distance * cos(temp[2]);
	temp[1] = View_distance * sin(temp[2]);
	View_from2[2] = View_to2[2] + temp[0];
	View_from2[0] = View_to2[0] + temp[1];
}
void initiation(){
	int i = 0;
	int j = 0;
	int k = 0;
	int h = 0;
	int s = 0;
	int trirem1[3];
	int trirem2[3];
	int trirem3[3];
	int trirem4[3];
	int trirem5[3];
	int tri_flag1;
	int tri_flag2;
	int tri_flag3;
	int tri_flag4;
	int tri_flag5;
	double tritemp_x;
	double tritemp_y;
	double natural_length_x;
	double natural_length_y;
	double natural_length = 1.0;
	int flag = 0;


	//static node node_surface[100][100][2];

	//for (s = 0; s < 1; s++){
	//	for (i = 0; i <= Hexa_h * 0.5; i++){
	//		for (j = 0; j <= i + Hexa_w - 1; j++){
	//			if (i == Hexa_h * 0.5){
	//					node_surface[i][j][s].pos.x[0] = (double)((((-Hexa_w - i + 1) * 0.5) + j) * natural_length);
	//					node_surface[i][j][s].pos.x[1] = (double)(i * sqrt(3) * 0.5);
	//					node_surface[i][j][s].pos.x[2] = (double)(0.0);
	//					//printf("pos = %f, %f, %d, %d\n", node_surface[i][j][s].pos.x[0], node_surface[i][j][s].pos.x[1], i, j);
	//			}
	//			else{
	//				if (i % 2 == 0){
	//					int I = Hexa_h - i;
	//					node_surface[i][j][s].pos.x[0] = (double)((j - ((i + Hexa_w - 1) * 0.5)) * natural_length);
	//					node_surface[I][j][s].pos.x[0] = (double)((j - ((i + Hexa_w - 1) * 0.5)) * natural_length);
	//					node_surface[i][j][s].pos.x[1] = (double)(i * sqrt(3) * 0.5);
	//					node_surface[I][j][s].pos.x[1] = (double)(I * sqrt(3) * 0.5);
	//					node_surface[i][j][s].pos.x[2] = (double)(0.0);
	//					node_surface[I][j][s].pos.x[2] = (double)(0.0);
	//					//printf("pos = %f, %f, %d, %d\n", node_surface[I][j][s].pos.x[0], node_surface[I][j][s].pos.x[1], I, j);
	//					//printf("pos = %f, %f, %d, %d\n", node_surface[i][j][s].pos.x[0], node_surface[i][j][s].pos.x[1], i, j);
	//				}
	//				else{
	//					int I = Hexa_h - i;
	//					node_surface[i][j][s].pos.x[0] = (double)((((-Hexa_w - i + 1) * 0.5) + j) * natural_length);
	//					node_surface[I][j][s].pos.x[0] = (double)((((-Hexa_w - i + 1) * 0.5) + j) * natural_length);
	//					node_surface[i][j][s].pos.x[1] = (double)(i * sqrt(3) * 0.5);
	//					node_surface[I][j][s].pos.x[1] = (double)(I * sqrt(3) * 0.5);
	//					node_surface[i][j][s].pos.x[2] = (double)(0.0);
	//					node_surface[I][j][s].pos.x[2] = (double)(0.0);
	//					//printf("pos = %f, %f, %d, %d\n", node_surface[i][j][s].pos.x[0], node_surface[i][j][s].pos.x[1], i, j);
	//					//printf("pos = %f, %f, %d, %d\n", node_surface[I][j][s].pos.x[0], node_surface[I][j][s].pos.x[1], I, j);
	//				}
	//			}
	//		}
	//	}
	//}

		for (s = 0; s < 2; s++){
			for (i = 0; i < Hexa_h * 0.5; i++){
				for (j = 0; j < 2 * Hexa_w - 1 - (Hexa_h * 0.5) + i; j++){
					//printf("i, j = %d, %d\n", i, j);
					if (i % 2 == 0){
						node_surface[i][j][s].pos.x[0] = (double)((-(Hexa_w - 1 - Hexa_h * 0.25 + i * 0.5) + j) * natural_length);
						node_surface[i][j][s].pos.x[1] = (double)(i * sqrt(3) * 0.5);
						node_surface[i][j][s].pos.x[2] = (double)(0.0);
					}
					else{
						node_surface[i][j][s].pos.x[0] = (double)(-(Hexa_w - 1 - Hexa_h * 0.25 + i * 0.5) + j) * natural_length;
						node_surface[i][j][s].pos.x[1] = (double)(i * sqrt(3) * 0.5);
						node_surface[i][j][s].pos.x[2] = (double)(0.0);
						//printf("pos = %d, %d, %f\n", i, j, node_surface[i][j][s].pos.x[0]);
					}
				}
			}
			for (i = Hexa_h * 0.5; i <= Hexa_h; i++){
				for (j = 0; j < 2 * Hexa_w - 1 - (i - Hexa_h * 0.5); j++){
					//printf("i, j = %d, %d\n", i, j);
					if (i % 2 == 0){
						node_surface[i][j][s].pos.x[0] = (double)((-(Hexa_w - 1 - (i - Hexa_h * 0.5) * 0.5) + j) * natural_length);
						node_surface[i][j][s].pos.x[1] = (double)(i * sqrt(3) * 0.5);
						node_surface[i][j][s].pos.x[2] = (double)(0.0);

					}
					else{
						node_surface[i][j][s].pos.x[0] = (double)((-(Hexa_w - 1 - (i - Hexa_h * 0.5) *0.5) + j) * natural_length);
						node_surface[i][j][s].pos.x[1] = (double)(i * sqrt(3) * 0.5);
						node_surface[i][j][s].pos.x[2] = (double)(0.0);
					}
				}
			}
		}
	


	/*for (s = 0; s < 1; s++){
		if (Hexa_w % 2 == 0){
			for (i = 0; i < Hexa_h; i++){
				if (i % 2 == 0){
					for (j = 0; j < Hexa_w - i; j++){
						node_surface[i][j][s].pos.x[0] = (double)(natural_length * 0.5 * (2 * j + 2 * i - (Hexa_w - 1)));
					}
				}
				else{
					for (j = 0; j < Hexa_w - i; j++){
						node_surface[i][j][s].pos.x[0] = (double)(natural_length *(j - (Hexa_w - 4)));
					}
				}
			}
		}
		else{
			for (i = 0; i < Hexa_h; i++){
				if (i % 2 == 0){
					for (j = 0; j < Hexa_w - i; j++){
						node_surface[i][j][s].pos.x[0] = (double)(natural_length * 0.5 * (2 * j + 2 * i - (Hexa_w + 1)));
					}
				}
				else{
					for (j = 0; j < Hexa_w - i; j++){
						node_surface[i][j][s].pos.x[0] = (double)(natural_length *(j - (Hexa_w - 3)));
					}
				}
			}
		}
		if (Hexa_w % 2 == 0){
			for (i = -1; i > -Hexa_h; i--){
				if (i % 2 == 0){
					for (j = 0; j < Hexa_w - i; j++){
						node_surface[i][j][s].pos.x[0] = (double)((0.0) + natural_length * 0.5 * (2 * j + 1));
					
					}
				for (j = -1; j > -Hexa_w + i; j--){
					node_surface[i][j][s].pos.x[0] = (double)((0.0) + natural_length * 0.5 * (2 * j + 1));
				}
			}
			else{
				for (j = 0; j < Hexa_w - i; j++){
					node_surface[i][j][s].pos.x[0] = (double)((0.0) + natural_length * j);
				}
				for (j = -1; j > -Hexa_w + i; j--){
					node_surface[i][j][s].pos.x[0] = (double)((0.0) + natural_length * j);
				}
			}
		}
	}*/
	for (s = 0; s < 2; s++){
		for (i = 0; i < Hexa_h * 0.5; i++){
			for (j = 0; j < 2 * Hexa_w - 1 - (Hexa_h * 0.5) + i; j++){
				node_surface[i][j][s].number = num_count;
				//printf("count = %d, %d, %d, %d\n", num_count, i, j, s);
				num_count++;
			}
		}
		for (i = Hexa_h * 0.5; i <= Hexa_h; i++){
			for (j = 0; j < 2 * Hexa_w - 1 - (i - Hexa_h * 0.5); j++){
				node_surface[i][j][s].number = num_count;
				num_count++;
				//printf("count1 = %d, %d, %d, %d\n", num_count, i, j, s);
			}
		}
	}

#if 0
	for (s = 0; s < 2; s++){
		for (i = 0; i < Hexa_h * 0.5; i++){
			for (j = 0; j <= i + Hexa_w - 1; j++){
				node_surface[i][j][s].number = num_count;
				num_count++;
				//printf("count = %d, %d, %d, %d\n", num_count, i, j, s);
			}
		}
			for (i = Hexa_h * 0.5; i <= Hexa_h; i++){
				for (j = 0; j < 2 * Hexa_w + Hexa_h * 0.5 - (i + 1); j++){
					node_surface[i][j][s].number = num_count;
					num_count++;
					//printf("count1 = %d, %d, %d, %d\n", num_count, i, j, s);
				}
			}
		}
#endif

	printf("num_count = %d\n", num_count);
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////mesh upside////////////////////////////////////////////////////////////////////////////////
#if 0
	for (s = 0; s < 2; s++){
		for (h = 0; h < num_count; h++){
			tri_flag3 = 0;
			for (i = 0; i <= Hexa_h * 0.5; i++){
				for (j = 0; j <= i + Hexa_w - 1; j++){
					if (node_surface[i][j][s].number == h){
						trirem1[0] = i;
						trirem1[1] = j;
						trirem1[2] = s;
						tri_flag3 = 1;
						//printf("trirem1 = %d, %d, %d\n", trirem1[0], trirem1[1], trirem1[2]);
					}
				}
			}
			tri_flag1 = 0;
			tri_flag2 = 0;
			tri_flag4 = 0;
			tri_flag5 = 0;
			for (i = 0; i <= Hexa_h * 0.5; i++){
				for (j = 0; j <= i + Hexa_w - 1; j++){
					tritemp_x = node_surface[trirem1[0]][trirem1[1]][trirem1[2]].pos.x[0] - node_surface[i][j][s].pos.x[0];
					tritemp_y = node_surface[trirem1[0]][trirem1[1]][trirem1[2]].pos.x[1] - node_surface[i][j][s].pos.x[1];
					if (fabs(tritemp_x + natural_length * 0.5) < judge && fabs(tritemp_y + (sqrt(3) * natural_length * 0.5)) < judge){
						trirem2[0] = i;
						trirem2[1] = j;
						trirem2[2] = s;
						tri_flag1 = 1;
						//printf("trirem2 = %d, %d, %d\n", trirem2[0], trirem2[1], trirem2[2]);
					}
					if (fabs(tritemp_x - natural_length * 0.5) < judge && fabs(tritemp_y + (sqrt(3) * natural_length * 0.5)) < judge){
						trirem3[0] = i;
						trirem3[1] = j;
						trirem3[2] = s;
						tri_flag2 = 1;
						//printf("trirem3 = %d, %d, %d\n", trirem3[0], trirem3[1], trirem3[2]);
					}
					if (fabs(tritemp_x + natural_length) < judge && fabs(tritemp_y) < judge){
						trirem4[0] = i;
						trirem4[1] = j;
						trirem4[2] = s;
						tri_flag4 = 1;
						//printf("trirem4 = %d, %d, %d\n", trirem4[0], trirem4[1], trirem4[2]);
					}
					if (fabs(tritemp_x + natural_length * 0.5) < judge && fabs(tritemp_y + (sqrt(3) * natural_length * 0.5)) < judge){
						trirem5[0] = i;
						trirem5[1] = j;
						trirem5[2] = s;
						tri_flag5 = 1;
						//printf("trirem5 = %d, %d, %d\n", trirem5[0], trirem5[1], trirem5[2]);
					}
				}
			}
			if (tri_flag1 == 1 && tri_flag2 == 1 && tri_flag3 == 1){
				if (tri_count == 0){
					triangle_data[tri_count].t[0] = h;
					triangle_data[tri_count].t[1] = node_surface[trirem2[0]][trirem2[1]][trirem2[2]].number;
					triangle_data[tri_count].t[2] = node_surface[trirem3[0]][trirem3[1]][trirem3[2]].number;
					triangle_data[tri_count].color = 1;
					tri_count++;
					//printf("tri_data = %d, %d, %d\n", triangle_data[tri_count].t[0], triangle_data[tri_count].t[1], triangle_data[tri_count].t[2]);
				}
				else{
					flag = 0;
					for (i = 0; i < tri_count; i++){
						if ((triangle_data[i].t[0] == h) && (triangle_data[i].t[1] == node_surface[trirem2[0]][trirem2[1]][trirem2[2]].number)
							&& (triangle_data[i].t[2] = node_surface[trirem3[0]][trirem3[1]][trirem3[2]].number)){
							flag = 1;
							//printf("tri_data1 = %d, %d, %d\n", triangle_data[i].t[0], triangle_data[i].t[1], triangle_data[i].t[2]);
						}
					}
					if (flag == 0){
						triangle_data[tri_count].t[0] = h;
						triangle_data[tri_count].t[1] = node_surface[trirem2[0]][trirem2[1]][trirem2[2]].number;
						triangle_data[tri_count].t[2] = node_surface[trirem3[0]][trirem3[1]][trirem3[2]].number;
						triangle_data[tri_count].color = 1;
						tri_count++;
						//printf("tri_data3 = %d, %d, %d\n", triangle_data[tri_count].t[0], triangle_data[tri_count].t[1], triangle_data[tri_count].t[2]);
					}
					else{
						flag = 0;
					}
				}
			}
			if (tri_flag3 == 1 && tri_flag5 == 1 && tri_flag4 == 1){
				if (tri_count == 0){
					triangle_data[tri_count].t[0] = h;
					triangle_data[tri_count].t[1] = node_surface[trirem4[0]][trirem4[1]][trirem4[2]].number;
					triangle_data[tri_count].t[2] = node_surface[trirem5[0]][trirem5[1]][trirem5[2]].number;
					triangle_data[tri_count].color = 1;
					tri_count++;
				}
				else{
					flag = 0;
					for (i = 0; i < tri_count; i++){
						if ((triangle_data[i].t[0] == h) && (triangle_data[i].t[1] == node_surface[trirem4[0]][trirem4[1]][trirem4[2]].number)
							&& (triangle_data[i].t[2] == node_surface[trirem5[0]][trirem5[1]][trirem5[2]].number)){
							flag = 1;
						}
					}
					if (flag == 0){
						triangle_data[tri_count].t[0] = h;
						triangle_data[tri_count].t[1] = node_surface[trirem4[0]][trirem4[1]][trirem4[2]].number;
						triangle_data[tri_count].t[2] = node_surface[trirem5[0]][trirem5[1]][trirem5[2]].number;
						triangle_data[tri_count].color = 1;
						tri_count++;
					}
					else{
						flag = 0;
					}
				}
			}
		}
		for (h = 0; h < num_count; h++){
			tri_flag3 = 0;
			for (i = Hexa_h * 0.5; i <= Hexa_h; i++){
				for (j = 0; j < 2 * Hexa_w + Hexa_h * 0.5 - (i + 1); j++){
					if (node_surface[i][j][s].number == h){
						trirem1[0] = i;
						trirem1[1] = j;
						trirem1[2] = s;
						tri_flag3 = 1;
						//printf("trirem1 = %d, %d, %d\n", trirem1[0], trirem1[1], trirem1[2]);
					}
				}
			}
			tri_flag1 = 0;
			tri_flag2 = 0;
			tri_flag4 = 0;
			tri_flag5 = 0;
			for (i = Hexa_h * 0.5; i <= Hexa_h; i++){
				for (j = 0; j < 2 * Hexa_w + Hexa_h * 0.5 - (i + 1); j++){
					tritemp_x = node_surface[trirem1[0]][trirem1[1]][trirem1[2]].pos.x[0] - node_surface[i][j][s].pos.x[0];
					tritemp_y = node_surface[trirem1[0]][trirem1[1]][trirem1[2]].pos.x[1] - node_surface[i][j][s].pos.x[1];
					if (fabs(tritemp_x - natural_length * 0.5) < judge && fabs(tritemp_y - (sqrt(3) * natural_length * 0.5)) < judge){
						trirem2[0] = i;
						trirem2[1] = j;
						trirem2[2] = s;
						tri_flag1 = 1;
						//printf("trirem2 = %d, %d, %d\n", trirem2[0], trirem2[1], trirem2[2]);
					}
					if (fabs(tritemp_x + natural_length * 0.5) < judge && fabs(tritemp_y - (sqrt(3) * natural_length * 0.5)) < judge){
						trirem3[0] = i;
						trirem3[1] = j;
						trirem3[2] = s;
						tri_flag2 = 1;
						//printf("trirem3 = %d, %d, %d\n", trirem3[0], trirem3[1], trirem3[2]);
					}
					if (fabs(tritemp_x + natural_length * 0.5) < judge && fabs(tritemp_y - (sqrt(3) * natural_length * 0.5)) < judge){
						trirem4[0] = i;
						trirem4[1] = j;
						trirem4[2] = s;
						tri_flag4 = 1;
					}
					if (fabs(tritemp_x + natural_length) < judge && fabs(tritemp_y) < judge){
						trirem5[0] = i;
						trirem5[1] = j;
						trirem5[2] = s;
						tri_flag5 = 1;
					}
				}
			}
			if (tri_flag1 == 1 && tri_flag2 == 1 && tri_flag3 == 1){
				flag = 0;
				for (i = 0; i < tri_count; i++){
					if ((triangle_data[i].t[0] == h) && (triangle_data[i].t[1] == node_surface[trirem2[0]][trirem2[1]][trirem2[2]].number)
						&& (triangle_data[i].t[2] = node_surface[trirem3[0]][trirem3[1]][trirem3[2]].number)){
						flag = 1;
					}
					//printf("a\n");
				}
				if (flag == 0){
					triangle_data[tri_count].t[0] = h;
					triangle_data[tri_count].t[1] = node_surface[trirem2[0]][trirem2[1]][trirem2[2]].number;
					triangle_data[tri_count].t[2] = node_surface[trirem3[0]][trirem3[1]][trirem3[2]].number;
					triangle_data[tri_count].color = 1;
					tri_count++;
				}
				else{
					flag = 0;
				}
			}
			if (tri_flag4 == 1 && tri_flag5 == 1 && tri_flag3 == 1){
				flag = 0;
				for (i = 0; i < tri_count; i++){
					if ((triangle_data[i].t[0] == h) && (triangle_data[i].t[1] == node_surface[trirem4[0]][trirem4[1]][trirem4[2]].number)
						&& (triangle_data[i].t[2] = node_surface[trirem5[0]][trirem5[1]][trirem5[2]].number)){
						flag = 1;
					}
					//printf("a\n");
				}
				if (flag == 0){
					triangle_data[tri_count].t[0] = h;
					triangle_data[tri_count].t[1] = node_surface[trirem4[0]][trirem4[1]][trirem4[2]].number;
					triangle_data[tri_count].t[2] = node_surface[trirem5[0]][trirem5[1]][trirem5[2]].number;
					triangle_data[tri_count].color = 1;
					tri_count++;
				}
				else{
					flag = 0;
				}
			}
		}
	}

#endif
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////mesh down////////////////////////////////////////////////////////////////////////////////
#if 1
	for (s = 0; s < 2; s++){
		for (h = 0; h < num_count; h++){
			tri_flag3 = 0;
			for (i = 0; i <= Hexa_h * 0.5; i++){
				for (j = 0; j < 2 * Hexa_w - 1 - (Hexa_h * 0.5) + i; j++){
					if (node_surface[i][j][s].number == h){
						trirem1[0] = i;
						trirem1[1] = j;
						trirem1[2] = s;
						tri_flag3 = 1;
						//printf("trirem1 = %d, %d, %d\n", trirem1[0], trirem1[1], trirem1[2]);
					}
				}
			}
			tri_flag1 = 0;
			tri_flag2 = 0;
			tri_flag4 = 0;
			tri_flag5 = 0;
			for (i = 0; i <= Hexa_h * 0.5; i++){
				for (j = 0; j < 2 * Hexa_w - 1 - (Hexa_h * 0.5) + i; j++){
					tritemp_x = node_surface[trirem1[0]][trirem1[1]][trirem1[2]].pos.x[0] - node_surface[i][j][s].pos.x[0];
					tritemp_y = node_surface[trirem1[0]][trirem1[1]][trirem1[2]].pos.x[1] - node_surface[i][j][s].pos.x[1];
					if (fabs(tritemp_x + natural_length * 0.5) < judge && fabs(tritemp_y + (sqrt(3) * natural_length * 0.5)) < judge){
						trirem2[0] = i;
						trirem2[1] = j;
						trirem2[2] = s;
						tri_flag1 = 1;
						//printf("trirem2 = %d, %d, %d\n", trirem2[0], trirem2[1], trirem2[2]);
					}
					if (fabs(tritemp_x - natural_length * 0.5) < judge && fabs(tritemp_y + (sqrt(3) * natural_length * 0.5)) < judge){
						trirem3[0] = i;
						trirem3[1] = j;
						trirem3[2] = s;
						tri_flag2 = 1;
						//printf("trirem3 = %d, %d, %d\n", trirem3[0], trirem3[1], trirem3[2]);
					}
					if (fabs(tritemp_x + natural_length) < judge && fabs(tritemp_y) < judge){
						trirem4[0] = i;
						trirem4[1] = j;
						trirem4[2] = s;
						tri_flag4 = 1;
						//printf("trirem4 = %d, %d, %d\n", trirem4[0], trirem4[1], trirem4[2]);
					}
					if (fabs(tritemp_x + natural_length * 0.5) < judge && fabs(tritemp_y + (sqrt(3) * natural_length * 0.5)) < judge){
						trirem5[0] = i;
						trirem5[1] = j;
						trirem5[2] = s;
						tri_flag5 = 1;
						//printf("trirem5 = %d, %d, %d\n", trirem5[0], trirem5[1], trirem5[2]);
					}
				}
			}
			if (tri_flag1 == 1 && tri_flag2 == 1 && tri_flag3 == 1){
				if (tri_count == 0){
					triangle_data[tri_count].t[0] = h;
					triangle_data[tri_count].t[1] = node_surface[trirem2[0]][trirem2[1]][trirem2[2]].number;
					triangle_data[tri_count].t[2] = node_surface[trirem3[0]][trirem3[1]][trirem3[2]].number;
					triangle_data[tri_count].color = 1;
					tri_count++;
					//printf("tri_data = %d, %d, %d\n", triangle_data[tri_count].t[0], triangle_data[tri_count].t[1], triangle_data[tri_count].t[2]);
				}
				else{
					flag = 0;
					for (i = 0; i < tri_count; i++){
						if ((triangle_data[i].t[0] == h) && (triangle_data[i].t[1] == node_surface[trirem2[0]][trirem2[1]][trirem2[2]].number)
							&& (triangle_data[i].t[2] = node_surface[trirem3[0]][trirem3[1]][trirem3[2]].number)){
							flag = 1;
							//printf("tri_data1 = %d, %d, %d\n", triangle_data[i].t[0], triangle_data[i].t[1], triangle_data[i].t[2]);
						}
					}
					if (flag == 0){
						triangle_data[tri_count].t[0] = h;
						triangle_data[tri_count].t[1] = node_surface[trirem2[0]][trirem2[1]][trirem2[2]].number;
						triangle_data[tri_count].t[2] = node_surface[trirem3[0]][trirem3[1]][trirem3[2]].number;
						triangle_data[tri_count].color = 1;
						tri_count++;
						//printf("tri_data3 = %d, %d, %d\n", triangle_data[tri_count].t[0], triangle_data[tri_count].t[1], triangle_data[tri_count].t[2]);
					}
					else{
						flag = 0;
					}
				}
			}
			if (tri_flag3 == 1 && tri_flag5 == 1 && tri_flag4 == 1){
				if (tri_count == 0){
					triangle_data[tri_count].t[0] = h;
					triangle_data[tri_count].t[1] = node_surface[trirem4[0]][trirem4[1]][trirem4[2]].number;
					triangle_data[tri_count].t[2] = node_surface[trirem5[0]][trirem5[1]][trirem5[2]].number;
					triangle_data[tri_count].color = 1;
					tri_count++;
				}
				else{
					flag = 0;
					for (i = 0; i < tri_count; i++){
						if ((triangle_data[i].t[0] == h) && (triangle_data[i].t[1] == node_surface[trirem4[0]][trirem4[1]][trirem4[2]].number)
							&& (triangle_data[i].t[2] == node_surface[trirem5[0]][trirem5[1]][trirem5[2]].number)){
							flag = 1;
						}
					}
					if (flag == 0){
						triangle_data[tri_count].t[0] = h;
						triangle_data[tri_count].t[1] = node_surface[trirem4[0]][trirem4[1]][trirem4[2]].number;
						triangle_data[tri_count].t[2] = node_surface[trirem5[0]][trirem5[1]][trirem5[2]].number;
						triangle_data[tri_count].color = 1;
						tri_count++;
					}
					else{
						flag = 0;
					}
				}
			}
		}
		for (h = 0; h < num_count; h++){
			tri_flag3 = 0;
			for (i = Hexa_h * 0.5; i <= Hexa_h; i++){
				for (j = 0; j < 2 * Hexa_w - 1 - (i - Hexa_h * 0.5); j++){
					if (node_surface[i][j][s].number == h){
						trirem1[0] = i;
						trirem1[1] = j;
						trirem1[2] = s;
						tri_flag3 = 1;
						//printf("trirem1 = %d, %d, %d\n", trirem1[0], trirem1[1], trirem1[2]);
					}
				}
			}
			tri_flag1 = 0;
			tri_flag2 = 0;
			tri_flag4 = 0;
			tri_flag5 = 0;
			for (i = Hexa_h * 0.5; i <= Hexa_h; i++){
				for (j = 0; j < 2 * Hexa_w - 1 - (i - Hexa_h * 0.5); j++){
					tritemp_x = node_surface[trirem1[0]][trirem1[1]][trirem1[2]].pos.x[0] - node_surface[i][j][s].pos.x[0];
					tritemp_y = node_surface[trirem1[0]][trirem1[1]][trirem1[2]].pos.x[1] - node_surface[i][j][s].pos.x[1];
					if (fabs(tritemp_x - natural_length * 0.5) < judge && fabs(tritemp_y - (sqrt(3) * natural_length * 0.5)) < judge){
						trirem2[0] = i;
						trirem2[1] = j;
						trirem2[2] = s;
						tri_flag1 = 1;
						//printf("trirem2 = %d, %d, %d\n", trirem2[0], trirem2[1], trirem2[2]);
					}
					if (fabs(tritemp_x + natural_length * 0.5) < judge && fabs(tritemp_y - (sqrt(3) * natural_length * 0.5)) < judge){
						trirem3[0] = i;
						trirem3[1] = j;
						trirem3[2] = s;
						tri_flag2 = 1;
						//printf("trirem3 = %d, %d, %d\n", trirem3[0], trirem3[1], trirem3[2]);
					}
					if (fabs(tritemp_x + natural_length * 0.5) < judge && fabs(tritemp_y - (sqrt(3) * natural_length * 0.5)) < judge){
						trirem4[0] = i;
						trirem4[1] = j;
						trirem4[2] = s;
						tri_flag4 = 1;
					}
					if (fabs(tritemp_x + natural_length) < judge && fabs(tritemp_y) < judge){
						trirem5[0] = i;
						trirem5[1] = j;
						trirem5[2] = s;
						tri_flag5 = 1;
					}
				}
			}
			if (tri_flag1 == 1 && tri_flag2 == 1 && tri_flag3 == 1){
				flag = 0;
				for (i = 0; i < tri_count; i++){
					if ((triangle_data[i].t[0] == h) && (triangle_data[i].t[1] == node_surface[trirem2[0]][trirem2[1]][trirem2[2]].number)
						&& (triangle_data[i].t[2] = node_surface[trirem3[0]][trirem3[1]][trirem3[2]].number)){
						flag = 1;
					}
					//printf("a\n");
				}
				if (flag == 0){
					triangle_data[tri_count].t[0] = h;
					triangle_data[tri_count].t[1] = node_surface[trirem2[0]][trirem2[1]][trirem2[2]].number;
					triangle_data[tri_count].t[2] = node_surface[trirem3[0]][trirem3[1]][trirem3[2]].number;
					triangle_data[tri_count].color = 1;
					tri_count++;
				}
				else{
					flag = 0;
				}
			}
			if (tri_flag4 == 1 && tri_flag5 == 1 && tri_flag3 == 1){
				flag = 0;
				for (i = 0; i < tri_count; i++){
					if ((triangle_data[i].t[0] == h) && (triangle_data[i].t[1] == node_surface[trirem4[0]][trirem4[1]][trirem4[2]].number)
						&& (triangle_data[i].t[2] = node_surface[trirem5[0]][trirem5[1]][trirem5[2]].number)){
						flag = 1;
					}
					//printf("a\n");
				}
				if (flag == 0){
					triangle_data[tri_count].t[0] = h;
					triangle_data[tri_count].t[1] = node_surface[trirem4[0]][trirem4[1]][trirem4[2]].number;
					triangle_data[tri_count].t[2] = node_surface[trirem5[0]][trirem5[1]][trirem5[2]].number;
					triangle_data[tri_count].color = 1;
					tri_count++;
				}
				else{
					flag = 0;
				}
			}
		}
	}

#endif

	printf("%d\n", tri_count);
	//node2 numbering nad edge(first trial)
#if 0
	for (s = 0; s < 2; s++){
		for (i = 0; i <= Hexa_h * 0.5; i++){
			for (j = 0; j <= i + Hexa_w - 1; j++){
				if (i == Hexa_h * 0.5){
					for (h = 0; h < 3; h++){
						node_surface2[node_surface[i][j][s].number].pos.x[h] = node_surface[i][j][s].pos.x[h];
					}
				}
				else{
					if (i % 2 == 0){
						for (h = 0; h < 3; h++){
							int I = Hexa_h - i;
							node_surface2[node_surface[i][j][s].number].pos.x[h] = node_surface[i][j][s].pos.x[h];
							node_surface2[node_surface[I][j][s].number].pos.x[h] = node_surface[I][j][s].pos.x[h];
						}
					}
					else{
						for (h = 0; h < 3; h++){
							int I = Hexa_h - i;
							node_surface2[node_surface[i][j][s].number].pos.x[h] = node_surface[i][j][s].pos.x[h];
							node_surface2[node_surface[I][j][s].number].pos.x[h] = node_surface[I][j][s].pos.x[h];
						}
					}
				}
			}
		}
	}

	for (s = 0; s < 2; s++){
		for (i = 0; i <= Hexa_h * 0.5; i++){
			for (j = 0; j <= i + Hexa_w - 1; j++){
				if (i == Hexa_h * 0.5){
					if (j == 0 || j == i + Hexa_w - 1){
						node_surface2[node_surface[i][j][s].number].edge_flag == 1;
						node_surface2[node_surface[i][j][s].number].none_flag == 0;
					}
					else{
						node_surface2[node_surface[i][j][s].number].edge_flag == 0;
						node_surface2[node_surface[i][j][s].number].none_flag == 1;
					}
				}
				else{
					if (i == 0 || j == 0 || j == i + Hexa_w - 1){
						int I = Hexa_h - i;
						node_surface2[node_surface[i][j][s].number].edge_flag == 1;
						node_surface2[node_surface[i][j][s].number].none_flag == 0;
						node_surface2[node_surface[I][j][s].number].edge_flag == 1;
						node_surface2[node_surface[I][j][s].number].none_flag == 0;
					}
					else{
						int I = Hexa_h - i;
						node_surface2[node_surface[i][j][s].number].edge_flag == 0;
						node_surface2[node_surface[i][j][s].number].none_flag == 1;
						node_surface2[node_surface[I][j][s].number].edge_flag == 0;
						node_surface2[node_surface[I][j][s].number].none_flag == 1;
					}
				}
			}
		}
	}
#endif
	//node2 numbering (second)
#if 0
	for (s = 0; s < 2; s++){
		for (i = 0; i < Hexa_h * 0.5; i++){
			for (j = 0; j <= i + Hexa_w - 1; j++){
				for (h = 0; h < 3; h++){
					node_surface2[node_surface[i][j][s].number].pos.x[h] = node_surface[i][j][s].pos.x[h];
					//printf("node2 = %d, %f, %f, %f\n", node_surface[i][j][s].number, node_surface2[node_surface[i][j][s].number].pos.x[0], node_surface2[node_surface[i][j][s].number].pos.x[1], node_surface2[node_surface[i][j][s].number].pos.x[2]);
				}
			}
		}
		for (i = Hexa_h * 0.5; i <= Hexa_h; i++){
			for (j = 0; j < 2 * Hexa_w + Hexa_h * 0.5 - (i + 1); j++){
				for (h = 0; h < 3; h++){
					node_surface2[node_surface[i][j][s].number].pos.x[h] = node_surface[i][j][s].pos.x[h];
					//printf("node3 = %d, %f, %f, %f\n", node_surface[i][j][s].number, node_surface2[node_surface[i][j][s].number].pos.x[0], node_surface2[node_surface[i][j][s].number].pos.x[1], node_surface2[node_surface[i][j][s].number].pos.x[2]);
				}
			}
		}
	}

	for (s = 0; s < 2; s++){
		for (i = 0; i < Hexa_h * 0.5; i++){
			for (j = 0; j <= i + Hexa_w - 1; j++){
				if (i == 0 || j == 0 || j == i + Hexa_w - 1){
					//printf("node num = %d\n", node_surface[i][j][s].number);
					node_surface2[node_surface[i][j][s].number].edge_flag = 1;
					node_surface2[node_surface[i][j][s].number].none_flag = 0;
				}
				else{
					node_surface2[node_surface[i][j][s].number].edge_flag = 0;
					node_surface2[node_surface[i][j][s].number].none_flag = 1;
				}
			}
		}
		for (i = Hexa_h * 0.5; i <= Hexa_h; i++){
			for (j = 0; j < 2 * Hexa_w + Hexa_h * 0.5 - (i + 1); j++){
				if (i == Hexa_h || j == 0 || j == 2 * Hexa_w + Hexa_h * 0.5 - (i + 1) - 1){
					node_surface2[node_surface[i][j][s].number].edge_flag = 1;
					node_surface2[node_surface[i][j][s].number].none_flag = 0;
				}
				else{
					node_surface2[node_surface[i][j][s].number].edge_flag = 0;
					node_surface2[node_surface[i][j][s].number].none_flag = 1;
				}
			}
		}
	}
	/*for (i = 0; i < num_count; i++){
		if (node_surface2[i].edge_flag == 1){
			printf("node = %d\n", i);
		}
	}*/
#endif
#if 1
	for (s = 0; s < 2; s++){
		for (i = 0; i < Hexa_h * 0.5; i++){
			for (j = 0; j < 2 * Hexa_w - 1 - (Hexa_h * 0.5) + i; j++){
				for (h = 0; h < 3; h++){
					node_surface2[node_surface[i][j][s].number].pos.x[h] = node_surface[i][j][s].pos.x[h];
					//printf("node2 = %d, %f, %f, %f\n", node_surface[i][j][s].number, node_surface2[node_surface[i][j][s].number].pos.x[0], node_surface2[node_surface[i][j][s].number].pos.x[1], node_surface2[node_surface[i][j][s].number].pos.x[2]);
				}
			}
		}
		for (i = Hexa_h * 0.5; i <= Hexa_h; i++){
			for (j = 0; j < 2 * Hexa_w - 1 - (i - Hexa_h * 0.5); j++){
				for (h = 0; h < 3; h++){
					node_surface2[node_surface[i][j][s].number].pos.x[h] = node_surface[i][j][s].pos.x[h];
					//printf("node3 = %d, %f, %f, %f\n", node_surface[i][j][s].number, node_surface2[node_surface[i][j][s].number].pos.x[0], node_surface2[node_surface[i][j][s].number].pos.x[1], node_surface2[node_surface[i][j][s].number].pos.x[2]);
				}
			}
		}
	}

	for (s = 0; s < 2; s++){
		for (i = 0; i < Hexa_h * 0.5; i++){
			for (j = 0; j < 2 * Hexa_w - 1 - (Hexa_h * 0.5) + i; j++){
				if (i == 0 || j == 0 || j == 2 * (Hexa_w - 1) - (Hexa_h * 0.5) + i){
					//printf("node num = %d\n", node_surface[i][j][s].number);
					node_surface2[node_surface[i][j][s].number].edge_flag = 1;
					node_surface2[node_surface[i][j][s].number].none_flag = 0;
				}
				else{
					node_surface2[node_surface[i][j][s].number].edge_flag = 0;
					node_surface2[node_surface[i][j][s].number].none_flag = 1;
				}
			}
		}
		for (i = Hexa_h * 0.5; i <= Hexa_h; i++){
			for (j = 0; j < 2 * Hexa_w - 1 - (i - Hexa_h * 0.5); j++){
				if (i == Hexa_h || j == 0 || j == 2 * (Hexa_w - 1) - (i - Hexa_h * 0.5)){
					node_surface2[node_surface[i][j][s].number].edge_flag = 1;
					node_surface2[node_surface[i][j][s].number].none_flag = 0;
				}
				else{
					node_surface2[node_surface[i][j][s].number].edge_flag = 0;
					node_surface2[node_surface[i][j][s].number].none_flag = 1;
				}
			}
		}
	}
#endif
	//scanning neigbhor
#if 1
	i = 0;
	j = 0;
	s = 0;

	while (i < num_count * 0.5){
		while (j < num_count * 0.5){
			if (natural_length * 3 / 2 > sqrt(pow((node_surface2[i].pos.x[0] - node_surface2[j].pos.x[0]), 2) + pow((node_surface2[i].pos.x[1] - node_surface2[j].pos.x[1]), 2)
				+ pow((node_surface2[i].pos.x[2] - node_surface2[j].pos.x[2]), 2))){
				if (i != j){
					if (node_surface2[i].pos.x[0] < node_surface2[j].pos.x[0]){
						node_surface2[i].N[s] = j;
						s++;
						j++;
						continue;
					}
					else{
						j++;
						continue;
					}
				}
				else{
					j++;
					continue;

				}
			}
			else{
				j++;
				continue;
			}
		}
		node_surface2[i].index_n = s;
		i++;
		j = 0;
		s = 0;
		continue;
	}


	i = num_count * 0.5 - 1;
	j = num_count * 0.5 - 1;
	s = 0;
	while (i >= 0){
		while (j >= 0){
			if (natural_length * 3 / 2 > sqrt(pow((node_surface2[i].pos.x[0] - node_surface2[j].pos.x[0]), 2) + pow((node_surface2[i].pos.x[1] - node_surface2[j].pos.x[1]), 2)
				+ pow((node_surface2[i].pos.x[2] - node_surface2[j].pos.x[2]), 2))){
				if (i != j){
					if (node_surface2[i].pos.x[0] > node_surface2[j].pos.x[0]){
						s = node_surface2[i].index_n;
						node_surface2[i].N[s] = j;
						//printf("s1 = %d\n", s);
						s++;
						node_surface2[i].index_n = s;
						j--;
						continue;
					}
					else{
						j--;
						continue;
					}
				}
				else{
					j--;
					continue;
				}
			}
			else{
				j--;
				continue;
			}
		}
		i--;
		j = num_count * 0.5 - 1;
		continue;
	}


	for (i = ((num_count + Hexa_h) / 2 - 2 * Hexa_w + 1); i < num_count * 0.5; i++){
		for (j = 0; j < 6; j++){
			node_surface2[i].N[j] = NULL;
			node_surface2[i].index_n = 0;
		}
	}
	//printf("%d\n", num_count);
	//printf("i = %d\n", ((num_count + Hexa_h) / 2 - 2 * Hexa_w + 1));
	for (i = ((num_count + Hexa_h) / 2 - 2 * Hexa_w + 1); i < num_count * 0.5; i++){
		for (j = 0; j < num_count * 0.5; j++){
			if (i != j){
				if (natural_length * 3 * 0.5 > sqrt(pow((node_surface2[i].pos.x[0] - node_surface2[j].pos.x[0]), 2.0) + pow((node_surface2[i].pos.x[1] - node_surface2[j].pos.x[1]), 2.0) + pow((node_surface2[i].pos.x[2] - node_surface2[j].pos.x[2]), 2.0))){
					if (node_surface2[i].pos.x[0] > node_surface2[j].pos.x[0] && node_surface2[i].pos.x[1] - node_surface2[j].pos.x[1] == 0){
						//printf("i = %d, j = %d\n", i, j);
						node_surface2[i].N[0] = j;
					}
					if (node_surface2[i].pos.x[0] > node_surface2[j].pos.x[0] && node_surface2[i].pos.x[1] > node_surface2[j].pos.x[1]){
						//printf("i = %d, j = %d\n", i, j);
						node_surface2[i].N[1] = j;
					}
					if (node_surface2[i].pos.x[0] < node_surface2[j].pos.x[0] && node_surface2[i].pos.x[1] > node_surface2[j].pos.x[1]){
						node_surface2[i].N[2] = j;
					}
					if (node_surface2[i].pos.x[0] < node_surface2[j].pos.x[0] && node_surface2[i].pos.x[1] - node_surface2[j].pos.x[1] == 0){
						//printf("i = %d, j = %d\n", i, j);
						node_surface2[i].N[3] = j;
					}
				}
			}
		}
		if (i == ((num_count + Hexa_h) / 2 - 2 * Hexa_w + 1) || i == num_count * 0.5 - 1){
			//printf("i = %d\n", i);
			node_surface2[i].index_n = 3;
			node_surface2[i].index_n = 3;
		}
		else{
			node_surface2[i].index_n = 4;
		}
	}

	for (i = 0; i < num_count * 0.5; i++){
		for (j = 0; j < 6; j++){
		//printf("i = %d, j = %d, N = %d, index = %d\n", i, j, node_surface2[i].N[j], node_surface2[i].index_n);
		}
	}

	i = 0;
	j = 0;
	s = 0;
	while (i < num_count * 0.5){
		while (j < tri_count * 0.5){
			if (triangle_data[j].t[0] == i || triangle_data[j].t[1] == i || triangle_data[j].t[2] == i){
				node_surface2[i].T[s] = j;
				s++;
				j++;
				continue;
			}
			else{
				j++;
				continue;
			}
		}
		node_surface2[i].index_t = s;
		i++;
		j = 0;
		s = 0;
		continue;
	}
	for (i = 0; i < num_count * 0.5; i++){
		for (j = 0; j < 6; j++){
			//printf("i = %d, j = %d, T = %d, index_t = %d\n", i, j, node_surface2[i].T[j], node_surface2[i].index_t);
		}
	}
#endif
	for (i = 0; i <= num_count - 1; i++){
		for (j = 0; j <= num_count - 1; j++){
			edge[i][j].torf = false;
		}
	}
	for (i = 0; i < tri_count; i++){
		edge[triangle_data[i].t[0]][triangle_data[i].t[1]].torf = true;
		edge[triangle_data[i].t[1]][triangle_data[i].t[0]].torf = true;
		edge[triangle_data[i].t[1]][triangle_data[i].t[2]].torf = true;
		edge[triangle_data[i].t[2]][triangle_data[i].t[1]].torf = true;
		edge[triangle_data[i].t[2]][triangle_data[i].t[0]].torf = true;
		edge[triangle_data[i].t[0]][triangle_data[i].t[2]].torf = true;
		edge[triangle_data[i].t[0]][triangle_data[i].t[1]].len = sqrt(
			pow((node_surface2[triangle_data[i].t[0]].pos.x[0] - node_surface2[triangle_data[i].t[1]].pos.x[0]), 2.0) +
			pow((node_surface2[triangle_data[i].t[0]].pos.x[1] - node_surface2[triangle_data[i].t[1]].pos.x[1]), 2.0) +
			pow((node_surface2[triangle_data[i].t[0]].pos.x[2] - node_surface2[triangle_data[i].t[1]].pos.x[2]), 2.0));
		edge[triangle_data[i].t[1]][triangle_data[i].t[0]].len = sqrt(
			pow((node_surface2[triangle_data[i].t[1]].pos.x[0] - node_surface2[triangle_data[i].t[0]].pos.x[0]), 2.0) +
			pow((node_surface2[triangle_data[i].t[1]].pos.x[1] - node_surface2[triangle_data[i].t[0]].pos.x[1]), 2.0) +
			pow((node_surface2[triangle_data[i].t[1]].pos.x[2] - node_surface2[triangle_data[i].t[0]].pos.x[2]), 2.0));
		edge[triangle_data[i].t[1]][triangle_data[i].t[2]].len = sqrt(
			pow((node_surface2[triangle_data[i].t[1]].pos.x[0] - node_surface2[triangle_data[i].t[2]].pos.x[0]), 2.0) +
			pow((node_surface2[triangle_data[i].t[1]].pos.x[1] - node_surface2[triangle_data[i].t[2]].pos.x[1]), 2.0) +
			pow((node_surface2[triangle_data[i].t[1]].pos.x[2] - node_surface2[triangle_data[i].t[2]].pos.x[2]), 2.0));
		edge[triangle_data[i].t[2]][triangle_data[i].t[1]].len = sqrt(
			pow((node_surface2[triangle_data[i].t[2]].pos.x[0] - node_surface2[triangle_data[i].t[1]].pos.x[0]), 2.0) +
			pow((node_surface2[triangle_data[i].t[2]].pos.x[1] - node_surface2[triangle_data[i].t[1]].pos.x[1]), 2.0) +
			pow((node_surface2[triangle_data[i].t[2]].pos.x[2] - node_surface2[triangle_data[i].t[1]].pos.x[2]), 2.0));
		edge[triangle_data[i].t[2]][triangle_data[i].t[0]].len = sqrt(
			pow((node_surface2[triangle_data[i].t[2]].pos.x[0] - node_surface2[triangle_data[i].t[0]].pos.x[0]), 2.0) +
			pow((node_surface2[triangle_data[i].t[2]].pos.x[1] - node_surface2[triangle_data[i].t[0]].pos.x[1]), 2.0) +
			pow((node_surface2[triangle_data[i].t[2]].pos.x[2] - node_surface2[triangle_data[i].t[0]].pos.x[2]), 2.0));
		edge[triangle_data[i].t[0]][triangle_data[i].t[2]].len = sqrt(
			pow((node_surface2[triangle_data[i].t[0]].pos.x[0] - node_surface2[triangle_data[i].t[2]].pos.x[0]), 2.0) +
			pow((node_surface2[triangle_data[i].t[0]].pos.x[1] - node_surface2[triangle_data[i].t[2]].pos.x[1]), 2.0) +
			pow((node_surface2[triangle_data[i].t[0]].pos.x[2] - node_surface2[triangle_data[i].t[2]].pos.x[2]), 2.0));
	}
}
void node_simulation(int view_con){
	if (first_count == 1){
		initiation();
		first_count--;
	}
	int i = 0;
	int j = 0;
	int k = 0;
	int l = 0;
	int h = 0;
	int s = 0;
	int E = i + num_count * 0.5;
	double min = 100000.0;
	double max = -10.0;
	double min_g = 0.00375;
	double max_g = 0.00875;
	double mass = 1000;
	double kizami = 0.01;
	double total_force[3] = { 0.0, 0.0, 0.0 };
	double temp_len;

	double normal_force[3];
	double normal_temp[9];
	double normal_temp3[3];

	for (i = 0; i <= num_count - 1; i++){
		for (j = 0; j < 3; j++){
			node_surface2[i].acc.x[j] = 0.0;
		}
	}

	for (i = 0; i <= num_count - 1; i++){
		node_surface2[i].color_grad = 0.0;
		for (j = 0; j < num_count; j++){
			if (edge[i][j].torf == 1){
				//printf("node_surface2[%d] = %f, %f, %f\n", i, node_surface2[i].pos.x[0], node_surface2[i].pos.x[1], node_surface2[i].pos.x[2]);
				temp_len = sqrt(pow((node_surface2[i].pos.x[0] - node_surface2[j].pos.x[0]), 2.0) + pow((node_surface2[i].pos.x[1] - node_surface2[j].pos.x[1]), 2.0) + pow((node_surface2[i].pos.x[2] - node_surface2[j].pos.x[2]), 2.0));
				//printf("%d %d %lf %lf\n", i, j, temp_len, edge[i][j].len);
				//printf("a");
				if (temp_len > edge[i][j].len){
				for (k = 0; k < 3; k++){
					node_surface2[i].color_grad += fabs(-1000.0 * damp_k * (node_surface2[i].pos.x[k] - node_surface2[j].pos.x[k]) * (temp_len - edge[i][j].len) / temp_len / mass);
					node_surface2[i].acc.x[k] += -1000.0 * damp_k * (node_surface2[i].pos.x[k] - node_surface2[j].pos.x[k]) * (temp_len - edge[i][j].len) / temp_len / mass;
					//printf("%lf\n", -1.0 * damp_k * (node_surface2[i].pos.x[k] - node_surface2[j].pos.x[k]) * (temp_len - edge[i][j].len) / temp_len / mass);
					}
				}
			}
		}
	}

		if (open_flag == true){
			wall_z += 0.01;
		}
		if (close_flag == true){
			wall_z -= 0.01;
		}
		for (i = 0; i < num_count; i++){
			for (k = 0; k < 3; k++){
				node_surface2[i].acc.x[k] += -dv * node_surface2[i].del_pos.x[k];
			}
			/*	else if (node_surface2[i].edge_flag == 1){
			for (k = 0; k < 2; k++){
			node_surface2[i].acc.x[k] += -dv * node_surface2[i].del_pos.x[k];
			}
			}*/
		}
		//printf("node_surface2 = %f, %f\n", node_surface2[1].del_pos.x[2], node_surface2[56].del_pos.x[2]);
		//pressure (external first noraml_temp3 will be b in a form of a x b)

#if 1

		for (i = 0; i < tri_count; i++){
			if (i < tri_count * 0.5){
				if (triangle_data[i].color == 1){
					//printf("a");
					for (j = 0; j < 3; j++){
						normal_temp3[j] = node_surface2[triangle_data[i].t[1]].pos.x[j] - node_surface2[triangle_data[i].t[0]].pos.x[j];
					}
					gaiseki_9_3(normal_temp, normal_temp3);
					for (j = 0; j < 3; j++){
						normal_temp3[j] = node_surface2[triangle_data[i].t[2]].pos.x[j] - node_surface2[triangle_data[i].t[0]].pos.x[j];
					}
					mult_matrix3x1(normal_force, normal_temp, normal_temp3);
				}
				else if (triangle_data[i].color == 2){
					for (j = 0; j < 3; j++){
						normal_temp3[j] = node_surface2[triangle_data[i].t[2]].pos.x[j] - node_surface2[triangle_data[i].t[0]].pos.x[j];
					}
					gaiseki_9_3(normal_temp, normal_temp3);
					for (j = 0; j < 3; j++){
						normal_temp3[j] = node_surface2[triangle_data[i].t[1]].pos.x[j] - node_surface2[triangle_data[i].t[0]].pos.x[j];
					}
					mult_matrix3x1(normal_force, normal_temp, normal_temp3);
				}
			}
			else{
				if (triangle_data[i].color == 1){
					//printf("a");
					for (j = 0; j < 3; j++){
						normal_temp3[j] = node_surface2[triangle_data[i].t[2]].pos.x[j] - node_surface2[triangle_data[i].t[0]].pos.x[j];
					}
					gaiseki_9_3(normal_temp, normal_temp3);
					for (j = 0; j < 3; j++){
						normal_temp3[j] = node_surface2[triangle_data[i].t[1]].pos.x[j] - node_surface2[triangle_data[i].t[0]].pos.x[j];
					}
					mult_matrix3x1(normal_force, normal_temp, normal_temp3);
				}
				else if (triangle_data[i].color == 2){
					for (j = 0; j < 3; j++){
						normal_temp3[j] = node_surface2[triangle_data[i].t[1]].pos.x[j] - node_surface2[triangle_data[i].t[0]].pos.x[j];
					}
					gaiseki_9_3(normal_temp, normal_temp3);
					for (j = 0; j < 3; j++){
						normal_temp3[j] = node_surface2[triangle_data[i].t[2]].pos.x[j] - node_surface2[triangle_data[i].t[0]].pos.x[j];
					}
					mult_matrix3x1(normal_force, normal_temp, normal_temp3);
				}
			}
			for (j = 0; j < 3; j++){
				triangle_data[i].normal[j] = normal_force[j] / sqrt(pow(normal_force[0], 2.0) + pow(normal_force[1], 2.0) + pow(normal_force[2], 2.0));// / normal_force[j];
				for (k = 0; k < 3; k++){
					node_surface2[triangle_data[i].t[j]].acc.x[k] += 0.1 * damp_k_normal * normal_force[k];
				}
			}
		}
#endif
#if 0
	for (i = 0; i < num_count; i++){
		for (k = 0; k < 3; k++){
			if (node_surface2[i].edge_flag == 1){
				for (k = 0; k < 3; k++){
					//printf("i = %d\n", i);
					if (i < num_count * 0.5){
						for (k = 0; k < 3; k++){
							int j = i + num_count * 0.5;
							node_surface2[i].pos.x[k] = node_surface2[j].pos.x[k];
							//	printf("i = %d,j = %d\n", i, j);
						}
						node_surface2[i].del_pos.x[k] = node_surface2[i].del_pos.x[k] + node_surface2[i].acc.x[k] * kizami;
						node_surface2[i].pos.x[k] = node_surface2[i].pos.x[k] + node_surface2[i].del_pos.x[k] * kizami + 1.0 / 2.0 * node_surface2[i].acc.x[k] * kizami * kizami;
					}
				}
			}
			else{
				node_surface2[i].del_pos.x[k] = node_surface2[i].del_pos.x[k] + node_surface2[i].acc.x[k] * kizami;
				node_surface2[i].pos.x[k] = node_surface2[i].pos.x[k] + node_surface2[i].del_pos.x[k] * kizami + 1.0 / 2.0 * node_surface2[i].acc.x[k] * kizami * kizami;
			}
		}
	}
#endif
#if 1
	//printf("node_surface_z = %f, %f\n", node_surface2[0].acc.x[1], node_surface2[55].acc.x[1]);
	for (i = 0; i < num_count; i++){
		if (node_surface2[i].none_flag == 1){
			for (k = 0; k < 3; k++){
				node_surface2[i].del_pos.x[k] = node_surface2[i].del_pos.x[k] + node_surface2[i].acc.x[k] * kizami;
				node_surface2[i].pos.x[k] = node_surface2[i].pos.x[k] + node_surface2[i].del_pos.x[k] * kizami + 1.0 / 2.0 * node_surface2[i].acc.x[k] * kizami * kizami;
			}
		}
	}
	for (i = 0; i < num_count * 0.5; i++){
		if (node_surface2[i].edge_flag == 1){
			for (k = 0; k < 3; k++){
				int j = i + num_count * 0.5;
				//printf("j = %d, i = %d\n", j, i);
				//node_surface2[i].pos.x[k] = node_surface2[j].pos.x[k];
				node_surface2[i].del_pos.x[k] = node_surface2[i].del_pos.x[k] + (node_surface2[i].acc.x[k] * 0.5 + node_surface2[j].acc.x[k] * 0.5) * kizami;
				node_surface2[j].del_pos.x[k] = node_surface2[j].del_pos.x[k] + (node_surface2[i].acc.x[k] * 0.5 + node_surface2[j].acc.x[k] * 0.5) * kizami;
				node_surface2[i].pos.x[k] = node_surface2[i].pos.x[k] + node_surface2[i].del_pos.x[k] * kizami + 1.0 / 4.0 * node_surface2[i].acc.x[k] * kizami * kizami + node_surface2[j].del_pos.x[k] * kizami + 1.0 / 4.0 * node_surface2[j].acc.x[k] * kizami * kizami;
				node_surface2[j].pos.x[k] = node_surface2[j].pos.x[k] + node_surface2[i].del_pos.x[k] * kizami + 1.0 / 4.0 * node_surface2[i].acc.x[k] * kizami * kizami + node_surface2[j].del_pos.x[k] * kizami + 1.0 / 4.0 * node_surface2[j].acc.x[k] * kizami * kizami;
			}
		}
	}
#endif

#if 1
	double P0P1[3];
	double P0P2[3];
	double P3P2[3];
	double P3P1[3];

	/////////////////////////////cosa cosb
	//index_n = 6
	for (i = 0; i < num_count * 0.5; i++){
		if (node_surface2[i].index_n == 6){
			for (j = 1; j < 5; j++){
				P0P1[0] = (node_surface2[node_surface2[i].N[j]].pos.x[0] - node_surface2[node_surface2[i].N[j - 1]].pos.x[0]);
				P0P1[1] = (node_surface2[node_surface2[i].N[j]].pos.x[1] - node_surface2[node_surface2[i].N[j - 1]].pos.x[1]);
				P0P1[2] = (node_surface2[node_surface2[i].N[j]].pos.x[2] - node_surface2[node_surface2[i].N[j - 1]].pos.x[2]);
				P0P2[0] = (node_surface2[i].pos.x[0] - node_surface2[node_surface2[i].N[j - 1]].pos.x[0]);
				P0P2[1] = (node_surface2[i].pos.x[1] - node_surface2[node_surface2[i].N[j - 1]].pos.x[1]);
				P0P2[2] = (node_surface2[i].pos.x[2] - node_surface2[node_surface2[i].N[j - 1]].pos.x[2]);
				P3P2[0] = (node_surface2[i].pos.x[0] - node_surface2[node_surface2[i].N[j + 1]].pos.x[0]);
				P3P2[1] = (node_surface2[i].pos.x[1] - node_surface2[node_surface2[i].N[j + 1]].pos.x[1]);
				P3P2[2] = (node_surface2[i].pos.x[2] - node_surface2[node_surface2[i].N[j + 1]].pos.x[2]);
				P3P1[0] = (node_surface2[node_surface2[i].N[j]].pos.x[0] - node_surface2[node_surface2[i].N[j + 1]].pos.x[0]);
				P3P1[1] = (node_surface2[node_surface2[i].N[j]].pos.x[1] - node_surface2[node_surface2[i].N[j + 1]].pos.x[1]);
				P3P1[2] = (node_surface2[node_surface2[i].N[j]].pos.x[2] - node_surface2[node_surface2[i].N[j + 1]].pos.x[2]);
				node_surface2[i].cosa[j] = ((P0P1[0] * P0P2[0]) + (P0P1[1] * P0P2[1]) + (P0P1[2] * P0P2[2])) / (sqrt(pow(P0P1[0], 2) + pow(P0P1[1], 2) + pow(P0P1[2], 2)) * sqrt(pow(P0P2[0], 2) + pow(P0P2[1], 2) + pow(P0P2[2], 2)));
				node_surface2[i].cosb[j] = ((P3P2[0] * P3P1[0]) + (P3P2[1] * P3P1[1]) + (P3P2[2] * P3P1[2])) / (sqrt(pow(P3P2[0], 2) + pow(P3P2[1], 2) + pow(P3P2[2], 2)) * sqrt(pow(P3P1[0], 2) + pow(P3P1[1], 2) + pow(P3P1[2], 2)));
				//printf("i = %d, j = %d, cosa = %lf, cosb = %lf\n", i, j, node_surface2[i].cosa[j], node_surface2[i].cosb[j]);
			}
			for (j = 0; j < 1; j++){
				P0P1[0] = (node_surface2[node_surface2[i].N[j]].pos.x[0] - node_surface2[node_surface2[i].N[5]].pos.x[0]);
				P0P1[1] = (node_surface2[node_surface2[i].N[j]].pos.x[1] - node_surface2[node_surface2[i].N[5]].pos.x[1]);
				P0P1[2] = (node_surface2[node_surface2[i].N[j]].pos.x[2] - node_surface2[node_surface2[i].N[5]].pos.x[2]);
				P0P2[0] = (node_surface2[i].pos.x[0] - node_surface2[node_surface2[i].N[5]].pos.x[0]);
				P0P2[1] = (node_surface2[i].pos.x[1] - node_surface2[node_surface2[i].N[5]].pos.x[1]);
				P0P2[2] = (node_surface2[i].pos.x[2] - node_surface2[node_surface2[i].N[5]].pos.x[2]);
				P3P2[0] = (node_surface2[i].pos.x[0] - node_surface2[node_surface2[i].N[j + 1]].pos.x[0]);
				P3P2[1] = (node_surface2[i].pos.x[1] - node_surface2[node_surface2[i].N[j + 1]].pos.x[1]);
				P3P2[2] = (node_surface2[i].pos.x[2] - node_surface2[node_surface2[i].N[j + 1]].pos.x[2]);
				P3P1[0] = (node_surface2[node_surface2[i].N[j]].pos.x[0] - node_surface2[node_surface2[i].N[j + 1]].pos.x[0]);
				P3P1[1] = (node_surface2[node_surface2[i].N[j]].pos.x[1] - node_surface2[node_surface2[i].N[j + 1]].pos.x[1]);
				P3P1[2] = (node_surface2[node_surface2[i].N[j]].pos.x[2] - node_surface2[node_surface2[i].N[j + 1]].pos.x[2]);
				node_surface2[i].cosa[j] = ((P0P1[0] * P0P2[0]) + (P0P1[1] * P0P2[1]) + (P0P1[2] * P0P2[2])) / (sqrt(pow(P0P1[0], 2) + pow(P0P1[1], 2) + pow(P0P1[2], 2)) * sqrt(pow(P0P2[0], 2) + pow(P0P2[1], 2) + pow(P0P2[2], 2)));
				node_surface2[i].cosb[j] = ((P3P2[0] * P3P1[0]) + (P3P2[1] * P3P1[1]) + (P3P2[2] * P3P1[2])) / (sqrt(pow(P3P2[0], 2) + pow(P3P2[1], 2) + pow(P3P2[2], 2)) * sqrt(pow(P3P1[0], 2) + pow(P3P1[1], 2) + pow(P3P1[2], 2)));
			}
			for (j = 5; j < 6; j++){
				for (k = 0; k < 3; k++){
					P0P1[k] = (node_surface2[node_surface2[i].N[j]].pos.x[k] - node_surface2[node_surface2[i].N[j - 1]].pos.x[k]);
					P0P2[k] = (node_surface2[i].pos.x[k] - node_surface2[node_surface2[i].N[j - 1]].pos.x[k]);
					P3P2[k] = (node_surface2[i].pos.x[k] - node_surface2[node_surface2[i].N[0]].pos.x[k]);
					P3P1[k] = (node_surface2[node_surface2[i].N[j]].pos.x[k] - node_surface2[node_surface2[i].N[0]].pos.x[k]);
					node_surface2[i].cosa[j] = ((P0P1[0] * P0P2[0]) + (P0P1[1] * P0P2[1]) + (P0P1[2] * P0P2[2])) / (sqrt(pow(P0P1[0], 2) + pow(P0P1[1], 2) + pow(P0P1[2], 2)) * sqrt(pow(P0P2[0], 2) + pow(P0P2[1], 2) + pow(P0P2[2], 2)));
					node_surface2[i].cosb[j] = ((P3P2[0] * P3P1[0]) + (P3P2[1] * P3P1[1]) + (P3P2[2] * P3P1[2])) / (sqrt(pow(P3P2[0], 2) + pow(P3P2[1], 2) + pow(P3P2[2], 2)) * sqrt(pow(P3P1[0], 2) + pow(P3P1[1], 2) + pow(P3P1[2], 2)));
				}
			}
		}
	}
	//index_n = 4
	for (i = 0; i < num_count * 0.5; i++){
		if (node_surface2[i].index_n == 4){
			for (j = 1; j < 3; j++){
				for (k = 0; k < 3; k++){
					P0P1[k] = (node_surface2[node_surface2[i].N[j]].pos.x[k] - node_surface2[node_surface2[i].N[j - 1]].pos.x[k]);
					P0P2[k] = (node_surface2[i].pos.x[k] - node_surface2[node_surface2[i].N[j - 1]].pos.x[k]);
					P3P2[k] = (node_surface2[i].pos.x[k] - node_surface2[node_surface2[i].N[j + 1]].pos.x[k]);
					P3P1[k] = (node_surface2[node_surface2[i].N[j]].pos.x[k] - node_surface2[node_surface2[i].N[j + 1]].pos.x[k]);
					node_surface2[i].cosa[j] = ((P0P1[0] * P0P2[0]) + (P0P1[1] * P0P2[1]) + (P0P1[2] * P0P2[2])) / (sqrt(pow(P0P1[0], 2) + pow(P0P1[1], 2) + pow(P0P1[2], 2)) * sqrt(pow(P0P2[0], 2) + pow(P0P2[1], 2) + pow(P0P2[2], 2)));
					node_surface2[i].cosb[j] = ((P3P2[0] * P3P1[0]) + (P3P2[1] * P3P1[1]) + (P3P2[2] * P3P1[2])) / (sqrt(pow(P3P2[0], 2) + pow(P3P2[1], 2) + pow(P3P2[2], 2)) * sqrt(pow(P3P1[0], 2) + pow(P3P1[1], 2) + pow(P3P1[2], 2)));
					//printf("i = %d, j = %d, cosa = %lf, cosb = %lf\n", i, j, node_surface2[i].cosa[j], node_surface2[i].cosb[j]);
				}
			}
		}
	}
	//inde_n = 3
	for (i = 0; i < num_count * 0.5; i++){
		if (i != ((num_count + Hexa_h) / 2 - 2 * Hexa_w + 1)){
			if (node_surface2[i].index_n == 3){
				for (j = 1; j < 2; j++){
					for (k = 0; k < 3; k++){
						P0P1[k] = (node_surface2[node_surface2[i].N[j]].pos.x[k] - node_surface2[node_surface2[i].N[j - 1]].pos.x[k]);
						P0P2[k] = (node_surface2[i].pos.x[k] - node_surface2[node_surface2[i].N[j - 1]].pos.x[k]);
						P3P2[k] = (node_surface2[i].pos.x[k] - node_surface2[node_surface2[i].N[j + 1]].pos.x[k]);
						P3P1[k] = (node_surface2[node_surface2[i].N[j]].pos.x[k] - node_surface2[node_surface2[i].N[j + 1]].pos.x[k]);
						node_surface2[i].cosa[j] = ((P0P1[0] * P0P2[0]) + (P0P1[1] * P0P2[1]) + (P0P1[2] * P0P2[2])) / (sqrt(pow(P0P1[0], 2) + pow(P0P1[1], 2) + pow(P0P1[2], 2)) * sqrt(pow(P0P2[0], 2) + pow(P0P2[1], 2) + pow(P0P2[2], 2)));
						node_surface2[i].cosb[j] = ((P3P2[0] * P3P1[0]) + (P3P2[1] * P3P1[1]) + (P3P2[2] * P3P1[2])) / (sqrt(pow(P3P2[0], 2) + pow(P3P2[1], 2) + pow(P3P2[2], 2)) * sqrt(pow(P3P1[0], 2) + pow(P3P1[1], 2) + pow(P3P1[2], 2)));
						//printf("i = %d, j = %d, cosa = %lf, cosb = %lf\n", i, j, node_surface2[i].cosa[j], node_surface2[i].cosb[j]);
					}
				}
			}
		}
		if (i == ((num_count + Hexa_h) / 2 - 2 * Hexa_w + 1)){
			for (j = 2; j < 3; j++){
				for (k = 0; k < 3; k++){
					P0P1[k] = (node_surface2[node_surface2[i].N[j]].pos.x[k] - node_surface2[node_surface2[i].N[j - 1]].pos.x[k]);
					P0P2[k] = (node_surface2[i].pos.x[k] - node_surface2[node_surface2[i].N[j - 1]].pos.x[k]);
					P3P2[k] = (node_surface2[i].pos.x[k] - node_surface2[node_surface2[i].N[j + 1]].pos.x[k]);
					P3P1[k] = (node_surface2[node_surface2[i].N[j]].pos.x[k] - node_surface2[node_surface2[i].N[j + 1]].pos.x[k]);
					node_surface2[i].cosa[j] = ((P0P1[0] * P0P2[0]) + (P0P1[1] * P0P2[1]) + (P0P1[2] * P0P2[2])) / (sqrt(pow(P0P1[0], 2) + pow(P0P1[1], 2) + pow(P0P1[2], 2)) * sqrt(pow(P0P2[0], 2) + pow(P0P2[1], 2) + pow(P0P2[2], 2)));
					node_surface2[i].cosb[j] = ((P3P2[0] * P3P1[0]) + (P3P2[1] * P3P1[1]) + (P3P2[2] * P3P1[2])) / (sqrt(pow(P3P2[0], 2) + pow(P3P2[1], 2) + pow(P3P2[2], 2)) * sqrt(pow(P3P1[0], 2) + pow(P3P1[1], 2) + pow(P3P1[2], 2)));
				}
			}
		}
	}
	//////////////////////////////cota cotb
	//inde_n = 6
	for (i = 0; i < num_count * 0.5; i++){
		if (node_surface2[i].index_n == 6){
			for (j = 0; j < 6; j++){
				node_surface2[i].cota[j] = node_surface2[i].cosa[j] / sqrt(1 - pow(node_surface2[i].cosa[j], 2));
				node_surface2[i].cotb[j] = node_surface2[i].cosb[j] / sqrt(1 - pow(node_surface2[i].cosb[j], 2));
			}
		}
		if (node_surface2[i].index_n == 4){
			for (j = 1; j < 3; j++){
				node_surface2[i].cota[j] = node_surface2[i].cosa[j] / sqrt(1 - pow(node_surface2[i].cosa[j], 2));
				node_surface2[i].cotb[j] = node_surface2[i].cosb[j] / sqrt(1 - pow(node_surface2[i].cosb[j], 2));
			}
		}
		if (node_surface2[i].index_n == 3){
			if (i != ((num_count + Hexa_h) / 2 - 2 * Hexa_w + 1)){
				for (j = 1; j < 2; j++){
					node_surface2[i].cota[j] = node_surface2[i].cosa[j] / sqrt(1 - pow(node_surface2[i].cosa[j], 2));
					node_surface2[i].cotb[j] = node_surface2[i].cosb[j] / sqrt(1 - pow(node_surface2[i].cosb[j], 2));
					//printf("i = %d, j = %d, cota = %lf, cotb = %lf\n", i, j, node_surface2[i].cota[j], node_surface2[i].cotb[j]);
				}
			}
			if (i == ((num_count + Hexa_h) / 2 - 2 * Hexa_w + 1)){
				for (j = 2; j < 3; j++){
					node_surface2[i].cota[j] = node_surface2[i].cosa[j] / sqrt(1 - pow(node_surface2[i].cosa[j], 2));
					node_surface2[i].cotb[j] = node_surface2[i].cosb[j] / sqrt(1 - pow(node_surface2[i].cosb[j], 2));
				}
			}
		}
	}
	///////////////////////////////////////////
	//normal
	///////////////////////////////////////////
	for (i = 0; i < num_count * 0.5; i++){
		if (node_surface2[i].index_t == 6){
			node_surface2[i].n_normal[0] = triangle_data[node_surface2[i].T[0]].normal[0] + triangle_data[node_surface2[i].T[1]].normal[0] + triangle_data[node_surface2[i].T[2]].normal[0] +
				triangle_data[node_surface2[i].T[3]].normal[0] + triangle_data[node_surface2[i].T[4]].normal[0] + triangle_data[node_surface2[i].T[5]].normal[0];
			node_surface2[i].n_normal[1] = triangle_data[node_surface2[i].T[0]].normal[1] + triangle_data[node_surface2[i].T[1]].normal[1] + triangle_data[node_surface2[i].T[2]].normal[1] +
				triangle_data[node_surface2[i].T[3]].normal[1] + triangle_data[node_surface2[i].T[4]].normal[1] + triangle_data[node_surface2[i].T[5]].normal[1];
			node_surface2[i].n_normal[2] = triangle_data[node_surface2[i].T[0]].normal[2] + triangle_data[node_surface2[i].T[1]].normal[2] + triangle_data[node_surface2[i].T[2]].normal[2] +
				triangle_data[node_surface2[i].T[3]].normal[2] + triangle_data[node_surface2[i].T[4]].normal[2] + triangle_data[node_surface2[i].T[5]].normal[2];
		}
		if (node_surface2[i].index_t == 3){
			node_surface2[i].n_normal[0] = triangle_data[node_surface2[i].T[0]].normal[0] + triangle_data[node_surface2[i].T[1]].normal[0] + triangle_data[node_surface2[i].T[2]].normal[0];
			node_surface2[i].n_normal[1] = triangle_data[node_surface2[i].T[0]].normal[1] + triangle_data[node_surface2[i].T[1]].normal[1] + triangle_data[node_surface2[i].T[2]].normal[1];
			node_surface2[i].n_normal[2] = triangle_data[node_surface2[i].T[0]].normal[2] + triangle_data[node_surface2[i].T[1]].normal[2] + triangle_data[node_surface2[i].T[2]].normal[2];
		}
		if (node_surface2[i].index_t == 2){
			if (i != ((num_count + Hexa_h) / 2 - 2 * Hexa_w + 1)){
				node_surface2[i].n_normal[0] = triangle_data[node_surface2[i].T[0]].normal[0] + triangle_data[node_surface2[i].T[1]].normal[0];
				node_surface2[i].n_normal[1] = triangle_data[node_surface2[i].T[0]].normal[1] + triangle_data[node_surface2[i].T[1]].normal[1];
				node_surface2[i].n_normal[2] = triangle_data[node_surface2[i].T[0]].normal[2] + triangle_data[node_surface2[i].T[1]].normal[2];

			}
			if (i == ((num_count + Hexa_h) / 2 - 2 * Hexa_w + 1)){
				node_surface2[i].n_normal[0] = triangle_data[node_surface2[i].T[1]].normal[0] + triangle_data[node_surface2[i].T[2]].normal[0];
				node_surface2[i].n_normal[1] = triangle_data[node_surface2[i].T[1]].normal[1] + triangle_data[node_surface2[i].T[2]].normal[1];
				node_surface2[i].n_normal[2] = triangle_data[node_surface2[i].T[1]].normal[2] + triangle_data[node_surface2[i].T[2]].normal[2];
			}
		}
	}

	for (i = 0; i < num_count * 0.5; i++){
			node_surface2[i].m_normal[0] = node_surface2[i].n_normal[0] / sqrt(pow(node_surface2[i].n_normal[0], 2) + pow(node_surface2[i].n_normal[1], 2) + pow(node_surface2[i].n_normal[2], 2));
			node_surface2[i].m_normal[1] = node_surface2[i].n_normal[1] / sqrt(pow(node_surface2[i].n_normal[0], 2) + pow(node_surface2[i].n_normal[1], 2) + pow(node_surface2[i].n_normal[2], 2));
			node_surface2[i].m_normal[2] = node_surface2[i].n_normal[2] / sqrt(pow(node_surface2[i].n_normal[0], 2) + pow(node_surface2[i].n_normal[1], 2) + pow(node_surface2[i].n_normal[2], 2));
	}
	///////////////////////////////////////////
	//K
	///////////////////////////////////////////
	
	for (i = 0; i < num_count * 0.5; i++){
		if (node_surface2[i].index_t == 6){
			node_surface2[i].K = (0.25 * ((node_surface2[i].cota[0] + node_surface2[i].cotb[0]) *((node_surface2[i].pos.x[0] - node_surface2[node_surface2[i].N[0]].pos.x[0]) * node_surface2[i].m_normal[0]
				+ (node_surface2[i].pos.x[1] - node_surface2[node_surface2[i].N[0]].pos.x[1]) * node_surface2[i].m_normal[1] + (node_surface2[i].pos.x[2] - node_surface2[node_surface2[i].N[0]].pos.x[2]) * node_surface2[i].m_normal[2]) 
				+ (node_surface2[i].cota[1] + node_surface2[i].cotb[1]) * ((node_surface2[i].pos.x[0] - node_surface2[node_surface2[i].N[1]].pos.x[0]) * node_surface2[i].m_normal[0]
				+ (node_surface2[i].pos.x[1] - node_surface2[node_surface2[i].N[1]].pos.x[1]) * node_surface2[i].m_normal[1] + (node_surface2[i].pos.x[2] - node_surface2[node_surface2[i].N[1]].pos.x[2]) * node_surface2[i].m_normal[2])
				+ (node_surface2[i].cota[2] + node_surface2[i].cotb[2]) * ((node_surface2[i].pos.x[0] - node_surface2[node_surface2[i].N[2]].pos.x[0]) * node_surface2[i].m_normal[0]
				+ (node_surface2[i].pos.x[1] - node_surface2[node_surface2[i].N[2]].pos.x[1]) * node_surface2[i].m_normal[1] + (node_surface2[i].pos.x[2] - node_surface2[node_surface2[i].N[2]].pos.x[2]) * node_surface2[i].m_normal[2])
				+ (node_surface2[i].cota[3] + node_surface2[i].cotb[3]) * ((node_surface2[i].pos.x[0] - node_surface2[node_surface2[i].N[3]].pos.x[0]) * node_surface2[i].m_normal[0]
				+ (node_surface2[i].pos.x[1] - node_surface2[node_surface2[i].N[3]].pos.x[1]) * node_surface2[i].m_normal[1] + (node_surface2[i].pos.x[2] - node_surface2[node_surface2[i].N[3]].pos.x[2]) * node_surface2[i].m_normal[2])
				+ (node_surface2[i].cota[4] + node_surface2[i].cotb[4]) * ((node_surface2[i].pos.x[0] - node_surface2[node_surface2[i].N[4]].pos.x[0]) * node_surface2[i].m_normal[0]
				+ (node_surface2[i].pos.x[1] - node_surface2[node_surface2[i].N[4]].pos.x[1]) * node_surface2[i].m_normal[1] + (node_surface2[i].pos.x[2] - node_surface2[node_surface2[i].N[4]].pos.x[2]) * node_surface2[i].m_normal[2])
				+ (node_surface2[i].cota[5] + node_surface2[i].cotb[5]) * ((node_surface2[i].pos.x[0] - node_surface2[node_surface2[i].N[5]].pos.x[0]) * node_surface2[i].m_normal[0]
				+ (node_surface2[i].pos.x[1] - node_surface2[node_surface2[i].N[5]].pos.x[1]) * node_surface2[i].m_normal[1] + (node_surface2[i].pos.x[2] - node_surface2[node_surface2[i].N[5]].pos.x[2]) * node_surface2[i].m_normal[2]))) / 6;
		
				//printf("i = %d, K = %lf\n", i, node_surface2[i].K);
		}
		if (node_surface2[i].index_t == 3){
			node_surface2[i].K = (0.25 * ((node_surface2[i].cota[1] + node_surface2[i].cotb[1]) *((node_surface2[i].pos.x[0] - node_surface2[node_surface2[i].N[1]].pos.x[0]) * node_surface2[i].m_normal[0]
				+ (node_surface2[i].pos.x[1] - node_surface2[node_surface2[i].N[1]].pos.x[1]) * node_surface2[i].m_normal[1] + (node_surface2[i].pos.x[2] - node_surface2[node_surface2[i].N[1]].pos.x[2]) * node_surface2[i].m_normal[2])
				+ (node_surface2[i].cota[2] + node_surface2[i].cotb[2]) * ((node_surface2[i].pos.x[0] - node_surface2[node_surface2[i].N[2]].pos.x[0]) * node_surface2[i].m_normal[0]
				+ (node_surface2[i].pos.x[1] - node_surface2[node_surface2[i].N[2]].pos.x[1]) * node_surface2[i].m_normal[1] + (node_surface2[i].pos.x[2] - node_surface2[node_surface2[i].N[2]].pos.x[2]) * node_surface2[i].m_normal[2])
				+ (node_surface2[i].cota[3] + node_surface2[i].cotb[3]) * ((node_surface2[i].pos.x[0] - node_surface2[node_surface2[i].N[3]].pos.x[0]) * node_surface2[i].m_normal[0]
				+ (node_surface2[i].pos.x[1] - node_surface2[node_surface2[i].N[3]].pos.x[1]) * node_surface2[i].m_normal[1] + (node_surface2[i].pos.x[2] - node_surface2[node_surface2[i].N[3]].pos.x[2]) * node_surface2[i].m_normal[2])))  / 3;
	}
		if (node_surface2[i].index_t == 2){
			if (i != ((num_count + Hexa_h) / 2 - 2 * Hexa_w + 1)){
				for (j = 1; j < 2; j++){
					node_surface2[i].K = (0.25 * ((node_surface2[i].cota[1] + node_surface2[i].cotb[1]) *((node_surface2[i].pos.x[0] - node_surface2[node_surface2[i].N[1]].pos.x[0]) * node_surface2[i].m_normal[0]
						+ (node_surface2[i].pos.x[1] - node_surface2[node_surface2[i].N[1]].pos.x[1]) * node_surface2[i].m_normal[1] + (node_surface2[i].pos.x[2] - node_surface2[node_surface2[i].N[1]].pos.x[2]) * node_surface2[i].m_normal[2])
						+ (node_surface2[i].cota[2] + node_surface2[i].cotb[2]) * ((node_surface2[i].pos.x[0] - node_surface2[node_surface2[i].N[2]].pos.x[0]) * node_surface2[i].m_normal[0]
						+ (node_surface2[i].pos.x[1] - node_surface2[node_surface2[i].N[2]].pos.x[1]) * node_surface2[i].m_normal[1] + (node_surface2[i].pos.x[2] - node_surface2[node_surface2[i].N[2]].pos.x[2]) * node_surface2[i].m_normal[2]))) / 2;
				}
			}
			if (i == ((num_count + Hexa_h) / 2 - 2 * Hexa_w + 1)){
				for (j = 2; j < 3; j++){
					node_surface2[i].K = (0.25 * ((node_surface2[i].cota[1] + node_surface2[i] .cotb[1]) *((node_surface2[i].pos.x[0] - node_surface2[node_surface2[i].N[1]].pos.x[0]) * node_surface2[i].m_normal[0]
						+ (node_surface2[i].pos.x[1] - node_surface2[node_surface2[i].N[1]].pos.x[1]) * node_surface2[i].m_normal[1] + (node_surface2[i].pos.x[2] - node_surface2[node_surface2[i].N[1]].pos.x[2]) * node_surface2[i].m_normal[2])
						+ (node_surface2[i].cota[2] + node_surface2[i].cotb[2]) * ((node_surface2[i].pos.x[0] - node_surface2[node_surface2[i].N[2]].pos.x[0]) * node_surface2[i].m_normal[0]
						+ (node_surface2[i].pos.x[1] - node_surface2[node_surface2[i].N[2]].pos.x[1]) * node_surface2[i].m_normal[1] + (node_surface2[i].pos.x[2] - node_surface2[node_surface2[i].N[2]].pos.x[2]) * node_surface2[i].m_normal[2]))) / 2;

				}
			}
		}
	}

	for (i = 0; i < num_count * 0.5; i++){
		//printf("i = %d, K = %lf\n", i, node_surface2[i].K);
	}
	
	///////////////////////////////////////////////////////////////////////////
	//TEXT out
#if 0
	FILE *fp;

	fopen_s(&fp, "curvature.txt", "w");
	for (i = 0; i < num_count * 0.5; i++){
		fprintf(fp, "%lf\n", node_surface2[i].K);
	}

	/*fopen_s(&fp, "coordinate.txt", "w");
	for (i = 0; i < num_points; i++){
	fprintf(fp, "%lf %lf %lf\n", point[i][0], point[i][1], point[i][2]);
	}*/
	fclose(fp);
#endif


#endif

	//for (i = 0; i < num_count; i++){
	//	glPushMatrix();
	//	glColor3d(0.0, 1.0, 1.0);
	//	/*printf("%lf, %lf, %lf\n", node_surface2[i].pos.x[0], node_surface2[i].pos.x[1], node_surface2[i].pos.x[2]);*/
	//	glTranslated((GLdouble)node_surface2[i].pos.x[0], (GLdouble)node_surface2[i].pos.x[2], (GLdouble)node_surface2[i].pos.x[1]);
	//	if (view_con == 1){
	//		glutSolidSphere(0.08, 10, 10); 
	//	}
	//	glPopMatrix();
	//}
	//glPopMatrix();
	GLfloat changing[] = { 0.5, 0.5, 1.0, 1.0 }; // Blue

	for (i = 0; i < num_count; i++){
		//node_surface2[i].color_grad = sqrt(pow(node_surface2[i].acc.x[0], 2.0) + pow(node_surface2[i].acc.x[1], 2.0) + pow(node_surface2[i].acc.x[2], 2.0));
		if (node_surface2[i].pos.x[2] > max){
			max = node_surface2[i].pos.x[2];
		}
		if (node_surface2[i].pos.x[2] < min){
			min = node_surface2[i].pos.x[2];
		}
		if (node_surface2[i].K > max_g){
			node_surface2[i].K = max_g;
		}
		if (node_surface2[i].K < min_g){
			node_surface2[i].K = min_g;
		}
	}
	printf("gap = %lf\n", max - min);
	//printf("node_surface2[0] = %f, %f, %f\n", node_surface2[0].pos.x[0], node_surface2[0].pos.x[1], node_surface2[0].pos.x[2]);
	/*printf("max = %lf min = %lf\n", max, min);*/
	for (i = 0; i < num_count; i++){
		glPushMatrix();
		glCullFace(GL_BACK);
		//changing[0] = (node_surface2[i].color_grad - min) / (max - min);
		//changing[1] = 0.0;
		//changing[2] = 1.0 - (node_surface2[i].color_grad) / (max - min);
		if (node_surface2[i].color_grad < (max - min) / 2.0){
			changing[0] = (node_surface2[i].color_grad - min) / ((max + min) / 2.0 - min);
			changing[1] = 0.1;
			changing[2] = ((max + min) / 2.0 - node_surface2[i].color_grad) / ((max + min) / 2.0 - min);
		}
		else{
			changing[0] = (max - node_surface2[i].color_grad) / ((max + min) / 2.0 - min);
			changing[1] = (node_surface2[i].color_grad - (max + min) / 2.0) / ((max + min) / 2.0 - min);
			changing[2] = 0.1;
			//changing[0] = 1.0;
			//changing[1] = 0.1;
			//changing[2] = 0.1;
			//printf("%f %f %f\n", changing[0], changing[1], changing[2]);
		}
		//glTranslated((GLdouble)node_surface2[i].pos.x[0] - 3, (GLdouble)node_surface2[i].pos.x[2], (GLdouble)node_surface2[i].pos.x[1] - 7);
		/*if (view_con == 1) glutSolidSphere(node_Radius, 10, 10);
		else if (view_con == 2){
			glMaterialfv(GL_FRONT, GL_DIFFUSE, changing);
			glutSolidCube(node_Radius * 4.0);*/
	//	}
		glPopMatrix();
	}
	/*glDisable(GL_LIGHTING);
	glBegin(GL_QUADS);
	glColor3d(1.0, 0.0, 0.0);
	glVertex3d(0, 30, 30);
	glVertex3d(0, 25, 30);
	glVertex3d(0, 30, 25);
	glVertex3d(0, 25, 25);
	glEnd();*/


	glPushMatrix();
	//glTranslated( - rec_x / 2.0,0.0,  - rec_y / 2.0);
	for (i = 0; i < tri_count; i++){
		if (i < tri_count * 0.5){
			if (triangle_data[i].color == 1){
				glCullFace(GL_FRONT);
			}
			else{
				glCullFace(GL_BACK);
			}
		}
		else if (i >= tri_count * 0.5){
			if (triangle_data[i].color == 1){
				glCullFace(GL_BACK);
			}
			else {
				glCullFace(GL_FRONT);
			}
		}
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		if (view_con == 2){
			glMaterialfv(GL_FRONT, GL_DIFFUSE, white);
		}
		else{
			glMaterialfv(GL_FRONT, GL_DIFFUSE, blue2);
		}
		glEnable(GL_LIGHTING);
		//glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
		//glEnable(GL_COLOR_MATERIAL);
		glBegin(GL_TRIANGLES);
		//printf("%f %f %f %f\n", blue[0], blue[1], blue[2], blue[3]);
		if (1){
			if (i < tri_count * 0.5){
				glNormal3d(triangle_data[i].normal[0], triangle_data[i].normal[2], triangle_data[i].normal[1]);
				for (j = 0; j < 3; j++){
					if (node_surface2[triangle_data[i].t[j]].K < (max_g - min_g) / 2.0){
						//glColor3d((node_surface2[triangle_data[i].t[j]].K - min_g) / ((max_g + min_g) / 2 - min_g), 0.1, ((max_g + min_g) / 2.0 - node_surface2[triangle_data[i].t[j]].K) / ((max_g + min_g) / 2.0 - min_g));
						//glColor3d(0.0, (node_surface2[triangle_data[i].t[j]].K - min_g) / ((max_g + min_g) / 2 - min_g), ((max_g + min_g) / 2.0 - node_surface2[triangle_data[i].t[j]].K) / ((max_g + min_g) / 2.0 - min_g));
						
						glColor3d(((max_g + min_g) / 2.0 - node_surface2[triangle_data[i].t[j]].K) / ((max_g + min_g) / 2.0 - min_g), 0.0, (node_surface2[triangle_data[i].t[j]].K - min_g) / ((max_g + min_g) / 2 - min_g));
						glVertex3d(node_surface2[triangle_data[i].t[j]].pos.x[0] - 3, node_surface2[triangle_data[i].t[j]].pos.x[2], node_surface2[triangle_data[i].t[j]].pos.x[1] - 7);
					}
					else{
						//glColor3d(((node_surface2[triangle_data[i].t[j]].K) - (max_g + min_g) / 2.0) / ((max_g + min_g) / 2 - min_g), (max_g - (node_surface2[triangle_data[i].t[j]].K)) / ((max_g + min_g) / 2.0 - min_g), 0.1);
						//glColor3d(((node_surface2[triangle_data[i].t[j]].K) - (max_g + min_g) / 2.0) / ((max_g + min_g) / 2 - min_g), (max_g - (node_surface2[triangle_data[i].t[j]].K)) / ((max_g + min_g) / 2.0 - min_g), 0.0);
						
						glColor3d(0.0, ((node_surface2[triangle_data[i].t[j]].K) - (max_g + min_g) / 2.0) / ((max_g + min_g) / 2 - min_g), (max_g - (node_surface2[triangle_data[i].t[j]].K)) / ((max_g + min_g) / 2.0 - min_g));
						glVertex3d(node_surface2[triangle_data[i].t[j]].pos.x[0] - 3, node_surface2[triangle_data[i].t[j]].pos.x[2], node_surface2[triangle_data[i].t[j]].pos.x[1] - 7);
					}
				}
			}
		}
		glEnd();
	}
	glPopMatrix();
}
void idle(void)
{
	glutPostRedisplay();
}
void display(void)
{

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();

	gluLookAt(View_from[0], View_from[1], View_from[2], View_to[0], View_to[1], View_to[2], 0.0, 1.0, 0.0);
	glViewport(0, 0, w_view * 2.0 / 3.0, h_view);
	glPushMatrix();
	node_simulation(1);
	glPopMatrix();

	if (MouseFlagLeft){
		if (View_point_flag) View_control(false);
		else View_control2(false);
	}
	else if (MouseFlagRight){
		if (View_point_flag) View_control(true);
		else View_control2(true);
	}
	if (up_flag) View_control_up_down(true);
	if (down_flag) View_control_up_down(false);
	//for (i = 0; i < num_count; i++){
	//		glPushMatrix();
	//		glColor3d(0.0, 1.0, 1.0);
	//		/*printf("%lf, %lf, %lf\n", node_surface2[i].pos.x[0], node_surface2[i].pos.x[1], node_surface2[i].pos.x[2]);*/
	//		glTranslated((GLdouble)node_surface2[i].pos.x[0], (GLdouble)node_surface2[i].pos.x[2], (GLdouble)node_surface2[i].pos.x[1]);
	//		glutSolidSphere(0.08, 10, 10);
	//		glPopMatrix();
	//	}
	//glPopMatrix();
	//glPushMatrix();

	//	for (i = 0; i < tri_count; i++){
	//		if (triangle_data[i].color == 1){
	//				glCullFace(GL_FRONT);
	//			}
	//			else{
	//				glCullFace(GL_BACK);
	//			}
	//			glEnable(GL_BLEND);
	//			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	//			glBegin(GL_TRIANGLES);
	//			//printf("%f %f %f %f\n", blue[0], blue[1], blue[2], blue[3]);
	//			//glNormal3d(-triangle_data[i].normal[0], -triangle_data[i].normal[2], -triangle_data[i].normal[1]);
	//			glNormal3d(triangle_data[i].normal[0], triangle_data[i].normal[2], triangle_data[i].normal[1]);
	//			//glNormal3d(-triangle_data[i].normal[0],- triangle_data[i].normal[1], -triangle_data[i].normal[2]);
	//			glVertex3d(node_surface2[triangle_data[i].t[0]].pos.x[0], node_surface2[triangle_data[i].t[0]].pos.x[2], node_surface2[triangle_data[i].t[0]].pos.x[1]);
	//			glVertex3d(node_surface2[triangle_data[i].t[1]].pos.x[0], node_surface2[triangle_data[i].t[1]].pos.x[2], node_surface2[triangle_data[i].t[1]].pos.x[1]);
	//			glVertex3d(node_surface2[triangle_data[i].t[2]].pos.x[0], node_surface2[triangle_data[i].t[2]].pos.x[2], node_surface2[triangle_data[i].t[2]].pos.x[1]);
	//			/*glVertex3d(node_surface2[triangle_data[i].t[0]].pos.x[0] - rec_x / 2.0, node_surface2[triangle_data[i].t[0]].pos.x[2], node_surface2[triangle_data[i].t[0]].pos.x[1] - rec_y / 2.0);
	//			glVertex3d(node_surface2[triangle_data[i].t[1]].pos.x[0] - rec_x / 2.0, node_surface2[triangle_data[i].t[1]].pos.x[2], node_surface2[triangle_data[i].t[1]].pos.x[1] - rec_y / 2.0);
	//			glVertex3d(node_surface2[triangle_data[i].t[2]].pos.x[0] - rec_x / 2.0, node_surface2[triangle_data[i].t[2]].pos.x[2], node_surface2[triangle_data[i].t[2]].pos.x[1] - rec_y / 2.0);*/
	//		}
	//	glEnd();
	//glPopMatrix();

	glutSwapBuffers();

}
void mouse(int button, int state, int x, int y)
{
	switch (button) {
	case GLUT_LEFT_BUTTON:
		switch (state) {
		case GLUT_UP:
			if (MouseFlagLeft){
				MouseFlagLeft = false;
			}
			break;
		case GLUT_DOWN:
			MouseFlagLeft = true;
			if (x < window_size_x * 2 / 3) View_point_flag = true;
			else View_point_flag = false;
			break;
		default:
			break;
		}
		break;
	case GLUT_MIDDLE_BUTTON:
		switch (state) {
		case GLUT_UP:
			if (MouseFlagRight) MouseFlagMiddle = false;
			break;
		case GLUT_DOWN:
			MouseFlagMiddle = true;
			break;
		default:
			break;
		}
		break;
	case GLUT_RIGHT_BUTTON:
		switch (state) {
		case GLUT_UP:
			if (MouseFlagRight) MouseFlagRight = false;
			break;
		case GLUT_DOWN:
			MouseFlagRight = true;
			if (x < window_size_x * 2 / 3) View_point_flag = true;
			else View_point_flag = false;
			break;
		default:
			break;
		}
		break;
	default:
		break;
	}
}
void resize(int w, int h)
{
	w_view = w;
	h_view = h;

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	gluPerspective(50.0, (double)w / (double)h * 2.0 / 3.0, 1.0, 1000.0);

	glMatrixMode(GL_MODELVIEW);

}
void keyboard(unsigned char key, int x, int y){
	switch (key){
	case 'r':
		close_flag = true;
		break;
	case 'e':
		close_flag = false;
		break;
	case 'w':
		open_flag = true;
		break;
	case 'q':
		open_flag = false;
		break;
	case 'y':
		close_flag_n = true;
		break;
	case 'u':
		close_flag_n = false;
		break;
	case 'i':
		open_flag_n = true;
		break;
	case 'o':
		open_flag_n = false;
		break;
	}
}
void init(){

	glClearColor(1.0, 1.0, 1.0, 1.0);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	GLfloat direction[] = { 0.0, 1.0, 0.0 };
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);//アルファの設定
	glEnable(GL_BLEND);//アルファのブレンド有効
	glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, direction);
	glLightfv(GL_LIGHTING, GL_SPOT_DIRECTION, direction);
	glDisable(GL_LIGHT2);
	glDisable(GL_LIGHT3);
	glDisable(GL_LIGHT4);

}
int main(int argc, char *argv[])
{
	//get_info();

	glutInit(&argc, argv);
	glutInitWindowPosition(0, 0);
	glutInitWindowSize(window_size_x, window_size_y);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
	glutCreateWindow(argv[0]);
	glutInitWindowPosition(0, 0);
	glutDisplayFunc(display);
	glutReshapeFunc(resize);
	glutIdleFunc(idle);
	glutMouseFunc(mouse);
	glutKeyboardFunc(keyboard);
	/*node_simulation();*/
	init();
	glutMainLoop();

	return 0;

}
