// air_bag_2014_10_21.cpp : コンソール アプリケーションのエントリ ポイントを定義します。
//

// finger_contact_model.cpp : コンソール アプリケーションのエントリ ポイントを定義します。
//

#include <stdio.h>
#include <stdlib.h>
#include <GL/glut.h>
#include <math.h>
#include <string.h>
#include <vector>

using namespace std;
#define node_Num_m 11
//#define node_Num_n 21
//#define node_Num 10002 //must be >node_Num_m*node_Num_n*2
#define rec_x 3.5
#define rec_y 7
#define judge 0.0001
#define MPI 3.14159265358979323846
double wall_z = 0.5;
double wall_n = 5.0;
int num_count = 0;
int tri_count = 0;
bool close_flag = false;
bool open_flag = false;
bool close_flag_n = false;
bool open_flag_n = false;

int w_view;
int h_view;
int windows_size_x = 1000;
int windows_size_y = windows_size_x;
int first_count = 1;
//GLfloat light0pos[] = { -300.0, 300.0, 500.0, 1.0 }; // 光源1
//GLfloat light0_ambient[] = { 0.1, 0.1, 0.1, 1.0 }; // 光源1
//GLfloat light0_diffuse[] = { 0.99, 0.99, 0.99, 1.0 }; // 光源1
//GLfloat light0_specular[] = { 0.1, 0.1, 0.1, 0.2 }; // 光源1
//GLfloat light0_shiniess[] = { 30 }; // 光源1
//GLfloat light1pos[] = { 560.0, 30.0, -500.0, 1.0 }; // 光源2
//GLfloat light1_ambient[] = { 0.01, 0.01, 0.01, 1.0 }; // 光源2
//GLfloat light1_diffuse[] = { 0.0, 0.0, 0.0, 1.0 }; // 光源2
//GLfloat light1_specular[] = { 0.01, 0.01, 0.01, 1.0 }; // 光源2
//GLfloat light1_shiniess[] = { 20 }; // 光源2
//GLfloat light2pos[] = { 0.0, 0.0, 100.0, 0.0 };//光源3
//GLfloat light3pos[] = { 0.0, -100.0, 0.0, 0.0 }; // 光源4
//GLfloat light4pos[] = { 0.0, -100.0, 100.0, 0.0 }; // 光源5
//GLfloat light5pos[] = { 0.0, -100.0, 0.0, 0.0 };//光源6
GLfloat orange[] = { 255.0 / 256.0, 153.0 / 256.0, 0.0 / 256.0, 0.9 };
GLfloat blue[] = { 0.0 / 256.0, 65.0 / 256.0, 255.0 / 256.0, 0.4 };
GLfloat white[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat blue2[] = { 102.0 / 256.0, 204.0 / 256.0, 255.0 / 256.0, 0.9 };
GLfloat blue_node[] = { 0.5, 0.5, 1.0, 1.0 }; // Blue
bool up_flag = false;
bool down_flag = false;
double damp_k = 1000.0;	//各粒子同士のばね定数
double damp_k_normal = 10;
double dv = 1.0;
double Finger_glut_Radius = 7.0;
double node_Radius = 0.04;
double View_from[3] = { 0.0, 8.0, 13.0 };
double View_to[3] = { 0.0, 0.0, 0.0 };
double View_from2[3] = { 0.0, 13.0, 0.01 };
double View_to2[3] = { 0.0, -10.0, 0.0 };
double View_from3[3] = { 0.0, 13.0, 0.01 };
double View_to3[3] = { 0.0, -10.0, 0.0 };
bool MouseFlagRight = false;
bool MouseFlagLeft = false;
bool MouseFlagMiddle = false;
bool View_point_flag = false;
//double t_sum = 0.0;
#define DEG_UPDATE 5
GLUnurbsObj *theNurb;
typedef struct{
	double x[3];
}position;
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
	position pos;		//位置
	position del_pos;	//速度
	position acc;		//加速度
	double color_grad;
}node2;
typedef struct{
	bool flag;
	position ini_position;
}ini_flag;
static node2 node_surface2[100];
//node *node_surface2 = NULL;
static edge_d edge[11][11];
//edge_d *edge = NULL;
static triangle_d triangle_data[node_Num_m * 11 * 4];
//triangle_d *triangle_data = NULL;
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
	// 球の描画
	//半透明表示
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, sph_col);
	GLUquadricObj *sphere; //オブジェクトポインタを準備
	sphere = gluNewQuadric();	//オブジェクトを生成 
	//オブジェクトの描画タイプを設定（省略可）
	gluQuadricDrawStyle(sphere, GLU_FILL);
	//球を描画 半径1.0，緯経それぞれ10.0分割
	gluSphere(sphere, R, precise, precise);
}
void View_control(bool vector_flag){
	double View_distance;
	double temp[5];
	temp[2] = View_from[2] - View_to[2];
	temp[1] = View_from[0] - View_to[0];
	temp[0] = pow(temp[1], 2.0) + pow(temp[2], 2.0);
	View_distance = pow(temp[0], 0.5);
	//	printf("%f\n", View_distance);
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


	static node node_surface[100][100][2];
	
	for (s = 0; s < 2; s++){
		for (i = 5; i < node_Num_m; i++){
			for (j = 0; j <= i; j++){
				if (i % 2 == 0){
					node_surface[i][j][s].pos.x[0] = (double)((0.0) + (j - (i / 2))*natural_length);
				}
				else{
					node_surface[i][j][s].pos.x[0] = (double)(((-i * 0.5) + j)*natural_length);

				}
				node_surface[i][j][s].pos.x[1] = (double)(i * sqrt(3) * 0.5);
				node_surface[i][j][s].pos.x[2] = (double)(0.0);
				//printf("node = %f\n", node_surface[0][0][0].pos.x[1]);
			}
		}
	}
	
	for (s = 0; s < 2; s++){
		for (i = 5; i < node_Num_m; i++){
			for (j = 0; j <= i; j++){
				node_surface[i][j][s].number = num_count;
				num_count++;
			}
		}
	}
	printf("num_count = %d\n", num_count);

	printf("%d\n", tri_count);

	for (s = 0; s < 2; s++){
		for (i = 5; i < node_Num_m; i++){
			for (j = 0; j <= i; j++){
				for (h = 0; h < 3; h++){
					node_surface2[node_surface[i][j][s].number].pos.x[h] = node_surface[i][j][s].pos.x[h];
				}
			}
		}
	}

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
	double min = 100000.0;
	double max = -10.0;
	double mass = 1000;
	double kizami = 0.01;
	double total_force[3] = { 0.0, 0.0, 0.0 };
	double temp_len;
	double normal_force[3];
	double normal_temp[9];
	double normal_temp3[3];
	double temp_wall_n;
	double temp_wall_n_2;
	for (i = 0; i <= num_count - 1; i++){
		for (j = 0; j < 3; j++){
			node_surface2[i].acc.x[j] = 0.0;
		}
	}
	for (i = 0; i <= num_count - 1; i++){
		node_surface2[i].color_grad = 0.0;
		for (j = 0; j <= num_count - 1; j++){
			if (edge[i][j].torf == 1){
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
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	if (open_flag == true){
		wall_z += 0.01;
	}
	if (close_flag == true){
		wall_z -= 0.01;
	}
	/*if (open_flag_n == true){
	wall_n += 0.01;
	}
	if (close_flag_n == true){
	wall_n -= 0.01;
	}
	for (i = 0; i <= num_count - 1; i++){
	if (node_surface2[i].pos.x[2] > wall_z){
	if (node_surface2[i].pos.x[1] - rec_y / 2.0 > -rec_x / 2.0 / 2.0 && node_surface2[i].pos.x[1] - rec_y / 2.0 < rec_x / 2.0 / 2.0){
	if (node_surface2[i].pos.x[0] - rec_x / 2.0 > -rec_x / 2.0 / 2.0 && node_surface2[i].pos.x[0] - rec_x / 2.0 < rec_x / 2.0 / 2.0){
	node_surface2[i].acc.x[2] -= 1000.0 * damp_k * (node_surface2[i].pos.x[2] - wall_z) / mass;
	node_surface2[i].color_grad += fabs(1000.0 * damp_k * (node_surface2[i].pos.x[2] - wall_z) / mass);
	}
	}
	}
	else if (node_surface2[i].pos.x[2] < (-1.0 * wall_z)){
	node_surface2[i].acc.x[2] -= 1000.0 * damp_k * (node_surface2[i].pos.x[2] + wall_z) / mass;
	node_surface2[i].color_grad += fabs(1000.0 * damp_k * (node_surface2[i].pos.x[2] + wall_z) / mass);
	}
	}
	printf("node_surface2 = %f, wall_z = %f\n", node_surface2[8].acc.x[2], wall_z);*/
	//for (i = 0; i <= num_count - 1; i++){
	//	temp_wall_n = fabs((node_surface2[i].pos.x[0] - rec_x / 2.0) + node_surface2[i].pos.x[2] - sqrt(2.0) * wall_n) / sqrt(2.0);
	//	temp_wall_n_2 = fabs((node_surface2[i].pos.x[0] - rec_x / 2.0) + node_surface2[i].pos.x[2] + sqrt(2.0) * wall_n) / sqrt(2.0);

	//	//		if(node_surface2[i].pos.x[2] > (wall_n * sqrt(2.0) - node_surface2[i].pos.x[0]) && node_surface2[i].pos.x[2] > -1.0 * node_surface2[i].pos.x[0]){
	//	if (node_surface2[i].pos.x[2] > (wall_n * sqrt(2.0) - (node_surface2[i].pos.x[0] - rec_x / 2.0))){
	//		//printf("a\n");
	//		//printf("%lf %lf\n",temp_wall_n, wall_n);
	//		node_surface2[i].acc.x[0] -= 1000.0 * damp_k * (temp_wall_n) / mass;
	//		node_surface2[i].acc.x[2] -= 1000.0 * damp_k * (temp_wall_n) / mass;
	//		node_surface2[i].color_grad += fabs(1000.0 * damp_k * (temp_wall_n) / mass * sqrt(2.0));
	//	}
	//	else if (node_surface2[i].pos.x[2] < -wall_n * sqrt(2.0) - (node_surface2[i].pos.x[0] - rec_x / 2.0)){
	//		node_surface2[i].acc.x[0] += 1000.0 * damp_k * (temp_wall_n_2) / mass / sqrt(2.0);
	//		node_surface2[i].acc.x[2] += 1000.0 * damp_k * (temp_wall_n_2) / mass / sqrt(2.0);
	//		node_surface2[i].color_grad += fabs(1000.0 * damp_k * (temp_wall_n_2) / mass * sqrt(2.0));
	//	}
	//}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////dc

	for (i = 0; i <= num_count - 1; i++){
		for (k = 0; k < 3; k++){
			node_surface2[i].acc.x[k] += -dv * node_surface2[i].del_pos.x[k];
		}
	}
	//pressure
	for (i = 0; i < tri_count; i++){
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
		for (j = 0; j < 3; j++){
			triangle_data[i].normal[j] = normal_force[j] / sqrt(pow(normal_force[0], 2.0) + pow(normal_force[1], 2.0) + pow(normal_force[2], 2.0));// / normal_force[j];
			for (k = 0; k < 3; k++){
				node_surface2[triangle_data[i].t[j]].acc.x[k] += 0.1 * damp_k_normal * normal_force[k];
			}
		}
	}
	//if (flag_loop == 1){
	//	for (i = 0; i <= num_count - 1; i++){
	//		for (k = 0; k < 3; k++){
	//			printf("node_surface2[%d].acc.x[%d] = %lf\n", i, k, node_surface2[i].acc.x[k]);
	//		}
	//	}
	//	flag_loop--;
	//}
	for (i = 0; i < num_count; i++){
		for (k = 0; k < 3; k++){
			node_surface2[i].del_pos.x[k] = node_surface2[i].del_pos.x[k] + node_surface2[i].acc.x[k] * kizami;
			node_surface2[i].pos.x[k] = node_surface2[i].pos.x[k] + node_surface2[i].del_pos.x[k] * kizami + 1.0 / 2.0 * node_surface2[i].acc.x[k] * kizami * kizami;
		}
	}
	//	t_sum += kizami;
	GLfloat changing[] = { 0.5, 0.5, 1.0, 1.0 }; // Blue

	for (i = 0; i < num_count; i++){
		//node_surface2[i].color_grad = sqrt(pow(node_surface2[i].acc.x[0], 2.0) + pow(node_surface2[i].acc.x[1], 2.0) + pow(node_surface2[i].acc.x[2], 2.0));
		if (node_surface2[i].color_grad > max){
			max = node_surface2[i].color_grad;
		}
		if (node_surface2[i].color_grad < min){
			min = node_surface2[i].color_grad;
		}
	}
	//printf("max = %lf min = %lf\n", max, min);
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
		glTranslated((GLdouble)node_surface2[i].pos.x[0] - rec_x / 2.0, (GLdouble)node_surface2[i].pos.x[2], (GLdouble)node_surface2[i].pos.x[1] - rec_y / 2.0);
		if (view_con == 1) sphere(node_Radius, 10.0, changing);
		else if (view_con == 2){
			glMaterialfv(GL_FRONT, GL_DIFFUSE, changing);
			glutSolidCube(node_Radius * 4.0);
		}
		glPopMatrix();
	}
	glPushMatrix();
	//glTranslated( - rec_x / 2.0,0.0,  - rec_y / 2.0);
	for (i = 0; i < tri_count; i++){
		if (triangle_data[i].color == 1){
			glCullFace(GL_FRONT);
		}
		else{
			glCullFace(GL_BACK);
		}
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		if (view_con == 2){
			glMaterialfv(GL_FRONT, GL_DIFFUSE, white);
		}
		else{
			glMaterialfv(GL_FRONT, GL_DIFFUSE, blue2);
		}
		glBegin(GL_TRIANGLES);
		//printf("%f %f %f %f\n", blue[0], blue[1], blue[2], blue[3]);
		if (1){
			//glNormal3d(-triangle_data[i].normal[0], -triangle_data[i].normal[2], -triangle_data[i].normal[1]);
			glNormal3d(triangle_data[i].normal[0], triangle_data[i].normal[2], triangle_data[i].normal[1]);
			//glNormal3d(-triangle_data[i].normal[0],- triangle_data[i].normal[1], -triangle_data[i].normal[2]);
			glVertex3d(node_surface2[triangle_data[i].t[0]].pos.x[0] - rec_x / 2.0, node_surface2[triangle_data[i].t[0]].pos.x[2], node_surface2[triangle_data[i].t[0]].pos.x[1] - rec_y / 2.0);
			glVertex3d(node_surface2[triangle_data[i].t[1]].pos.x[0] - rec_x / 2.0, node_surface2[triangle_data[i].t[1]].pos.x[2], node_surface2[triangle_data[i].t[1]].pos.x[1] - rec_y / 2.0);
			glVertex3d(node_surface2[triangle_data[i].t[2]].pos.x[0] - rec_x / 2.0, node_surface2[triangle_data[i].t[2]].pos.x[2], node_surface2[triangle_data[i].t[2]].pos.x[1] - rec_y / 2.0);
		}
		glEnd();
	}
	glPopMatrix();

	//glPushMatrix();
	//glEnable(GL_BLEND);
	//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	//glMaterialfv(GL_FRONT, GL_DIFFUSE, blue);
	//glBegin(GL_TRIANGLES);
	//glNormal3d(0.0, 1.0, 0.0);
	//glVertex3d(-rec_x / 2.0/ 2.0, wall_z, -rec_x / 2.0 / 2.0);
	//glVertex3d(rec_x / 2.0 / 2.0, wall_z, -rec_x / 2.0 / 2.0);
	//glVertex3d(rec_x / 2.0 / 2.0, wall_z, rec_x / 2.0 / 2.0);
	//glEnd();
	//glPopMatrix();
	//glPushMatrix();
	//glBegin(GL_TRIANGLES);
	//glNormal3d(0.0, 1.0, 0.0);
	//glVertex3d(-rec_x / 2.0 / 2.0, wall_z, rec_x / 2.0 / 2.0);
	//glVertex3d(-rec_x / 2.0 / 2.0, wall_z, -rec_x / 2.0 / 2.0);
	//glVertex3d(rec_x / 2.0 / 2.0, wall_z, rec_x / 2.0 / 2.0);
	//glEnd();
	//glPopMatrix();
	//glPushMatrix();
	//glBegin(GL_TRIANGLES);
	//glNormal3d(0.0, 1.0, 0.0);
	//glVertex3d(-rec_x / 2.0, -1.0 * wall_z, -rec_y / 2.0);
	//glVertex3d(rec_x / 2.0, -1.0 * wall_z, -rec_y / 2.0);
	//glVertex3d(rec_x / 2.0, -1.0 * wall_z, rec_y / 2.0);
	//glEnd();
	//glPopMatrix();
	//glPushMatrix();
	//glBegin(GL_TRIANGLES);
	//glNormal3d(0.0, 1.0, 0.0);
	//glVertex3d(-rec_x / 2.0, -1.0 * wall_z, rec_y / 2.0);
	//glVertex3d(-rec_x / 2.0, -1.0 * wall_z, -rec_y / 2.0);
	//glVertex3d(rec_x / 2.0, -1.0 * wall_z, rec_y / 2.0);
	//glEnd();

	//glEnable(GL_BLEND);
	//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	//glMaterialfv(GL_FRONT, GL_DIFFUSE, blue);
	//glBegin(GL_TRIANGLES);
	//glNormal3d(0.0,1.0,0.0);
	//glVertex3d(- rec_x / 2.0, wall_z, -rec_y/2.0);
	//glVertex3d(rec_x / 2.0,wall_z, -rec_y/2.0);
	//glVertex3d(rec_x / 2.0,wall_z, rec_y/2.0);
	//glEnd();
	//glPopMatrix();
	//glPushMatrix();
	//glBegin(GL_TRIANGLES);
	//glNormal3d(0.0, 1.0, 0.0);
	//glVertex3d(- rec_x / 2.0,wall_z, rec_y/2.0);
	//glVertex3d(-rec_x / 2.0,wall_z, -rec_y/2.0);
	//glVertex3d(rec_x / 2.0,wall_z, rec_y/2.0);
	//glEnd();
	//glPopMatrix();
	//glPushMatrix();
	//glBegin(GL_TRIANGLES);
	//glNormal3d(0.0, 1.0, 0.0);
	//glVertex3d(-rec_x / 2.0, -1.0 * wall_z, -rec_y / 2.0);
	//glVertex3d(rec_x / 2.0,-1.0 * wall_z, -rec_y/2.0);
	//glVertex3d(rec_x / 2.0,-1.0 * wall_z, rec_y/2.0);
	//glEnd();
	//glPopMatrix();
	//glPushMatrix();
	//glBegin(GL_TRIANGLES);
	//glNormal3d(0.0, 1.0, 0.0);
	//glVertex3d(-rec_x / 2.0, -1.0 * wall_z, rec_y / 2.0);
	//glVertex3d(-rec_x / 2.0,-1.0 * wall_z, -rec_y/2.0);
	//glVertex3d(rec_x / 2.0,-1.0 * wall_z, rec_y/2.0);
	//glEnd();
	glPopMatrix();


	glPushMatrix();
	//glEnable(GL_BLEND);
	//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	//glMaterialfv(GL_FRONT, GL_DIFFUSE, blue);
	//glBegin(GL_TRIANGLES);
	//glVertex3d(wall_n * sqrt(2.0) + 1.0 * (rec_x / 2.0 - wall_n)/ 2.0 , -1.0 * (rec_x / 2.0 - wall_n)/ 2.0, -rec_y/2.0);
	//glVertex3d(wall_n * sqrt(2.0) + 1.0 * (rec_x / 2.0 - wall_n)/ 2.0 , -1.0 * (rec_x / 2.0 - wall_n)/ 2.0, rec_y/2.0);
	//glVertex3d(- 1.0 * (rec_x / 2.0 - wall_n)/ 2.0 , wall_n * sqrt(2.0) + 1.0 * (rec_x / 2.0 - wall_n)/ 2.0, rec_y/2.0);
	//glEnd();
	//glPopMatrix();
	//glPushMatrix();
	//glBegin(GL_TRIANGLES);
	//glVertex3d(- 1.0 * (rec_x / 2.0 - wall_n)/ 2.0 , wall_n * sqrt(2.0) + 1.0 * (rec_x / 2.0 - wall_n)/ 2.0, rec_y/2.0);
	//glVertex3d(- 1.0 * (rec_x / 2.0 - wall_n)/ 2.0 , wall_n * sqrt(2.0) + 1.0 * (rec_x / 2.0 - wall_n)/ 2.0, - rec_y/2.0);
	//glVertex3d(wall_n * sqrt(2.0) + 1.0 * (rec_x / 2.0 - wall_n)/ 2.0 , -1.0 * (rec_x / 2.0 - wall_n)/ 2.0, -rec_y/2.0);
	//glEnd();
	//glPopMatrix();
	//glPushMatrix();
	//glBegin(GL_TRIANGLES);
	//glVertex3d(-1.0 * (wall_n * sqrt(2.0) + 1.0 * (rec_x / 2.0 - wall_n)/ 2.0) , 1.0 * (rec_x / 2.0 - wall_n)/ 2.0, rec_y/2.0);
	//glVertex3d(-1.0 * (wall_n * sqrt(2.0) + 1.0 * (rec_x / 2.0 - wall_n)/ 2.0) , 1.0 * (rec_x / 2.0 - wall_n)/ 2.0, - rec_y/2.0);
	//glVertex3d( 1.0 * (rec_x / 2.0 - wall_n)/ 2.0 , -1.0 * (wall_n * sqrt(2.0) + 1.0 * (rec_x / 2.0 - wall_n)/ 2.0), - rec_y/2.0);
	//glEnd();
	//glPushMatrix();
	//glPopMatrix();
	//glBegin(GL_TRIANGLES);
	//glVertex3d( 1.0 * (rec_x / 2.0 - wall_n)/ 2.0 , -1.0 * (wall_n * sqrt(2.0) + 1.0 * (rec_x / 2.0 - wall_n)/ 2.0), - rec_y/2.0);
	//glVertex3d( 1.0 * (rec_x / 2.0 - wall_n)/ 2.0 , -1.0 * (wall_n * sqrt(2.0) + 1.0 * (rec_x / 2.0 - wall_n)/ 2.0), rec_y/2.0);
	//glVertex3d(-1.0 * (wall_n * sqrt(2.0) + 1.0 * (rec_x / 2.0 - wall_n)/ 2.0) , 1.0 * (rec_x / 2.0 - wall_n)/ 2.0,  rec_y/2.0);
	//glEnd();
	glPopMatrix();
}
void idle(void)
{
	glutPostRedisplay(); // 1コマ進める関数
}
void display(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();
	/* 光源の位置設定 */

	//glLightfv(GL_LIGHT1, GL_POSITION, light2pos); // 光源２
	//glLightfv(GL_LIGHT2, GL_POSITION, light3pos); // 光源２
	//glLightfv(GL_LIGHT3, GL_POSITION, light4pos); // 光源２
	//glLightfv(GL_LIGHT4, GL_POSITION, light5pos); // 光源２
	//	glLightfv(GL_LIGHT0, GL_POSITION, light0pos); // 光源２
	//	glLightfv(GL_LIGHT1, GL_POSITION, light1pos); // 光源２
	/* 視点位置と視線方向 */
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
	//斜視
	gluLookAt(View_from[0], View_from[1], View_from[2], View_to[0], View_to[1], View_to[2], 0.0, 1.0, 0.0);//正面
	glViewport(0, 0, w_view * 2.0 / 3.0, h_view);
	glScissor(0, 0, w_view * 2.0 / 3.0, h_view);
	glPushMatrix();
	node_simulation(1);
	glPopMatrix();
	glLoadIdentity();
	//上方向
	gluLookAt(View_from2[0], View_from2[1], View_from2[2], View_to2[0], View_to2[1], View_to2[2], 0.0, 1.0, 0.0);//正面
	glViewport(w_view * 2.0 / 3.0, 0, w_view / 3.0, h_view / 2.0);
	glScissor(w_view * 2.0 / 3.0, 0, w_view / 3.0, h_view / 2.0);
	glPushMatrix();
	node_simulation(2);
	glPopMatrix();
	glLoadIdentity();
	//寄り
	gluLookAt(View_from3[0], View_from3[1], View_from3[2], View_to3[0], View_to3[1], View_to3[2], 0.0, 1.0, 0.0);//正面
	glViewport(w_view * 2.0 / 3.0, h_view / 2.0, w_view / 3.0, h_view / 2.0);
	glScissor(w_view * 2.0 / 3.0, h_view / 2.0, w_view / 3.0, h_view / 2.0);
	glPushMatrix();
	node_simulation(3);
	glPopMatrix();
	glLoadIdentity();
	glutSwapBuffers(); // ちらつき防止（1秒間に60フレーム 60Hz）
	//glutPostRedisplay();
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
			if (x < windows_size_x * 2 / 3) View_point_flag = true;
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
			if (x < windows_size_x * 2 / 3) View_point_flag = true;
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
	/* 透視変換行列の設定 */
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	/* スクリーン上の座標系をマウスの座標系に一致させる */
	// glOrtho(-0.5, (GLdouble)w - 0.5, (GLdouble)h - 0.5, -0.5, -1.0, 1.0);

	gluPerspective(50.0, (double)w / (double)h * 2.0 / 3.0, 1.0, 1000.0);
	//gluPerspective(50.0, (double)w / (double)h, 1.0, 100.0);
	//glViewport(0,0,w/2.0,h);
	//gluLookAt(View_from[0], View_from[1], View_from[2], View_to[0], View_to[1], View_to[2], 0.0, 1.0, 0.0);//正面


	//gluPerspective(50.0, (double)w / (double)h / 2.0, 1.0, 100.0);
	//gluLookAt(View_from2[0], View_from2[1], View_from2[2], View_to2[0], View_to2[1], View_to2[2], 0.0, 1.0, 0.0);//正面

	//gluPerspective(50.0, (double)w / (double)h, 1.0, 100.0);

	/* モデルビュー変換行列の設定 */
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
void init(void)
{
	//	theSpere = glGenLits(1);
	glClearColor(1.0, 1.0, 1.0, 0.0); // 背景色
	glEnable(GL_DEPTH_TEST);
	// 陰面消去用
	glEnable(GL_CULL_FACE); // 隠れている面を計算しない
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
	//glLightfv(GL_LIGHT1, GL_AMBIENT, light1_ambient);
	//glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_diffuse);
	//glLightfv(GL_LIGHT1, GL_SPECULAR, light1_specular);
	//glLightfv(GL_LIGHT1, GL_SHININESS, light1_shiniess);
}
int main(int argc, char *argv[])
{
	//initiation();
	glutInitWindowPosition(0, 0);
	glutInitWindowSize(windows_size_x, windows_size_y);
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
	glutCreateWindow(argv[0]);
	glutInitWindowPosition(0, 0);
	glutDisplayFunc(display);
	glutReshapeFunc(resize);
	glutIdleFunc(idle);
	glutMouseFunc(mouse);
	glutKeyboardFunc(keyboard);
	init();
	glutMainLoop();
	return 0;
}