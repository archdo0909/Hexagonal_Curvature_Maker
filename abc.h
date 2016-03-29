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
	int edge_flag;			//�m�[�h�ԍ�
	int none_flag;
	position pos;		//�ʒu
}node;
typedef struct{
	int number;
	int edge_flag;			//�m�[�h�ԍ�
	int none_flag;
	position pos;		//�ʒu
	position del_pos;	//���x
	position acc;		//�����x
	double color_grad;
}node2;
typedef struct{
	bool flag;
	position ini_position;
}ini_flag;

static node2 node_surface2[50000];
static node node_surface[500][500][3];
static edge_d edge[700][700];
static triangle_d triangle_data[50000];