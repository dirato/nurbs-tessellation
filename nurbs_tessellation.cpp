#include <SFML/Graphics.hpp>
#include<iostream>
#include<cmath>

using std::cout, std::cin, std::endl;

float TA[30][30], T[30][30], P[30][30];

bool ponto_fora_tela = false;

struct Point2D {
	int x;
	int y;
};
typedef struct Point2D Point2D;

Point2D* newPoint2D(int x, int y){
	Point2D* p2D = new Point2D;

	p2D->x = x;
	p2D->y = y;

	return p2D;
}

struct Point {
	float x;
	float y;
	float z;
};
typedef struct Point Point;

Point* newPoint(float x, float y, float z){
	Point* p = new Point;

	p->x = x;
	p->y = y;
	p->z = z;

	return p;
}

struct Vector {
	float x;
	float y;
	float z;
};
typedef struct Vector Vector;

Vector* newVector(float x, float y, float z){
	Vector* v = new Vector;

	v->x = x;
	v->y = y;
	v->z = z;

	return v;
}

struct Camera {
	Point* C;
	Vector* U;
	Vector* V;
	Vector* N;
};
typedef struct Camera Camera;

Camera* newCamera(Point* C, Vector* U, Vector* V, Vector* N){
	Camera* cam = new Camera;

	cam->C = C;
	cam->U = U;
	cam->V = V;
	cam->N = N;

	return cam;
}


float norma(Vector* v){

	return sqrt( pow((v->x), 2.0) + pow((v->y), 2.0) + pow((v->z), 2.0) );
}

Vector* sub_point_point(Point* p1, Point* p2){
// Retorna p2 - p1
	Vector* v = new Vector;

	v->x = p2->x - p1->x;
	v->y = p2->y - p1->y;
	v->z = p2->z - p1->z;

	return v;
}

Vector* sub_vector_vector(Vector* v1, Vector* v2){
// Retorna v2 - v1
	Vector* v = new Vector;

	v->x = v2->x - v1->x;
	v->y = v2->y - v1->y;
	v->z = v2->z - v1->z;

	return v;
}

Point* add_point_vector(Point* p, Vector* v){

	Point* _p = new Point;

	_p->x = p->x + v->x;
	_p->y = p->y + v->y;
	_p->z = p->z + v->z;

	return _p;
}

Point* add_point_point(Point* p1, Point* p2){

	Point* _p = new Point;

	_p->x = p1->x + p2->x;
	_p->y = p1->y + p2->y;
	_p->z = p1->z + p2->z;

	return _p;
}

Vector* prod_vector_escal(Vector* v, float t){
	Vector* vec_temp = new Vector;

	vec_temp->x = (v->x)*t;
	vec_temp->y = (v->y)*t;
	vec_temp->z = (v->z)*t;

	return vec_temp;
}

Vector* div_vector_escal(Vector* v, float t){
	Vector* vec_temp = new Vector;

	vec_temp->x = (float)((v->x)/t);
	vec_temp->y = (float)((v->y)/t);
	vec_temp->z = (float)((v->z)/t);

	return vec_temp;
}

Point* div_point_escal(Point* p, float t){
	Point* point_temp = new Point;

	point_temp->x = (float)((p->x)/t);
	point_temp->y = (float)((p->y)/t);
	point_temp->z = (float)((p->z)/t);

	return point_temp;
}

Point* prod_point_escal(Point* p, float t){
	Point* p_temp = new Point;

	p_temp->x = (p->x)*t;
	p_temp->y = (p->y)*t;
	p_temp->z = (p->z)*t;

	return p_temp;
}

Point* interpolacao(Point* p1, Point* p2, float t){

	Vector* vec_temp = new Vector;
	vec_temp->x = (p2->x - p1->x)*t;
	vec_temp->y = (p2->y - p1->y)*t;
	vec_temp->z = (p2->z - p1->z)*t;

	Point* p_t = new Point;
	p_t->x = p1->x + vec_temp->x;
	p_t->y = p1->y + vec_temp->y;
	p_t->z = p1->z + vec_temp->z;

	delete vec_temp;

	return p_t;
}

float prod_interno(Vector* v, Vector* w){
	return v->x * w->x + v->y * w->y + v->z * w->z;
}

Vector* prod_vectorial(Vector* v1, Vector* v2){
	Vector* r = new Vector;
	r->x = v1->y*v2->z - v1->z*v2->y;
	r->y = v1->z*v2->x - v1->x*v2->z;
	r->z = v1->x*v2->y - v1->y*v2->x;

	return r;
}

Camera* normalizar_camera(Camera* cam){
// Algoritmo de Gram-Schmmidt de Ortogonalização e normalização:

	Vector* _N = cam->N;
	_N = div_vector_escal(_N, norma(_N)); // N normalizado

	Vector* _V = cam->V;
	_V = div_vector_escal(_V, norma(_V));// V normalizado
	_V = sub_vector_vector(prod_vector_escal(_N, prod_interno(_V, _N)), _V);

	Vector* _U = prod_vectorial(_N, _V);

	return newCamera(cam->C, _U, _V, _N);
}

float **mult_pointer(float A[][30], int rows_A, int columns_A, float B[][30], int rows_B, int columns_B){
/*
Retorna A*B ou NULL, se as matrizes não são compatíveis com
a multiplicação.
*/
    float **R = NULL;

    if(columns_A == rows_B) {
        //As matrizes são compatíveis com a multiplicação:
        R = new float*[rows_A]; // Não esquecer de deletar: delete[] R;

        for(int i = 0; i < rows_A; i++) {
            R[i] = new float[columns_B];// Não esquecer de deletar: delete[] R[i];
        }

        for (int i = 0; i < rows_A; i++) {
            for (int j = 0; j < columns_B; j++) {
                float temp = 0.0;

                for(int k = 0; k < columns_A; k++) { // ou k < rows_B
                    temp = temp + A[i][k] * B[k][j];
                }

                R[i][j] = temp;

            }
        }
    } else {
        // As matrizes ão são compatíveis com a multiplicação
        // R = NULL;
    }

    return R;
}

Point* mudanca_coord_mundial_vista(Point* p, Camera* cam){
	Camera* cam_norm = new Camera;
	cam_norm = normalizar_camera(cam);

	//Vector* v_p = new Vector;
	Vector* v_p = sub_point_point(cam_norm->C, p);
	P[0][0] = v_p->x;
	P[1][0] = v_p->y;
	P[2][0] = v_p->z;


	//float** T = NULL;
	T[0][0] = cam_norm->U->x;
	T[0][1] = cam_norm->U->y;
	T[0][2] = cam_norm->U->z;

	T[1][0] = cam_norm->V->x;
	T[1][1] = cam_norm->V->y;
	T[1][2] = cam_norm->V->z;

	T[2][0] = cam_norm->N->x;
	T[2][1] = cam_norm->N->y;
	T[2][2] = cam_norm->N->z;

	float** T_A;
	T_A = mult_pointer(T, 3, 3, P, 3, 1);

	delete v_p;
	delete cam_norm;

	return newPoint(T_A[0][0], T_A[1][0], T_A[2][0]);
}

Point* projecao_tela(Point* p, float d, float hx, float hy){
/*
	hx e hy é o tamanho da tela.
	d é a distância entre a tela e a camera.
	p é um ponto em coordenadas de vista.
	Retorna um ponto projetado na tela normalizada.
*/
	//Projetando os ptos na tela
	float px = (float)(p->x*d)/p->z;
	float py = (float)(p->y*d)/p->z;

	if(px <= hx && px >= -hx && py <= hy && py >= -hy){
		ponto_fora_tela = false;
	}
	else {
		cout << "Ponto não pode ser plotado, pois excede o tamanho da tela." << endl;
		ponto_fora_tela = true;
	}

	//Retorna ptos projetados e normalizados na tela
	return newPoint(px/hx, py/hy, d);

}

Point2D* coord_tela_pixels(int Rx, int Ry, Point* p){
/*
	p é um ponto projetado na tela normalizada.
	Rx e Ry são as resoluções da tela (em PIXELS).
	Retorna as coordenadas de tela em pixels (arredondadas para inteiro).
*/

	return newPoint2D( (p->x + 1.0f)*Rx/2.0f, (1.0f - p->y)*Ry/2.0f );
}

Point* bezier_quadratic_bernstein(Point* b1, Point* b2, Point* b3, float t){

	return add_point_point(add_point_point(prod_point_escal(b1, pow((1 - t), 2)), prod_point_escal(b2, 2*(1 - t)*t )), prod_point_escal(b3, pow(t, 2)));
}

float N0_i(int i, float u, float* u_list){
/*
i = 1, 2, 3, ... n.
u é variálvel real.
u_list = {u0, u1, ..., u_k_r_u}.
*/
	if (i <= 0) return 0.0f;
	if (u >= u_list[i-1] && u < u_list[i]) return 1.0f;
	else return 0.0f;
}

float Nk_l(int k, int l, float u, float* u_list){
/*
retorna Nk_l em função de N(k-1)_l e N(k-1)_(l+1)
k = 1, 2, ...
l = 1, 2, ...
*/
	if (k == 1) return (u - u_list[l-1])/(u_list[l+k-1] - u_list[l-1]) * N0_i( l, u, u_list) + (u_list[l+k] - u)/(u_list[l+k] - u_list[l]) * N0_i(l+1, u, u_list);
	return (u - u_list[l-1])/(u_list[l+k-1] - u_list[l-1]) * Nk_l(k - 1, l, u, u_list) + (u_list[l+k] - u)/(u_list[l+k] - u_list[l]) * Nk_l(k - 1, l+1, u, u_list);
}



Point* nurbs(float u, float v, int n_u, int n_v, int c_u, int c_v, Point** d, float** w, float* u_list, float* v_list){
	Point* x_u_v = newPoint(0.0f, 0.0f, 0.0f);
	float coeficientes_d_i_j = 0;
	float soma_coeficientes_d_i_j = 0;

	for(int row = 0; row < c_v - 1; row++){
		for(int col = 0; col < c_u - 1; col++){
			coeficientes_d_i_j = w[row][col]*Nk_l(n_u, col, u, u_list)*Nk_l(n_v, row, v, v_list);
			soma_coeficientes_d_i_j = soma_coeficientes_d_i_j + coeficientes_d_i_j;

			Point* temp = add_point_point(x_u_v, prod_point_escal(&d[row][col], coeficientes_d_i_j));
			delete x_u_v;
			x_u_v = temp;
			temp = NULL;
		}
	}

	Point* temp = div_point_escal(x_u_v, soma_coeficientes_d_i_j);
	delete x_u_v;
	x_u_v = temp;
	temp = NULL;
	return x_u_v;
}

int main()
{
	// Parametros de entrada BSpline:
	int n_u = 3;
	int n_v = 3;
	int k_u;
	int k_v;

	int k_r_u = 5;
	int k_r_v = 5;
	//float u_list[k_r_u];
	//float v_list[k_r_v];

	int r_ui;
	int r_vi;
	// d_ij = (x_ij, y_ij, z_ij);
	// float w_ij;
	int c_u = 3;
	int c_v = 3;

	int p_u;
	int p_v;
	//P_ij = (x_ij, y_ij, z_ij); Ptos da superfície indexados.

	//u_list = {u0, u1, ..., u_k_r_u}.
	float* u_list = new float[k_r_u];
	u_list[0] = 1.2f;
	u_list[1] = 8.3f;
	u_list[2] = 9.7f;
	u_list[3] = 12.5f;
	u_list[4] = 14.0f;
	//v_list = {u0, u1, ..., u_k_r_v}.
	float* v_list = new float[k_r_v];
	v_list[0] = 0.21f;
	v_list[1] = 3.87f;
	v_list[2] = 7.34f;
	v_list[3] = 9.6f;
	v_list[4] = 11.4f;

	//d são os pontos de controle da malha:
	Point** d = NULL;
	d = new Point*[c_v]; // linhas
	for (int i = 0; i < c_v; i++){
		d[i] = new Point[c_u]; // colunas
	}

	d[0][0].x = 0;
	d[0][0].y = 0;
	d[0][0].z = 0;

	d[0][1].x = 0;
	d[0][1].y = 1;
	d[0][1].z = 1;

	d[0][2].x = 0;
	d[0][2].y = 2;
	d[0][2].z = 0;

	d[1][0].x = 1;
	d[1][0].y = 1;
	d[1][0].z = 2;

	d[1][1].x = 1;
	d[1][1].y = 2;
	d[1][1].z = -1;

	d[1][2].x = 1;
	d[1][2].y = 3;
	d[1][2].z = 1;

	d[2][0].x = -1;
	d[2][0].y = 0;
	d[2][0].z = 0;

	d[2][1].x = -1;
	d[2][1].y = 1;
	d[2][1].z = 1.5f;

	d[2][2].x = -1;
	d[2][2].y = 2;
	d[2][2].z = -1;



	//w são os pesos associados aos pontos de controle da malha:
	float** w = NULL;
	w = new float*[c_v]; // linhas
	for (int i = 0; i < c_v; i++){
		w[i] = new float[c_u]; // colunas
	}

	for(int row = 0; row < c_v; row++){
		for(int col = 0; col < c_u; col++){
			w[row][col] = 1;
		}
	}
	//cout << N0_i(0, 8.3f, u_list) << endl;
	//cout << Nk_l(3, 1, 8.35f, u_list) << endl;
//-----------------------------------------------------------------------------
	// Resoluções Rx e Ry (em pixels):
	int Rx = 512;
	int Ry = 512;

	sf::RenderWindow window(sf::VideoMode(Rx, Ry), "Bezier curve", sf::Style::Close | sf::Style::Resize);

	// Plotando ponto 2D
	//sf::Vertex point(sf::Vector2f(10, 10), sf::Color::White);

	float hx = 1.0f;
	float hy = 1.0f;
	float t = 0.0f;

	float u = u_list[n_u-1];
	float v = v_list[n_v-1];

//-*------------------------------------------------------------------------------
	while(window.isOpen())
	{
		// Curva de Bezier
		//Point* b2_0 = bezier_quadratic_bernstein(newPoint(1, 0, 1), newPoint(0,0,0), newPoint(1, 1, 1), t);
		//Point* b2_0 = bezier_quadratic_bernstein(newPoint(0, 0, 2), newPoint(0,0.5f,0), newPoint(0, 1, 2), t);
		//cout << "b2_0(" << t << "): " << b2_0->x << ", " << b2_0->y << ", " << b2_0->z << endl;

		Point* x_u_v = nurbs(u, v, n_u, n_v, c_u, c_v, d, w, u_list, v_list);

		//Procedimento para o plot
		//Camera* C = newCamera(newPoint(3.0f, 3.0f, 1.25f), NULL, newVector(0, 0, 1), newVector(-sqrt(2)/2, -sqrt(2)/2, 0));
		//Camera* C = newCamera(newPoint(2.0f, 0.0f, 0.0f), NULL, newVector(0, 0, 1), newVector(-1, 0, 0));
		//NURBS camera:
		Camera* C = newCamera(newPoint(2.0f, 2.0f, -1.0f), NULL, newVector(0, 0, 1), newVector(-1, 0, 0));
		//Camera* C = newCamera(newPoint(0.5f, 0.5f, 2.0f), NULL, newVector(0, 1, 0), newVector(0, 0, -1));
		C = normalizar_camera(C);

		//b2_0 = mudanca_coord_mundial_vista(b2_0, C);
		//b2_0 = projecao_tela(b2_0, 1, hx, hy);
		x_u_v = mudanca_coord_mundial_vista(x_u_v, C);
		x_u_v = projecao_tela(x_u_v, 1, hx, hy);

		// Plotando ponto 2D
		//Point2D* ponto = coord_tela_pixels(Rx, Ry, b2_0);
		//cout << "Proj_b2_0(" << t << "): " << ponto->x << ", " << ponto->y << endl;
		//sf::Vertex point(sf::Vector2f(ponto->x, ponto->y), sf::Color::White);

		Point2D* ponto = coord_tela_pixels(Rx, Ry, x_u_v);
		//cout << "Proj_b2_0(" << t << "): " << ponto->x << ", " << ponto->y << endl;
		sf::Vertex point(sf::Vector2f(ponto->x, ponto->y), sf::Color::White);


		sf::Event evnt;
		while(window.pollEvent(evnt))
		{
			switch(evnt.type)
			{
				case sf::Event::Closed:
					window.close();
					break;
				case sf::Event::Resized:
					//cout << "Window width: " << evnt.size.width << ", Window height: " << evnt.size.height << endl;
					break;
				case sf::Event::TextEntered:
					// Pegando entrada de texto do usuario (player)
					if(evnt.text.unicode <128){
						cout << evnt.text.unicode << endl;
					}
			}
		}

		if(sf::Keyboard::isKeyPressed(sf::Keyboard::Key::A)) {
			//t += 0.0001f;
			u += 0.0001f;
		}
		if(sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Z)) {
			//t -= 0.0001f;
			u -= 0.0001f;
		}
		if(sf::Keyboard::isKeyPressed(sf::Keyboard::Key::S)) {
			//t += 0.0001f;
			v += 0.0001f;
		}
		if(sf::Keyboard::isKeyPressed(sf::Keyboard::Key::X)) {
			//t -= 0.0001f;
			v -= 0.0001f;
		}

		//window.clear();
		window.draw(&point, 1, sf::Points);
		window.display();
		//window.close();
	}

	// Desalocando **d e **w:
	for (int i = 0; i < c_v; i++) {
			delete[] d[i];
			delete[] w[i];
	}
	delete[] d;
	delete[] w;

	// Desalocando u_list e v_list:
	delete[] u_list;
	delete[] v_list;

	return 0;
}
