#include "Render.h"
using namespace std;

#include <sstream>
#include <iostream>

#include <windows.h>
#include <GL\GL.h>
#include <GL\GLU.h>

#include "MyOGL.h"

#include "Camera.h"
#include "Light.h"
#include "Primitives.h"
#include <math.h>
#include <chrono>
#include <vector>

#include "GUItextRectangle.h"
#define PI 3.14159265358979323846

class Point { //класс, помогающий работать с матаном и облегчающий мне жизнь
public:
	double x;
	double y;
	double z;
	Point(double x, double y, double z) {
		this->x = x;
		this->y = y;
		this->z = z;
	}

	static Point SearchVector(Point A, Point B)
	{
		return Point(B.x - A.x, B.y - A.y, B.z - A.z);
	}

	//ищу длину вектора
	static double SearchVectorLength(Point vec)
	{
		double length = sqrt(pow(vec.x, 2) + pow(vec.y, 2) + pow(vec.z, 2));
		return length;
	}

	//нормализация вектора
	static Point VectorNormal(Point vec)
	{
		double length = SearchVectorLength(vec);

		vec.x = vec.x / length;
		vec.y = vec.y / length;
		vec.z = vec.z / length;

		return vec;
	}

	//векторное произведение
	static Point VectorProduct(Point vecA, Point vecB)
	{
		Point result(0, 0, 0);
		result.x = vecA.y * vecB.z - vecB.y * vecA.z;
		result.y = -1 * vecA.x * vecB.z + vecB.x * vecA.z;
		result.z = vecA.x * vecB.y - vecB.x * vecA.y;
		return result;
	}

	//скалярное произведение
	static double ScalarProduct(Point A, Point B)
	{
		return A.x * B.x + A.y * B.y + A.z * B.z;
	}
};

vector<vector<Point>> massiv =
{
{Point(0, 0, 0), Point(0, 1, 1), Point(0, 2, 1), Point(0, 3, 1), Point(0, 4, 0)},
{Point(1, 0, 1), Point(1, 1, 1), Point(1, 2, 1), Point(1, 3, 1), Point(1, 4, 1)},
{Point(2, 0, 1), Point(2, 1, 1), Point(2, 2, 1), Point(2, 3, 1), Point(2, 4, 1)},
{Point(3, 0, 1), Point(3, 1, 1), Point(3, 2, 1), Point(3, 3, 0), Point(3, 4, 1)},
{Point(4, 0, 0), Point(4, 1, 1), Point(4, 2, 1), Point(4, 3, 1), Point(4, 4, 0)}

};

bool textureMode = true;
bool lightMode = true;

//класс для настройки камеры
class CustomCamera : public Camera
{
public:
	//дистанция камеры
	double camDist;
	//углы поворота камеры
	double fi1, fi2;

	
	//значния масеры по умолчанию
	CustomCamera()
	{
		camDist = 15;
		fi1 = 1;
		fi2 = 1;
	}

	
	//считает позицию камеры, исходя из углов поворота, вызывается движком
	void SetUpCamera()
	{
		//отвечает за поворот камеры мышкой
		lookPoint.setCoords(0, 0, 0);

		pos.setCoords(camDist*cos(fi2)*cos(fi1),
			camDist*cos(fi2)*sin(fi1),
			camDist*sin(fi2));

		if (cos(fi2) <= 0)
			normal.setCoords(0, 0, -1);
		else
			normal.setCoords(0, 0, 1);

		LookAt();
	}

	void CustomCamera::LookAt()
	{
		//функция настройки камеры
		gluLookAt(pos.X(), pos.Y(), pos.Z(), lookPoint.X(), lookPoint.Y(), lookPoint.Z(), normal.X(), normal.Y(), normal.Z());
	}



}  camera;   //создаем объект камеры


//Класс для настройки света
class CustomLight : public Light
{
public:
	CustomLight()
	{
		//начальная позиция света
		pos = Vector3(1, 1, 3);
	}

	
	//рисует сферу и линии под источником света, вызывается движком
	void  DrawLightGhismo()
	{
		glDisable(GL_LIGHTING);

		
		glColor3d(0.9, 0.8, 0);
		Sphere s;
		s.pos = pos;
		s.scale = s.scale*0.08;
		s.Show();
		
		if (OpenGL::isKeyPressed('G'))
		{
			glColor3d(0, 0, 0);
			//линия от источника света до окружности
			glBegin(GL_LINES);
			glVertex3d(pos.X(), pos.Y(), pos.Z());
			glVertex3d(pos.X(), pos.Y(), 0);
			glEnd();

			//рисуем окруность
			Circle c;
			c.pos.setCoords(pos.X(), pos.Y(), 0);
			c.scale = c.scale*1.5;
			c.Show();
		}

	}

	void SetUpLight()
	{
		GLfloat amb[] = { 0.2, 0.2, 0.2, 0 };
		GLfloat dif[] = { 1.0, 1.0, 1.0, 0 };
		GLfloat spec[] = { .7, .7, .7, 0 };
		GLfloat position[] = { pos.X(), pos.Y(), pos.Z(), 1. };

		// параметры источника света
		glLightfv(GL_LIGHT0, GL_POSITION, position);
		// характеристики излучаемого света
		// фоновое освещение (рассеянный свет)
		glLightfv(GL_LIGHT0, GL_AMBIENT, amb);
		// диффузная составляющая света
		glLightfv(GL_LIGHT0, GL_DIFFUSE, dif);
		// зеркально отражаемая составляющая света
		glLightfv(GL_LIGHT0, GL_SPECULAR, spec);

		glEnable(GL_LIGHT0);
	}


} light;  //создаем источник света




//старые координаты мыши
int mouseX = 0, mouseY = 0;

void mouseEvent(OpenGL *ogl, int mX, int mY)
{
	int dx = mouseX - mX;
	int dy = mouseY - mY;
	mouseX = mX;
	mouseY = mY;

	//меняем углы камеры при нажатой левой кнопке мыши
	if (OpenGL::isKeyPressed(VK_RBUTTON))
	{
		camera.fi1 += 0.01*dx;
		camera.fi2 += -0.01*dy;
	}

	
	//двигаем свет по плоскости, в точку где мышь
	if (OpenGL::isKeyPressed('G') && !OpenGL::isKeyPressed(VK_LBUTTON))
	{
		LPPOINT POINT = new tagPOINT();
		GetCursorPos(POINT);
		ScreenToClient(ogl->getHwnd(), POINT);
		POINT->y = ogl->getHeight() - POINT->y;

		Ray r = camera.getLookRay(POINT->x, POINT->y);

		double z = light.pos.Z();

		double k = 0, x = 0, y = 0;
		if (r.direction.Z() == 0)
			k = 0;
		else
			k = (z - r.origin.Z()) / r.direction.Z();

		x = k*r.direction.X() + r.origin.X();
		y = k*r.direction.Y() + r.origin.Y();

		light.pos = Vector3(x, y, z);
	}

	if (OpenGL::isKeyPressed('G') && OpenGL::isKeyPressed(VK_LBUTTON))
	{
		light.pos = light.pos + Vector3(0, 0, 0.02*dy);
	}

	if (!OpenGL::isKeyPressed('G') && OpenGL::isKeyPressed(VK_LBUTTON))
	{
		LPPOINT point = new tagPOINT();
		GetCursorPos(point);
		ScreenToClient(ogl->getHwnd(), point);
		point->y = ogl->getHeight() - point->y;

		GLint viewport[5];
		GLdouble projection[25];
		GLdouble modelview[25];

		glGetIntegerv(GL_VIEWPORT, viewport);
		glGetDoublev(GL_PROJECTION_MATRIX, projection);
		glGetDoublev(GL_MODELVIEW_MATRIX, modelview);

		double delta = 10;

		for (auto& v : massiv)
		{
			for (auto& elem : v)
			{
				double Tpoint[3];
				gluProject(elem.x, elem.y, elem.z, modelview, projection, viewport, &Tpoint[0], &Tpoint[1], &Tpoint[2]);
				if (Tpoint[0] > point->x - delta && Tpoint[0]<point->x + delta && Tpoint[1]>point->y - delta && Tpoint[1] < point->y + delta)
				{
					Tpoint[0] -= dx;
					Tpoint[1] += dy;

					gluUnProject(Tpoint[0], Tpoint[1], Tpoint[2], modelview, projection, viewport, &elem.x, &elem.y, &elem.z);
				}
			}
		}
	}


	
}

void mouseWheelEvent(OpenGL *ogl, int delta)
{

	if (delta < 0 && camera.camDist <= 1)
		return;
	if (delta > 0 && camera.camDist >= 100)
		return;

	camera.camDist += 0.01*delta;

}

void keyDownEvent(OpenGL *ogl, int key)
{
	if (key == 'L')
	{
		lightMode = !lightMode;
	}

	if (key == 'T')
	{
		textureMode = !textureMode;
	}

	if (key == 'R')
	{
		camera.fi1 = 1;
		camera.fi2 = 1;
		camera.camDist = 15;

		light.pos = Vector3(1, 1, 3);
	}

	if (key == 'F')
	{
		light.pos = camera.pos;
	}
}

void keyUpEvent(OpenGL *ogl, int key)
{
	
}



std::string NameTexture = "Sneg.bmp";
GLuint texId;

//выполняется перед первым рендером
void initRender(OpenGL *ogl)
{
	//настройка текстур

	//4 байта на хранение пикселя
	glPixelStorei(GL_UNPACK_ALIGNMENT, 4);

	//настройка режима наложения текстур
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	//включаем текстуры
	glEnable(GL_TEXTURE_2D);
	

	//массив трехбайтных элементов  (R G B)
	RGBTRIPLE *texarray;

	//массив символов, (высота*ширина*4      4, потомучто   выше, мы указали использовать по 4 байта на пиксель текстуры - R G B A)
	char *texCharArray;
	int texW, texH;
	const char* name = NameTexture.c_str();
	OpenGL::LoadBMP(name, &texW, &texH, &texarray);
	OpenGL::RGBtoChar(texarray, texW, texH, &texCharArray);

	
	
	//генерируем ИД для текстуры
	glGenTextures(1, &texId);
	//биндим айдишник, все что будет происходить с текстурой, будте происходить по этому ИД
	glBindTexture(GL_TEXTURE_2D, texId);

	//загружаем текстуру в видеопямять, в оперативке нам больше  она не нужна
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, texW, texH, 0, GL_RGBA, GL_UNSIGNED_BYTE, texCharArray);

	//отчистка памяти
	free(texCharArray);
	free(texarray);

	//наводим шмон
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST); 
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);


	//камеру и свет привязываем к "движку"
	ogl->mainCamera = &camera;
	ogl->mainLight = &light;

	// нормализация нормалей : их длины будет равна 1
	glEnable(GL_NORMALIZE);

	// устранение ступенчатости для линий
	glEnable(GL_LINE_SMOOTH); 


	//   задать параметры освещения
	//  параметр GL_LIGHT_MODEL_TWO_SIDE - 
	//                0 -  лицевые и изнаночные рисуются одинаково(по умолчанию), 
	//                1 - лицевые и изнаночные обрабатываются разными режимами       
	//                соответственно лицевым и изнаночным свойствам материалов.    
	//  параметр GL_LIGHT_MODEL_AMBIENT - задать фоновое освещение, 
	//                не зависящее от сточников
	// по умолчанию (0.2, 0.2, 0.2, 1.0)

	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 0);

	camera.fi1 = -1.3;
	camera.fi2 = 0.8;
}




//формула Базье
double funBazie(double p0, double p1, double p2, double p3, double t)
{
	return pow((1 - t), 3) * p0 + 3 * t * pow((1 - t), 2) * p1 + 3 * pow(t, 2) * (1 - t) * p2 + pow(t, 3) * p3;
}

//функция Эрмиту
double funErmit(double p1, double p4, double r1, double r4, double t)
{
	return p1 * (2 * pow(t, 3) - 3 * pow(t, 2) + 1) + p4 * (3 * pow(t, 2) - 2 * pow(t, 3)) + r1 * (pow(t, 3) - 2 * pow(t, 2) + t) + r4 * (pow(t, 3) - pow(t, 2));
}

//расчёт вектора
void CalculationVec(double x1[3], double x2[3], double x3[3], double x4[3], double* v1, double* v4)
{
	v1[0] = 3 * (x2[0] - x1[0]);
	v1[1] = 3 * (x2[1] - x1[1]);
	v1[2] = 3 * (x2[2] - x1[2]);

	v4[0] = 3 * (x4[0] - x3[0]);
	v4[1] = 3 * (x4[1] - x3[1]);
	v4[2] = 3 * (x4[2] - x3[2]);
}

//рисование кривой Эрмита
void DrawErmit(double* P1, double* P2, double* P3, double* P4)
{
	double R1[3], R4[3];
	CalculationVec(P1, P2, P3, P4, R1, R4);

	glLineWidth(3); //ширина линии

	glColor3d(0, 0, 1);
	glBegin(GL_LINE_STRIP);
	for (double t = 0; t <= 1.0001; t += 0.01)
	{
		double P[3];

		P[0] = funErmit(P1[0], P4[0], R1[0], R4[0], t);
		P[1] = funErmit(P1[1], P4[1], R1[1], R4[1], t);
		P[2] = funErmit(P1[2], P4[2], R1[2], R4[2], t);
		glVertex3dv(P);
	}
	glEnd();
	glColor3d(1, 0, 1);
	glLineWidth(1);

	glBegin(GL_LINES); //построим касательные вектора в точках, когда t=0 и t=1
	glVertex3dv(P1);
	glVertex3dv(R1);

	glVertex3dv(P4);
	glVertex3dv(R4);
	glEnd();
}

//отрисовка кривой Базье
void DrawBazie(double* P0, double* P1, double* P2, double* P3)
{
	glLineWidth(3); //ширина линии

	glColor3d(0, 0, 0);
	glBegin(GL_LINE_STRIP);
	for (double t = 0; t <= 1.0001; t += 0.01)
	{
		double P[3];

		P[0] = funBazie(P0[0], P1[0], P2[0], P3[0], t);
		P[1] = funBazie(P0[1], P1[1], P2[1], P3[1], t);
		P[2] = funBazie(P0[2], P1[2], P2[2], P3[2], t);
		glVertex3dv(P);
	}
	glEnd();
	glColor3d(1, 0, 1);

	glLineWidth(1);

	glBegin(GL_LINES); //построим отрезки P0P1 и P1P2 P2P3
	glVertex3dv(P0);
	glVertex3dv(P1);

	glVertex3dv(P1);
	glVertex3dv(P2);

	glVertex3dv(P2);
	glVertex3dv(P3);
	glEnd();

}

//вращение фигуры
void MoveHouse(Point point, Point next_point)
{
	Point dir = Point::SearchVector(point, next_point); //ищу вектор по двум точкам
	dir = Point::VectorNormal(dir); //нормализирую вектор: нахожу длину и делю каждую координату на длину вектора

	Point orig(1, 0, 0);
	Point rotX(dir.x, dir.y, 0);
	rotX = Point::VectorNormal(rotX);
	double cosU = Point::ScalarProduct(orig, rotX); //ищу скалярное произведение
	Point vecpr = Point::VectorProduct(orig, rotX); //ищу векторное произведение

	double sinSign = vecpr.z / abs(vecpr.z);
	double U = acos(cosU) * 180.0 / PI * sinSign;
	double ZU = acos(dir.z) * 180.0 / PI - 90;

	//выполняю поворот осей
	glRotated(U, 0, 0, 1);
	glRotated(ZU, 0, 1, 0);

}

void SearchStep(double t_max, double next, double* P0, double* P1, double* P2, double* P3)
{
	double P[3];
	P[0] = funBazie(P0[0], P1[0], P2[0], P3[0], t_max);
	P[1] = funBazie(P0[1], P1[1], P2[1], P3[1], t_max);
	P[2] = funBazie(P0[2], P1[2], P2[2], P3[2], t_max);

	glTranslated(P[0], P[1], P[2]);

	double Next[3];
	Next[0] = funBazie(P0[0], P1[0], P2[0], P3[0], next);
	Next[1] = funBazie(P0[1], P1[1], P2[1], P3[1], next);
	Next[2] = funBazie(P0[2], P1[2], P2[2], P3[2], next);



	MoveHouse(Point(P[0], P[1], P[2]), Point(Next[0], Next[1], Next[2])); //вращаю мой домик
}

//факториал
int f(int M)
{
	if (M == 0)
		return 1;
	if (M < 0)
		return 0;
	else return M * f(M - 1);
}

//базис Бернштейна
double BasisBernstein(int m, int i, double t)
{
	return (f(m) / (f(m - i) * f(i))) * pow(t, i) * pow(1 - t, m - i);
}

//нахожу радиус вектора
Point RadiusOfVector(vector<vector<Point>> massiv, double u, double v)
{
	Point Rad(0, 0, 0);
	int Row = (int)massiv.size();
	for (int i = 0; i < Row; i++)
	{
		int Column = (int)massiv[i].size();

		for (int j = 0; j < Column; j++) {
			double Bi = BasisBernstein(Row - 1, i, u);
			double Bj = BasisBernstein(Column - 1, j, v);
			Rad.x += Bj * Bi * massiv[i][j].x;
			Rad.y += Bj * Bi * massiv[i][j].y;
			Rad.z += Bj * Bi * massiv[i][j].z;
		}
	}
	return Rad;
}

void Normal(Point first, Point second, Point third)
{
	//ищу два вектора, кторые исходят из одной точки
	Point vector_first = Point::SearchVector(second, first);
	Point vector_second = Point::SearchVector(second, third);

	Point MyNormal = Point::VectorProduct(vector_first, vector_second);
	//double Lenght = Point::SearchVectorLength(MyNormal);

	glNormal3d(MyNormal.x, MyNormal.y, MyNormal.z);
}

//поверхность Базье
void DrawBazieSurface(vector<vector<Point>> massiv)
{

	double t_max = 1.0001;
	vector<Point> points;
	vector<vector<Point>> MassivePoints;
	for (double u = 0; u <= t_max; u += 0.01)
	{
		for (double v = 0; v <= t_max; v += 0.01)
		{
			points.push_back(RadiusOfVector(massiv, u, v));
		}
		MassivePoints.push_back(points);
		points.clear();
	}

	glBindTexture(GL_TEXTURE_2D, texId);

	glBegin(GL_TEXTURE_2D);

	glBegin(GL_QUADS);

	for (int i = 0; i < MassivePoints.size() - 1; i++)
	{
		for (int j = 0; j < MassivePoints[i].size() - 1; j++)
		{
			Point A = MassivePoints[i][j];
			Point B = MassivePoints[i][j + 1];
			Point D = MassivePoints[i + 1][j];
			Point C = MassivePoints[i + 1][j+1];
		

			double R = 1./MassivePoints.size();
			double Col = 1./MassivePoints[i].size(); 

			double A_Sneg[2] = { R*i, Col*j };
			double B_Sneg[2] = { R * i, Col * (j + 1) };
			double C_Sneg[2] = { R * (i + 1), Col * (j + 1) };
			double D_Sneg[2] = { R * (i + 1), Col * j };

			Normal(A, B, C); //ищу нормаль к треугольнику

			glColor3d(0, 1, 0);
			glTexCoord2dv(A_Sneg);
			glVertex3d(A.x, A.y, A.z);
			glTexCoord2dv(B_Sneg);
			glVertex3d(B.x, B.y, B.z);
			glTexCoord2dv(C_Sneg);
			glVertex3d(C.x, C.y, C.z);
			glTexCoord2dv(D_Sneg);
			glVertex3d(D.x, D.y, D.z);
		}
	}
	glEnd();
	glDisable(GL_TEXTURE_2D);

}

//поверхность Базье из точек
void DrawBazieSurfaceWithPoints(vector<vector<Point>> massiv)
{
	glPointSize(12);
	glColor3d(0, 1, 1);

	glBegin(GL_POINTS);
	for (int i = 0; i < (int)massiv.size(); i++)
	{
		for (int j = 0; j < (int)massiv[i].size(); j++)
		{
			glVertex3d(massiv[i][j].x, massiv[i][j].y, massiv[i][j].z);
		}

	}
	glEnd();

}

//поверхность Базье из линий
void DrawBazieSurfaceWithLine(vector<vector<Point>> massiv)
{
	glLineWidth(2);
	glColor3d(0, 0, 0);

	glBegin(GL_LINES);
	for (int i = 0; i < (int)massiv.size(); i++)
	{
		for (int j = 0; j < (int)massiv[i].size() - 1; j++)
		{
			glVertex3d(massiv[i][j].x, massiv[i][j].y, massiv[i][j].z);
			glVertex3d(massiv[i][j + 1].x, massiv[i][j + 1].y, massiv[i][j + 1].z);
		}
	}
	for (int i = 0; i < (int)massiv.size() - 1; i++)
	{
		for (int j = 0; j < (int)massiv[i].size(); j++)
		{
			glVertex3d(massiv[i][j].x, massiv[i][j].y, massiv[i][j].z);
			glVertex3d(massiv[i + 1][j].x, massiv[i + 1][j].y, massiv[i + 1][j].z);
		}
	}

	glEnd();
	glLineWidth(1);
}

//рисую мы милый и ,специально, не кривой, а несимметричный домик
void DrawHouse()
{
	double A[] = { 3, 0.5, 1 };
	double B[] = { 0, 0.5, 2 };
	double C[] = { -2, 0.5, 1 };
	double D[] = { -1, 0.5, 1 };
	double E[] = { -1, 0.5, -1 };
	double F[] = { 2, 0.5, -1 };
	double G[] = { 2, 0.5, 1 };

	double A1[] = { 3, -0.7, 1 };
	double B1[] = { 0, -0.7, 2 };
	double C1[] = { -2, -0.7, 1 };
	double D1[] = { -1, -0.7, 1 };
	double E1[] = { -1, -0.7, -1 };
	double F1[] = { 2, -0.7, -1 };
	double G1[] = { 2, -0.7, 1 };

	glBegin(GL_TRIANGLES);
	glColor3d(0, 1, 1);
	glVertex3dv(A);
	glVertex3dv(B);
	glVertex3dv(G);

	glVertex3dv(B);
	glVertex3dv(G);
	glVertex3dv(D);

	glVertex3dv(B);
	glVertex3dv(D);
	glVertex3dv(C);

	glColor3d(0, 1, 0);
	glVertex3dv(G);
	glVertex3dv(F);
	glVertex3dv(D);

	glVertex3dv(D);
	glVertex3dv(F);
	glVertex3dv(E);
	glEnd();

	glBegin(GL_TRIANGLES);
	glColor3d(0, 1, 1);
	glVertex3dv(A1);
	glVertex3dv(B1);
	glVertex3dv(G1);

	glVertex3dv(B1);
	glVertex3dv(G1);
	glVertex3dv(D1);

	glVertex3dv(B1);
	glVertex3dv(D1);
	glVertex3dv(C1);

	glColor3d(0, 1, 0);
	glVertex3dv(G1);
	glVertex3dv(F1);
	glVertex3dv(D1);

	glVertex3dv(D1);
	glVertex3dv(F1);
	glVertex3dv(E1);
	glEnd();

	glBegin(GL_QUADS);
	glColor3d(0, 1, 1);
	glVertex3dv(A);
	glVertex3dv(A1);
	glVertex3dv(B1);
	glVertex3dv(B);

	glVertex3dv(B);
	glVertex3dv(B1);
	glVertex3dv(C1);
	glVertex3dv(C);

	glVertex3dv(C);
	glVertex3dv(C1);
	glVertex3dv(D1);
	glVertex3dv(D);

	glVertex3dv(A);
	glVertex3dv(A1);
	glVertex3dv(G1);
	glVertex3dv(G);

	glColor3d(0, 1, 0);
	glVertex3dv(G);
	glVertex3dv(G1);
	glVertex3dv(F1);
	glVertex3dv(F);

	glVertex3dv(F);
	glVertex3dv(F1);
	glVertex3dv(E1);
	glVertex3dv(E);

	glVertex3dv(E);
	glVertex3dv(E1);
	glVertex3dv(D1);
	glVertex3dv(D);
	glEnd();
}


double DeltaTime()
{
	static auto end = std::chrono::steady_clock::now();
	auto time = std::chrono::steady_clock::now();
	auto deltatime = time - end;
	double delta = 1.0 * std::chrono::duration_cast<std::chrono::microseconds>(deltatime).count() / 1000000;
	end = time;
	return delta;
}

double t_max = 0;
	double next_t_max = 0;
	bool flag_t_max = true;

void Render(OpenGL *ogl)
{


	glDisable(GL_TEXTURE_2D);
	glDisable(GL_LIGHTING);

	glEnable(GL_DEPTH_TEST);
	if (textureMode)
		glEnable(GL_TEXTURE_2D);

	if (lightMode)
		glEnable(GL_LIGHTING);


	//альфаналожение
	//glEnable(GL_BLEND);
	//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


	//настройка материала
	GLfloat amb[] = { 0.2, 0.2, 0.1, 1. };
	GLfloat dif[] = { 0.4, 0.65, 0.5, 1. };
	GLfloat spec[] = { 0.9, 0.8, 0.3, 1. };
	GLfloat sh = 0.1f * 256;


	//фоновая
	glMaterialfv(GL_FRONT, GL_AMBIENT, amb);
	//дифузная
	glMaterialfv(GL_FRONT, GL_DIFFUSE, dif);
	//зеркальная
	glMaterialfv(GL_FRONT, GL_SPECULAR, spec); \
		//размер блика
		glMaterialf(GL_FRONT, GL_SHININESS, sh);

	//чтоб было красиво, без квадратиков (сглаживание освещения)
	glShadeModel(GL_SMOOTH);

	

	//настраиваю время
	double delta_time = DeltaTime();
	double step = delta_time / 5; //t_max становится = 1 за 5 секунд
	//t_max сама по себе изменяется от 0 до 1 постепенно от кадра к кадру
	if (flag_t_max) {
		t_max += step;
		next_t_max = t_max + step;
		if (t_max > 1) {
			t_max = 1;
			flag_t_max = false;
		}
		if (next_t_max > 1) {
			next_t_max = 1;
		}
	}
	else {
		t_max -= step;
		next_t_max = t_max - step;
		if (t_max < 0) {
			t_max = 0;
			flag_t_max = true;
		}
		if (next_t_max < 0) {
			next_t_max = 0;
		}
	}

	double P1_1[] = { 0,0,0 };
	double P2_1[] = { 10,10,15 };
	double P3_1[] = { -5,-7, 20 };
	double P4_1[] = { 9,0,0 };

	//DrawErmit(P1_1, P2_1, P3_1, P4_1); //Кривая Эрмита 1

	double P1_2[] = { 0,0,0 };
	double P2_2[] = { 8,10,7 };
	double P3_2[] = { 15,9,8 };
	double P4_2[] = { 9,9,0 };

	//DrawErmit(P1_2, P2_2, P3_2, P4_2); //Кривая Эрмита 2

	double P0_3[] = { 0, 0, 0 };
	double P1_3[] = { -5,-10, -2 };
	double P2_3[] = { -8, 10, 2 };
	double P3_3[] = { 15 , 15 ,0 };

	//DrawBazie(P0_3, P1_3, P2_3, P3_3); //Кривая Базье 1

	double P0_4[] = { 0, 0, 0 };
	double P1_4[] = { 7,4,1 };
	double P2_4[] = { -15,10,12 };
	double P3_4[] = { 3,-1,1 };

    //DrawBazie(P0_4, P1_4, P2_4, P3_4); //Кривая Базье 2

	glPushMatrix();
	//SearchStep(t_max, next_t_max, P0_4, P1_4, P2_4, P3_4);
	//DrawHouse();
	glPopMatrix();

	DrawBazieSurfaceWithPoints(massiv); //ура, рисую поверхность Базье из точек
	DrawBazieSurfaceWithLine(massiv); //ура, рисую поверхность Базье из линий
	DrawBazieSurface(massiv); //ура, рисую поверхность Базье


	


   //Сообщение вверху экрана

	
	glMatrixMode(GL_PROJECTION);	//Делаем активной матрицу проекций. 
	                                //(всек матричные операции, будут ее видоизменять.)
	glPushMatrix();   //сохраняем текущую матрицу проецирования (которая описывает перспективную проекцию) в стек 				    
	glLoadIdentity();	  //Загружаем единичную матрицу
	glOrtho(0, ogl->getWidth(), 0, ogl->getHeight(), 0, 1);	 //врубаем режим ортогональной проекции

	glMatrixMode(GL_MODELVIEW);		//переключаемся на модел-вью матрицу
	glPushMatrix();			  //сохраняем текущую матрицу в стек (положение камеры, фактически)
	glLoadIdentity();		  //сбрасываем ее в дефолт

	glDisable(GL_LIGHTING);



	GuiTextRectangle rec;		   //классик моего авторства для удобной работы с рендером текста.
	rec.setSize(300, 150);
	rec.setPosition(10, ogl->getHeight() - 150 - 10);


	std::stringstream ss;
	ss << "T - вкл/выкл текстур" << std::endl;
	ss << "L - вкл/выкл освещение" << std::endl;
	ss << "F - Свет из камеры" << std::endl;
	ss << "G - двигать свет по горизонтали" << std::endl;
	ss << "G+ЛКМ двигать свет по вертекали" << std::endl;
	ss << "Коорд. света: (" << light.pos.X() << ", " << light.pos.Y() << ", " << light.pos.Z() << ")" << std::endl;
	ss << "Коорд. камеры: (" << camera.pos.X() << ", " << camera.pos.Y() << ", " << camera.pos.Z() << ")" << std::endl;
	ss << "Параметры камеры: R="  << camera.camDist << ", fi1=" << camera.fi1 << ", fi2=" << camera.fi2 << std::endl;
	
	rec.setText(ss.str().c_str());
	rec.Draw();

	glMatrixMode(GL_PROJECTION);	  //восстанавливаем матрицы проекции и модел-вью обратьно из стека.
	glPopMatrix();


	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();

}