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

class Point { //�����, ���������� �������� � ������� � ����������� ��� �����
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

	//��� ����� �������
	static double SearchVectorLength(Point vec)
	{
		double length = sqrt(pow(vec.x, 2) + pow(vec.y, 2) + pow(vec.z, 2));
		return length;
	}

	//������������ �������
	static Point VectorNormal(Point vec)
	{
		double length = SearchVectorLength(vec);

		vec.x = vec.x / length;
		vec.y = vec.y / length;
		vec.z = vec.z / length;

		return vec;
	}

	//��������� ������������
	static Point VectorProduct(Point vecA, Point vecB)
	{
		Point result(0, 0, 0);
		result.x = vecA.y * vecB.z - vecB.y * vecA.z;
		result.y = -1 * vecA.x * vecB.z + vecB.x * vecA.z;
		result.z = vecA.x * vecB.y - vecB.x * vecA.y;
		return result;
	}

	//��������� ������������
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

//����� ��� ��������� ������
class CustomCamera : public Camera
{
public:
	//��������� ������
	double camDist;
	//���� �������� ������
	double fi1, fi2;

	
	//������� ������ �� ���������
	CustomCamera()
	{
		camDist = 15;
		fi1 = 1;
		fi2 = 1;
	}

	
	//������� ������� ������, ������ �� ����� ��������, ���������� �������
	void SetUpCamera()
	{
		//�������� �� ������� ������ ������
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
		//������� ��������� ������
		gluLookAt(pos.X(), pos.Y(), pos.Z(), lookPoint.X(), lookPoint.Y(), lookPoint.Z(), normal.X(), normal.Y(), normal.Z());
	}



}  camera;   //������� ������ ������


//����� ��� ��������� �����
class CustomLight : public Light
{
public:
	CustomLight()
	{
		//��������� ������� �����
		pos = Vector3(1, 1, 3);
	}

	
	//������ ����� � ����� ��� ���������� �����, ���������� �������
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
			//����� �� ��������� ����� �� ����������
			glBegin(GL_LINES);
			glVertex3d(pos.X(), pos.Y(), pos.Z());
			glVertex3d(pos.X(), pos.Y(), 0);
			glEnd();

			//������ ���������
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

		// ��������� ��������� �����
		glLightfv(GL_LIGHT0, GL_POSITION, position);
		// �������������� ����������� �����
		// ������� ��������� (���������� ����)
		glLightfv(GL_LIGHT0, GL_AMBIENT, amb);
		// ��������� ������������ �����
		glLightfv(GL_LIGHT0, GL_DIFFUSE, dif);
		// ��������� ���������� ������������ �����
		glLightfv(GL_LIGHT0, GL_SPECULAR, spec);

		glEnable(GL_LIGHT0);
	}


} light;  //������� �������� �����




//������ ���������� ����
int mouseX = 0, mouseY = 0;

void mouseEvent(OpenGL *ogl, int mX, int mY)
{
	int dx = mouseX - mX;
	int dy = mouseY - mY;
	mouseX = mX;
	mouseY = mY;

	//������ ���� ������ ��� ������� ����� ������ ����
	if (OpenGL::isKeyPressed(VK_RBUTTON))
	{
		camera.fi1 += 0.01*dx;
		camera.fi2 += -0.01*dy;
	}

	
	//������� ���� �� ���������, � ����� ��� ����
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

//����������� ����� ������ ��������
void initRender(OpenGL *ogl)
{
	//��������� �������

	//4 ����� �� �������� �������
	glPixelStorei(GL_UNPACK_ALIGNMENT, 4);

	//��������� ������ ��������� �������
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	//�������� ��������
	glEnable(GL_TEXTURE_2D);
	

	//������ ����������� ���������  (R G B)
	RGBTRIPLE *texarray;

	//������ ��������, (������*������*4      4, ���������   ����, �� ������� ������������ �� 4 ����� �� ������� �������� - R G B A)
	char *texCharArray;
	int texW, texH;
	const char* name = NameTexture.c_str();
	OpenGL::LoadBMP(name, &texW, &texH, &texarray);
	OpenGL::RGBtoChar(texarray, texW, texH, &texCharArray);

	
	
	//���������� �� ��� ��������
	glGenTextures(1, &texId);
	//������ ��������, ��� ��� ����� ����������� � ���������, ����� ����������� �� ����� ��
	glBindTexture(GL_TEXTURE_2D, texId);

	//��������� �������� � �����������, � ���������� ��� ������  ��� �� �����
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, texW, texH, 0, GL_RGBA, GL_UNSIGNED_BYTE, texCharArray);

	//�������� ������
	free(texCharArray);
	free(texarray);

	//������� ����
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST); 
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);


	//������ � ���� ����������� � "������"
	ogl->mainCamera = &camera;
	ogl->mainLight = &light;

	// ������������ �������� : �� ����� ����� ����� 1
	glEnable(GL_NORMALIZE);

	// ���������� ������������� ��� �����
	glEnable(GL_LINE_SMOOTH); 


	//   ������ ��������� ���������
	//  �������� GL_LIGHT_MODEL_TWO_SIDE - 
	//                0 -  ������� � ���������� �������� ���������(�� ���������), 
	//                1 - ������� � ���������� �������������� ������� ��������       
	//                �������������� ������� � ���������� ��������� ����������.    
	//  �������� GL_LIGHT_MODEL_AMBIENT - ������ ������� ���������, 
	//                �� ��������� �� ���������
	// �� ��������� (0.2, 0.2, 0.2, 1.0)

	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 0);

	camera.fi1 = -1.3;
	camera.fi2 = 0.8;
}




//������� �����
double funBazie(double p0, double p1, double p2, double p3, double t)
{
	return pow((1 - t), 3) * p0 + 3 * t * pow((1 - t), 2) * p1 + 3 * pow(t, 2) * (1 - t) * p2 + pow(t, 3) * p3;
}

//������� ������
double funErmit(double p1, double p4, double r1, double r4, double t)
{
	return p1 * (2 * pow(t, 3) - 3 * pow(t, 2) + 1) + p4 * (3 * pow(t, 2) - 2 * pow(t, 3)) + r1 * (pow(t, 3) - 2 * pow(t, 2) + t) + r4 * (pow(t, 3) - pow(t, 2));
}

//������ �������
void CalculationVec(double x1[3], double x2[3], double x3[3], double x4[3], double* v1, double* v4)
{
	v1[0] = 3 * (x2[0] - x1[0]);
	v1[1] = 3 * (x2[1] - x1[1]);
	v1[2] = 3 * (x2[2] - x1[2]);

	v4[0] = 3 * (x4[0] - x3[0]);
	v4[1] = 3 * (x4[1] - x3[1]);
	v4[2] = 3 * (x4[2] - x3[2]);
}

//��������� ������ ������
void DrawErmit(double* P1, double* P2, double* P3, double* P4)
{
	double R1[3], R4[3];
	CalculationVec(P1, P2, P3, P4, R1, R4);

	glLineWidth(3); //������ �����

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

	glBegin(GL_LINES); //�������� ����������� ������� � ������, ����� t=0 � t=1
	glVertex3dv(P1);
	glVertex3dv(R1);

	glVertex3dv(P4);
	glVertex3dv(R4);
	glEnd();
}

//��������� ������ �����
void DrawBazie(double* P0, double* P1, double* P2, double* P3)
{
	glLineWidth(3); //������ �����

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

	glBegin(GL_LINES); //�������� ������� P0P1 � P1P2 P2P3
	glVertex3dv(P0);
	glVertex3dv(P1);

	glVertex3dv(P1);
	glVertex3dv(P2);

	glVertex3dv(P2);
	glVertex3dv(P3);
	glEnd();

}

//�������� ������
void MoveHouse(Point point, Point next_point)
{
	Point dir = Point::SearchVector(point, next_point); //��� ������ �� ���� ������
	dir = Point::VectorNormal(dir); //������������ ������: ������ ����� � ���� ������ ���������� �� ����� �������

	Point orig(1, 0, 0);
	Point rotX(dir.x, dir.y, 0);
	rotX = Point::VectorNormal(rotX);
	double cosU = Point::ScalarProduct(orig, rotX); //��� ��������� ������������
	Point vecpr = Point::VectorProduct(orig, rotX); //��� ��������� ������������

	double sinSign = vecpr.z / abs(vecpr.z);
	double U = acos(cosU) * 180.0 / PI * sinSign;
	double ZU = acos(dir.z) * 180.0 / PI - 90;

	//�������� ������� ����
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



	MoveHouse(Point(P[0], P[1], P[2]), Point(Next[0], Next[1], Next[2])); //������ ��� �����
}

//���������
int f(int M)
{
	if (M == 0)
		return 1;
	if (M < 0)
		return 0;
	else return M * f(M - 1);
}

//����� ����������
double BasisBernstein(int m, int i, double t)
{
	return (f(m) / (f(m - i) * f(i))) * pow(t, i) * pow(1 - t, m - i);
}

//������ ������ �������
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
	//��� ��� �������, ������ ������� �� ����� �����
	Point vector_first = Point::SearchVector(second, first);
	Point vector_second = Point::SearchVector(second, third);

	Point MyNormal = Point::VectorProduct(vector_first, vector_second);
	//double Lenght = Point::SearchVectorLength(MyNormal);

	glNormal3d(MyNormal.x, MyNormal.y, MyNormal.z);
}

//����������� �����
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

			Normal(A, B, C); //��� ������� � ������������

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

//����������� ����� �� �����
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

//����������� ����� �� �����
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

//����� �� ����� � ,����������, �� ������, � �������������� �����
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


	//��������������
	//glEnable(GL_BLEND);
	//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


	//��������� ���������
	GLfloat amb[] = { 0.2, 0.2, 0.1, 1. };
	GLfloat dif[] = { 0.4, 0.65, 0.5, 1. };
	GLfloat spec[] = { 0.9, 0.8, 0.3, 1. };
	GLfloat sh = 0.1f * 256;


	//�������
	glMaterialfv(GL_FRONT, GL_AMBIENT, amb);
	//��������
	glMaterialfv(GL_FRONT, GL_DIFFUSE, dif);
	//����������
	glMaterialfv(GL_FRONT, GL_SPECULAR, spec); \
		//������ �����
		glMaterialf(GL_FRONT, GL_SHININESS, sh);

	//���� ���� �������, ��� ����������� (����������� ���������)
	glShadeModel(GL_SMOOTH);

	

	//���������� �����
	double delta_time = DeltaTime();
	double step = delta_time / 5; //t_max ���������� = 1 �� 5 ������
	//t_max ���� �� ���� ���������� �� 0 �� 1 ���������� �� ����� � �����
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

	//DrawErmit(P1_1, P2_1, P3_1, P4_1); //������ ������ 1

	double P1_2[] = { 0,0,0 };
	double P2_2[] = { 8,10,7 };
	double P3_2[] = { 15,9,8 };
	double P4_2[] = { 9,9,0 };

	//DrawErmit(P1_2, P2_2, P3_2, P4_2); //������ ������ 2

	double P0_3[] = { 0, 0, 0 };
	double P1_3[] = { -5,-10, -2 };
	double P2_3[] = { -8, 10, 2 };
	double P3_3[] = { 15 , 15 ,0 };

	//DrawBazie(P0_3, P1_3, P2_3, P3_3); //������ ����� 1

	double P0_4[] = { 0, 0, 0 };
	double P1_4[] = { 7,4,1 };
	double P2_4[] = { -15,10,12 };
	double P3_4[] = { 3,-1,1 };

    //DrawBazie(P0_4, P1_4, P2_4, P3_4); //������ ����� 2

	glPushMatrix();
	//SearchStep(t_max, next_t_max, P0_4, P1_4, P2_4, P3_4);
	//DrawHouse();
	glPopMatrix();

	DrawBazieSurfaceWithPoints(massiv); //���, ����� ����������� ����� �� �����
	DrawBazieSurfaceWithLine(massiv); //���, ����� ����������� ����� �� �����
	DrawBazieSurface(massiv); //���, ����� ����������� �����


	


   //��������� ������ ������

	
	glMatrixMode(GL_PROJECTION);	//������ �������� ������� ��������. 
	                                //(���� ��������� ��������, ����� �� ������������.)
	glPushMatrix();   //��������� ������� ������� ������������� (������� ��������� ������������� ��������) � ���� 				    
	glLoadIdentity();	  //��������� ��������� �������
	glOrtho(0, ogl->getWidth(), 0, ogl->getHeight(), 0, 1);	 //������� ����� ������������� ��������

	glMatrixMode(GL_MODELVIEW);		//������������� �� �����-��� �������
	glPushMatrix();			  //��������� ������� ������� � ���� (��������� ������, ����������)
	glLoadIdentity();		  //���������� �� � ������

	glDisable(GL_LIGHTING);



	GuiTextRectangle rec;		   //������� ����� ��������� ��� ������� ������ � �������� ������.
	rec.setSize(300, 150);
	rec.setPosition(10, ogl->getHeight() - 150 - 10);


	std::stringstream ss;
	ss << "T - ���/���� �������" << std::endl;
	ss << "L - ���/���� ���������" << std::endl;
	ss << "F - ���� �� ������" << std::endl;
	ss << "G - ������� ���� �� �����������" << std::endl;
	ss << "G+��� ������� ���� �� ���������" << std::endl;
	ss << "�����. �����: (" << light.pos.X() << ", " << light.pos.Y() << ", " << light.pos.Z() << ")" << std::endl;
	ss << "�����. ������: (" << camera.pos.X() << ", " << camera.pos.Y() << ", " << camera.pos.Z() << ")" << std::endl;
	ss << "��������� ������: R="  << camera.camDist << ", fi1=" << camera.fi1 << ", fi2=" << camera.fi2 << std::endl;
	
	rec.setText(ss.str().c_str());
	rec.Draw();

	glMatrixMode(GL_PROJECTION);	  //��������������� ������� �������� � �����-��� �������� �� �����.
	glPopMatrix();


	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();

}