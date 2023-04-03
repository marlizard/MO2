#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>

using namespace std;

typedef double real;

//Тип экстремума(1 - минимум, 2 - максимум)
#define extr_type 2

//Номер целевой функции (Основная - 1, Квадратичная - 2, Розенброка - 3)
#define function 1

//Одномерный метод поиска лямбды (Фибоначчи - 1, Парабол - 2)
#define lambda_max 2

//Номер метода (Ньютон - 1, Бройден - 2)
#define method 2

class Method {
public:
	//Вектор х(k)
	vector<real> x;
	//Вектор х(k+1)
	vector<real> xn;
	//Напраление поиска
	vector<real> s;
	//Точность для приближений
	vector<real> ex = { 10e-7, 10e-7 };
	//Точность для функций
	real ef = 10e-7;
	//Матрица Гессена
	vector<vector<real>> Hessen;
	//Вектор градиента обновленный
	vector<real> Gradient;
	//Вектор градиента устаревший
	vector<real> grad;
	//Количество обращений к функции
	int used_f;

	//Целевая функция
	real f(real x, real y);

	//Первая производная целевой функции по х
	real f_x(real x, real y);

	//Первая производная целевой функции по у
	real f_y(real x, real y);

	//Вторая производная целевой функции по х и х
	real f_x2(real x, real y);

	//Вторая производная целевой функции по х и у
	real f_xy(real x, real y);

	//Вторая производная целевой функции по у и у
	real f_y2(real x, real y);

	//Ввод данных
	void input();

	//Пересчет матрицы Гессена
	void Hessen_matrix();

	//Пересчет вектора градиента
	void Gradient_vector();

	//Угол между векторами х и s
	real angle();

	//Пересчет направления поиска
	void search_direction();

	//Функция от лямбды
	real lambda_f(real lambda_number);

	//Нахождение промежутка
	void interval(real& a, real& b);

	//Нахождение лямбда методом Фибоначчиа
	real Fibonacci_method();

	//Число Фибоначчи под номером n
	real fibonacci(int n);

	//Нахождение лямбда методом парабол
	real Parabolic_method();

	//Реализация метода Ньютона
	void Newton_method();

	//Пересчет добавки к аппроксимированной матрице Гессена
	void Recount_Hessen_approx(real lambda_num);

	//Подсчет произведения векторов
	real vector_mult(vector<real> c, vector<real> p);
	
	//Нахождение нормы вектора
	real norm(vector<real> vect);

	//Реализация метода Бройдена
	void Broyden_method();
};


real Method::f(real x, real y) {
	switch (function) {
	case 1:
		return 1 / (1 + pow(((x - 2) / 3), 2) + pow(((y - 2) / 3), 2))
			+ 3 / (1 + pow((x - 1), 2) + pow(((y - 1) / 2), 2));
		break;
	case 2:
		return 100 * pow((y - x), 2) + pow((1 - x), 2);
		break;
	case 3:
		return 100 * pow((y - pow(x, 2)), 2) + pow((1 - x), 2);
		break;
	}
}


real Method::f_x(real x, real y) {
	switch (function) {
	case 1:
		return -(2 * x - 4) / 9 / pow((1 + pow(((x - 2) / 3), 2) + pow(((y - 2) / 3), 2)), 2)
		-(6 * x - 6) / pow((1 + pow((x - 1), 2) + pow(((y - 1) / 2), 2)), 2);
		break;
	case 2:
		return 202 * x - 200 * y - 2;
		break;
	case 3:
		return 400 * pow(x, 3) + (2 - 400 * y) * x - 2;
		break;
	}
}


real Method::f_y(real x, real y) {
	switch (function) {
	case 1:
		return -(3 * y - 3) / 2 / pow((1 + pow((x - 1), 2) + pow(((y - 1) / 2), 2)), 2)
		-(2 * y - 4) / 9 / pow((1 + pow(((x - 2) / 3), 2) + pow(((y - 2) / 3), 2)), 2);
		break;
	case 2:
		return -200 * x + 200 * y;
		break;
	case 3:
		return 200 * y - 200 * pow(x, 2);
		break;
	}
}


real Method::f_x2(real x, real y) {
	switch (function) {
	case 1:
		return -6 * (pow(1 + pow((x - 1), 2) + pow(((y - 1) / 2), 2), 2) - 4 * (1 + pow((x - 1), 2) + pow(((y - 1) / 2), 2)) * pow((x - 1), 2)) / pow(1 + pow((x - 1), 2) + pow(((y - 1) / 2), 2), 4)
		-2 * (pow(1 + pow(((x - 2) / 3), 2) + pow(((y - 2) / 3), 2), 2) - 4 / 9 * (1 + pow(((x - 2) / 3), 2) + pow(((y - 2) / 3), 2)) * pow(x - 2, 2)) / 9 / pow(1 + pow(((x - 2) / 3), 2) + pow(((y - 2) / 3), 2), 4);
		break;
	case 2:
		return 202;
		break;
	case 3:
		return 1200 * pow(x, 2) - 400 * y + 2;
		break;
	}
}


real Method::f_xy(real x, real y) {
	switch (function) {
	case 1:
		return 6 * (y - 1) * (x - 1) / pow((1 + pow((x - 1), 2) + pow(((y - 1) / 2), 2)), 3)
		+ 8 * (y + 2) * (x - 2) / 81 / pow((1 + pow(((x - 2) / 3), 2) + pow(((y - 2) / 3), 2)), 3);
		break;
	case 2:
		return -200;
		break;
	case 3:
		return -400 * x;
		break;
	}
}


real Method::f_y2(real x, real y) {
	switch (function) {
	case 1:
		return -3 * (pow((1 + pow((x - 1), 2) + pow(((y - 1) / 2), 2)), 2) - (1 + pow((x - 1), 2) + pow(((y - 1) / 2), 2)) * pow((y - 1), 2)) / 2 / pow((1 + pow((x - 1), 2) + pow(((y - 1) / 2), 2)), 4)
		-2 * (pow((1 + pow(((x - 2) / 3), 2) + pow(((y - 2) / 3), 2)), 2) - 4 / 9 * (1 + pow(((x - 2) / 3), 2) + pow(((y - 2) / 3), 2)) * pow((y - 2), 2)) / 9 / pow((1 + pow(((x - 2) / 3), 2) + pow(((y - 2) / 3), 2)), 4);
		break;
	case 2:
		return 200;
		break;
	case 3:
		return 200;
		break;
	}
}


void Method::input() {
	used_f = 0;

	x.resize(2);
	xn.resize(2);
	s.resize(2);
	Gradient.resize(2);
	grad.resize(2);
	Hessen.resize(2);
	for (int i = 0; i < Hessen.size(); i++) Hessen[i].resize(2);

	ifstream init("initial_approximation.txt");
	for (int i = 0; i < x.size(); i++) {
		init >> xn[i];
		x[i] = 0;
	}
	init.close();
}


void Method::Hessen_matrix() {
	if (method == 1) {
		Hessen[0][0] = f_x2(x[0], x[1]);
		Hessen[0][1] = f_xy(x[0], x[1]);
		Hessen[1][0] = f_xy(x[0], x[1]);
		Hessen[1][1] = f_y2(x[0], x[1]);

		used_f += 4;
	}
	else {
		Hessen[0][0] = 1;
		Hessen[0][1] = 0;
		Hessen[1][0] = 0;
		Hessen[1][1] = 1;
	}
}


void Method::Gradient_vector() {
	Gradient[0] = f_x(xn[0], xn[1]);
	Gradient[1] = f_y(xn[0], xn[1]);

	used_f += 2;
}


real Method::angle() {
	return acos((x[0] * s[0] + x[1] * s[1]) / sqrt((pow(x[0], 2) + pow(x[1], 2)) * (pow(x[0], 2) + pow(x[1], 2))));
}


void Method::search_direction() {
	Hessen_matrix();

	if (method == 1) {
		Gradient_vector();

		s[0] = (Hessen[1][1] * Gradient[0] - Hessen[0][1] * Gradient[1]) / (Hessen[0][0] * Hessen[1][1] - Hessen[1][0] * Hessen[0][1]);
		s[1] = (Hessen[0][0] * Gradient[1] - Hessen[1][0] * Gradient[0]) / (Hessen[0][0] * Hessen[1][1] - Hessen[1][0] * Hessen[0][1]);
	}
	else {
		grad[0] = f_x(x[0], x[1]);
		grad[1] = f_y(x[0], x[1]);
		used_f += 2;

		s[0] = -(Hessen[1][1] * grad[0] - Hessen[0][1] * grad[1]) / (Hessen[0][0] * Hessen[1][1] - Hessen[0][1] * Hessen[1][0]);
		s[1] = -(-Hessen[1][0] * grad[0] + Hessen[0][0] * grad[1]) / (Hessen[0][0] * Hessen[1][1] - Hessen[0][1] * Hessen[1][0]);
	}
}


real Method::lambda_f(real lambda_number) {
	vector<real> z(2);
	z[0] = x[0] - s[0] * lambda_number;
	z[1] = x[1] - s[1] * lambda_number;
	return f(z[0], z[1]);
}


void Method::interval(real& a, real& b) {
	real lambda = 0;
	real h = 10e-5;
	real prev = 0, now, next;

	if (extr_type == 1) {
		if (lambda_f(lambda) < lambda_f(lambda + h)) h *= -1;
	}
	else if (lambda_f(lambda) > lambda_f(lambda + h)) h *= -1;
	used_f += 2;

	now = lambda;
	next = lambda + h;

	bool mark = true;

	while (mark) {
		if (extr_type == 1) {
			if (lambda_f(now) > lambda_f(next)) {
				prev = now;
				now = next;
				h *= 2;
				next = now + h;
				used_f++;
			}
			else {
				mark = false;
				if (h < 0) {
					a = next;
					b = prev;
				}
				else {
					a = prev;
					b = next;
				}
			}
		}
		else {
			if (lambda_f(now) < lambda_f(next)) {
				prev = now;
				now = next;
				h *= 2;
				next = now + h;
				used_f++;
			}
			else {
				mark = false;
				if (h < 0) {
					a = next;
					b = prev;
				}
				else {
					a = prev;
					b = next;
				}
			}
		}
		used_f++;
	}
}


real Method::fibonacci(int n) {
	return (pow(1 + sqrt(5), n) / pow(2, n) - pow(1 - sqrt(5), n) / pow(2, n)) / sqrt(5);
}


real Method::Fibonacci_method() {
	real a = 0, b = 0;
	interval(a, b);

	real len = abs(b - a);
	int n = 0;
	real e = 10e-15;

	while (len / e >= fibonacci(n + 2)) n++;

	real x1 = a + fibonacci(n) / fibonacci(n + 2) * len;
	real x2 = a + fibonacci(n + 1) / fibonacci(n + 2) * len;
	real f1 = lambda_f(x1);
	real f2 = lambda_f(x2);

	for (int k = 1; k <= n; k++) {
		if (extr_type == 1) {
			if (f1 > f2) {
				a = x1;
				x1 = x2;
				f1 = f2;
				x2 = a + fibonacci(n - k + 2) / fibonacci(n + 2) * len;
				f2 = lambda_f(x2);
			}
			else {
				b = x2;
				x2 = x1;
				f2 = f1;
				x1 = a + fibonacci(n - k + 1) / fibonacci(n + 2) * len;
				f1 = lambda_f(x1);
			}
		}
		else {
			if (f1 < f2) {
				a = x1;
				x1 = x2;
				f1 = f2;
				x2 = a + fibonacci(n - k + 2) / fibonacci(n + 2) * len;
				f2 = lambda_f(x2);
			}
			else {
				b = x2;
				x2 = x1;
				f2 = f1;
				x1 = a + fibonacci(n - k + 1) / fibonacci(n + 2) * len;
				f1 = lambda_f(x1);
			}
		}
		used_f++;
	}

	if (extr_type == 1) {
		f1 <= lambda_f(x1 + e) ? b = x1 : a = x1;
	}
	else {
		f1 >= lambda_f(x1 + e) ? b = x1 : a = x1;
	}
	used_f++;

	return (a + b) / 2;
}

real Method::Parabolic_method() {
	real a = 0, b = 0;
	interval(a, b);

	real mid = (a + b) / 2;
	vector<vector<real>> A(3);
	A[0].resize(4);
	A[1].resize(4);
	A[2].resize(4);
	vector<real> res(3, 0);

	A[0][0] = pow(a, 2);
	A[0][1] = a;
	A[0][2] = 1;
	A[0][3] = lambda_f(a);
	used_f++;

	A[1][0] = pow(mid, 2);
	A[1][1] = mid;
	A[1][2] = 1;
	A[1][3] = lambda_f(mid);
	used_f++;

	A[2][0] = pow(b, 2);
	A[2][1] = b;
	A[2][2] = 1;
	A[2][3] = lambda_f(b);
	used_f++;

	real tmp;

	for (int i = 0; i < res.size(); i++)
	{
		tmp = A[i][i];
		for (int j = res.size(); j >= i; j--)
			A[i][j] /= tmp;
		for (int j = i + 1; j < res.size(); j++)
		{
			tmp = A[j][i];
			for (int k = res.size(); k >= i; k--)
				A[j][k] -= tmp * A[i][k];
		}
	}

	res[2] = A[2][3];
	for (int i = res.size() - 2; i >= 0; i--)
	{
		res[i] = A[i][3];
		for (int j = i + 1; j < res.size(); j++) res[i] -= A[i][j] * res[j];
	}

	return -res[1] / 2 / res[0];
}


void Method::Newton_method() {
	int iteration = 0;
	real lambda_num;
	ofstream out("Newton.txt");

	used_f++;
	while (abs(f(x[0], x[1]) - f(xn[0], xn[1])) >= ef && abs(x[0] - xn[0]) >= ex[0] && abs(x[1] - xn[1]) >= ex[1]) {
		used_f++;

		x[0] = xn[0];
		x[1] = xn[1];
		search_direction();

		switch (lambda_max) {
		case 1:
			lambda_num = Fibonacci_method();
			break;
		case 2:
			lambda_num = Parabolic_method();
			break;
		}

		xn[0] = x[0] - s[0] * lambda_num;
		xn[1] = x[1] - s[1] * lambda_num;

		iteration++;

		out << iteration << " | ";
		out << fixed;
		out << "(" << setprecision(14) << xn[0] << ", " << xn[1] << ") | ";
		out << f(xn[0], xn[1]) << " | ";
		out << "(" << setprecision(14) << s[0] << ", " << s[1] << ") | ";
		out << setprecision(14) << lambda_num << " | ";
		out << "(" << setprecision(7) << scientific << abs(x[0] - xn[0]) << ", " << abs(x[1] - xn[1]) << ", "
			<< abs(f(x[0], x[1]) - f(xn[0], xn[1])) << ") | ";
		out << setprecision(14) << fixed << angle() << " | ";
		out << "(" << setprecision(14) << Gradient[0] << ", " << Gradient[1] << ") | ";
		out << "(" << setprecision(14) << Hessen[0][0] << ", " << Hessen[0][1] << ", "
			<< Hessen[1][0] << ", " << Hessen[1][1] << ")" << endl;
	}
	out << "Number of iterations = " << iteration << endl;
	out << "Number of function call = " << used_f << endl;
	out.close();
}


real Method::vector_mult(vector<real> c, vector<real> p) {
	real sum = 0;
	for (int i = 0; i < c.size(); i++) sum += c[i] * p[i];
	return sum;
}


real Method::norm(vector<real> vect) {
	real sum = 0;
	for (int i = 0; i < vect.size(); i++) sum += pow(vect[i], 2);
	return sqrt(sum);
}


void Method::Recount_Hessen_approx(real lambda_num) {
	vector<real> dx(2, 0), g(2, 0), tmp(2, 0);

	dx[0] = xn[0] - x[0];
	dx[1] = xn[1] - x[1];

	Gradient_vector();

	g[0] = Gradient[0] - grad[0];
	g[1] = Gradient[1] - grad[1];

	s[0] = -(Hessen[1][1] * g[0] - Hessen[0][1] * g[1]) / (Hessen[0][0] * Hessen[1][1] - Hessen[0][1] * Hessen[1][0]);
	s[1] = -(-Hessen[1][0] * g[0] + Hessen[0][0] * g[1]) / (Hessen[0][0] * Hessen[1][1] - Hessen[0][1] * Hessen[1][0]);

	tmp[0] = dx[0] + lambda_num * s[0];
	tmp[1] = dx[1] + lambda_num * s[1];

	if (norm(tmp) > ef)
	{
		double det = tmp[0] * g[0] + tmp[1] * g[1];

		Hessen[0][0] += tmp[0] * tmp[0] / det;
		Hessen[0][1] += tmp[0] * tmp[1] / det;
		Hessen[1][0] += tmp[0] * tmp[1] / det;
		Hessen[1][1] += tmp[1] * tmp[1] / det;
	}
}


void Method::Broyden_method() {
	int iteration = 0;
	real lambda_num;
	ofstream out("Broyden.txt");

	Gradient_vector();

	used_f++;
	while (abs(f(x[0], x[1]) - f(xn[0], xn[1])) >= ef && abs(x[0] - xn[0]) >= ex[0] && abs(x[1] - xn[1]) >= ex[1] && norm(Gradient) >= ef) {
		used_f++;
		x[0] = xn[0];
		x[1] = xn[1];

		search_direction();

		switch (lambda_max) {
		case 1:
			lambda_num = Fibonacci_method();
			break;
		case 2:
			lambda_num = Parabolic_method();
			break;
		}

		xn[0] = x[0] - s[0] * lambda_num;
		xn[1] = x[1] - s[1] * lambda_num;
		iteration++;

		Recount_Hessen_approx(lambda_num);

		out << iteration << " | ";
		out << fixed;
		out << "(" << setprecision(14) << xn[0] << ", " << xn[1] << ") | ";
		out << f(xn[0], xn[1]) << " | ";
		out << "(" << setprecision(14) << s[0] << ", " << s[1] << ") | ";
		out << setprecision(14) << lambda_num << " | ";
		out << "(" << setprecision(7) << scientific << abs(x[0] - xn[0]) << ", " << abs(x[1] - xn[1]) << ", "
			<< abs(f(x[0], x[1]) - f(xn[0], xn[1])) << ") | ";
		out << setprecision(14) << fixed << angle() << " | ";
		out << "(" << setprecision(14) << Gradient[0] << ", " << Gradient[1] << ") | ";
		out << "(" << setprecision(14) << Hessen[0][0] << ", " << Hessen[0][1] << ", "
			<< Hessen[1][0] << ", " << Hessen[1][1] << ")" << endl;
	}
	out << "Number of iterations = " << iteration << endl;
	out << "Number of function call = " << used_f << endl;
	out.close();
}


int main() {
	Method m;
	m.input();
	switch (method) {
	case 1:
		m.Newton_method();
		break;
	case 2:
		m.Broyden_method();
		break;
	default:
		cout << "Choose correct number of method";
	}
	return 0;
}
