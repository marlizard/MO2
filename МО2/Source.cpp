#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>

using namespace std;

typedef double real;

real lambda_f(real x) {
	return (x - 4) * (x + 4);
}


real fibonacci(int n) {
	return (pow(1 + sqrt(5), n) / pow(2, n) - pow(1 - sqrt(5), n) / pow(2, n)) / sqrt(5);
}

int main000() {
	real a = -10, b = 10;

	real len = abs(b - a);
	int n = 0;
	real e = 10e-15;

	while (len / e >= fibonacci(n + 2)) n++;

	real x1 = a + fibonacci(n) / fibonacci(n + 2) * len;
	real x2 = a + fibonacci(n + 1) / fibonacci(n + 2) * len;
	real f1 = lambda_f(x1);
	real f2 = lambda_f(x2);

	for (int k = 1; k <= n; k++) {
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

	f1 <= lambda_f(x1 + e) ? b = x1 : a = x1;

	return (a + b) / 2;
}