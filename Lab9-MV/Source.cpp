#include <iostream>
#include "windows.h"
#include <vector>
#include "string"
#include <cmath>

using namespace std;

class Matrix {
private:
	vector<vector<double>> matrix;
	vector<int> width_form;
	int m = 0, n = 0;

	void Format() {
		width_form.clear();
		int width, buf_width, width_null;
		string str_width;
		width = 1;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				str_width = to_string(matrix[j][i]);
				for (int i = str_width.size() - 1; i >= 0; i--) {
					if (str_width[i] == '0' || str_width[i] == ',') {
						str_width.erase(i);
					}
					else break;
				}
				buf_width = str_width.size();
				if (width < buf_width) width = buf_width;
			}
			width_form.push_back(width);
		}
	}
	void Size() {
		m = matrix.size();
		n = matrix[0].size();
		for (int i = 1; i < m; i++) {
			if (matrix[i].size() != n) {
				n = 0;
				break;
			}
		}
	}
public:
	void Input_Matrix(int m, int n)
	{
		double number;
		vector<double> str_matrix;
		cout << "Enter the [" << m << "," << n << "] matrix:" << endl;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				cin >> number;
				str_matrix.push_back(number);
			}
			matrix.push_back(str_matrix);
			str_matrix.clear();
		}
	}
	void Output_Matrix() {
		Size();
		if (n == 0 || m == 0) {
			return;
		}
		Format();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				cout.width(width_form[j]);
				cout << matrix[i][j] << " ";
			}
			cout << endl;
		}
		cout << endl;
	}
	/*void Output_SLAU() {
		Size();
		if (n == 0 || m == 0) {
			return;
		}
		Format();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				cout.width(width_form[j]);
				cout << matrix[i][j] << " ";
				if (j == n - 2) cout << "| ";
			}
			cout << endl;
		}
		cout << endl;
	}*/
	vector<vector<double>> get_Matrix() {
		return matrix;
	}
	void set_Matrix(vector<vector<double>> matrix_) {
		matrix = matrix_;
	}
	void Matrix_dot_B(vector<vector<double>> matrix_b) {
		m = matrix.size();
		n = matrix_b[0].size();
		int n_ = matrix_b.size();
		if (matrix[0].size() != matrix_b.size()) return;
		vector<vector<double>> dot_matrix;
		for (int i = 0; i < m; i++) {
			vector<double> buf_vector;
			for (int j = 0; j < n; j++) {
				double buf_element = 0;
				for (int k = 0; k < n_; k++) {
					buf_element += matrix[i][k] * matrix_b[k][j];
				}
				buf_vector.push_back(buf_element);
			}
			dot_matrix.push_back(buf_vector);
		}
		matrix = dot_matrix;
	}
};

class Determinant {
private:
	vector<vector<double>> matrix_det;
	double det = 0;
	Matrix M;
	double Det_2x(vector<vector<double>> matrix_2x) {
		double right = matrix_2x[0][0] * matrix_2x[1][1];
		double left = matrix_2x[0][1] * matrix_2x[1][0];
		return right - left;
	}
	double Det_3x(vector<vector<double>> matrix_3x) {
		double right = 0;
		double buf = 1, buf1 = 1;
		for (int i = 0; i < 3; i++) {
			buf *= matrix_3x[i][i];
		}
		right += buf;
		buf = 1;
		for (int i = 0; i < 3; i++) {
			if (i + 1 < 3) {
				buf *= matrix_3x[i][i + 1];
				buf1 *= matrix_3x[i + 1][i];
			}
			else {
				buf *= matrix_3x[i][0];
				buf1 *= matrix_3x[0][i];
			}
		}
		right += buf;
		right += buf1;
		double left = 0;
		buf = matrix_3x[0][2] * matrix_3x[2][0] * matrix_3x[1][1];
		left += buf;
		buf = matrix_3x[0][1] * matrix_3x[1][0] * matrix_3x[2][2];
		left += buf;
		buf = matrix_3x[1][2] * matrix_3x[2][1] * matrix_3x[0][0];
		left += buf;
		return right - left;
	}
	double Det_mx(vector<vector<double>> matrix_mx) {
		if (matrix_mx.size() == 2) return Det_2x(matrix_mx);
		else if (matrix_mx.size() == 3) return Det_3x(matrix_mx);
		else {
			double det_mx = 0;
			for (int i = 0; i < matrix_mx.size(); i++) {
				vector<vector<double>> buf_m = matrix_mx;
				buf_m.erase(buf_m.begin());
				int buf_size = buf_m.size();
				for (int buf_index = 0; buf_index < buf_size; buf_index++) {
					buf_m[buf_index].erase(buf_m[buf_index].begin() + i);
				}
				if (i % 2 == 0) {
					det_mx += matrix_mx[0][i] * Det_mx(buf_m);
				}
				else det_mx -= matrix_mx[0][i] * Det_mx(buf_m);
			}
			return det_mx;
		}
	}
public:
	void Input_Det_Matrix(vector<vector<double>> matrix) {
		matrix_det = matrix;
	}
	void Decision() {
		det = Det_mx(matrix_det);
		cout << "det(A) = " << det << endl;
	}
	double get_det() {
		return det;
	}
};

class Danilevskiy {
private:
	vector<vector<double>> matrix_A, matrix_B, matrix_E;
	Matrix A;
	Matrix B;
	Matrix Buf;
	int m;
	void Accur_str(int str_index) {
		double buf;
		vector<double> str_slau = matrix_A[str_index];
		for (int i = 0; i < str_slau.size(); i++) {
			buf = abs(str_slau[i]);
			if (buf != 0) {
				if (fmod(buf, 0.00001) < 1) {
					str_slau[i] = round(str_slau[i] * 100000) / 100000;
				}
			}
		}
		matrix_A[str_index] = str_slau;
	}
	void Output_A() {
		A.set_Matrix(matrix_A);
		A.Output_Matrix();
	}
	void Output_Buf(vector<vector<double>> matrix_b) {
		Buf.set_Matrix(matrix_b);
		Buf.Output_Matrix();
	}
	vector<vector<double>> Create_B_index(int index_i) {
		vector<vector<double>> B_index = matrix_E;
		int index = m - 2 - index_i;
		for (int i = 0; i < m; i++) {
			if(i == index) B_index[index][i] = (1 / matrix_A[index + 1][index]);
			else B_index[index][i] = -(matrix_A[index + 1][i] / matrix_A[index + 1][index]);
		
		}
		cout << "B" << index_i+1 << " =" << endl;
		Output_Buf(B_index);
		return B_index;
	}
	vector<vector<double>> Create_B_1_index(int index_i) {
		vector<vector<double>> B_1_index = matrix_E;
		int index = m - 2 - index_i;
		for (int i = 0; i < m; i++) {
			B_1_index[index][i] = matrix_A[index + 1][i];
		}
		cout << "B" << index_i+1 << "^(-1) =" << endl;
		Output_Buf(B_1_index);
		return B_1_index;
	}
	void Create_Frobenius() {
		vector<vector<double>> buf_B;
		vector<vector<double>> buf_B_;
		buf_B = Create_B_index(0);
		buf_B_ = Create_B_1_index(0);
		B.set_Matrix(buf_B);
		A.set_Matrix(buf_B_);
		A.Matrix_dot_B(matrix_A);
		A.Matrix_dot_B(buf_B);
		matrix_A = A.get_Matrix();
		cout << "D" << 1 << " =" << endl;
		Output_A();
		for (int i = 1; i < m - 1; i++) {
			buf_B = Create_B_index(i);
			buf_B_ = Create_B_1_index(i);
			B.Matrix_dot_B(buf_B);
			A.set_Matrix(buf_B_);
			A.Matrix_dot_B(matrix_A);
			A.Matrix_dot_B(buf_B);
			matrix_A = A.get_Matrix();
			cout << "D" << i + 1 << " =" << endl;
			if(i != m - 2) Output_A();
		}
		for (int i = 1; i < m - 1; i++) {
			Accur_str(i);
		}
		Output_A();
		cout << "B =" << endl;
		B.Output_Matrix();
	}
public:
	bool Input(vector<vector<double>> matrix_a) {
		if (matrix_a.size() != matrix_a[0].size()) {
			return false;
		}
		matrix_A = matrix_a;
		m = matrix_A.size();
		Output_A();
		for (int i = 0; i < m; i++) {
			vector<double> buf_vector;
			for (int j = 0; j < m; j++) {
				if(i != j) buf_vector.push_back(0);
				else buf_vector.push_back(1);
			}
			matrix_E.push_back(buf_vector);
		}
		return true;
	}
	void Decision() {
		Create_Frobenius();
	}
	void Rezult() {
		cout << "Polynomial: ";
		cout << "k^" << m << " + (";
		for (int i = 0; i < m - 1; i++) {
			cout << -matrix_A[0][i] << ")*" << " k^" << m - (i + 1) << " + (";
		}
		cout << -matrix_A[0][m - 1] << ") = 0";
		cout << endl;
	}

};

class Task {
private:
	Matrix M;
	Matrix M1;
	Danilevskiy D;
	Determinant Det;
	int m;
	void Create() {
		M.Input_Matrix(m, m);
		Det.Input_Det_Matrix(M.get_Matrix());
	}
public:
	void Task_1() {
		cout << "Input the size of the matrix." << endl << "m = ";
		cin >> m;
		if (m < 3 || m > 49) {
			cout << "Incorrect size of the matrix." << endl;
			return;
		}
		Create();
		Det.Decision();
		if (Det.get_det() == 0) {
			cout << "Error! The matrix det = 0." << endl;
			return;
		}
		D.Input(M.get_Matrix());
		D.Decision();
		D.Rezult();
	}
};

int main() {
	setlocale(LC_ALL, "Rus");
	SetConsoleCP(1251);
	int ans, exit = 1;
	while (exit == 1) {
		Task T;
		cout << "1.Danilevskiy" << endl << "2.Exit" << endl << "Choose a way:" << endl;
		cin >> ans;
		switch (ans)
		{
		case 1:
			T.Task_1();
			break;
		case 2:
			exit = 0;
			break;
		default:
			cout << "This task does not exist" << endl;
			break;
		}
	}
	system("pause");
	return 0;
}
/*
Example:
2.2 1 0.5 2
1 1.3 2 1
0.5 2 0.5 1.6
2 1 1.6 2

V9
2 1.5 3.5 4.5
1.5 2 2 1.6
3.5 2 2 1.7
4.5 1.6 1.7 2


*/