#pragma once

#include "includes.h"
#include "Array.h"

template <typename T> class Matrix
{
public:
	T** array;
	int n, m;

	Matrix();
	Matrix(int n, int m);
	Matrix(int n) : Matrix(n, n) {}
	Matrix(int n, int m, T** src);
	Matrix(Matrix<T> const &src);
	Matrix(Array<T> const &src);

	~Matrix();

	T* operator[](int i);
	T& operator()(int i);
	T& operator()(int i, int j);
	Matrix<T>& operator=(Matrix<T> const &src);
	Matrix<T> operator-(Matrix<T> a);
	Matrix<T> operator*(Matrix<T> a);
	Matrix<T> operator/(double a);

	Matrix<T> trans();

	template <typename T> friend std::ostream &operator<<(std::ostream &os, Matrix<T> &mat);

private:
	Matrix<T> subtract(Matrix<T> &a);
	Matrix<T> multiply(Matrix<T> &a);
	Matrix<T> divide(double a);
};



template <typename T>
Matrix<T>::Matrix() {
	this->n = 0;
	this->m = 0;
	this->array = NULL;
}

template <typename T>
Matrix<T>::Matrix(int n, int m) {
	this->n = n;
	this->m = m;
	this->array = new T*[n];
	for (int i = 0; i < n; i++) {
		this->array[i] = new T[m]();
	}
}

template <typename T>
Matrix<T>::Matrix(int n, int m, T** src) : Matrix(n, m) {
	for (int i = 0; i < n; i++) {
		memcpy(this->array[i], src[i], m * sizeof(T));
	}
}

template <typename T>
Matrix<T>::Matrix(Matrix<T> const &src) {
	this->n = src.n;
	this->m = src.m;
	this->array = new T*[this->n];
	for (int i = 0; i < this->n; i++) {
		this->array[i] = new T[m]();
		memcpy(this->array[i], src.array[i], this->m * sizeof(T));
	}
}

template <typename T>
Matrix<T>::Matrix(Array<T> const &src) {
	n = 1;
	m = src.size();
	array = new T*[1];
	array[0] = new T[m];
	for (int i = 0; i < m; i++) {
		array[0][i] = src[i];
	}
}

template <typename T>
Matrix<T>::~Matrix() {
	for (int i = 0; i < this->n; i++) {
		delete[] this->array[i];
	}
	delete[] this->array;
	this->array = NULL;
	this->n = 0;
	this->m = 0;
}

template <typename T>
T* Matrix<T>::operator[](int i) {
	return array[i];
}

template <typename T>
T& Matrix<T>::operator()(int i) {
	assert(this->n == 1 || this->m == 1);
	assert(i >= 0 && ((this->n == 1) ? i < this->m : i < this->n));
	return (this->n == 1) ? this->array[0][i] : this->array[i][0];
}

template <typename T>
T& Matrix<T>::operator()(int i, int j) {
	assert((i >= 0 && i < this->n) && (j >= 0 && j < this->m));
	return this->array[i][j];
}

template <typename T>
Matrix<T>& Matrix<T>::operator=(Matrix<T> const &src) {
	for (int i = 0; i < this->n; i++) {
		delete[] this->array[i];
	}
	delete[] this->array;
	this->array = NULL;

	this->n = src.n;
	this->m = src.m;
	this->array = new T*[this->n];
	for (int i = 0; i < this->n; i++) {
		this->array[i] = new T[m]();
		memcpy(this->array[i], src.array[i], this->m * sizeof(T));
	}

	return *this;
}

template <typename T>
Matrix<T> Matrix<T>::trans() {
	Matrix<T> rslt(this->m, this->n);
	for (int i = 0; i < rslt.n; i++) {
		for (int j = 0; j < rslt.m; j++) {
			rslt.array[i][j] = this->array[j][i];
		}
	}

	return rslt;
}

template <typename T>
Matrix<T> Matrix<T>::operator-(Matrix<T> a) {
	return subtract(a);
}

template <typename T>
Matrix<T> Matrix<T>::operator*(Matrix<T> a) {
	return multiply(a);
}

template <typename T>
Matrix<T> Matrix<T>::operator/(double a) {
	return divide(a);
}

template <typename T>
std::ostream &operator<<(std::ostream &os, Matrix<T> &mat) {
	for (int i = 0; i < mat.n; i++) {
		for (int j = 0; j < mat.m; j++) {
			//os << std::setw(10) << std::left << mat.array[i][j] << " ";
			os << mat.array[i][j] << " ";
		}
		if (i != mat.n - 1)
			os << std::endl;
	}
	return os;
}

template <typename T>
Matrix<T> Matrix<T>::subtract(Matrix<T> &a) {
	assert((this->n == a.n) && (this->m == a.m));
	Matrix<T> rslt(this->n, this->m);
	for (int i = 0; i < this->n; i++) {
		for (int j = 0; j < this->m; j++) {
			rslt[i][j] = this->array[i][j] - a.array[i][j];
		}
	}
	return rslt;
}

template <typename T>
Matrix<T> Matrix<T>::multiply(Matrix<T> &a) {
	assert(this->m == a.n);
	Matrix<T> rslt(this->n, a.m);
	for (int i = 0; i < this->n; i++) {
		for (int j = 0; j < a.m; j++) {
			rslt[i][j] = 0;
			for (int k = 0; k < this->m; k++) {
				rslt[i][j] += this->array[i][k] * a.array[k][j];
			}
		}
	}
	return rslt;
}

template <typename T>
Matrix<T> Matrix<T>::divide(double a) {
	Matrix<T> rslt(this->n, this->m);
	for (int i = 0; i < this->n; i++) {
		for (int j = 0; j < this->m; j++) {
			rslt[i][j] = this->array[i][j] / a;
		}
	}
	return rslt;
}