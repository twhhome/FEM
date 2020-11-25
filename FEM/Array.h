#pragma once

#include "includes.h"

template <typename T> class Array {
public:
	Array();
	Array(int n);
	Array(Array<T> const &src);
	~Array();

	T& operator[](int i);
	T& operator[](int i) const;
	Array<T>& operator=(Array<T> const &src);
	int size() const;

private:
	T *data;
	int n;
};

template <typename T>
Array<T>::Array() {
	n = 0;
	data = NULL;
}

template <typename T>
Array<T>::Array(int n) {
	this->n = n;
	data = new T[this->n]();
}

template <typename T>
Array<T>::Array(Array<T> const &src) {
	delete[] data;
	n = src.n;
	data = new T[n];
	memcpy(data, src.data, n * sizeof(T));
}

template <typename T>
Array<T>::~Array() {
	delete[] data;
	data = NULL;
	n = 0;
}

template <typename T>
T& Array<T>::operator[](int i) {
	return data[i];
}

template <typename T>
T& Array<T>::operator[](int i) const {
	return data[i];
}

template <typename T>
Array<T>& Array<T>::operator=(Array<T> const &src) {
	delete[] data;
	n = src.n;
	data = new T[n];
	memcpy(data, src.data, n * sizeof(T));
	return *this;
}

template <typename T>
int Array<T>::size() const {
	return n;
}