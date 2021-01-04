#ifndef HEADER_gslWrapper_
#define HEADER_gslWrapper_

#include <vector>
#include <iostream>
#include <gsl/gsl_vector.h>


class gslVector {
public:
    
    gslVector() 
	: base(nullptr) {}
    
    gslVector(int len) {
	base = gsl_vector_alloc(len);
    }

    ~gslVector() {
	if (base != nullptr)
	    gsl_vector_free(base);
    }

    // Copy constructor
    gslVector(const gslVector& obj)
	: base(nullptr) {
	// Default GSL error handler aborts program if out of memory
	if (obj.gsl() != nullptr) {
	    base = gsl_vector_alloc(obj.size());
	    gsl_vector_memcpy(base, obj.gsl());
	}
    }

    // Copy from base gsl_vector
    gslVector(const gsl_vector* obj)
	: base(nullptr) {
	if (obj != nullptr) {
	    base = gsl_vector_alloc(obj->size);
	    gsl_vector_memcpy(base, obj);
	}
    }
    
    // Copy assignment
    gslVector& operator=(const gslVector& obj) {
	if (this != &obj) {
	    if (base != nullptr)
		gsl_vector_free(base);
	    if (obj.gsl() != nullptr) {
		base = gsl_vector_alloc(obj.size());
		gsl_vector_memcpy(base, obj.gsl());
	    } else {
		base = nullptr;
	    }
	}
	return *this;
    }

    // Assign from base gsl_vector
    gslVector& operator=(const gsl_vector* obj) {
	if (base != nullptr)
	    gsl_vector_free(base);
	if (obj != nullptr) {
	    base = gsl_vector_alloc(obj->size);
	    gsl_vector_memcpy(base, obj);
	} else {
	    base = nullptr;
	}
	return *this;
    }

    // Move constructor
    gslVector(gslVector&& obj) {
	base = obj.base;
	obj.base = nullptr;
    }

    // Move assignment
    gslVector& operator=(gslVector&& obj) {
	if (this != &obj) {
	    if (base != nullptr)
		gsl_vector_free(base);
	    base = obj.base;
	    obj.base = nullptr;
	}
	return *this;
    }
    
    void create(int len) {
	if (base != nullptr)
	    std::cerr << "gslVector error: create() called on a non-null vector\n";
	base = gsl_vector_alloc(len);
    }
    
    void create(std::vector<double> in) {
	if (base != NULL)
	    std::cerr << "gslVector error: create() called on non-null vector\n";
	base = gsl_vector_alloc(in.size());
	for (unsigned i = 0; i < in.size(); i++) 
	    gsl_vector_set(base, i, in[i]);
    }

    void eraseAndResize(int len) {
	if (base != nullptr)
	    gsl_vector_free(base);
	base = gsl_vector_alloc(len);
    }	

    // Dereference operator
    gsl_vector* operator*() const {
	return base;
    }
    
    gsl_vector* gsl() const {
	return base;
    }
    
    const double operator[](int idx) const {
	return gsl_vector_get(base, idx);
    }
    
    double& operator[](int idx) {
	// To correctly serve as lvalue, have to point to the data directly
	// Skips gsl range checking
	return base->data[idx * base->stride];
    }
    
    int size() const {
	return base->size;
    }

    friend std::ostream& operator<< (std::ostream& stream, const gslVector& vec) {
	for (int i = 0; i < vec.size(); i++)
	    stream << vec[i] << " ";
	return stream;
    }
	    
 private:
    gsl_vector *base;
};



class gslVectorInt {
public:

    gslVectorInt()
	: base(nullptr) {}

    gslVectorInt(int len) {
	base = gsl_vector_int_alloc(len);
    }

    ~gslVectorInt() {
	if (base != nullptr)
	    gsl_vector_int_free(base);
    }

    gslVectorInt(const gslVectorInt& obj) {
	base = gsl_vector_int_alloc(obj.size());
	gsl_vector_int_memcpy(base, obj.gsl());
    }

    gslVectorInt(const gsl_vector_int* obj) {
	base = gsl_vector_int_alloc(obj->size);
	gsl_vector_int_memcpy(base, obj);
    }

    // Copy assignment
    gslVectorInt& operator=(const gslVectorInt& obj) {
	if (this != &obj) {
	    if (base != nullptr)
		gsl_vector_int_free(base);
	    base = gsl_vector_int_alloc(obj.size());
	    gsl_vector_int_memcpy(base, obj.gsl());
	}
	return *this;
    }

    // Assign from base gsl_vector_int
    gslVectorInt& operator=(const gsl_vector_int* obj) {
	if (base != nullptr)
	    gsl_vector_int_free(base);
	base = gsl_vector_int_alloc(obj->size);
	gsl_vector_int_memcpy(base, obj);
	return *this;
    }
    
    // Move constructor
    gslVectorInt(gslVectorInt&& obj) {
	base = obj.base;
	obj.base = nullptr;
    }

    // Move assignment
    gslVectorInt& operator=(gslVectorInt&& obj) {
	if (this != &obj) {
	    if (base != nullptr)
		gsl_vector_int_free(base);
	    base = obj.base;
	    obj.base = nullptr;
	}
	return *this;
    }
    
    void create(int len) {
	if (base != NULL)
	    std::cerr << "gslVector error: create() called on a non-null vector\n";
	base = gsl_vector_int_alloc(len);
    }
    
    void create(std::vector<int> in) {
	if (base != nullptr)
	    std::cerr << "gslVector error: create() called on non-null vector\n";
	base = gsl_vector_int_alloc(in.size());
	for (unsigned i = 0; i < in.size(); i++) 
	    gsl_vector_int_set(base, i, in[i]);
    }

    void eraseAndResize(int len) {
	if (base != nullptr)
	    gsl_vector_int_free(base);
	base = gsl_vector_int_alloc(len);
    }

    // Dereference operator
    gsl_vector_int* operator*() const {
	return base;
    }

    gsl_vector_int* gsl() const {
	return base;
    }
    
    const int operator[](int idx) const {
	return gsl_vector_int_get(base, idx);
    }
    
    int& operator[](int idx) {
	// To correctly serve as lvalue, have to point to the data directly
	// Skips gsl range checking
	return base->data[idx * base->stride];
    }
    
    int size() const {
	return base->size;
    }
	    
 private:
    gsl_vector_int *base;
};


void gsl_vector_print(const gsl_vector* vec, int size = -1);


#endif // gslWrapper_h_
