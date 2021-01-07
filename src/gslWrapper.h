#ifndef HEADER_gslWrapper_
#define HEADER_gslWrapper_

#include <vector>
#include <iostream>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


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


class gslMatrix {
 public:
    gslMatrix() 
	: base(NULL) {}       

    gslMatrix(int rows, int cols) {
	base = gsl_matrix_alloc(rows, cols);
    } 

    ~gslMatrix() {
	if (base != NULL)
	    gsl_matrix_free(base);
    }

    void create(int rows, int cols) {
	if (base != NULL)
	    std::cerr << "gslMatrix error: create() called on a non-null matrix\n";
	base = gsl_matrix_alloc(rows, cols);
    }

    gslMatrix(const gslMatrix &obj) {
	if (base != NULL)
	    std::cerr << "copy constructor called on non-null base\n";
	base = gsl_matrix_alloc(obj.nrows(), obj.ncols());
	gsl_matrix_memcpy(base, obj.gsl());
    }

    gslMatrix &operator= (const gslMatrix &obj) {
	if (base != NULL)
	    gsl_matrix_free(base);
	base = gsl_matrix_alloc(obj.nrows(), obj.ncols());
	gsl_matrix_memcpy(base, obj.gsl());
	return *this;
    }

    int nrows() const {
	return base->size1;
    }
    int ncols() const {
	return base->size2;
    }
    gsl_matrix *gsl() const {
	return base;
    }

    void print() {
	std::cout << "{\n";
	for (unsigned i = 0; i < nrows(); i++) {
	    for (unsigned j = 0; j < ncols(); j++)
		std::cout << (*this)[i][j] << " ";
	    std::cout << std::endl;
	}
	std::cout << "}\n";
    }
    friend std::ostream& operator<<(std::ostream &stream, gslMatrix mat) {
	stream << "{\n";
	for (unsigned i = 0; i < mat.nrows(); i++) {
	    for (unsigned j = 0; j < mat.ncols(); j++)
		stream << mat[i][j] << " ";
	    stream << std::endl;
	}
	stream << "}\n";
	return stream;
    }

    // TODO: This is potentially dangerous, as the view is in scope only as long as the matrix it 
    // came from is in scope. In theory, a user could create an object of class gslMatrixRow, and
    // keep it in scope. Figure out how to stop this.

    class gslMatrixRow {
	friend class gslMatrix;
    public:
	const double operator[](int idx) const {
	    return gsl_vector_get(&(rowBase.vector), idx);
	}

	// TODO: Decide whether to keep bounds checking here (no gsl check because of direct access to data)
	double& operator[](int idx) {
	    if (idx >= rowBase.vector.size)
	      std::cerr << "Invalid col access: " << idx << " " << rowBase.vector.size << std::endl;
	    return rowBase.vector.data[idx * rowBase.vector.stride];
	}
	
	// WARNING: Potentially very unsafe
	gsl_vector_view *asView() {
	    return &rowBase;
	}

    private:
	
        gslMatrixRow(gsl_vector_view row) 
	    : rowBase(row) {}

	gsl_vector_view rowBase;
    };

    // These cannot be references, as it would return a reference
    // to the gslMatrixRow object. But as long as we have the view, it should be OK, I think.
    const gslMatrixRow operator[](int idx) const {
	return gslMatrixRow(gsl_matrix_row(base, idx));
    }
    gslMatrixRow operator[](int idx) {
	return gslMatrixRow(gsl_matrix_row(base, idx));
    }

 private:
    gsl_matrix *base;
};


#endif // gslWrapper_h_
