#ifndef HEADER_gslWrapper_
#define HEADER_gslWrapper_

#include <cassert>
#include <vector>
#include <initializer_list>
#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

class gslVector {
public:
    
  gslVector()
    : base(nullptr) {}
    
  gslVector(size_t len) {
    base = gsl_vector_alloc(len);
  }

  ~gslVector() {
    if (base != nullptr)
      gsl_vector_free(base);
  }

  // Initializer list constructor
  gslVector(std::initializer_list<double> list) {
    base = gsl_vector_alloc(list.size());
    int i = 0;
    for (double element : list)
      (*this)[i++] = element;
  }

  // Initializer list assignment
  gslVector& operator=(std::initializer_list<double> list) {
    if (base != nullptr)
      gsl_vector_free(base);
    base = gsl_vector_alloc(list.size());
    int i = 0;
    for (double element : list)
      (*this)[i++] = element;
    return *this;
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

  // Swap two vectors without copying contents
  void swapWith(gslVector& other) {
    std::swap(base, other.base);
  }
    
  void alloc(size_t len) {
    assert(base == nullptr);
    base = gsl_vector_alloc(len);
  }
    
  void alloc(std::vector<double> in) {
    assert(base == nullptr);
    base = gsl_vector_alloc(in.size());
    for (unsigned i = 0; i < in.size(); i++) 
      gsl_vector_set(base, i, in[i]);
  }

  void allocZero(size_t len) {
    assert(base == nullptr);
    base = gsl_vector_calloc(len);
  }
  
  void eraseRealloc(size_t len) {
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

  double operator[](size_t idx) const {
    return gsl_vector_get(base, idx);
  }
  double& operator[](size_t idx) {
    // To correctly serve as lvalue, have to point to the data directly
#ifndef GSL_RANGE_CHECK_OFF
    if (idx >= base->size) {
      gsl_error("index out of range", __FILE__, __LINE__, GSL_EINVAL);
      throw std::invalid_argument("gslVector operator[]: Index out of range");
    }
#endif
    return base->data[idx * base->stride];
  }

  // Return a pointer to the i-th element of the vector
  // Example use: Writing a binary value to file with ofstream::write()
  const double* ptr(size_t idx) const {
#ifndef GSL_RANGE_CHECK_OFF
    if (idx >= base->size) {
      gsl_error("index out of range", __FILE__, __LINE__, GSL_EINVAL);
      throw std::invalid_argument("gslVector ptr(): Index out of range");
    }
#endif
    return base->data + (idx * base->stride);
  }
    
  
  // Vector addition
  gslVector& operator+=(const gslVector& rhs) {
    gsl_vector_add(this->gsl(), rhs.gsl());
    return *this;
  }
  // lhs passed by value makes copy, so we don't modify original object
  friend gslVector operator+(gslVector lhs, const gslVector& rhs) {
    lhs += rhs;
    return lhs;
  }
  
  // Vector subtraction
  gslVector& operator-=(const gslVector& rhs) {
    gsl_vector_sub(this->gsl(), rhs.gsl());
    return *this;
  }
  friend gslVector operator-(gslVector lhs, const gslVector& rhs) {
    lhs -= rhs;
    return lhs;
  }

  // Scalar multiplication
  gslVector& operator*=(const double scalar) {
    gsl_vector_scale(this->gsl(), scalar);
    return *this;
  }
  friend gslVector operator*(gslVector vec, const double scalar) {
    vec *= scalar;
    return vec;
  }
  friend gslVector operator*(const double scalar, gslVector vec) {
    return vec * scalar;
  }

  size_t size() const {
    if (base == nullptr)
      return 0;
    return base->size;
  }

  friend std::ostream& operator<< (std::ostream& stream, const gslVector& vec) {
    for (size_t i = 0; i < vec.size(); i++)
      stream << vec[i] << " ";
    return stream;
  }

  void print(std::ostream& stream = std::cout) const {
    print(0, size(), stream);
  }
  void print(size_t num, std::ostream& stream = std::cout) const {
    print(0, num, stream);
  }
  void print(size_t start, size_t num, std::ostream& stream = std::cout) const {
    size_t end = std::min(start+num, size());
    for (size_t i = start; i < end; i++)
      stream << gsl_vector_get(base, i) << " ";
    stream << std::endl;
  }
  
private:
  gsl_vector *base;
};


class gslVectorInt {
public:

  gslVectorInt()
    : base(nullptr) {}

  gslVectorInt(size_t len) {
    base = gsl_vector_int_alloc(len);
  }

  ~gslVectorInt() {
    if (base != nullptr)
      gsl_vector_int_free(base);
  }

  gslVectorInt(const gslVectorInt& obj)
    : base(nullptr) {
    if (obj.gsl() != nullptr) {
      base = gsl_vector_int_alloc(obj.size());
      gsl_vector_int_memcpy(base, obj.gsl());
    }
  }

  gslVectorInt(const gsl_vector_int* obj)
    : base(nullptr) {
    if (obj != nullptr) {
      base = gsl_vector_int_alloc(obj->size);
      gsl_vector_int_memcpy(base, obj);
    }
  }

  // Initializer list constructor
  gslVectorInt(std::initializer_list<int> list) {
    base = gsl_vector_int_alloc(list.size());
    int i = 0;
    for (double element : list)
      (*this)[i++] = element;
  }

  // Initializer list assignment
  gslVectorInt& operator=(std::initializer_list<int> list) {
    if (base != nullptr)
      gsl_vector_int_free(base);
    base = gsl_vector_int_alloc(list.size());
    int i = 0;
    for (double element : list)
      (*this)[i++] = element;
    return *this;
  }

  // Copy assignment
  gslVectorInt& operator=(const gslVectorInt& obj) {
    if (this != &obj) {
      if (base != nullptr)
	gsl_vector_int_free(base);
      if (obj.gsl() != nullptr) {
	base = gsl_vector_int_alloc(obj.size());
	gsl_vector_int_memcpy(base, obj.gsl());
      } else {
	base = nullptr;
      }
    }
    return *this;
  }

  // Assign from base gsl_vector_int
  gslVectorInt& operator=(const gsl_vector_int* obj) {
    if (base != nullptr)
      gsl_vector_int_free(base);
    if (obj != nullptr) {
      base = gsl_vector_int_alloc(obj->size);
      gsl_vector_int_memcpy(base, obj);
    } else {
      base = nullptr;
    }
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
    
  void alloc(size_t len) {
    assert(base == nullptr);
    base = gsl_vector_int_alloc(len);
  }
    
  void alloc(std::vector<int> in) {
    assert(base == nullptr);
    base = gsl_vector_int_alloc(in.size());
    for (unsigned i = 0; i < in.size(); i++) 
      gsl_vector_int_set(base, i, in[i]);
  }

  void allocZero(size_t len) {
    assert(base == nullptr);
    base = gsl_vector_int_calloc(len);
  }
  
  void eraseRealloc(size_t len) {
    if (base != nullptr)
      gsl_vector_int_free(base);
    base = gsl_vector_int_alloc(len);
  }

  void setAll(int in) {
    assert(base != nullptr);
    gsl_vector_int_set_all(base, in);
  }

  // Dereference operator
  gsl_vector_int* operator*() const {
    return base;
  }
  gsl_vector_int* gsl() const {
    return base;
  }
  
  int operator[](size_t idx) const {
    return gsl_vector_int_get(base, idx);
  }
  int& operator[](size_t idx) {
    // To correctly serve as lvalue, have to point to the data directly
#ifndef GSL_RANGE_CHECK_OFF
    if (idx >= base->size) {
      gsl_error("index out of range", __FILE__, __LINE__, GSL_EINVAL);
      throw std::invalid_argument("gslVectorInt: Index out of range");
    }
#endif
    return base->data[idx * base->stride];
  }

  // Scalar multiplication
  gslVectorInt& operator*=(const int scalar) {
    gsl_vector_int_scale(this->gsl(), scalar);
    return *this;
  }
  friend gslVectorInt operator*(gslVectorInt vec, const int scalar) {
    vec *= scalar;
    return vec;
  }
  friend gslVectorInt operator*(const int scalar, gslVectorInt vec) {
    return vec * scalar;
  }
    
  size_t size() const {
    if (base == nullptr)
      return 0;
    return base->size;
  }

  friend std::ostream& operator<< (std::ostream& stream, const gslVectorInt& vec) {
    for (size_t i = 0; i < vec.size(); i++)
      stream << vec[i] << " ";
    return stream;
  }
  void print(std::ostream& stream = std::cout) const {
    print(0, size(), stream);
  }
  void print(size_t num, std::ostream& stream = std::cout) const {
    print(0, num, stream);
  }
  void print(size_t start, size_t num, std::ostream& stream = std::cout) const {
    size_t end = std::min(start+num, size());
    for (size_t i = start; i < end; i++)
      stream << gsl_vector_int_get(base, i) << " ";
    stream << std::endl;
  }
  
private:
  gsl_vector_int *base;
};


class gslMatrix {
public:

  gslMatrix() 
    : base(nullptr) {}       
  
  gslMatrix(size_t rows, size_t cols) {
    base = gsl_matrix_alloc(rows, cols);
  } 
  
  ~gslMatrix() {
    if (base != nullptr)
      gsl_matrix_free(base);
  }
  
  gslMatrix(const gslMatrix &obj)
    : base(nullptr) {
    if (obj.gsl() != nullptr) {
      base = gsl_matrix_alloc(obj.nrows(), obj.ncols());
      gsl_matrix_memcpy(base, obj.gsl());
    }
  }
  
  gslMatrix& operator=(const gslMatrix &obj) {
    if (this != &obj) {
      if (base != nullptr)
	gsl_matrix_free(base);
      if (obj.gsl() != nullptr) {
	base = gsl_matrix_alloc(obj.nrows(), obj.ncols());
	gsl_matrix_memcpy(base, obj.gsl());
      } else {
	base = nullptr;
      }
    }
    return *this;
  }
  gslMatrix& operator=(const gsl_matrix* obj) {
    if (base != nullptr)
      gsl_matrix_free(base);
    if (obj != nullptr) {
      base = gsl_matrix_alloc(obj->size1, obj->size2);
      gsl_matrix_memcpy(base, obj);
    } else {
      base = nullptr;
    }
    return *this;
  }
  
  // Move semantics
  gslMatrix(gslMatrix&& obj) {
    base = obj.base;
    obj.base = nullptr;
  }
  
  gslMatrix &operator=(gslMatrix&& obj) {
    if (this != &obj) {
      if (base != nullptr)
	gsl_matrix_free(base);
      base = obj.base;
      obj.base = nullptr;
    }
    return *this;
  }

  void alloc(size_t rows, size_t cols) {
    assert(base == nullptr);
    base = gsl_matrix_alloc(rows, cols);
  }

  void allocZero(size_t rows, size_t cols) {
    assert(base == nullptr);
    base = gsl_matrix_calloc(rows, cols);
  }
  
  void eraseRealloc(size_t rows, size_t cols) {
    if (base != nullptr)
      gsl_matrix_free(base);
    base = gsl_matrix_alloc(rows, cols);
  }
  
  size_t nrows() const {
    if (base == nullptr)
      return 0;
    return base->size1;
  }
  size_t ncols() const {
    if (base == nullptr)
      return 0;
    return base->size2;
  }
  
  gsl_matrix* gsl() const {
    return base;
  }
  gsl_matrix* operator*() const {
    return base;
  }

  // Matrix addition
  gslMatrix& operator+=(const gslMatrix& rhs) {
    gsl_matrix_add(this->gsl(), rhs.gsl());
    return *this;
  }
  friend gslMatrix operator+(gslMatrix lhs, const gslMatrix& rhs) {
    lhs += rhs;
    return lhs;
  }
  
  // Scalar multiplication. Overload for mat * scalar and scalar * mat
  gslMatrix& operator*=(const double mult) {
    gsl_matrix_scale(this->gsl(), mult);
    return *this;
  }
  friend gslMatrix operator*(gslMatrix mat, const double mult) {
    mat *= mult;
    return mat;
  }
  friend gslMatrix operator*(const double mult, gslMatrix mat) {
    return mat * mult;
  }

  // Proxy class to enable matrix[i][j] syntax
  class gslMatrixRow {
  private:
    friend class gslMatrix;
    gslMatrixRow(gsl_matrix* mat_, size_t row_)
      : mat(mat_), row(row_) { }
    
  public:
    double operator[](size_t col) const {
      return gsl_matrix_get(mat, row, col);
    }
    double& operator[](size_t col) {
#ifndef GSL_RANGE_CHECK_OFF
      if (col >= mat->size2) {
	gsl_error("index out of range", __FILE__, __LINE__, GSL_EINVAL);
	throw std::invalid_argument("gslMatrix: Column index out of range");
      }
#endif
      return mat->data[row * mat->tda + col];
    }
    
  private:
    gsl_matrix* mat;
    size_t row;
  };

  gslMatrixRow operator[](size_t row) {
#ifndef GSL_RANGE_CHECK_OFF
    if (row >= base->size1) {
      gsl_error("index out of range", __FILE__, __LINE__, GSL_EINVAL);
      throw std::invalid_argument("gslMatrix: Row index out of range");
    }
#endif
    return gslMatrixRow(base, row);
  }

  const gslMatrixRow operator[](size_t row) const {
#ifndef GSL_RANGE_CHECK_OFF
    if (row >= base->size1) {
      gsl_error("index out of range", __FILE__, __LINE__, GSL_EINVAL);
      throw std::invalid_argument("gslMatrix: Row index out of range");
    }
#endif
    return gslMatrixRow(base, row);
  }
  
  void print(std::ostream& stream = std::cout)  {
    print(0, 0, nrows(), ncols(), stream);
  }
  void print(size_t numRows, size_t numCols, std::ostream& stream = std::cout)  {
    print(0, 0, numRows, numCols, stream);
  }
  void print(size_t rowStart, size_t numRows, size_t colStart, size_t numCols, std::ostream& stream = std::cout)  {
    size_t rowEnd = std::min(rowStart + numRows, nrows());
    size_t colEnd = std::min(colStart + numCols, ncols());
    std::cout << "{\n";
    for (size_t i = rowStart; i < rowEnd; i++) {
      for (size_t j = colStart; j < colEnd; j++)
	std::cout << (*this)[i][j] << " ";
      std::cout << std::endl;
    }
    std::cout << "}\n";
  }
  friend std::ostream& operator<<(std::ostream &stream, gslMatrix mat) {
    stream << "{\n";
    for (size_t i = 0; i < mat.nrows(); i++) {
      for (size_t j = 0; j < mat.ncols(); j++)
	stream << mat[i][j] << " ";
      stream << std::endl;
    }
    stream << "}\n";
    return stream;
  }

private:
  gsl_matrix *base;
};


#endif // gslWrapper_h_
