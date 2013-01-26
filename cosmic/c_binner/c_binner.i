%module c_binner

 //%include "std_shared_ptr.i";
 //%include "std_string.i"
 //%include "std_vector.i"
 //%include "std_complex.i"
%include "stdint.i"

%{
#define SWIG_FILE_WITH_INIT
#include "Python.h"

#include <vector>
#include <complex>
#include <iostream>
%}

// Get the NumPy typemaps
//must come after def of SWIG_FILE_WITH_INIT as the numpy def needs it. (swig parses in-order)
%include "numpy.i"

%init %{
  import_array();
  //if we need ufunc support from numpy
  //import_ufuncs();
%}

//%template(VectorFloat) std::vector<float>;
//%template(VectorDouble) std::vector<double>;
//%template(VectorCDouble) std::vector<std::complex<double> >;

//save for a rainy day
//%typemap(out) std::vector<double> {
//    std::cout << "test" << std::endl;
//    int length = $1.size();
//    $result = PyArray_FromDims(1, &length, PyArray_DOUBLE);
//    memcpy(PyArray_DATA($result), $1.data(), sizeof(double) * length);
//}

%{
void binner(uint64_t *timestamps, unsigned timestamps_len, int64_t *bins, unsigned bins_len)
{
  //std::cout << "bin that crap" << std::endl;
  //std::cout << "size of timestamps is " << timestamps_len << std::endl;
  //std::cout << "size of bins is " << bins_len << std::endl;
  //code goes here
  int64_t ts;
  for (int i=0; i<timestamps_len; i++) {
    ts = timestamps[i];
    if (ts >= 0 && ts < bins_len) bins[ts]++;
  }
}
%}

%typemap(in) (uint64_t *timestamps, unsigned timestamps_len) {
    /////
    //should check numpy array type (make sure is uint)
  if (!(PyArray_Check($input))) {
      PyErr_SetString(PyExc_RuntimeError, "not  NUMPY array");
    }
  int type = PyArray_TYPE($input);

  $1 = (uint64_t*) PyArray_GETPTR1($input, 0);   //get array ptr
  $2 = PyArray_DIM($input, 0);   //get num elements
}

// PyErr_SetString(PyExc_RuntimeError, "abc");

%typemap(in) (int64_t *bins, unsigned bins_len) {
    /////
    //should check numpy array type (make sure is int64)
  $1 = (int64_t*) PyArray_GETPTR1($input, 0);   //get array ptr
  $2 = PyArray_DIM($input, 0);   //get num elements
}
void binner(uint64_t *timestamps, unsigned timestamps_len, int64_t *bins, unsigned bins_len);

