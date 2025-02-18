// EXAMPLE EXTERNAL FUNCTION 

#include <math.h>
#include <stdlib.h>
#include <string.h>
 
#ifdef _MSC_VER
#define EXPORT __declspec(dllexport)
#else
#define EXPORT
#endif
 
static const char *error = NULL;
 
EXPORT int init(const char *str) {
  return 1;
}
 
EXPORT const char * getLastError() {
  return error;
}
 
EXPORT int eval(const char *func,
                              int nArgs,
                              const double **inReal,
                              const double **inImag,
                              int blockSize,
                              double *outReal,
                              double *outImag) {
  int i;
 
  if (strcmp("extsinc", func) == 0) {
    if (nArgs != 1) {
      error = "One argument expected";
      return 0;
    }
    for (i = 0; i < blockSize; i++) {
      double x = inReal[0][i];
      outReal[i] = (x == 0) ? 1 : sin(x) / x;
    }
    return 1;
  }
  else {
    error = "Unknown function";
    return 0;
  }
}