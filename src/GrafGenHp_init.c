
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void C_main(void *, void *, void *);

static const R_CMethodDef CEntries[] = {
  {"C_main", (DL_FUNC) &C_main, 3},
  {NULL, NULL, 0}
};

void R_init_GrafGen(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
