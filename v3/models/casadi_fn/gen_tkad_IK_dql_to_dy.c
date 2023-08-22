/* This file was automatically generated by CasADi 3.6.2.
 *  It consists of: 
 *   1) content generated by CasADi runtime: not copyrighted
 *   2) template code copied from CasADi source: permissively licensed (MIT-0)
 *   3) user code: owned by the user
 *
 */
#ifdef __cplusplus
extern "C" {
#endif

/* How to prefix internal symbols */
#ifdef CASADI_CODEGEN_PREFIX
  #define CASADI_NAMESPACE_CONCAT(NS, ID) _CASADI_NAMESPACE_CONCAT(NS, ID)
  #define _CASADI_NAMESPACE_CONCAT(NS, ID) NS ## ID
  #define CASADI_PREFIX(ID) CASADI_NAMESPACE_CONCAT(CODEGEN_PREFIX, ID)
#else
  #define CASADI_PREFIX(ID) gen_tkad_IK_dql_to_dy_ ## ID
#endif

#include <math.h>
#include <string.h>
#ifdef MATLAB_MEX_FILE
#include <mex.h>
#endif

#ifndef casadi_real
#define casadi_real double
#endif

#ifndef casadi_int
#define casadi_int long long int
#endif

/* Add prefix to internal symbols */
#define casadi_f0 CASADI_PREFIX(f0)
#define casadi_fill CASADI_PREFIX(fill)
#define casadi_from_mex CASADI_PREFIX(from_mex)
#define casadi_s0 CASADI_PREFIX(s0)
#define casadi_sq CASADI_PREFIX(sq)
#define casadi_to_mex CASADI_PREFIX(to_mex)

/* Symbol visibility in DLLs */
#ifndef CASADI_SYMBOL_EXPORT
  #if defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__)
    #if defined(STATIC_LINKED)
      #define CASADI_SYMBOL_EXPORT
    #else
      #define CASADI_SYMBOL_EXPORT __declspec(dllexport)
    #endif
  #elif defined(__GNUC__) && defined(GCC_HASCLASSVISIBILITY)
    #define CASADI_SYMBOL_EXPORT __attribute__ ((visibility ("default")))
  #else
    #define CASADI_SYMBOL_EXPORT
  #endif
#endif

casadi_real casadi_sq(casadi_real x) { return x*x;}

void casadi_fill(casadi_real* x, casadi_int n, casadi_real alpha) {
  casadi_int i;
  if (x) {
    for (i=0; i<n; ++i) *x++ = alpha;
  }
}

#ifdef MATLAB_MEX_FILE
casadi_real* casadi_from_mex(const mxArray* p, casadi_real* y, const casadi_int* sp, casadi_real* w) {
  casadi_int nrow, ncol, is_sparse, c, k, p_nrow, p_ncol;
  const casadi_int *colind, *row;
  mwIndex *Jc, *Ir;
  const double* p_data;
  if (!mxIsDouble(p) || mxGetNumberOfDimensions(p)!=2)
    mexErrMsgIdAndTxt("Casadi:RuntimeError",
      "\"from_mex\" failed: Not a two-dimensional matrix of double precision.");
  nrow = *sp++;
  ncol = *sp++;
  colind = sp;
  row = sp+ncol+1;
  p_nrow = mxGetM(p);
  p_ncol = mxGetN(p);
  is_sparse = mxIsSparse(p);
  Jc = 0;
  Ir = 0;
  if (is_sparse) {
    Jc = mxGetJc(p);
    Ir = mxGetIr(p);
  }
  p_data = (const double*)mxGetData(p);
  if (p_nrow==1 && p_ncol==1) {
    casadi_int nnz;
    double v = is_sparse && Jc[1]==0 ? 0 : *p_data;
    nnz = sp[ncol];
    casadi_fill(y, nnz, v);
  } else {
    casadi_int tr = 0;
    if (nrow!=p_nrow || ncol!=p_ncol) {
      tr = nrow==p_ncol && ncol==p_nrow && (nrow==1 || ncol==1);
      if (!tr) mexErrMsgIdAndTxt("Casadi:RuntimeError",
                                 "\"from_mex\" failed: Dimension mismatch. "
                                 "Expected %d-by-%d, got %d-by-%d instead.",
                                 nrow, ncol, p_nrow, p_ncol);
    }
    if (is_sparse) {
      if (tr) {
        for (c=0; c<ncol; ++c)
          for (k=colind[c]; k<colind[c+1]; ++k) w[row[k]+c*nrow]=0;
        for (c=0; c<p_ncol; ++c)
          for (k=Jc[c]; k<(casadi_int) Jc[c+1]; ++k) w[c+Ir[k]*p_ncol] = p_data[k];
        for (c=0; c<ncol; ++c)
          for (k=colind[c]; k<colind[c+1]; ++k) y[k] = w[row[k]+c*nrow];
      } else {
        for (c=0; c<ncol; ++c) {
          for (k=colind[c]; k<colind[c+1]; ++k) w[row[k]]=0;
          for (k=Jc[c]; k<(casadi_int) Jc[c+1]; ++k) w[Ir[k]]=p_data[k];
          for (k=colind[c]; k<colind[c+1]; ++k) y[k]=w[row[k]];
        }
      }
    } else {
      for (c=0; c<ncol; ++c) {
        for (k=colind[c]; k<colind[c+1]; ++k) {
          y[k] = p_data[row[k]+c*nrow];
        }
      }
    }
  }
  return y;
}

#endif

#define casadi_to_double(x) ((double) x)

#ifdef MATLAB_MEX_FILE
mxArray* casadi_to_mex(const casadi_int* sp, const casadi_real* x) {
  casadi_int nrow, ncol, c, k;
#ifndef CASADI_MEX_NO_SPARSE
  casadi_int nnz;
#endif
  const casadi_int *colind, *row;
  mxArray *p;
  double *d;
#ifndef CASADI_MEX_NO_SPARSE
  casadi_int i;
  mwIndex *j;
#endif /* CASADI_MEX_NO_SPARSE */
  nrow = *sp++;
  ncol = *sp++;
  colind = sp;
  row = sp+ncol+1;
#ifndef CASADI_MEX_NO_SPARSE
  nnz = sp[ncol];
  if (nnz!=nrow*ncol) {
    p = mxCreateSparse(nrow, ncol, nnz, mxREAL);
    for (i=0, j=mxGetJc(p); i<=ncol; ++i) *j++ = *colind++;
    for (i=0, j=mxGetIr(p); i<nnz; ++i) *j++ = *row++;
    if (x) {
      d = (double*)mxGetData(p);
      for (i=0; i<nnz; ++i) *d++ = casadi_to_double(*x++);
    }
    return p;
  }
#endif /* CASADI_MEX_NO_SPARSE */
  p = mxCreateDoubleMatrix(nrow, ncol, mxREAL);
  if (x) {
    d = (double*)mxGetData(p);
    for (c=0; c<ncol; ++c) {
      for (k=colind[c]; k<colind[c+1]; ++k) {
        d[row[k]+c*nrow] = casadi_to_double(*x++);
      }
    }
  }
  return p;
}

#endif

static const casadi_int casadi_s0[6] = {2, 1, 0, 2, 0, 1};

/* tkad_IK_vel:(i0[2],i1[2])->(o0[2]) */
static int casadi_f0(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem) {
  casadi_real a0, a1, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a2, a20, a3, a4, a5, a6, a7, a8, a9;
  a0=arg[1]? arg[1][0] : 0;
  a1=2.0799999999999999e-02;
  a2=7.;
  a3=arg[0]? arg[0][1] : 0;
  a4=1.9896753472735356e+00;
  a3=(a3+a4);
  a4=cos(a3);
  a4=(a2*a4);
  a5=2500.;
  a4=(a4/a5);
  a1=(a1-a4);
  a4=3.1573672058406521e-03;
  a1=(a1-a4);
  a4=2.0522886837964237e-02;
  a6=91.;
  a7=cos(a3);
  a7=(a6*a7);
  a8=5000.;
  a7=(a7/a8);
  a4=(a4+a7);
  a7=1.3813481525552854e+02;
  a8=cos(a3);
  a8=(a7*a8);
  a9=50000.;
  a8=(a8/a9);
  a4=(a4-a8);
  a8=5.0276961068873298e+01;
  a10=sin(a3);
  a10=(a8*a10);
  a10=(a10/a9);
  a4=(a4-a10);
  a10=2.6135840000000000e-02;
  a4=(a4-a10);
  a10=(a1*a4);
  a9=-1.1491876815742470e-03;
  a11=sin(a3);
  a11=(a2*a11);
  a11=(a11/a5);
  a9=(a9-a11);
  a11=casadi_sq(a9);
  a5=casadi_sq(a1);
  a11=(a11+a5);
  a5=casadi_sq(a4);
  a11=(a11-a5);
  a11=sqrt(a11);
  a5=(a9*a11);
  a10=(a10-a5);
  a5=(a1*a11);
  a12=(a9*a4);
  a5=(a5+a12);
  a12=casadi_sq(a5);
  a13=casadi_sq(a10);
  a12=(a12+a13);
  a10=(a10/a12);
  a13=4.0000000000000002e-04;
  a14=sin(a3);
  a15=arg[1]? arg[1][1] : 0;
  a14=(a14*a15);
  a14=(a2*a14);
  a14=(a13*a14);
  a16=(a11*a14);
  a17=(a1+a1);
  a17=(a17*a14);
  a18=(a9+a9);
  a19=cos(a3);
  a19=(a19*a15);
  a2=(a2*a19);
  a13=(a13*a2);
  a18=(a18*a13);
  a17=(a17-a18);
  a18=(a4+a4);
  a2=2.0000000000000002e-05;
  a19=sin(a3);
  a19=(a19*a15);
  a7=(a7*a19);
  a7=(a2*a7);
  a19=2.0000000000000001e-04;
  a20=sin(a3);
  a20=(a20*a15);
  a6=(a6*a20);
  a19=(a19*a6);
  a7=(a7-a19);
  a3=cos(a3);
  a3=(a3*a15);
  a8=(a8*a3);
  a2=(a2*a8);
  a7=(a7-a2);
  a18=(a18*a7);
  a17=(a17-a18);
  a18=(a11+a11);
  a17=(a17/a18);
  a18=(a1*a17);
  a16=(a16+a18);
  a18=(a9*a7);
  a2=(a4*a13);
  a18=(a18-a2);
  a16=(a16+a18);
  a10=(a10*a16);
  a5=(a5/a12);
  a4=(a4*a14);
  a1=(a1*a7);
  a4=(a4+a1);
  a9=(a9*a17);
  a11=(a11*a13);
  a9=(a9-a11);
  a4=(a4-a9);
  a5=(a5*a4);
  a10=(a10-a5);
  a5=(a0-a10);
  if (res[0]!=0) res[0][0]=a5;
  a0=(a0+a10);
  if (res[0]!=0) res[0][1]=a0;
  return 0;
}

CASADI_SYMBOL_EXPORT int tkad_IK_vel(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem){
  return casadi_f0(arg, res, iw, w, mem);
}

CASADI_SYMBOL_EXPORT int tkad_IK_vel_alloc_mem(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT int tkad_IK_vel_init_mem(int mem) {
  return 0;
}

CASADI_SYMBOL_EXPORT void tkad_IK_vel_free_mem(int mem) {
}

CASADI_SYMBOL_EXPORT int tkad_IK_vel_checkout(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT void tkad_IK_vel_release(int mem) {
}

CASADI_SYMBOL_EXPORT void tkad_IK_vel_incref(void) {
}

CASADI_SYMBOL_EXPORT void tkad_IK_vel_decref(void) {
}

CASADI_SYMBOL_EXPORT casadi_int tkad_IK_vel_n_in(void) { return 2;}

CASADI_SYMBOL_EXPORT casadi_int tkad_IK_vel_n_out(void) { return 1;}

CASADI_SYMBOL_EXPORT casadi_real tkad_IK_vel_default_in(casadi_int i) {
  switch (i) {
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* tkad_IK_vel_name_in(casadi_int i) {
  switch (i) {
    case 0: return "i0";
    case 1: return "i1";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* tkad_IK_vel_name_out(casadi_int i) {
  switch (i) {
    case 0: return "o0";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* tkad_IK_vel_sparsity_in(casadi_int i) {
  switch (i) {
    case 0: return casadi_s0;
    case 1: return casadi_s0;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* tkad_IK_vel_sparsity_out(casadi_int i) {
  switch (i) {
    case 0: return casadi_s0;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT int tkad_IK_vel_work(casadi_int *sz_arg, casadi_int* sz_res, casadi_int *sz_iw, casadi_int *sz_w) {
  if (sz_arg) *sz_arg = 2;
  if (sz_res) *sz_res = 1;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 0;
  return 0;
}

#ifdef MATLAB_MEX_FILE
void mex_tkad_IK_vel(int resc, mxArray *resv[], int argc, const mxArray *argv[]) {
  casadi_int i;
  int mem;
  casadi_real w[27];
  casadi_int *iw = 0;
  const casadi_real* arg[2] = {0};
  casadi_real* res[1] = {0};
  if (argc>2) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"tkad_IK_vel\" failed. Too many input arguments (%d, max 2)", argc);
  if (resc>1) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"tkad_IK_vel\" failed. Too many output arguments (%d, max 1)", resc);
  if (--argc>=0) arg[0] = casadi_from_mex(argv[0], w, casadi_s0, w+6);
  if (--argc>=0) arg[1] = casadi_from_mex(argv[1], w+2, casadi_s0, w+6);
  --resc;
  res[0] = w+4;
  tkad_IK_vel_incref();
  mem = tkad_IK_vel_checkout();
  i = tkad_IK_vel(arg, res, iw, w+6, mem);
  if (i) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"tkad_IK_vel\" failed.");
  tkad_IK_vel_release(mem);
  tkad_IK_vel_decref();
  if (res[0]) resv[0] = casadi_to_mex(casadi_s0, res[0]);
}
#endif


#ifdef MATLAB_MEX_FILE
void mexFunction(int resc, mxArray *resv[], int argc, const mxArray *argv[]) {
  char buf[12];
  int buf_ok = argc > 0 && !mxGetString(*argv, buf, sizeof(buf));
  if (!buf_ok) {
    mex_tkad_IK_vel(resc, resv, argc, argv);
    return;
  } else if (strcmp(buf, "tkad_IK_vel")==0) {
    mex_tkad_IK_vel(resc, resv, argc-1, argv+1);
    return;
  }
  mexErrMsgTxt("First input should be a command string. Possible values: 'tkad_IK_vel'");
}
#endif
#ifdef __cplusplus
} /* extern "C" */
#endif
