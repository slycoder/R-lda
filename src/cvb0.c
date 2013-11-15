#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>

#define CHECK(VAL, TYPE) if (!is##TYPE(VAL)) { \
    error(#VAL " must be a(n) " #TYPE "."); \
  }

#define CHECKLEN(VAL, TYPE, LEN) if (!is##TYPE(VAL) || length(VAL) != LEN) { \
  error(#VAL " must be a length " #LEN " " #TYPE "."); \
  }

#define CHECKMATROW(VAL, TYPE, NROW) if (!isMatrix(VAL) || !is##TYPE(VAL) || NUMROWS(VAL) != NROW) { \
    error(#VAL " must be a matrix with " #NROW " rows of type " #TYPE "."); \
  }

#define NUMROWS(MAT) (INTEGER(GET_DIM(MAT))[0])
#define NUMCOLS(MAT) (INTEGER(GET_DIM(MAT))[1])

SEXP cvb0(
    SEXP documents,
    SEXP K_,
    SEXP V_,
    SEXP N_,
    SEXP alpha_,
    SEXP eta_,
    SEXP trace_) {
  GetRNGstate();
  int dd;
  int ii;
  int kk;

  CHECKLEN(alpha_, Real, 1);
  double alpha = REAL(alpha_)[0];

  CHECKLEN(eta_, Real, 1);
  double eta = REAL(eta_)[0];

  CHECKLEN(K_, Integer, 1);
  int K = INTEGER(K_)[0];

  CHECKLEN(N_, Integer, 1);
  int N = INTEGER(N_)[0];

  CHECKLEN(trace_, Integer, 1);
  int trace = INTEGER(trace_)[0];

  CHECKLEN(V_, Integer, 1);
  int V = INTEGER(V_)[0];

  CHECK(documents, NewList);
  int nd = length(documents);

  SEXP retval;
  PROTECT(retval = allocVector(VECSXP, 4));

  SEXP assignments;
  SEXP topics = NULL;
  SEXP topic_sums = NULL;
  SEXP document_sums = NULL;

  SET_VECTOR_ELT(retval, 0, assignments = allocVector(VECSXP, nd));
  SET_VECTOR_ELT(retval, 1, topics = allocMatrix(REALSXP, K, V));
  SET_VECTOR_ELT(retval, 2, topic_sums = allocMatrix(REALSXP, K, length(V_)));
  SET_VECTOR_ELT(retval, 3, document_sums = allocMatrix(REALSXP, K, nd));

  for (ii = 0; ii < K * nd; ++ii) {
    REAL(document_sums)[ii] = 0;
  }

  for (dd = 0; dd < nd; ++dd) {
    int ww;
    SEXP document = VECTOR_ELT(documents, dd);

    CHECKMATROW(document, Integer, 2);

    int nw = INTEGER(GET_DIM(document))[1];
    SET_VECTOR_ELT(assignments, dd, allocMatrix(REALSXP, K, nw));
    SEXP zs = VECTOR_ELT(assignments, dd);

    for (ww = 0; ww < nw; ++ww) {
      int word = INTEGER(document)[ww * 2];
      int count = INTEGER(document)[ww * 2 + 1];
      if (count < 0) {
        error("Count must be positive.");
      }
      if (word >= V || word < 0) {
        error("Word (%d) must be positive and less than or "
            "equal to the number of words (%d).", word, V);
      }
    }
  }

  double* p = (double *)R_alloc(K, sizeof(double));

  int initialized = 0;
  int iteration;
  for (iteration = 0; iteration < N; ++iteration) {
    if (trace >= 1) {
      REprintf("Iteration %d\n", iteration);
    }
    for (dd = 0; dd < nd; ++dd) {
      R_CheckUserInterrupt();
      SEXP zs = VECTOR_ELT(assignments, dd);
      SEXP document = VECTOR_ELT(documents, dd);
      int ww;
      int nw = INTEGER(GET_DIM(document))[1];

      for (ww = 0; ww < nw; ++ww) {
        double* z = &REAL(zs)[ww * K];
        int word = -1;
        int count = 1;
        double* topic_wk;
        double* topic_k;
        double* document_k;

        word = INTEGER(document)[ww * 2];
        count = INTEGER(document)[ww * 2 + 1];

        topic_wk = &REAL(topics)[K * word];
        topic_k = REAL(topic_sums);
        document_k = &REAL(document_sums)[K * dd];
        if (initialized) {
          for (kk = 0; kk < K; ++kk) {
            topic_wk[kk] -= count * z[kk];
            topic_k[kk] -= count * z[kk];
            document_k[kk] -= count * z[kk];
          }
        }

        double p_sum = 0.0;
        for (kk = 0; kk < K; ++kk) {
          if (!initialized) {
            p[kk] = rgamma(1.0, 1.0);
          } else {
            p[kk] = (REAL(document_sums)[K * dd + kk] + alpha);
            p[kk] *= (REAL(topics)[kk + K * word] + eta);
            p[kk] /= (REAL(topic_sums)[kk] + V * eta);
          }
          p_sum += p[kk];
        }

        if (p_sum <= 0.0) {
          for (kk = 0; kk < K; ++kk) {
            p[kk] = 1.0;
          }
          p_sum = K;
          REprintf("Warning:Sum of probabilities is zero, assigning equal probabilities.\n");
        }

        for (kk = 0; kk < K; ++kk) {
          document_k[kk] += p[kk] / p_sum * count;
          topic_k[kk] += p[kk] / p_sum * count;
          topic_wk[kk] += p[kk] / p_sum * count;
          z[kk] = p[kk] / p_sum;
        }
      }
    }
    initialized = 1;
  }

  PutRNGstate();
  UNPROTECT(1);
  return retval;
}
