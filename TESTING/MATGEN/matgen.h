#include "slu_cdefs.h"
#include "slu_ddefs.h"
#include "slu_sdefs.h"
#include "slu_zdefs.h"

int clatms_slu(int *m, int *n, char *dist, int * iseed,
               char *sym, float *d, int *mode, float *cond, float *dmax__,
               int *kl, int *ku, char *pack, singlecomplex *a, int *lda,
               singlecomplex *work, int *info);

int clatb4_slu(char *path, int *imat, int *m, int * n,
               char *type, int *kl, int *ku, float *anorm, int *mode,
               float *cndnum, char *dist);

int dlatms_slu(int *m, int *n, char *dist, int * iseed,
               char *sym, double *d, int *mode, double *cond, double *dmax__,
               int *kl, int *ku, char *pack, double *a, int *lda,
               double *work, int *info);

int dlatb4_slu(char *path, int *imat, int *m, int * n,
               char *type, int *kl, int *ku, double *anorm, int * mode,
               double *cndnum, char *dist);

int slatb4_slu(char *path, int *imat, int *m, int *n,
               char *type, int *kl, int *ku, float *anorm, int *mode,
               float *cndnum, char *dist);

int slatms_slu(int *m, int *n, char *dist, int *iseed,
               char *sym, float *d, int *mode, float *cond, float *dmax__,
               int *kl, int *ku, char *pack, float *a, int *lda,
               float * work, int *info);

int zlatb4_slu(char *path, int *imat, int *m, int *n,
               char *type, int *kl, int *ku, double *anorm, int *mode,
               double *cndnum, char *dist);

int zlatms_slu(int *m, int *n, char *dist, int *iseed,
               char *sym, double *d, int *mode, double *cond, double *dmax__,
               int *kl, int *ku, char *pack, doublecomplex *a, int *lda,
               doublecomplex *work, int *info);
