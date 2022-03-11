#include "definition.h"
#include "vec_mtx.h"
#include "vec_mtx.c"
#include "jacobi.c"

/* ---------------------------------------
 * computing the error of A*v - lamdba*v. v, lamdba is eigenvector and eigenvalue respectively.
 *      ori_A: a matrix N by N
 *      A: a diagonal matrix of multiple eigenvalues
 *      P: a matrix of multiple eigenvectors
 *      n: dimension      
 */
void error(double** ori_A, double** A,double** P, int n){

    int i, j, k;
    double *t, *v;
    double lambda;
    t = alloc_vec(n);
    v = alloc_vec(n);
    printf("-------------------------\nerror: ||A*v-lambda*v||\n");
    for(j = 0; j < n; j++){
        lambda = A[j][j];
        for(i = 0; i < n; i++){
            v[i] = P[i][j];
        }
        mtx_vec_mult(t, ori_A, v, n);
        for(i = 0; i < n; i++){
            t[i] = t[i] - lambda*v[i];
        }
        printf("j=%d: %.10lf\n", j, vec_norm(t, n));
    }
    printf("------------------------------------------\n");

    free(t);
    free(v);
}


/* ---------------------------------------------------------
 * computing inner_product every two eigenvector
 *      P: a matrix of multiple eigenvectors
 *      n: dimension       
 */
void inner_product_forall(double** P, int n){
    double *a, *b;
    int i, j, k;
    a = alloc_vec(n);
    b = alloc_vec(n);
    printf("------------------------------------------\n");
    for(i = 0; i < n; i++){
        for(j = i; j < n; j++){
            for(k = 0; k < n; k++){
                if(i == j) a[k] = P[i][k];
                b[k] = P[j][k];
            }
            printf("< v_%d, v_%d > : %.10lf\n", i, j, inner_product(a, b, n));
        }
    }
    printf("----------------------------------------\n");
    free(a);
    free(b);
}


int main(){

    int N = 4;
    double **A, **P, **ori_A;
    int i, j, k;
    FILE* fp;
    char filename[30] = "iterations.txt";
    fp = fopen(filename, "w");

    for(N; N <=20; N++){

        ori_A = alloc_mtx(N);
        gen_sym_mtx(ori_A, N);
        A = alloc_mtx(N);
        gen_sym_mtx(A, N);
    

        P = alloc_mtx(N);
        make_identity_mtx(P, N);
        //print_mtx(A, N);
        //print_mtx(P, N);

        k = jacobian_method(A, P, N);
        //print_mtx(P, N);
        //error(ori_A, A, P, N);
        //inner_product_forall(P, N);

        fprintf(fp, "%d %d\n", N, k);

        free_mtx(A, N);
        free_mtx(P, N);
        free_mtx(ori_A, N);
    }
    fclose(fp);
    return 0;
}
