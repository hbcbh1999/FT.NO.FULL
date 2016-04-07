#include<stdio.h>
#include<assert.h>
#include<math.h>

int FindInterpPoint(double x0, double* x, int dim, int* index);
double BilinearInterp(double** field, int li, int lj, double* x_field, double* y_field, double x, double y);
double TrilinearInterp(double*** field, int li, int lj, int lk, double* x_field, double* y_field, double* z_field, double x, double y, double z);
double TrilinearInterp_test(double (*f)(double,double,double), int li, int lj, int lk, double* x_field, double* y_field, double* z_field, double x, double y, double z);
double SimpleInterp_MAC(double*** field, int li, int lj, int lk, int dirc);
double ExactSolution(double x, double y, double z);

//Refer to JCP 101 pp334-348(1992), by Bell & Marcus
int main()
{
    FILE* file1 = fopen("10.dat", "r");
    FILE* file2 = fopen("20.dat", "r");
    FILE* fileout = fopen("norm10to20.dat","w");
/*    FILE* file1 = fopen("20.dat", "r");
    FILE* file2 = fopen("40.dat", "r");
    FILE* fileout = fopen("norm20to40.dat","w");
*/
    double ldomain[3]={0, 0, 0}, rdomain[3]={1, 0.4, 1};

    int i,j,k;
    int ind_max[3];

    //for first file

    int gridsize1[3];
    fscanf(file1, "%d  %d  %d", &gridsize1[0], &gridsize1[1], &gridsize1[2]);

    double U1[gridsize1[0]][gridsize1[1]][gridsize1[2]]; //u
    double V1[gridsize1[0]][gridsize1[1]][gridsize1[2]]; //v
    double W1[gridsize1[0]][gridsize1[1]][gridsize1[2]]; //w
    double p1[gridsize1[0]][gridsize1[1]][gridsize1[2]]; //density
    double c1[gridsize1[0]][gridsize1[1]][gridsize1[2]]; //concentration

    for (k = 0; k<gridsize1[2]; k++)
    for (j = 0; j<gridsize1[1]; j++)
    for (i = 0; i<gridsize1[0]; i++)
    {
        fscanf(file1, "%lf %lf %lf %lf %lf", &U1[i][j][k], &V1[i][j][k], &W1[i][j][k], &p1[i][j][k], &c1[i][j][k]);
    }
    double cellArea1;
    fscanf(file1, "%lf", &cellArea1);
    printf("cellArea1 = %12.8e\n", cellArea1);
    fclose(file1);

    //for second file

    int gridsize2[3];
    fscanf(file2, "%d  %d  %d", &gridsize2[0], &gridsize2[1], &gridsize2[2]);

    double*** U2;
    double*** V2;
    double*** W2;
    double*** p2;
    double*** c2;

    U2 = new double**[gridsize2[0]];
    V2 = new double**[gridsize2[0]];
    W2 = new double**[gridsize2[0]];
    p2 = new double**[gridsize2[0]];
    c2 = new double**[gridsize2[0]];

    for (i = 0; i<gridsize2[0]; i++)
    {
        U2[i] = new double*[gridsize2[1]];
        V2[i] = new double*[gridsize2[1]];
        W2[i] = new double*[gridsize2[1]];
        p2[i] = new double*[gridsize2[1]];
        c2[i] = new double*[gridsize2[1]];
        for (j = 0; j<gridsize2[1]; j++)
        {
            U2[i][j] = new double[gridsize2[2]];
            V2[i][j] = new double[gridsize2[2]];
            W2[i][j] = new double[gridsize2[2]];
            p2[i][j] = new double[gridsize2[2]];
            c2[i][j] = new double[gridsize2[2]];
        }
    }

    for (k = 0; k<gridsize2[2]; k++)
    for (j = 0; j<gridsize2[1]; j++)
    for (i = 0; i<gridsize2[0]; i++)
    {
        fscanf(file2, "%lf %lf %lf %lf %lf", &U2[i][j][k], &V2[i][j][k], &W2[i][j][k], &p2[i][j][k], &c2[i][j][k]);
    }
    double cellArea2;
    fscanf(file2, "%lf", &cellArea2);
    printf("cellArea2 = %12.8e\n", cellArea2);
    fclose(file2);

    //compute position for each point
    double x1[gridsize1[0]];
    double y1[gridsize1[1]];
    double z1[gridsize1[2]];
    double x2[gridsize2[0]];
    double y2[gridsize2[1]];
    double z2[gridsize2[2]];

    //file1
    double dx = (rdomain[0] - ldomain[0])/gridsize1[0];
    double dy = (rdomain[1] - ldomain[1])/gridsize1[1];
    double dz = (rdomain[2] - ldomain[2])/gridsize1[2];
    for (i=0; i<gridsize1[0]; i++)
        x1[i] = ldomain[0] + dx/2 + i*dx;

    for (j=0; j<gridsize1[1]; j++)
        y1[j] = ldomain[1] + dy/2 + j*dy;

    for (k=0; k<gridsize1[2]; k++)
        z1[k] = ldomain[2] + dz/2 + k*dz;

    //file2
    dx = (rdomain[0] - ldomain[0])/gridsize2[0];
    dy = (rdomain[1] - ldomain[1])/gridsize2[1];
    dz = (rdomain[2] - ldomain[2])/gridsize2[2];
    for (i=0; i<gridsize2[0]; i++)
        x2[i] = ldomain[0] + dx/2 + i*dx;

    for (j=0; j<gridsize2[1]; j++)
        y2[j] = ldomain[1] + dy/2 + j*dy;

    for (k=0; k<gridsize2[2]; k++)
        z2[k] = ldomain[2] + dz/2 + k*dz;

    //assume file2 to be the fine grid, to interpolate the points in file1
    double U1_new[gridsize1[0]][gridsize1[1]][gridsize1[2]];
    double V1_new[gridsize1[0]][gridsize1[1]][gridsize1[2]];
    double W1_new[gridsize1[0]][gridsize1[1]][gridsize1[2]];
    double p1_new[gridsize1[0]][gridsize1[1]][gridsize1[2]];
    double c1_new[gridsize1[0]][gridsize1[1]][gridsize1[2]];
    double error_U[gridsize1[0]][gridsize1[1]][gridsize1[2]];
    double error_V[gridsize1[0]][gridsize1[1]][gridsize1[2]];
    double error_W[gridsize1[0]][gridsize1[1]][gridsize1[2]];
    double error_p[gridsize1[0]][gridsize1[1]][gridsize1[2]];
    double error_c[gridsize1[0]][gridsize1[1]][gridsize1[2]];
    int li,lj,lk; //the lower left point index in file2 used to interpolate

    //error norm
    double U_Linf=0, V_Linf=0, W_Linf=0, p_Linf=0, c_Linf=0;
    double U_L1=0, V_L1=0, W_L1=0, p_L1=0, c_L1=0;
    double U_L2=0, V_L2=0, W_L2=0, p_L2=0, c_L2=0;

/*
    //debug lines
    double error;
    double (*p)(double, double, double);
    i = 9;
    j = 1;
    k = 0;
    p = ExactSolution;

    FindInterpPoint(x1[i], x2, gridsize2[0], &li);
    if (li==gridsize2[0])   li--;
    FindInterpPoint(y1[j], y2, gridsize2[1], &lj);
    if (lj==gridsize2[1])   lj--;
    FindInterpPoint(z1[k], z2, gridsize2[2], &lk);
    if (lk==gridsize2[2])   lk--;
    printf("(i, j, k) = (%d,%d,%d) in coarse grid\n", i, j, k);
    printf("(li, lj, lk) = (%d,%d,%d) in fine grid\n", li, lj, lk);

    error = TrilinearInterp_test(p, li, lj, lk, x2, y2, z2, x1[i], y1[j], z1[k]);
    printf("error from TrilinearInterp is %lf\n",error);
*/

    for (i=0; i<gridsize1[0]; i++)
    for (j=0; j<gridsize1[1]; j++)
    for (k=0; k<gridsize1[2]; k++)
    {
        FindInterpPoint(x1[i], x2, gridsize2[0], &li);
        if (li==gridsize2[0])   li--;

        FindInterpPoint(y1[j], y2, gridsize2[1], &lj);
        if (lj==gridsize2[1])   lj--;

        FindInterpPoint(z1[k], z2, gridsize2[2], &lk);
        if (lk==gridsize2[2])   lk--;

        //for MAC-grid, with staggered U and density
        U1_new[i][j][k] = SimpleInterp_MAC(U2, li, lj, lk, 0);
        V1_new[i][j][k] = SimpleInterp_MAC(V2, li, lj, lk, 1);
        W1_new[i][j][k] = SimpleInterp_MAC(W2, li, lj, lk, 2);
        p1_new[i][j][k] = TrilinearInterp(p2, li, lj, lk, x2, y2, z2, x1[i], y1[j], z1[k]);
        c1_new[i][j][k] = TrilinearInterp(c2, li, lj, lk, x2, y2, z2, x1[i], y1[j], z1[k]);

        error_U[i][j][k] = fabs(U1[i][j][k] - U1_new[i][j][k]);
        error_V[i][j][k] = fabs(V1[i][j][k] - V1_new[i][j][k]);
        error_W[i][j][k] = fabs(W1[i][j][k] - W1_new[i][j][k]);
        error_p[i][j][k] = fabs(p1[i][j][k] - p1_new[i][j][k]);
        error_c[i][j][k] = fabs(c1[i][j][k] - c1_new[i][j][k]);

        if (error_U[i][j][k]>U_Linf)
            U_Linf = error_U[i][j][k];
        if (error_V[i][j][k]>V_Linf)
            V_Linf = error_V[i][j][k];
        if (error_W[i][j][k]>W_Linf)
        {
            W_Linf = error_W[i][j][k];
            ind_max[0] = i;
            ind_max[1] = j;
            ind_max[2] = k;
        }
        if (error_p[i][j][k]>p_Linf)
            p_Linf = error_p[i][j][k];
        if (error_c[i][j][k]>c_Linf)
            c_Linf = error_c[i][j][k];

        U_L1 += error_U[i][j][k]*cellArea1;
        V_L1 += error_V[i][j][k]*cellArea1;
        W_L1 += error_W[i][j][k]*cellArea1;
        p_L1 += error_p[i][j][k]*cellArea1;
        c_L1 += error_c[i][j][k]*cellArea1;

        U_L2 += error_U[i][j][k]*error_U[i][j][k]*cellArea1;
        V_L2 += error_V[i][j][k]*error_V[i][j][k]*cellArea1;
        W_L2 += error_W[i][j][k]*error_W[i][j][k]*cellArea1;
        p_L2 += error_p[i][j][k]*error_p[i][j][k]*cellArea1;
        c_L2 += error_c[i][j][k]*error_c[i][j][k]*cellArea1;
    }

    U_L2 = sqrt(U_L2);
    V_L2 = sqrt(V_L2);
    W_L2 = sqrt(W_L2);
    p_L2 = sqrt(p_L2);
    c_L2 = sqrt(c_L2);

    fprintf(fileout,"u_Linf = %20.16g\n", U_Linf);
    fprintf(fileout,"v_Linf = %20.16g\n", V_Linf);
    fprintf(fileout,"w_Linf = %20.16g\n", W_Linf);
    fprintf(fileout,"rho_Linf = %20.16g\n", p_Linf);
    fprintf(fileout,"conc_Linf = %20.16g\n", c_Linf);

    fprintf(fileout,"u_L1 = %20.16g\n", U_L1);
    fprintf(fileout,"v_L1 = %20.16g\n", V_L1);
    fprintf(fileout,"w_L1 = %20.16g\n", W_L1);
    fprintf(fileout,"rho_L1 = %20.16g\n", p_L1);
    fprintf(fileout,"conc_L1 = %20.16g\n", c_L1);

    fprintf(fileout,"u_L2 = %20.16g\n", U_L2);
    fprintf(fileout,"v_L2 = %20.16g\n", V_L2);
    fprintf(fileout,"w_L2 = %20.16g\n", W_L2);
    fprintf(fileout,"rho_L2 = %20.16g\n", p_L2);
    fprintf(fileout,"conc_L2 = %20.16g\n", c_L2);
    fclose(fileout);

    printf("ind_max = (%d,%d,%d) of maximum W_Linf\n", ind_max[0], ind_max[1], ind_max[2]);
    return 0;
}

//find index, s.t. x[index] <= x0 < x[index+1]
int FindInterpPoint(double x0, double* x, int dim, int* index_addr)
{
    if (dim<=2)
    {
        (*index_addr) = 0;
        return 0;
    }

    int i,j;
    if (x[dim/2]>x0)
    {
        i = 0;
        j = dim/2;
    }
    else
    {
        i = dim/2;
        j = dim-1;
    }
    while ((j-i)>1)
    {
        if (x[i+(j-i)/2]>x0)
            j = i+(j-i)/2;
        else
            i = i + (j-i)/2;
    }
    (*index_addr) = i;

    if ((j-i)!=1)
        assert(false);

    return 0;
}

double BilinearInterp(double** field, int li, int lj,
        double* x_field, double* y_field,
        double x, double y)
{

    double value1, value2, value;
    value1 = field[li][lj] * (x_field[li+1] - x) / (x_field[li+1] - x_field[li])
         + field[li+1][lj] * (x - x_field[li]) / (x_field[li+1] - x_field[li]);
    value2 = field[li][lj+1] * (x_field[li+1] - x) / (x_field[li+1] - x_field[li])
         + field[li+1][lj+1] * (x - x_field[li]) / (x_field[li+1] - x_field[li]);
    value = value1 * (y_field[lj+1] - y) / (y_field[lj+1] - y_field[lj])
         + value2 * (y - y_field[lj]) / (y_field[lj+1] - y_field[lj]);

    return value;
}

double TrilinearInterp(double*** field, int li, int lj, int lk,
        double* x_field, double* y_field, double* z_field,
        double x, double y, double z)
{

    double c00, c10, c01, c11, c0, c1, c;

    c00 = field[li][lj][lk] * (x_field[li+1]-x) / (x_field[li+1]-x_field[li]) +
          field[li+1][lj][lk] * (x-x_field[li]) / (x_field[li+1]-x_field[li]);
    c10 = field[li][lj+1][lk] * (x_field[li+1]-x) / (x_field[li+1]-x_field[li]) +
          field[li+1][lj+1][lk] * (x-x_field[li]) / (x_field[li+1]-x_field[li]);
    c01 = field[li][lj][lk+1] * (x_field[li+1]-x) / (x_field[li+1]-x_field[li]) +
          field[li+1][lj][lk+1] * (x-x_field[li]) / (x_field[li+1]-x_field[li]);
    c11 = field[li][lj+1][lk+1] * (x_field[li+1]-x) / (x_field[li+1]-x_field[li]) +
          field[li+1][lj+1][lk+1] * (x-x_field[li]) / (x_field[li+1]-x_field[li]);

    c0 = c00 * (y_field[lj+1] - y) / (y_field[lj+1] - y_field[lj]) +
         c10 * (y - y_field[lj]) / (y_field[lj+1] - y_field[lj]);
    c1 = c01 * (y_field[lj+1] - y) / (y_field[lj+1] - y_field[lj]) +
         c11 * (y - y_field[lj]) / (y_field[lj+1] - y_field[lj]);

    c = c0 * (z_field[lk+1] - z) / (z_field[lk+1] - z_field[lk]) +
        c1 * (z - z_field[lk]) / (z_field[lk+1] - z_field[lk]);

    return c;
}

double TrilinearInterp_test(double (*f)(double,double,double), int li, int lj, int lk,
        double* x_field, double* y_field, double* z_field,
        double x, double y, double z)
{
    double c00, c10, c01, c11, c0, c1, c;
    double field[2][2][2];

    field[li][lj][lk] = (*f)(x_field[li],y_field[lj],z_field[lk]);
    field[li+1][lj][lk] = (*f)(x_field[li+1],y_field[lj],z_field[lk]);
    field[li][lj+1][lk] = (*f)(x_field[li],y_field[lj+1],z_field[lk]);
    field[li][lj][lk+1] = (*f)(x_field[li],y_field[lj],z_field[lk+1]);
    field[li+1][lj+1][lk] = (*f)(x_field[li+1],y_field[lj+1],z_field[lk]);
    field[li][lj+1][lk+1] = (*f)(x_field[li],y_field[lj+1],z_field[lk+1]);
    field[li+1][lj][lk+1] = (*f)(x_field[li+1],y_field[lj],z_field[lk+1]);
    field[li+1][lj+1][lk+1] = (*f)(x_field[li+1],y_field[lj+1],z_field[lk+1]);

    c00 = field[li][lj][lk] * (x_field[li+1]-x) / (x_field[li+1]-x_field[li]) +
          field[li+1][lj][lk] * (x-x_field[li]) / (x_field[li+1]-x_field[li]);
    c10 = field[li][lj+1][lk] * (x_field[li+1]-x) / (x_field[li+1]-x_field[li]) +
          field[li+1][lj+1][lk] * (x-x_field[li]) / (x_field[li+1]-x_field[li]);
    c01 = field[li][lj][lk+1] * (x_field[li+1]-x) / (x_field[li+1]-x_field[li]) +
          field[li+1][lj][lk+1] * (x-x_field[li]) / (x_field[li+1]-x_field[li]);
    c11 = field[li][lj+1][lk+1] * (x_field[li+1]-x) / (x_field[li+1]-x_field[li]) +
          field[li+1][lj+1][lk+1] * (x-x_field[li]) / (x_field[li+1]-x_field[li]);

    c0 = c00 * (y_field[lj+1] - y) / (y_field[lj+1] - y_field[lj]) +
         c10 * (y - y_field[lj]) / (y_field[lj+1] - y_field[lj]);
    c1 = c01 * (y_field[lj+1] - y) / (y_field[lj+1] - y_field[lj]) +
         c11 * (y - y_field[lj]) / (y_field[lj+1] - y_field[lj]);

    c = c0 * (z_field[lk+1] - z) / (z_field[lk+1] - z_field[lk]) +
        c1 * (z - z_field[lk]) / (z_field[lk+1] - z_field[lk]);

    c = fabs(c-(*f)(x,y,z));
    return c;
}

double ExactSolution(double x, double y, double z)
{
    double value;
    value = 1.0 + 0.5*(x+y+z) + 0.25*(x*y+y*z+z*x) + 0.125*x*y*z;
    return value;
}

double SimpleInterp_MAC(double*** field, int li, int lj, int lk, int dirc)
{
    double value;

    if (dirc == 0)
        value = (field[li+1][lj][lk] + field[li+1][lj][lk+1] +
                 field[li+1][lj+1][lk+1] + field[li+1][lj+1][lk])/4.0;
    else if (dirc == 1)
        value = 0.0;
    else if (dirc == 2)
        value = (field[li][lj][lk+1] + field[li+1][lj][lk+1] +
                 field[li+1][lj+1][lk+1] + field[li][lj+1][lk+1])/4.0;
    else assert(false);

    return value;
}
