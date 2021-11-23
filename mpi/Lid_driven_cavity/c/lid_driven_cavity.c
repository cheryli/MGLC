/* this is a lattice boltzmann method code to solve a 2D lid driven cavity problem
   D2Q9 with multi-relaxation time and BGK model
   author: Btli
*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define NX 200
#define NY 200

const int SRT = 1, MRT = 2;
const int model = 2;


const double pi = 3.141592653589793;

// flow variables and distribution functions
double rho[NX][NY], u[NX][NY], v[NX][NY];
double f[NX][NY][9], f_post[NX][NY][9];
double up[NX][NY], vp[NX][NY];

// D2Q9 model, consts
const double ex[9] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
const double ey[9] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
const double omega[9] = {4.0/9.0, \
                        1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, \
                        1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
double tau, s_e, s_epsilon, s_q, s_nu;


// study parameters
const double Re = 1000;  // Reynolds number
const double height = NX;     // length and height of cavity
const double u_zero = 0.1;  // upper velocity
const double rho_zero = 1.0;
double nu;


int initial();
int collision();
int streaming();
int boundary();
int macro();
double check(int itc);
int output_tecplot(int itc);
int output_binary();

int main()
{
    const double eps = 1e-6;
    const int itc_max = 5e7;

    double error = 1.0;
    int itc = 0;
    
    initial();
    output_tecplot(itc);
    while (error >= eps && itc++ <= itc_max)
    {
        collision();
        streaming();
        boundary();
        macro();
        if (itc % 2000 == 0)
        {
            error = check(itc);
        }
        if(itc % 20000 == 0)
        {
            output_tecplot(itc);
        }
    }
    
    output_tecplot(itc);
    output_binary();

    return 0;
}

int initial()
{
    double u2, ue[9];

    // // initial relaxation parameters
    // nu = u_zero * height / Re;
    // tau = 3.0 * nu + 0.5;
    // s_nu = 1.0 / tau;
    // s_e = s_nu;     // 1.05
    // s_epsilon = s_nu;       // 1.1
    // s_q = 8.0 * (2.0 * tau - 1) / (8.0 * tau - 1);      // 1.25


    // initial relaxation parameters
    nu = u_zero * height / Re;
    tau = 3.0 * nu + 0.5;
    s_nu = 1.0 / tau;
    s_e = s_nu;     // 1.05
    s_epsilon = s_nu;       // 1.1
    s_q = 8.0 * (2.0 * tau - 1) / (8.0 * tau - 1);      // 1.25

    printf("Lid driven cavity flow.\n");
    printf("Flow parameters: \n");
    printf("Re = %f\nHeight = Length = %f\n", Re, height);
    printf("Upper velocity = %f\n", u_zero);
    printf("Flow viscosity = %f\n", nu);
    if(model == SRT)
    {
        printf("We are using SRT/BGK model...\n");
        printf("tau = %f\n", tau);
    } 

    if(model == MRT) 
    {
        printf("We are using MRT model...\n");
        printf("s_e = %f,    s_epsilon = %f\ns_q = %f,    s_nu = %f\n", s_e, s_epsilon, s_q, s_nu);
    }
    printf("------------Start--------------\n");
    printf("itc      error\n");


    for (int i = 0; i < NX; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            u[i][j] = 0.0;
            v[i][j] = 0.0;
            rho[i][j] =  rho_zero;
        }
    }

    for (int i = 0; i < NX; i++)
    {
        u[i][NY-1] = u_zero;
    }

    for (int i = 0; i < NX; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            u2 = u[i][j] * u[i][j]  + v[i][j] * v[i][j];
            for (int alpha = 0; alpha < 9; alpha++)
            {
                ue[alpha] = u[i][j] * ex[alpha] + v[i][j] * ey[alpha];
                f[i][j][alpha] = rho[i][j] * omega[alpha] * (1.0 + 3.0 * ue[alpha] \
                        + 4.5 * ue[alpha] * ue[alpha] - 1.5 * u2);
            }
        }
    }

    return 0;
}

int collision()
{
    double u2, ue[9], feq[9];    
    double m[9], m_eq[9], m_post[9], s[9];

    if (model == SRT)
    {
        for (int i = 0; i < NX; i++)
        {
            for (int j = 0; j < NY; j++)
            {
                /* collision --------- SRT/BGK model --------------*/
                u2 = u[i][j] * u[i][j]  + v[i][j] * v[i][j];
                for (int alpha = 0; alpha < 9; alpha++)
                {
                    ue[alpha] = u[i][j] * ex[alpha] + v[i][j] * ey[alpha];
                    feq[alpha] = rho[i][j] * omega[alpha] * (1.0 + 3.0 * ue[alpha] \
                            + 4.5 * ue[alpha] * ue[alpha] - 1.5 * u2);
                    f_post[i][j][alpha] = f[i][j][alpha] - 1.0 / tau * (f[i][j][alpha] - feq[alpha]);
                }  
            }
        }
    }else if (model == MRT)
    {
        for (int i = 0; i < NX; i++)
        {
            for (int j = 0; j < NY; j++)
            { 
                /* collision --------- MRT model ------------------*/
                // compute the momentum
                m[0] = f[i][j][0] + f[i][j][1] + f[i][j][2] + f[i][j][3] + f[i][j][4] \
                    + f[i][j][5] + f[i][j][6] + f[i][j][7] + f[i][j][8];
                m[1] = -4*f[i][j][0] - f[i][j][1] - f[i][j][2] - f[i][j][3] - f[i][j][4] \
                    + 2*f[i][j][5] + 2*f[i][j][6] + 2*f[i][j][7] + 2*f[i][j][8];
                m[2] = 4*f[i][j][0] - 2*f[i][j][1] - 2*f[i][j][2] - 2*f[i][j][3] - 2*f[i][j][4] \
                    + f[i][j][5] + f[i][j][6] + f[i][j][7] + f[i][j][8];
                m[3] = f[i][j][1] - f[i][j][3] \
                    + f[i][j][5] - f[i][j][6] - f[i][j][7] + f[i][j][8];
                m[4] = -2*f[i][j][1] + 2*f[i][j][3] \
                    + f[i][j][5] - f[i][j][6] - f[i][j][7] + f[i][j][8];
                m[5] = f[i][j][2] - f[i][j][4] \
                    + f[i][j][5] + f[i][j][6] - f[i][j][7] - f[i][j][8];
                m[6] = -2*f[i][j][2] + 2*f[i][j][4] \
                    + f[i][j][5] + f[i][j][6] - f[i][j][7] - f[i][j][8];
                m[7] = f[i][j][1] - f[i][j][2] + f[i][j][3] - f[i][j][4];
                m[8] = f[i][j][5] - f[i][j][6] + f[i][j][7] - f[i][j][8];

                // compute the equilibrium momentum
                m_eq[0] = rho[i][j];
                m_eq[1] = rho[i][j] * (-2.0 + 3.0 * (u[i][j] * u[i][j] + v[i][j] * v[i][j]));
                m_eq[2] = rho[i][j] * (1.0 - 3.0 * (u[i][j] * u[i][j] + v[i][j] * v[i][j]));
                m_eq[3] = rho[i][j] * u[i][j];
                m_eq[4] = -1.0 * rho[i][j] * u[i][j];
                m_eq[5] = rho[i][j] * v[i][j];
                m_eq[6] = -1.0 * rho[i][j] * v[i][j];
                m_eq[7] = rho[i][j] * (u[i][j] * u[i][j] - v[i][j] * v[i][j]);
                m_eq[8] = u[i][j] * v[i][j];

                // // set the relaxation matrix
                // s[0] = 0.0;     //s_{\rho}  
                // s[1] = s_e;     //s_{e}
                // s[2] = s_epsilon;    //s_{\epsilon}
                // s[3] = 0.0;     //s_{j}
                // s[4] = s_q;     //s_{q}
                // s[5] = 0.0;     //s_{j}
                // s[6] = s_q;     //s_{q}
                // s[7] = s_nu;    //s_{\nu}
                // s[8] = s_nu;    //s_{\nu}

                // set the relaxation matrix
                s[0] = 0.0;     //s_{\rho}  
                s[1] = s_nu;     //s_{e}
                s[2] = s_nu;    //s_{\epsilon}
                s[3] = 0.0;     //s_{j}
                s[4] = s_q;     //s_{q}
                s[5] = 0.0;     //s_{j}
                s[6] = s_q;     //s_{q}
                s[7] = s_nu;    //s_{\nu}
                s[8] = s_nu;    //s_{\nu}

                // collision
                for (int alpha = 0; alpha < 9; alpha++)
                {
                    m_post[alpha] = m[alpha] - s[alpha] * (m[alpha] - m_eq[alpha]);
                }

                // compute the f
                f_post[i][j][0] = (m_post[0] - m_post[1] + m_post[2]) / 9.0;
                f_post[i][j][1] = (4.0*m_post[0] - m_post[1] - 2.0*m_post[2] + 6.0*m_post[3] - 6.0*m_post[4] + 9.0*m_post[7]) / 36.0;
                f_post[i][j][2] = (4.0*m_post[0] - m_post[1] - 2.0*m_post[2] + 6.0*m_post[5] - 6.0*m_post[6] - 9.0*m_post[7]) / 36.0;
                f_post[i][j][3] = (4.0*m_post[0] - m_post[1] - 2.0*m_post[2] - 6.0*m_post[3] + 6.0*m_post[4] + 9.0*m_post[7]) / 36.0;
                f_post[i][j][4] = (4.0*m_post[0] - m_post[1] - 2.0*m_post[2] - 6.0*m_post[5] + 6.0*m_post[6] - 9.0*m_post[7]) / 36.0;
                f_post[i][j][5] = (4.0*m_post[0] + 2.0*m_post[1] + m_post[2] + 6.0*m_post[3] + 3.0*m_post[4] \
                        + 6.0*m_post[5] + 3.0*m_post[6] + 9.0*m_post[8]) / 36.0;
                f_post[i][j][6] = (4.0*m_post[0] + 2.0*m_post[1] + m_post[2] - 6.0*m_post[3] - 3.0*m_post[4] \
                        + 6.0*m_post[5] + 3.0*m_post[6] - 9.0*m_post[8]) / 36.0;
                f_post[i][j][7] = (4.0*m_post[0] + 2.0*m_post[1] + m_post[2] - 6.0*m_post[3] - 3.0*m_post[4] \
                        - 6.0*m_post[5] - 3.0*m_post[6] + 9.0*m_post[8]) / 36.0;
                f_post[i][j][8] = (4.0*m_post[0] + 2.0*m_post[1] + m_post[2] + 6.0*m_post[3] + 3.0*m_post[4] \
                        - 6.0*m_post[5] - 3.0*m_post[6] - 9.0*m_post[8]) / 36.0;
            }
        }
    }

    return 0;
}

int streaming()
{
    int ip, jp;

    for (int i = 0; i < NX; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            for (int alpha = 0; alpha < 9; alpha++)
            {
                ip = i + ex[alpha];
                jp = j + ey[alpha];
                // period
                if(ip == NX) ip = 0;
                if(ip == -1) ip = NX-1;
                if(jp == NY) jp = 0;
                if(jp == -1) jp = NY-1;

                f[ip][jp][alpha] = f_post[i][j][alpha];
            }
        }   
    }
    return 0;
}

int boundary()
{
    // half way bounce back
    for (int j = 0; j < NY; j++)
    {
        // left surface
        f[0][j][1] = f_post[0][j][3];
        f[0][j][5] = f_post[0][j][7];
        f[0][j][8] = f_post[0][j][6];

        // right surface
        f[NX-1][j][3] = f_post[NX-1][j][1];
        f[NX-1][j][6] = f_post[NX-1][j][8];
        f[NX-1][j][7] = f_post[NX-1][j][5];
    }

    for (int i = 0; i < NX; i++)
    {
        // lower surface
        f[i][0][2] = f_post[i][0][4];
        f[i][0][6] = f_post[i][0][8];
        f[i][0][5] = f_post[i][0][7];

        // upper surface
        f[i][NY-1][4] = f_post[i][NY-1][2];
        f[i][NY-1][8] = f_post[i][NY-1][6] - rho[i][NY-1] * (-1.0 * u_zero) / 6.0;
        f[i][NY-1][7] = f_post[i][NY-1][5] - rho[i][NY-1] * (1.0 * u_zero) / 6.0;;
    }


    
    return 0;
}

int macro()
{
    for (int i = 0; i < NX; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            u[i][j] = 0;
            v[i][j] = 0;
            rho[i][j] = 0;
            for (int alpha = 0; alpha < 9; alpha++)
            {
                rho[i][j] += f[i][j][alpha];
                u[i][j] += f[i][j][alpha] * ex[alpha];
                v[i][j] += f[i][j][alpha] * ey[alpha];
            }
            u[i][j] = u[i][j] / rho[i][j];
            v[i][j] = v[i][j] / rho[i][j];
        }       
    }
    
    return 0;
}

double check(int itc)
{
    FILE *fp;
    double error, error_abs = 0, sum = 0;
    for (int i = 0; i < NX; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            error_abs += pow(u[i][j] - up[i][j], 2) + pow(v[i][j] - vp[i][j], 2);
            sum += pow(u[i][j], 2) + pow(v[i][j], 2);

            up[i][j] = u[i][j];
            vp[i][j] = v[i][j];
        }      
    }  
    error = pow(error_abs / sum, 0.5);
    error = pow(error_abs, 0.5) / pow(sum, 0.5);
    printf("%d, %.12f\n", itc, error);

    fp = fopen("log.dat", "a+");
    fprintf(fp, "%d  %.12f\n", itc, error);
    fclose(fp);    

    return error;
}

int output_tecplot(int itc)
{
    char s[40];
    FILE *fp = NULL, *fp2 = NULL;
    double x = 0.0, delta_x = height / (NX-1);
    double y = 0.0, delta_y = height / (NY-1);
    double vot, vot_max = 0.0;
    int x_max, y_max;

    sprintf(s, "Result\\lid_driven__%d.dat", itc);
    fp = fopen(s, "w+");

    fprintf(fp, "TITLE = \"Lid driven cavity.\"\n");
    fprintf(fp, "VARIABLES = \"X\", \"Y\", \"P\", \"U\", \"V\", \"Vorticity\"\n");
    fprintf(fp, "ZONE\n");
    fprintf(fp, "I = %d, J = %d\n", NX, NY);

    for (int i = 0; i < NX; i++)
    {
        y = 0.0;
        for (int j = 0; j < NY; j++)
        {
            if(i <= 1 || i >= NX-2 || j <= 1 || j >= NY-2)
            {
                vot = 0;
            }
            else
            {
                vot = (v[i+1][j] - v[i-1][j]) / (3.0 * delta_x) \
                    + (v[i+1][j+1] - v[i-1][j+1]) / (12.0 * delta_x) \
                    + (v[i+1][j-1] - v[i-1][j-1]) / (12.0 * delta_x) \
                    - (u[i][j+1] - u[i][j-1]) / (3.0 * delta_y) \
                    - (u[i+1][j+1] - u[i+1][j-1]) / (12.0 * delta_y) \
                    - (u[i-1][j+1] - u[i-1][j-1]) / (12.0 * delta_y); 
            }
            if(fabs(vot) > fabs(vot_max))
            {
                vot_max = vot;
                x_max = x;
                y_max = y;
            }

            fprintf(fp, "%.12f  %.12f  %.12f  %.12f  %.12f  %.12f\n", x, y, rho[i][j] / 3.0, u[i][j], v[i][j], vot);
            y += delta_y;
        }
        x += delta_x;
    }
    
    fclose(fp);

    fp2 = fopen("vorticity core.dat", "a+");
    fprintf(fp2, "%f  %f  %.12f\n", x_max, y_max, vot_max);
    fclose(fp2);

    return 0;
}


int output_binary()
{
    FILE *fp;
    double *x, *y, delta_x = height / (NX-1), delta_y = height / (NY-1);

    x = (double *)malloc(NX * NY * sizeof(double));
    y = (double *)malloc(NX * NY * sizeof(double));

    // mesh generation
    for (int i = 0; i < NX; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            x[i * NX + j] = i * delta_x;
            y[i * NX + j] = j * delta_y;
        }
    }

    fp = fopen("flow_binary", "wb+");
    fwrite(x, sizeof(double), NX*NY, fp);
    fwrite(y, sizeof(double), NX*NY, fp);
    fwrite(rho, sizeof(double), NX*NY, fp);
    fwrite(u, sizeof(double), NX*NY, fp);
    fwrite(v, sizeof(double), NX*NY, fp);
    fclose(fp);
    free(x);
    free(y);

    return 0;
}
