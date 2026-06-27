#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* -------------------- 常量与宏 -------------------- */
#define PI 3.14159265358979323846
#define DEG2RAD (PI / 180.0)
#define RAD2DEG (180.0 / PI)
#define EPS 1e-12

/* ---------- 辅助内存分配（与参考代码一致） ---------- */
static double* double_alloc(int n) {
    return (double*)malloc(n * sizeof(double));
}
static double* double_calloc(int n) {
    return (double*)calloc(n, sizeof(double));
}
static double** double_alloc2(int rows, int cols) {
    double **m = (double**)malloc(rows * sizeof(double*));
    for (int i = 0; i < rows; i++)
        m[i] = (double*)malloc(cols * sizeof(double));
    return m;
}
static void free_double2(double **m, int rows) {
    for (int i = 0; i < rows; i++) free(m[i]);
    free(m);
}

/* ---------- 雅可比特征分解（参考代码） ---------- */
#define ROTATE(a,i,j,k,l) do { \
    double g = a[i][j], h = a[k][l]; \
    a[i][j] = g - s*(h + g*tau); \
    a[k][l] = h + s*(g - h*tau); \
} while(0)

void jacobi(double **a, int n, int np, double *d, double **v, int *nrot) {
    int NMAX = 500, i, ip, iq, j;
    double c, g, h, s, sm, t, tau, theta, tresh;
    double *b = double_alloc(NMAX);
    double *z = double_calloc(NMAX);

    for (ip = 0; ip < n; ip++) {
        for (iq = 0; iq < n; iq++) v[ip][iq] = 0.0;
        v[ip][ip] = 1.0;
    }
    for (ip = 0; ip < n; ip++) {
        b[ip] = d[ip] = a[ip][ip];
        z[ip] = 0.0;
    }
    *nrot = 0;
    for (i = 1; i <= 50; i++) {
        sm = 0.0;
        for (ip = 0; ip < n-1; ip++)
            for (iq = ip+1; iq < n; iq++)
                sm += fabs(a[ip][iq]);
        if (sm == 0.0) {
            free(z); free(b); return;
        }
        tresh = (i < 4) ? 0.2 * sm / (n*n) : 0.0;
        for (ip = 0; ip < n-1; ip++) {
            for (iq = ip+1; iq < n; iq++) {
                g = 100.0 * fabs(a[ip][iq]);
                if (i > 4 && (fabs(d[ip])+g == fabs(d[ip])) &&
                    (fabs(d[iq])+g == fabs(d[iq]))) {
                    a[ip][iq] = 0.0;
                } else if (fabs(a[ip][iq]) > tresh) {
                    h = d[iq] - d[ip];
                    if (fabs(h)+g == fabs(h))
                        t = a[ip][iq] / h;
                    else {
                        theta = 0.5 * h / a[ip][iq];
                        t = 1.0 / (fabs(theta) + sqrt(1.0 + theta*theta));
                        if (theta < 0.0) t = -t;
                    }
                    c = 1.0 / sqrt(1.0 + t*t);
                    s = t * c;
                    tau = s / (1.0 + c);
                    h = t * a[ip][iq];
                    z[ip] -= h;
                    z[iq] += h;
                    d[ip] -= h;
                    d[iq] += h;
                    a[ip][iq] = 0.0;
                    for (j = 0; j < ip; j++) ROTATE(a, j, ip, j, iq);
                    for (j = ip+1; j < iq; j++) ROTATE(a, ip, j, j, iq);
                    for (j = iq+1; j < n; j++) ROTATE(a, ip, j, iq, j);
                    for (j = 0; j < n; j++) ROTATE(v, j, ip, j, iq);
                    ++(*nrot);
                }
            }
        }
        for (ip = 0; ip < n; ip++) {
            b[ip] += z[ip];
            d[ip] = b[ip];
            z[ip] = 0.0;
        }
    }
    fprintf(stderr, "too many iterations in jacobi\n");
    free(z); free(b);
}

void eigsrt(double *d, double **v, int n) {
    int i, j, k;
    double p;
    for (i = 0; i < n-1; i++) {
        p = d[k = i];
        for (j = i+1; j < n; j++)
            if (d[j] >= p) p = d[k = j];
        if (k != i) {
            d[k] = d[i];
            d[i] = p;
            for (j = 0; j < n; j++) {
                p = v[j][i];
                v[j][i] = v[j][k];
                v[j][k] = p;
            }
        }
    }
}

/* ---------- 将 vm[6] 转换为矩阵（北、东、下） ---------- */
static void set_mt(double *vm, double **TM) {
    // vm 顺序: Mrr, Mtt, Mpp, Mrt, Mrp, Mtp
    TM[0][0] =  vm[1];   // Mxx = Mtt
    TM[1][1] =  vm[2];   // Myy = Mpp
    TM[2][2] =  vm[0];   // Mzz = Mrr
    TM[0][1] = -vm[5];   // Mxy = -Mtp
    TM[0][2] =  vm[3];   // Mxz = Mrt
    TM[1][2] = -vm[4];   // Myz = -Mrp
    TM[1][0] = TM[0][1];
    TM[2][0] = TM[0][2];
    TM[2][1] = TM[1][2];
}

/* ---------- 计算两个节面参数（返回度数） ---------- */
void moment_tensor_to_planes(double *vm,
                             double *s1, double *d1, double *r1,
                             double *s2, double *d2, double *r2) {
    double **TM = double_alloc2(3,3);
    double **evec3 = double_alloc2(3,3);
    double *eval3 = double_alloc(3);
    double *vn = double_alloc(3);
    double *vs = double_alloc(3);
    int nrot;

    set_mt(vm, TM);

    jacobi(TM, 3, 3, eval3, evec3, &nrot);
    eigsrt(eval3, evec3, 3);  // 降序：T轴最大，P轴最小

    // 确保特征向量向上（与参考代码一致）
    if (evec3[2][0] < 0.0)
        for (int i=0; i<3; i++) evec3[i][0] = -evec3[i][0];   // P轴
    if (evec3[2][2] < 0.0)
        for (int i=0; i<3; i++) evec3[i][2] = -evec3[i][2];   // T轴

    // 叉积计算中间轴（null轴） evec3[:,1] = evec3[:,0] × evec3[:,2]
    evec3[0][1] = evec3[1][0]*evec3[2][2] - evec3[1][2]*evec3[2][0];
    evec3[1][1] = evec3[2][0]*evec3[0][2] - evec3[2][2]*evec3[0][0];
    evec3[2][1] = evec3[0][0]*evec3[1][2] - evec3[0][2]*evec3[1][0];

    // ---- 第一个节面 ----
    for (int i=0; i<3; i++) {
        vn[i] = (evec3[i][0] + evec3[i][2]) / sqrt(2.0);
        vs[i] = (evec3[i][0] - evec3[i][2]) / sqrt(2.0);
    }
    if (vn[2] > 0.0) {   // 确保法向量向下（z负）
        for (int i=0; i<3; i++) { vn[i] = -vn[i]; vs[i] = -vs[i]; }
    }
    double strike1 = atan2(-vn[0], vn[1]);
    double dip1    = acos(-vn[2]);
    double si = sin(strike1), co = cos(strike1);
    double rake1   = atan2((vs[0]*si - vs[1]*co), -(vs[0]*co + vs[1]*si)*vn[2]);

    // ---- 第二个节面 ----
    for (int i=0; i<3; i++) {
        vn[i] = (evec3[i][0] - evec3[i][2]) / sqrt(2.0);
        vs[i] = (evec3[i][0] + evec3[i][2]) / sqrt(2.0);
    }
    if (vn[2] > 0.0) {
        for (int i=0; i<3; i++) { vn[i] = -vn[i]; vs[i] = -vs[i]; }
    }
    double strike2 = atan2(-vn[0], vn[1]);
    double dip2    = acos(-vn[2]);
    si = sin(strike2); co = cos(strike2);
    double rake2   = atan2((vs[0]*si - vs[1]*co), -(vs[0]*co + vs[1]*si)*vn[2]);

    // 转换为度并规范化
    *s1 = RAD2DEG * strike1;
    if (*s1 < 0.0) *s1 += 360.0;
    *d1 = RAD2DEG * dip1;
    *r1 = RAD2DEG * rake1;

    *s2 = RAD2DEG * strike2;
    if (*s2 < 0.0) *s2 += 360.0;
    *d2 = RAD2DEG * dip2;
    *r2 = RAD2DEG * rake2;

    free(eval3); free(vn); free(vs);
    free_double2(TM, 3);
    free_double2(evec3, 3);
}

/* ---------- 以下为绘图相关函数 ---------- */

/* 生成节面曲线点集（等面积投影），输入角度制 */
static void pnodal(double strike, double dip, double rake,
                   double **px, double **py, int *npts) {
    (void)rake; // 消除未使用警告
    strike *= DEG2RAD;
    dip    *= DEG2RAD;
    int n = 180;
    double *xarr = (double*)malloc(n * sizeof(double));
    double *yarr = (double*)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        double ii = (double)i * PI / 180.0;
        double x = cos(strike)*cos(ii) - sin(strike)*sin(ii)*cos(dip);
        double y = sin(strike)*cos(ii) + cos(strike)*sin(ii)*cos(dip);
        double z = sin(ii)*sin(dip);
        double azimuth = atan2(y, x);
        double ain = atan2(sqrt(1.0 - z*z), z);
        double r = sqrt(2.0) * sin(ain/2.0);
        xarr[i] = r * sin(azimuth);
        yarr[i] = r * cos(azimuth);
    }
    *px = xarr;
    *py = yarr;
    *npts = n;
}

/* 计算辐射花样值 r(θ,φ) = r^T·M·r (M为3x3矩阵) */
static double get_mt_value(double x, double y, double mt[3][3]) {
    double rad = sqrt(x*x + y*y);
    if (rad > 1.0) return 0.0;
    double azi = atan2(x, y);
    double ain = 2.0 * asin(rad / sqrt(2.0));
    double r[3] = { sin(ain)*cos(azi), sin(ain)*sin(azi), cos(ain) };
    double val = 0.0;
    for (int i=0;i<3;i++)
        for (int j=0;j<3;j++)
            val += r[i] * mt[i][j] * r[j];
    return val;
}

/* Marching Squares 提取零值线 (返回线段列表) */
typedef struct { double x1,y1,x2,y2; } Segment;
static Segment* marching_squares(double mt[3][3], int res, int *nseg) {
    double *x_range = (double*)malloc(res * sizeof(double));
    double *y_range = (double*)malloc(res * sizeof(double));
    for (int i=0;i<res;i++) {
        x_range[i] = -1.0 + 2.0*i/(res-1);
        y_range[i] = -1.0 + 2.0*i/(res-1);
    }
    double **grid = (double**)malloc(res * sizeof(double*));
    for (int i=0;i<res;i++) grid[i] = (double*)malloc(res * sizeof(double));
    for (int i=0;i<res;i++)
        for (int j=0;j<res;j++)
            grid[i][j] = get_mt_value(x_range[i], y_range[j], mt);

    Segment *segs = NULL;
    int cap = 0, count = 0;
    for (int i=0;i<res-1;i++) {
        for (int j=0;j<res-1;j++) {
            double v00 = grid[i][j], v10 = grid[i+1][j];
            double v01 = grid[i][j+1], v11 = grid[i+1][j+1];
            double x0=x_range[i], x1=x_range[i+1], y0=y_range[j], y1=y_range[j+1];
            if (x0*x0+y0*y0 > 1.05) continue;
            double edges[4][2]; int ec=0;
            if (v00*v10 < 0) { double t=-v00/(v10-v00); edges[ec][0]=x0+t*(x1-x0); edges[ec][1]=y0; ec++; }
            if (v10*v11 < 0) { double t=-v10/(v11-v10); edges[ec][0]=x1; edges[ec][1]=y0+t*(y1-y0); ec++; }
            if (v01*v11 < 0) { double t=-v01/(v11-v01); edges[ec][0]=x0+t*(x1-x0); edges[ec][1]=y1; ec++; }
            if (v00*v01 < 0) { double t=-v00/(v01-v00); edges[ec][0]=x0; edges[ec][1]=y0+t*(y1-y0); ec++; }
            if (ec>=2) {
                if (count+1 > cap) {
                    cap = cap?cap*2:100;
                    segs = (Segment*)realloc(segs, cap * sizeof(Segment));
                }
                segs[count].x1 = edges[0][0]; segs[count].y1 = edges[0][1];
                segs[count].x2 = edges[1][0]; segs[count].y2 = edges[1][1];
                count++;
            }
        }
    }
    for (int i=0;i<res;i++) free(grid[i]);
    free(grid); free(x_range); free(y_range);
    *nseg = count;
    return segs;
}

/* 连接线段形成连续轮廓，返回三维数组：lines[nlines][2][len] */
static double*** connect_lines(Segment *segs, int nseg, int *nlines, int **line_lens) {
    if (nseg==0) { *nlines=0; *line_lens=NULL; return NULL; }
    int *used = (int*)calloc(nseg, sizeof(int));
    double ***lines = NULL;
    int line_count = 0, line_cap = 0;
    int *lens = NULL;
    for (int i=0;i<nseg;i++) {
        if (used[i]) continue;
        double *xs = (double*)malloc(2*nseg*sizeof(double));
        double *ys = (double*)malloc(2*nseg*sizeof(double));
        xs[0]=segs[i].x1; ys[0]=segs[i].y1;
        xs[1]=segs[i].x2; ys[1]=segs[i].y2;
        int len = 2;
        used[i] = 1;
        int changed = 1;
        while (changed) {
            changed = 0;
            for (int j=0;j<nseg;j++) {
                if (used[j]) continue;
                double ends[4][2] = {{segs[j].x1,segs[j].y1},{segs[j].x2,segs[j].y2},
                                     {segs[j].x2,segs[j].y2},{segs[j].x1,segs[j].y1}};
                for (int k=0;k<2;k++) {
                    if (fabs(xs[0]-ends[k][0])<1e-8 && fabs(ys[0]-ends[k][1])<1e-8) {
                        double nx = ends[1-k][0], ny = ends[1-k][1];
                        for (int m=len; m>0; m--) { xs[m]=xs[m-1]; ys[m]=ys[m-1]; }
                        xs[0]=nx; ys[0]=ny; len++;
                        used[j]=1; changed=1; break;
                    }
                    if (fabs(xs[len-1]-ends[k][0])<1e-8 && fabs(ys[len-1]-ends[k][1])<1e-8) {
                        double nx = ends[1-k][0], ny = ends[1-k][1];
                        xs[len]=nx; ys[len]=ny; len++;
                        used[j]=1; changed=1; break;
                    }
                }
                if (changed) break;
            }
        }
        if (line_count+1 > line_cap) {
            line_cap = line_cap?line_cap*2:10;
            lines = (double***)realloc(lines, line_cap * sizeof(double**));
            lens = (int*)realloc(lens, line_cap * sizeof(int));
        }
        double **line = (double**)malloc(2 * sizeof(double*));
        line[0] = xs; line[1] = ys;
        lines[line_count] = line;
        lens[line_count] = len;
        line_count++;
    }
    free(used);
    *nlines = line_count;
    *line_lens = lens;
    return lines;
}

/* 绘制一条节面曲线 (输入角度制) */
static void draw_curve(FILE *ps, double strike, double dip, double rake,
                       double r, double g, double b, double lw) {
    double *px, *py; int n;
    pnodal(strike, dip, rake, &px, &py, &n);
    fprintf(ps, "%g %g %g setrgbcolor\n", r, g, b);
    fprintf(ps, "%g setlinewidth\n", lw);
    fprintf(ps, "newpath\n");
    for (int i=0;i<n;i++) {
        double xp = 72.0 + (px[i]+1.0)/2.0 * 432.0; // margin=72, plot_size=432
        double yp = 72.0 + (py[i]+1.0)/2.0 * 432.0;
        if (i==0) fprintf(ps, "%g %g moveto\n", xp, yp);
        else fprintf(ps, "%g %g lineto\n", xp, yp);
    }
    fprintf(ps, "stroke\n");
    fprintf(ps, "0 setlinewidth\n");
    fprintf(ps, "0 0 0 setrgbcolor\n");
    free(px); free(py);
}

/* ---------- 主绘图函数 ---------- */
void plot_beachball(double mt[3][3], FILE *ps) {
    // 1. 从 mt 提取 vm[6] (Mrr, Mtt, Mpp, Mrt, Mrp, Mtp)
    double vm[6];
    vm[0] = mt[2][2];    // Mrr = Mzz
    vm[1] = mt[0][0];    // Mtt = Mxx
    vm[2] = mt[1][1];    // Mpp = Myy
    vm[3] = mt[0][2];    // Mrt = Mxz
    vm[4] = -mt[1][2];   // Mrp = -Myz
    vm[5] = -mt[0][1];   // Mtp = -Mxy

    // 2. 计算两个节面参数 (角度制)
    double s1, d1, r1, s2, d2, r2;
    moment_tensor_to_planes(vm, &s1, &d1, &r1, &s2, &d2, &r2);

    // 3. 计算辐射花样点阵 (用于绘制正圆点)
    double *points_x = NULL, *points_y = NULL, *values = NULL;
    int npts = 0;
    double aP = mt[0][0], aT = mt[0][0];
    double jxP=0, jyP=0, jxT=0, jyT=0;
    for (int i=-20; i<20; i++) {
        double x = i/20.0;
        double ymax = sqrt(1.0 - x*x);
        int m = (int)(20.0 * ymax);
        for (int j=-m; j<=m; j++) {
            double y = j/20.0;
            double val = get_mt_value(x, y, mt);
            if (fabs(val) < 1e-15) continue;
            points_x = (double*)realloc(points_x, (npts+1)*sizeof(double));
            points_y = (double*)realloc(points_y, (npts+1)*sizeof(double));
            values = (double*)realloc(values, (npts+1)*sizeof(double));
            points_x[npts] = x; points_y[npts] = y; values[npts] = val; npts++;
            if (val < aP) { aP = val; jxP=x; jyP=y; }
            if (val > aT) { aT = val; jxT=x; jyT=y; }
        }
    }

    double max_abs = 0.0;
    for (int i=0;i<npts;i++) if (fabs(values[i]) > max_abs) max_abs = fabs(values[i]);
    if (max_abs < 1e-15) max_abs = 1.0;
    double base_size = 70.0;

    // 4. 计算零值等高线
    int nseg;
    Segment *segs = marching_squares(mt, 300, &nseg);
    int nlines, *line_lens;
    double ***lines = connect_lines(segs, nseg, &nlines, &line_lens);
    free(segs);

    // 5. 生成 PostScript 文件
    const double margin = 72.0, plot_size = 432.0;
    double ps_cx = margin + plot_size/2.0;
    double ps_cy = margin + plot_size/2.0;
    double radius = plot_size/2.0;

    #define TO_PS_X(x) (margin + ((x)+1.0)/2.0 * plot_size)
    #define TO_PS_Y(y) (margin + ((y)+1.0)/2.0 * plot_size)

    fprintf(ps, "%%!PS-Adobe-3.0\n");
    fprintf(ps, "%%%%BoundingBox: 0 0 612 792\n");
    fprintf(ps, "%%%%Page: 1 1\n");
    fprintf(ps, "<</PageSize [612 792]>> setpagedevice\n");
    fprintf(ps, "/Helvetica findfont 12 scalefont setfont\n");
    fprintf(ps, "0 setlinewidth\n");
    fprintf(ps, "0 0 0 setrgbcolor\n");

    // 绘制正值圆点 (黑色填充，白边)
    for (int i=0;i<npts;i++) {
        if (values[i] > 0.0) {
            double s_norm = (values[i] / max_abs) * base_size;
            double rp = sqrt(s_norm / PI);
            if (rp < 0.1) continue;
            double xp = TO_PS_X(points_x[i]);
            double yp = TO_PS_Y(points_y[i]);
            fprintf(ps, "newpath %g %g %g 0 360 arc fill\n", xp, yp, rp);
            fprintf(ps, "gsave\n");
            fprintf(ps, "1 1 1 setrgbcolor\n");
            fprintf(ps, "0.5 setlinewidth\n");
            fprintf(ps, "newpath %g %g %g 0 360 arc stroke\n", xp, yp, rp);
            fprintf(ps, "grestore\n");
            fprintf(ps, "0 0 0 setrgbcolor\n");
        }
    }

    // P/T 标签 (带白色背景)
    void draw_label(double x, double y, const char *label) {
        double xp = TO_PS_X(x), yp = TO_PS_Y(y);
        fprintf(ps, "gsave\n");
        fprintf(ps, "/x_ps %g def /y_ps %g def\n", xp, yp);
        fprintf(ps, "(%s) dup stringwidth pop /tw exch def\n", label);
        fprintf(ps, "/th 12 def /pad 2 def\n");
        fprintf(ps, "1 1 1 setrgbcolor\n");
        fprintf(ps, "newpath x_ps tw 2 div sub pad sub y_ps th 2 div sub pad sub moveto\n");
        fprintf(ps, "x_ps tw 2 div add pad add y_ps th 2 div sub pad sub lineto\n");
        fprintf(ps, "x_ps tw 2 div add pad add y_ps th 2 div add pad add lineto\n");
        fprintf(ps, "x_ps tw 2 div sub pad sub y_ps th 2 div add pad add lineto closepath fill\n");
        fprintf(ps, "1 0 0 setrgbcolor\n");
        fprintf(ps, "(%s) x_ps tw 2 div sub y_ps th 2 div sub moveto show\n", label);
        fprintf(ps, "grestore\n");
        fprintf(ps, "0 0 0 setrgbcolor\n");
    }
    draw_label(jxP, jyP, "P");
    draw_label(jxT, jyT, "T");

    // 边界大圆
    fprintf(ps, "newpath %g %g %g 0 360 arc stroke\n", ps_cx, ps_cy, radius);

    // 绘制两条节面曲线（角度制）
    draw_curve(ps, s1, d1, r1, 0.5, 0.0, 0.5, 1.5);  // 紫色
    draw_curve(ps, s2, d2, r2, 0.0, 0.0, 1.0, 1.5);  // 蓝色

    // 零值等高线 (交替绿红)
    const double cols[2][3] = {{0,1,0},{1,0,0}};
    for (int i=0;i<nlines;i++) {
        double *xs = lines[i][0], *ys = lines[i][1];
        int len = line_lens[i];
        fprintf(ps, "%g %g %g setrgbcolor\n", cols[i%2][0], cols[i%2][1], cols[i%2][2]);
        fprintf(ps, "1.5 setlinewidth\n");
        fprintf(ps, "newpath\n");
        for (int j=0;j<len;j++) {
            double xp = TO_PS_X(xs[j]);
            double yp = TO_PS_Y(ys[j]);
            if (j==0) fprintf(ps, "%g %g moveto\n", xp, yp);
            else fprintf(ps, "%g %g lineto\n", xp, yp);
        }
        fprintf(ps, "stroke\n");
        fprintf(ps, "0 setlinewidth\n");
        fprintf(ps, "0 0 0 setrgbcolor\n");
        free(xs); free(ys); free(lines[i]);
    }
    free(lines); free(line_lens);

    // 标题 (显示第一条节面)
    fprintf(ps, "gsave\n");
    fprintf(ps, "0 0 0 setrgbcolor\n");
    fprintf(ps, "/Helvetica findfont 14 scalefont setfont\n");
    char title[100];
    sprintf(title, "strike=%.1f, dip=%.1f, rake=%.1f", s1, d1, r1);
    fprintf(ps, "(%s) dup stringwidth pop /tw exch def\n", title);
    fprintf(ps, "/cx %g def /cy %g def\n", ps_cx, margin+plot_size+10);
    fprintf(ps, "cx tw 2 div sub cy moveto\n");
    fprintf(ps, "(%s) show\n", title);
    fprintf(ps, "grestore\n");

    // 释放点阵内存
    free(points_x); free(points_y); free(values);
}

/* ---------- 测试主程序 ---------- */
int main() {
    // 示例矩张量（对应 Python 示例 cmt = [-0.003, 0.133, -0.130, 0.612, 0.084, 4.330]）
    // 注意：此处的 mt 矩阵按 "南(x)-东(y)-上(z)" 坐标系排列
    // Mxx = Mtt = 0.133, Myy = Mpp = -0.130, Mzz = Mrr = -0.003
    // Mxy = -Mtp = -4.330, Mxz = Mrt = 0.612, Myz = -Mrp = -0.084
    double mt[3][3] = {
        { 0.133,  -4.330,  0.612 },
        { -4.330, -0.130, -0.084 },
        { 0.612,  -0.084, -0.003 }
    };

    FILE *f = fopen("beachball.ps", "w");
    if (!f) {
        fprintf(stderr, "Cannot open output file\n");
        return 1;
    }
    plot_beachball(mt, f);
    fclose(f);
    printf("PostScript file 'beachball.ps' generated successfully.\n");
    return 0;
}