import math

import matplotlib.pyplot as plt
import numpy as np

# beachball plot of earthquake focal mechanism
# the beachball is circle of radius 1


### calculate moment tensor
# use formula from SRL(1989),60(2) paper "A student guide and reviewer of moment tensors"
def dislocationtomt(strike, dip, rake, m0):

    m = np.empty((3, 3))
    phi = (strike * np.pi) / 180.0
    twophi = phi * 2.0
    delta = (dip * np.pi) / 180.0
    twodelta = delta * 2.0
    lamda = (rake * np.pi) / 180.0

    mxx = -(np.sin(delta)*np.cos(lamda)*np.sin(twophi)
           +np.sin(twodelta)*np.sin(lamda)*np.sin(phi)*np.sin(phi))

    myy =  (np.sin(delta)*np.cos(lamda)*np.sin(twophi)
           -np.sin(twodelta)*np.sin(lamda)*np.cos(phi)*np.cos(phi))

    mzz =  (np.sin(twodelta)*np.sin(lamda))

    mxy =  (np.sin(delta)*np.cos(lamda)*np.cos(twophi)
           +0.5*np.sin(twodelta)*np.sin(lamda)*np.sin(twophi))

    mxz = -(np.cos(delta)*np.cos(lamda)*np.cos(phi)
           +np.cos(twodelta)*np.sin(lamda)*np.sin(phi))

    myz = -(np.cos(delta)*np.cos(lamda)*np.sin(phi)
           -np.cos(twodelta)*np.sin(lamda)*np.cos(phi))

    mrr = mzz * m0
    mrt = mxz * m0
    mrp = -myz * m0
    mtt = mxx * m0
    mtp = -mxy * m0
    mpp = myy * m0

    m[0][0] = mxx
    m[0][1] = mxy
    m[0][2] = mxz
    m[1][1] = myy
    m[1][2] = myz
    m[2][2] = mzz
    m[1][0] = m[0][1]
    m[2][0] = m[0][2]
    m[2][1] = m[1][2]

    return m
    # return mrr,mtt,mpp,mrt,mrp,mtp

    # using focal mechansim parameter strike dip rake to
    # get its stereographic equal area projction plane


def pnodal(strike, dip, rake):
    strike = (strike * np.pi) / 180.0
    dip = (dip * np.pi) / 180.0
    rake = (rake * np.pi) / 180.0
    px = []
    py = []
    # px.append(np.cos(strike))
    # py.append(np.sin(strike))
    for i in range(0, 180):
        ii = i * np.pi / 180
        x = np.cos(strike) * np.cos(ii) - np.sin(strike) * np.sin(ii) * np.cos(dip)
        y = np.sin(strike) * np.cos(ii) + np.cos(strike) * np.sin(ii) * np.cos(dip)
        z = np.sin(ii) * np.sin(dip)
        azimuth = np.atan2(y, x)
        ain = np.atan2(np.sqrt(1 - z**2), z)
        r = np.sqrt(2) * np.sin(ain / 2.0)
        px.append(r * np.sin(azimuth))
        py.append(r * np.cos(azimuth))

    return px, py


# GETAUX returns auxilary fault plane strike, dip & rake,
# given strike,dip,rake of main fault plane.
def getaux(strike1, dip1, rake1):

    degrad = 180.0 / 3.1415927
    s1 = strike1 / degrad
    d1 = dip1 / degrad
    r1 = rake1 / degrad

    d2 = np.arccos(np.sin(r1) * np.sin(d1))

    sr2 = np.cos(d1) / np.sin(d2)
    cr2 = -np.sin(d1) * np.cos(r1) / np.sin(d2)
    r2 = np.arctan2(sr2, cr2)

    s12 = np.cos(r1) / np.sin(d2)
    c12 = -1.0 / (np.tan(d1) * np.tan(d2))
    s2 = s1 - np.arctan2(s12, c12)

    strike2 = s2 * degrad
    dip2 = d2 * degrad
    rake2 = r2 * degrad

    if dip2 > 90.0:
        strike2 = strike2 + 180.0
        dip2 = 180.0 - dip2
        rake2 = 360.0 - rake2
    if strike2 > 360.0:
        strike2 = strike2 - 360.0

    return strike2, dip2, rake2


# ------- find contour =0 ---------------
# ...............................
# 核心计算函数：直接根据物理公式计算值，不需要插值查找
def get_mt_value_at_coord(x, y, mt):
    """直接根据坐标计算矩张量投影值，替代缓慢的插值过程"""
    rad = np.sqrt(x**2 + y**2)
    if rad > 1.0:
       return 0.0

    azi = np.atan2(x, y)
    # 等面积投影公式逆运算
    ain = 2.0 * np.asin(rad / np.sqrt(2.0))

    # 构建单位矢量 r
    r = np.zeros(3)
    r[0] = np.sin(ain) * np.cos(azi)
    r[1] = np.sin(ain) * np.sin(azi)
    r[2] = np.cos(ain)

    # 计算 rpp = r^T * M * r
    val = 0.0
    for i in range(3):
        for j in range(3):
            val += r[i] * mt[i][j] * r[j]
    return val


def marching_squares_optimized(mt, resolution=300):
    """优化版的 Marching Squares，直接生成网格并提取线段"""
    print(f"执行优化算法: 分辨率 {resolution}x{resolution}...")

    x_range = np.linspace(-1, 1, resolution)
    y_range = np.linspace(-1, 1, resolution)

    # 1. 预计算所有网格点的值 (向量化计算)
    grid_v = np.zeros((resolution, resolution))
    for i in range(resolution):
        for j in range(resolution):
            grid_v[i, j] = get_mt_value_at_coord(x_range[i], y_range[j], mt)

    lines = []
    # 2. 遍历单元格
    for i in range(resolution - 1):
        for j in range(resolution - 1):
            # 四个角的值
            v00, v10 = grid_v[i, j], grid_v[i + 1, j]
            v01, v11 = grid_v[i, j + 1], grid_v[i + 1, j + 1]

            p00, p10 = (x_range[i], y_range[j]), (x_range[i + 1], y_range[j])
            p01, p11 = (x_range[i], y_range[j + 1]), (x_range[i + 1], y_range[j + 1])

            # 仅处理圆内的单元格（简单判断）
            if p00[0] ** 2 + p00[1] ** 2 > 1.05:
                continue

            # 查找穿过 0 值的边
            edges = []
            # 检查四条边
            if (v00 * v10) < 0:  # 下边
                t = -v00 / (v10 - v00)
                edges.append((p00[0] + t * (p10[0] - p00[0]), p00[1]))
            if (v10 * v11) < 0:  # 右边
                t = -v10 / (v11 - v10)
                edges.append((p10[0], p10[1] + t * (p11[1] - p10[1])))
            if (v01 * v11) < 0:  # 上边
                t = -v01 / (v11 - v01)
                edges.append((p01[0] + t * (p11[0] - p01[0]), p01[1]))
            if (v00 * v01) < 0:  # 左边
                t = -v00 / (v01 - v00)
                edges.append((p00[0], p00[1] + t * (p01[1] - p00[1])))

            if len(edges) >= 2:
                lines.append([edges[0], edges[1]])
    return lines


def connect_contour_lines(lines):
    """
    将线段连接成连续的轮廓线
    """
    if not lines:
        return []

    remaining = lines.copy()
    contour_lines = []

    while remaining:
        current_line = [remaining[0][0], remaining[0][1]]
        remaining.pop(0)

        changed = True
        while changed:
            changed = False
            for i, (p1, p2) in enumerate(remaining):
                if (
                    abs(p1[0] - current_line[0][0]) < 1e-8
                    and abs(p1[1] - current_line[0][1]) < 1e-8
                ):
                    current_line.insert(0, p2)
                    remaining.pop(i)
                    changed = True
                    break
                elif (
                    abs(p2[0] - current_line[0][0]) < 1e-8
                    and abs(p2[1] - current_line[0][1]) < 1e-8
                ):
                    current_line.insert(0, p1)
                    remaining.pop(i)
                    changed = True
                    break
                elif (
                    abs(p1[0] - current_line[-1][0]) < 1e-8
                    and abs(p1[1] - current_line[-1][1]) < 1e-8
                ):
                    current_line.append(p2)
                    remaining.pop(i)
                    changed = True
                    break
                elif (
                    abs(p2[0] - current_line[-1][0]) < 1e-8
                    and abs(p2[1] - current_line[-1][1]) < 1e-8
                ):
                    current_line.append(p1)
                    remaining.pop(i)
                    changed = True
                    break

        contour_lines.append(current_line)

    return contour_lines


def calculate_contour_zero(mt, resolution=80, num_points_per_segment=30):
    """封装后的快速计算函数"""
    # 1. 快速提取线段
    lines = marching_squares_optimized(mt, resolution)

    # 2. 连接线段
    print("正在连接轮廓线...")
    contour_lines = connect_contour_lines(lines)
    return contour_lines


# ...............................


# find fault plane's strike, dip, rake angles
# using fault normal (vn) and slip vector (vs)
def vn2sdr(vn, vs):

    epsi = 0.001

    if vn[0] < 0.0:  # Upwards normal
        vn = -1.0 * vn
        vs = -1.0 * vs

    if vn[0] > 1.0 - epsi:  # Horizontal plane
        strike = 0
        dip = 0
        rake = np.rad2deg(np.arctan2(-vs[2], -vs[1]))
    elif vn[0] < epsi:  # Vertical plane
        strike = np.rad2deg(np.arctan2(vn[1], vn[2]))
        dip = np.rad2deg(np.pi / 2.0)
        rake = np.rad2deg(np.arctan2(vs[0], -vs[1] * vn[2] + vs[2] * vn[1]))
    else:  # Oblique plane
        strike = np.rad2deg(np.arctan2(vn[1], vn[2]))
        dip = np.rad2deg(np.arccos(vn[0]))
        rake = np.rad2deg(np.arctan2((-vs[1]*vn[1]-vs[2]*vn[2]),(-vs[1]*vn[2]+vs[2]*vn[1])*vn[0]))

    if strike < 0:
        strike = strike + 360.0

    if rake < -180:
        rake = rake + 360.0
    if rake > 180:
        rake = rake - 360.0

    return strike, dip, rake


# find strike, dip, rake from gCMT moment tensor
# cmt=[-1.350,5.410,-4.060,-3.210,-3.580,-0.736]
#        Mrr    Mtt   Mpp   Mrt    Mrp    Mtp
def mt2sdr(cmt):

    mt = np.array(
        [[cmt[0], cmt[3], cmt[4]],
        [cmt[3], cmt[1], cmt[5]],
        [cmt[4], cmt[5], cmt[2]]])

    m, v = np.linalg.eig(mt)
    # idx = m.argsort()[::-1]  # sort m in descending order
    idx = np.argsort(m)  # sort eigen value m in ascending order
    m = m[idx]  # sort eigen value m in ascending order (m1<m2<m3)
    # After sort m in ascening order,
    # the first element of m correspond to P axis,
    # the second element correspond to null axis,
    # the third element correspond to T axis
    v = v[:, idx]  # sort eigen vector according to eigen value

    # Now first colum of matrix v correspond to vector of P axis, the second colum of matrix v
    # correspond to vector of null axis and the third colum of v correspond to the vector of T axis
    # calculate P, null and T axis plunge and azimuth angle using eigen vector
    for i in range(3):
        azm = np.rad2deg(np.arctan2(v[2][i], -v[1][i]))
        scale = v[1][i] * v[1][i] + v[2][i] * v[2][i]
        plg = np.rad2deg(np.arctan2(-v[0][i], np.sqrt(scale)))
        if plg < 0.0:
            plg = -1.0 * plg
            azm = azm + 180.0

        azm = math.fmod(azm, 360.0)

        if azm < 0:
            azm = azm + 360.0

        print("eigval={0:.2f} plunge={1:.1f} azimuth={2:.1f}".format(m[i], plg, azm))

    # construct eigen vector correspond P, null and T axis
    p = v[:, 0]
    n = v[:, 1]
    t = v[:, 2]

    # calculate fault normal and slip vector
    v1 = (t + p) / np.sqrt(2.0)
    v2 = (t - p) / np.sqrt(2.0)

    # focal parameter for fault plane I
    strike, dip, rake = vn2sdr(v1, v2)
    print("strike={0:.1f} dip={1:.1f} rake={2:.1f}".format(strike, dip, rake))

    # focal parameter for fault plane II
    # strike,dip,rake=vn2sdr(v2,v1)

    return strike, dip, rake


######################## main program ####################

points = []
values = []


# ------- using strike, dip, rake to plot beachball --------------
# Focal mechanism represents with strike dpi, rake, m0
# strike, dip, rake = 315.0, 45.0, -12.0
# strike1, dip1, rake1 = getaux(strike, dip, rake)
# mt = dislocationtomt(strike,dip,rake,1)
# ----------------------------------------------------------------

# ++++++++ using gCMT moment tensor to plot beachball ++++++++++++
# “南(x) - 东(y) - 上(z)”系统
#     0(zz) 1(xx)  2(yy)  3(xz)  4(-yz) 5(-xy)
# cmt=[0.801, 0.026, -0.827, -1.010, -0.349, -0.598]
#      Mrr    Mtt   Mpp   Mrt    Mrp    Mtp
# cmt=[-1.770, 0.855, 0.919, -0.484, 0.235, -0.613]
# cmt=[-0.367, 2.770, -2.400, 0.238, 0.943, -0.351]
cmt = [0.817, -0.238, -0.579, 0.213, 0.334, -0.413]
# cmt=[0.854, -1.070, 0.212, -1.520, -4.270, -1.520]
# cmt=[-0.003, 0.133, -0.130, 0.612, 0.084, 4.330]

# Mxx = Mtt
# Myy = Mpp
# Mzz = Mrr
# Mxz = Mrt
# Myz = -Mrp
# Mxy = -Mtp

# 特别注意两个坐标系转换时对应关系
mt = np.array(
    [[cmt[1], -cmt[5], cmt[3]],
    [-cmt[5], cmt[2], -cmt[4]],
    [cmt[3], -cmt[4], cmt[0]]]
)

strike, dip, rake = mt2sdr(cmt)
strike1, dip1, rake1 = getaux(strike, dip, rake)
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

plt.figure(figsize=(6, 6))
plt.axis("equal")
r = [0.0, 0.0, 0.0]
aP = mt[0][0]
aT = mt[0][0]
jxP, jyP = 0.0, 0.0
jxT, jyT = 0.0, 0.0

# 第一步：计算值辐射花样的数值，不绘图
points = []
values = []
for i in range(-20, 20):
    x = i / 20.0
    ymax = np.sqrt(1 - x**2)
    m = int(20.0 * ymax)
    # print(m)
    if m > 0:
        # print(m)
        for j in range(-m, m):
            y = j / 20.0
            rad = np.sqrt(x**2 + y**2)
            azi = np.atan2(x, y)
            ain = 2.0 * np.asin(rad / np.sqrt(2.0))
            r[0] = np.sin(ain) * np.cos(azi)
            r[1] = np.sin(ain) * np.sin(azi)
            r[2] = np.cos(ain)
            rpp = 0.0
            for ii in range(3):
                for jj in range(3):
                    rpp += r[ii] * mt[ii][jj] * r[jj]

            points.append([x, y])
            values.append(rpp)

            # if rpp > 0.0:
            #    plt.scatter(x, y, s=rpp*35, c='black', edgecolors='white', linewidths=1)
            # plt.text(x, y, str(rpp), fontsize=8, ha='center', va='center')
            # if rpp < 0.0:
            # plt.scatter(x, y, s=abs(rpp*55), c='blue', marker='_', edgecolors='white', linewidths=1)
            # plt.text(x, y, str(rpp), fontsize=8, ha='center', va='center')
            # if rpp > -0.033 and rpp < 0.033:
            # plt.scatter(x, y, s=45, c='black')
            if rpp < aP:
                aP = rpp
                jxP = x
                jyP = y
            if rpp > aT:
                aT = rpp
                jxT = x
                jyT = y

# plt.scatter(jxP, jyP, s=abs(aP*55), marker='^', c='blue')
# plt.text(jxP, jyP, str(aP), fontsize=8, ha='center', va='center')
# plt.scatter(jxT, jyT, s=abs(aT*55), marker='^', c='green')
# plt.text(jxT, jyT, str(aT), fontsize=8, ha='center', va='center')

# 第二步：找出辐射花样数值的最大绝对值，用于归一化绘图中小圆的大小
max_abs_rpp = max([abs(v) for v in values]) if values else 1.0

# 第三步：用不同颜色及不同半径大小的圆点表示辐射花样的大小及正负
base_size = 70  # 你认为最合适的最大点大小
for p, rpp in zip(points, values):
    if rpp > 0.0:
        # 归一化处理：(rpp / max_abs_rpp) 的范围永远在 0 到 1 之间
        # 这样无论 CMT 值多大，点的大小都会锁定在 0 到 base_size 之间
        s_normalized = (rpp / max_abs_rpp) * base_size
        plt.scatter(p[0], p[1], s=s_normalized, c="black", edgecolors="white", linewidths=0.5)
    # if rpp < 0.0:
    #    s_normalized = (abs(rpp) / max_abs_rpp) * base_size
    #    plt.scatter(p[0], p[1], s=s_normalized, c='blue', edgecolors='white', linewidths=0.5)

circle = plt.Circle((0, 0), 1, color="black", fill=False)
plt.gca().add_patch(circle)

x, y = pnodal(strike, dip, rake)
plt.plot(x, y, color="purple", linewidth=1.5)

x, y = pnodal(strike1, dip1, rake1)
plt.plot(x, y, color="blue", linewidth=1.5)

plt.xlim(-1, 1)
plt.ylim(-1, 1)
title = "strike=%.1f,dip=%.1f,rake=%.1f" % (strike, dip, rake)
plt.title(title)

# *********** Extract contour = 0 lines, matplotlib and scipy function is not working right
#            write my ownself code

# 提取到辐射花样数值为零点的坐标，并且把线段用不同颜色画出来
smoothed_lines = calculate_contour_zero(mt, resolution=300, num_points_per_segment=30)

i = 0
if smoothed_lines:
    for line in smoothed_lines:
        line = np.array(line)
        if i == 0:
            plt.plot(line[:, 0], line[:, 1], color="green", linewidth=1.5)
        elif i == 1:
            plt.plot(line[:, 0], line[:, 1], color="red", linewidth=1.5)

        print(i)
        i = i + 1


# 显示网格
# plt.grid(True)
plt.axis("off")

plt.show()
