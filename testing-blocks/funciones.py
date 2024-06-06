# this module will be imported in the into your flowgraph



# function to calculate the d(xv,yv)
def dv(x1, y1, x2, y2, fv):


	return np.abs((y1 - y2) * fv[0] - (x1 - x2) * fv[1] - x2 * y1 + x1 * y2) / np.sqrt((y1 - y2) ** 2 + (x1 - x2) ** 2)


# function to calculate the s(xv,yv)
def sv(x1, y1, z1, x2, y2, z2, fv):
    if z1 <= z2:
        s_v = (((y1 - y2) ** 2 + (x1 - x2) ** 2 + (fv[0] - x1) ** 2 + (fv[1] - y1) ** 2 -
                ((fv[0] - x2) ** 2 + (fv[1] - y2) ** 2)) / (2 * np.sqrt((y1 - y2) ** 2 + (x1 - x2) ** 2))) + z1
    else:
        s_v = (((y1 - y2) ** 2 + (x2 - x1) ** 2 + (fv[0] - x2) ** 2 + (fv[1] - y2) ** 2 -
                ((fv[0] - x1) ** 2 + (fv[1] - y1) ** 2)) / (2 * np.sqrt((y1 - y2) ** 2 + (x2 - x1) ** 2))) + z2
    return s_v


# %funcion para calcular vector y distancia entre dos puntos
def point_to_vector(x1, y1, z1, x2, y2, z2):
    Vec = np.array([x2 - x1, y2 - y1, z2 - z1])
    length = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)
    
    return Vec, length


# Calculo del vector normal para transmisores
def norm_vec_trans(alpha, beta):
    n = np.array([np.cos(np.deg2rad(alpha)) * np.sin(np.deg2rad(beta)),
                  np.sin(np.deg2rad(alpha)) * np.sin(np.deg2rad(beta)),
                  -np.cos(np.deg2rad(beta))])
    return n

# funcion para calcular el producto punto
def dot_product(a, b):
    f = np.dot(a, b)
    return f

# Calculo del vector normal para receptores
def norm_vec_receiver(alpha, beta):
    m = np.array([np.cos(np.deg2rad(alpha)) * np.sin(np.deg2rad(beta)),
                  np.sin(np.deg2rad(alpha)) * np.sin(np.deg2rad(beta)),
                  np.cos(np.deg2rad(beta))])
    return m

# funcion para calcular la ganancia
def gain(eta, incide, incide_r, fov):
    if 0 <= incide_r <= 2 * fov:
        g = (eta ** 2) / (np.sin(np.deg2rad(incide)) ** 2)
    else:
        g = 0
    return g

# Function to calculate the  Pij for shadowing model
def P_expt(gv, fv, W, H, X, Y, t, es, d_v, s_v):
    p = 0.1
         
    if (gv[0] >= 2 * d_v) and (gv[1] >= s_v):
        w_int = gv[0] * W / 2
        h_int = gv[1] * H / 2
        A = [w_int, h_int]
        x_int = fv[0] * X / 2
        y_int = fv[1] * Y / 2
        B = [x_int, y_int]
        exp_value = np.dot(A, B)
        f = p * t
        est = -es * exp_value
        d = np.exp(est)
        Pij = d
    else:
        Pij = 0
    
    return Pij


# funcion para calcular el HLoS
def HLoS_direct(x_i, y_i, z_i, x_j, y_j, z_j, Ap, eta, alpha_i, alpha_j, beta_i, beta_j, incidencia, incidencia_r, m, fov, gv, fv, W, H, X, Y, t, es, c):
    dv_ij = dv(x_i, y_i, x_j, y_j, fv)
    sv_ij = sv(x_i, y_i, z_i, x_j, y_j, z_j, fv)
    Pij = P_expt(gv, fv, W, H, X, Y, t, es, dv_ij, sv_ij)
    v1, d1 = point_to_vector(x_i, y_i, z_i, x_j, y_j, z_j)
    Nnorm1 = norm_vec_trans(alpha_i, beta_i)
    p1 = dot_product(v1, Nnorm1)
    v2, d2 = point_to_vector(x_j, y_j, z_j, x_i, y_i, z_i)
    Nnorm2 = norm_vec_receiver(alpha_j, beta_j)
    p2 = dot_product(v2, Nnorm2)
    g = gain(eta, incidencia, incidencia_r, fov)
    dm = d1 / c

    if 0 <= incidencia <= 2 * fov:
        m_HLoS = abs(((m + 1) * Ap / (2 * np.pi * d1 ** 2)) * (p1 ** m / d1) * (p2 / d2) * g * Pij)
    else:
        m_HLoS = 0

    return m_HLoS, dm

