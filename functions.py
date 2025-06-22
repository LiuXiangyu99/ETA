from autograd import grad  
# import statement for autograd wrapped numpy
import autograd.numpy as np   
# for finding inverse function
from scipy.optimize import root
from scipy.stats import norm
from itertools import combinations
from scipy.special import comb
import matplotlib.pyplot as plt




def find_indices(data, k):
    n = len(data)
    indexed_data = sorted([(val, idx) for idx, val in enumerate(data)], reverse=True)
    J = set()
    R = dict()
    r = 0
    for (_, idx) in indexed_data:
        if r > k+1:
            break
        r += 1
        R[idx] = r 
        J.add(idx)
    return J, R


def find_indices_boot(data, multi, k):
    n = len(data)
    indexed_data = sorted([(val, idx) for idx, val in enumerate(data)], reverse=True)
    J = set()
    R = dict()

    multi_av = np.mean(multi)

    r = 0
    for (_, idx) in indexed_data:
        if r > k:
            break
        r += multi[idx]/multi_av
        R[idx] = r 
        J.add(idx)

    return J, R


def rank(data, J):
    R = dict()
    n = len(data)

    for idx in J:
        r = 0
        for i in range(n):
            if data[i] >= data[idx]:
                r += 1
        R[idx] = r
    return R


def rank_boot(data, multi, J):
    R = dict()
    n = len(data)
    multi_av = np.mean(multi)
    for idx in J:
        r = 0
        for i in range(n):
            if data[i] >= data[idx]:
                r += multi[i]/multi_av
        R[idx] = r 
    return R



# computing eta
def eta(xs, ys, k):
    # sorted_ys and corresponded xs
    corr_xs = np.zeros_like(xs)
    # indices
    indices = sorted(range(len(ys)), key=lambda k: ys[k], reverse=True)
    corr_xs = xs[indices]
    # compute rank
    rs = np.zeros_like(xs)
    for i in range(len(corr_xs)):
        r = 0
        for j in range(len(corr_xs)):
            if corr_xs[j] >= corr_xs[i]:
                r += 1
        rs[i] = r
    # compute
    estimator = 0.0
    for i in range(k):
        for j in range(k):
            estimator += 3/k**3 * max((k+1 - max(rs[i], rs[j])), 0)
    return estimator


# computing bootstrap eta
def eta_boot(xs, ys, multi, k):
    multi_av = np.mean(multi)
    J_boot, _ = find_indices_boot(ys, multi, k)
    R_boot = rank_boot(xs, multi, J_boot)

    estimator = 0.0
    for i in J_boot:
        for j in J_boot:
            estimator += multi[j]/multi_av * max(0, k+multi[j]/multi_av - max(R_boot[i], R_boot[j]))
    return 3/k**3 * estimator


# computing delta
def delta(xs, ys, k):
    return eta(xs, ys, k) - eta(ys, xs, k)


def delta_boot(xs, ys, multi, k):
    return eta_boot(xs, ys, multi, k) - eta_boot(ys, xs, multi, k)



def inverse_cdf(f):
    def bisearch(f, p, eps=1e-3):
        xl = 0
        xr = 1
        x = (xl + xr)/2
        while abs(f(x) - p) > eps:
            if f(x) > p:
                xr = x
                xl = xl
                x = (xl + xr)/2
            else:
                xl = x
                xr = xr
                x = (xl + xr)/2
        return x
    return (lambda p: bisearch(f, p))

def sample_one(c): 
    """    
    pdf of c is ctns bivariate copula
    define helper function
    """ 
    cdf_du = grad(c, 0)
    u, t = np.random.random(2)
    cond_cdf = lambda v: cdf_du(u, v)
    inv_cond_cdf = inverse_cdf(cond_cdf)
    v = inv_cond_cdf(t)
    return u, v

def sampling(c, n):
    data = np.zeros([n, 2])
    for i in range(n):
        data[i, :] = sample_one(c)
    return data


def gumbel(u, v, delta=3):
    expr = -((-np.log(u))**delta + (-np.log(v))**delta) ** (1/delta)
    return np.exp(expr)


def asymmetrization(c, a, b):
    return (lambda u, v: u**a * v**b * c(u**(1-a), v**(1-b)))


# Sample from copula C(u, v; t) = min(u, tv + (1-t)(u + v - 1)_{+})

def sample_C_min(n, t=2/3, rng=None):
    """
    Draw n samples from the copula
        C(u,v; t) = min(u, t v + (1-t) (u+v-1)_+),  0<=u,v<=1, 0<=t<=1.

    Parameters
    ----------
    n   : int
        Number of samples required.
    t   : float in [0,1]
        Copula parameter.  t=0 -> independence, t=1 -> comonotone (u=v).
    rng : np.random.Generator, optional
        Numpy random number generator (for reproducibility).

    Returns
    -------
    uv  : ndarray, shape (n, 2)
        Samples (U_i, V_i) from the copula.
    """
    if rng is None:
        rng = np.random.default_rng()

    # ----- 1. determine which points fall on the line and which in the triangle
    on_line     = rng.random(n) < t          # Bernoulli(t)
    n_line      = on_line.sum()
    n_triangle  = n - n_line

    # ----- 2. samples on the line  u = t v
    v_line = rng.random(n_line)              # V ~ Unif(0,1)
    u_line = t * v_line                      # U = t V

    # ----- 3. samples in the triangle  u+v>1  (uniform over that region)
    #
    #  • Draw (X,Y) uniform in the unit square
    #  • Keep only those with X+Y > 1  (triangular pdf = 2 inside region)
    #    but we can do it without rejection: take independent U',W and map
    #    (U = 1 - W + W*U',  V = 1 - W + W*(1-U'))  for W,U' ~ Unif(0,1).
    #
    if n_triangle > 0:
        W  = rng.random(n_triangle)          # controls distance from (1,1)
        U_ = rng.random(n_triangle)          # barycentric coordinate

        u_tri = 1 - W + W * U_
        v_tri = 1 - W + W * (1 - U_)
    else:
        u_tri = v_tri = np.empty(0)

    # ----- 4. stack the two parts together
    u = np.empty(n)
    v = np.empty(n)

    u[ on_line] = u_line
    v[ on_line] = v_line
    u[~on_line] = u_tri
    v[~on_line] = v_tri

    return np.column_stack((u, v))



def sample_max_model(n, k, rng=None):
    if rng is None:
        rng = np.random.default_rng()

    # Step 1: generate k independent vectors xs_i of shape (k, n)
    xs_all = rng.uniform(0, 1, size=(k, n))
    # Step 2: define xs = xs_1
    xs = xs_all[0, :]
    # Step 3: coordinate-wise max
    ys = np.max(xs_all, axis=0)
    # Step 4: stack into data array
    data = np.column_stack((xs, ys))
    return data



def asymmetric_tail_kendall_tau(xy_array, q=0.90, direction='XtoY'):
    """
    Estimate asymmetric tail Kendall's tau.
    
    Parameters:
    - xy_array: np.ndarray of shape (n, 2), each row is (X_i, Y_i)
    - q: quantile level to define threshold in X (or Y)
    - direction: 'XtoY' or 'YtoX' → condition on X or Y tail
    
    Returns:
    - tau_hat: estimated asymmetric tail Kendall's tau
    """
    X = xy_array[:, 0]
    Y = xy_array[:, 1]
    n = len(X)
    
    if direction == 'XtoY':
        threshold = np.quantile(X, q)
        mask = X > threshold
    elif direction == 'YtoX':
        threshold = np.quantile(Y, q)
        mask = Y > threshold
    else:
        raise ValueError("direction must be 'XtoY' or 'YtoX'")
    
    # Indices of points in the tail
    tail_indices = np.where(mask)[0]
    k = len(tail_indices)
    # print(k)
    
    if k < 2:
        # Not enough points to compute Kendall's tau
        return np.nan
    
    # Loop over all pairs
    count = 0
    total = comb(k, 2)
    
    for i, j in combinations(tail_indices, 2):
        sign_X = np.sign(X[i] - X[j])
        sign_Y = np.sign(Y[i] - Y[j])
        count += sign_X * sign_Y
    
    tau_hat = count / total
    return tau_hat



def edm_estimator(data, k, b=(lambda t:t)):
    """
    assume the margin is standard Pareto
    """
    n = data.shape[0]
    numerator = 0
    denumerator = 0

    for i in range(n):
        x, y = data[i]
        if (x**2 + y**2)**0.5 > b(n/k):
            theta_i = np.arctan(y/x)
            numerator += (theta_i - np.pi/4)**2
            denumerator += 1
    return 1 - (numerator / denumerator)/(np.pi/4)**2



def margin_transform(data):
    """
    Transform data with uniform margins (copula data) to CDF F(x) = 1 - 1/(x+1).

    Parameters:
    - data: np.ndarray, shape (n, d), assumed to have uniform[0,1] entries

    Returns:
    - transformed_data: np.ndarray, shape (n, d), with standard Pareto margins
    """
    epsilon = 1e-10
    data = np.clip(data, 0, 1 - epsilon)  # avoid division by zero
    return 1.0 / (1.0 - data) - 1.0