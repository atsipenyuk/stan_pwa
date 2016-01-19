import numpy as np
import matplotlib.pyplot as plt

def integral(func, bounds, N=10000):
    """
    Compute a definite Monte Carlo integral.

    Integrate func from bounds[i][0] to bounds[i][1] for i-th variable using 
    uniform random sampling.

    Parameters
    ----------
    func : function
        A Python function or method to integrate.
    bounds : list of 2-elemental lists
        List of elements of the type [a_i,b_i], a_i and b_i are floats, lower
        and upper limits of integration of i-th variable.
    N : int, optional
        2*N is the number of drawn samples.

    Returns
    -------
    res : float
        The integral of func from bounds[:][0] to bounds[:][1]
    abserr : float
        A crude estimate of the absolute error of the integration
    """
    n_vars = len(bounds)
    random_point = np.zeros(n_vars, dtype=float)
    total_volume = np.product([(bounds[i][1] - bounds[i][0]) for i in range(n_vars)])

    res1 =  sum(__eval_integrate(func, bounds, N, n_vars)) \
           / float(N) * total_volume
    res2 =  sum(__eval_integrate(func, bounds, N, n_vars)) \
           / float(N) * total_volume

    return [(res1+res2)/2., np.abs(res1-res2)/2]


def __eval_integrate(func, m, N_points, n_vars):
    # Generator yielding func evaluated at N_points random points.
    num = 0
    while num < N_points:
        yield func([np.random.uniform(m[i][0], m[i][1]) for i in range(n_vars)])
        num += 1


def integral_of_tensor_product(func, bounds, N=1000000):
    """
    Similar function, more oriented to our purposes.

    For the function func that takes a list of floats and
    returns a complex vector A, this function returns the matrix
    I[i,j] = \int A*[i] A[j].
    """
    n_vars = len(bounds)
    random_point = np.zeros(n_vars, dtype=float)
    total_volume = np.product([(bounds[i][1] - bounds[i][0]) for i in range(n_vars)])

    res1 =  sum(__eval_integrate_A(func, bounds, N, n_vars)) \
           / float(N) * total_volume
    res2 =  sum(__eval_integrate_A(func, bounds, N, n_vars)) \
           / float(N) * total_volume

    return [(res1+res2)/2., np.abs(res1-res2)/2]

def __eval_integrate_A(func, m, N_points, n_vars):
    # Generator yielding I for func evaluated at N_points random points,
    # as described above.
    num = 0
    while num < N_points:
        A = func([np.random.uniform(m[i][0], m[i][1]) for i in range(n_vars)])
        n_res = A.shape[0]
        yield np.asarray([[np.conj(A[i1]) * A[i2] for i1 in range(n_res)] for i2 in range(n_res)])
        num += 1


# Overload the functions above for a tensor product of two different functions
def integral_of_tensor_product_2func(func1, func2, bounds, N=1000000):
    """
    For functions func1, func2 that take a list of floats and
    return complex vectors A1, A2 this function returns the matrix
    I[i,j] = \int conj(A1)[i] *  A2[j].
    """
    n_vars = len(bounds)
    random_point = np.zeros(n_vars, dtype=float)
    total_volume = np.product([(bounds[i][1] - bounds[i][0]) for i in range(n_vars)])

    res1 =  sum(__eval_integrate_A(func1, func2, bounds, N, n_vars)) \
           / float(N) * total_volume
    res2 =  sum(__eval_integrate_A(func1, func2, bounds, N, n_vars)) \
           / float(N) * total_volume

    return [(res1+res2)/2., np.abs(res1-res2)/2]

def __eval_integrate_A_2func(func1, func2, m, N_points, n_vars):
    # Generator yielding I for func evaluated at N_points random points,
    # as described above.
    num = 0
    while num < N_points:
        rand_pt = [np.random.uniform(m[i][0], m[i][1]) for i in range(n_vars)]
        A1 = func1(rand_pt)
        A2 = func2(rand_pt)
        n1 = A1.shape[0]
        n2 = A2.shape[0]
        yield np.asarray([[np.conj(A[i1]) * A[i2] for i1 in range(n1)] for i2 in range(n2)])
        num += 1


# Overload the functions above for a tensor product of a function with
# a binned function (i.e., this latter function takes a data point and 
# returns an integer (bin number) to which this data point belongs.
def integral_of_tensor_product_w_bins(func, bin_fct, bounds, N=1000000):
    """
    Consider 
        func: [x_1..x_n_vars] |-> (A_1 .. A_R);
        bin_fct: [x_1..x_n_vars] |-> i;
    A_j are complex, i is an integer between 1 and NUM_BINS. 
    Then, this function returns a RxNUM_BINS matrix with following
    entries:
        I[r,b] = \int_{bin b} A_r(y) dy
    """
    n_vars = len(bounds)
    random_point = np.zeros(n_vars, dtype=float)
    total_volume = np.product([(bounds[i][1] - bounds[i][0]) for i in range(n_vars)])

    res1 =  sum(__eval_integrate_A(func, bin_fct, bounds, N, n_vars)) \
           / float(N) * total_volume
    res2 =  sum(__eval_integrate_A(func, bin_fct, bounds, N, n_vars)) \
           / float(N) * total_volume

    return [(res1+res2)/2., np.abs(res1-res2)/2]


def __eval_integrate_A_w_bins(func1, func2, m, N_points, n_vars):
    # Generator yielding I for func evaluated at N_points random points,
    # as described above.
    num = 0
    while num < N_points:
        rand_pt = [np.random.uniform(m[i][0], m[i][1]) for i in range(n_vars)]
        A1 = func1(rand_pt)
        A2 = func2(rand_pt)
        n1 = A1.shape[0]
        n2 = A2.shape[0]
        yield np.asarray([[np.conj(A[i1]) * A[i2] for i1 in range(n1)] for i2 in range(n2)])
        num += 1



#### Same functions - except random points are not generated within bounds,
#### but passed explicitely as y_data argument.


def integral_w_pts(func, y_data, n_pts, volume):
    """
    Compute a definite Monte Carlo integral.

    Parameters
    ----------
    func : function
        A Python function or method to integrate.
    pts_gen : point generator. Generates elements 
        like [x_1_i,x_2_i, ..., x_var_i] - with other words, generates
        points, at which func can be evaluated.
    n_pts : integer, number of points that will be generated

    Returns
    -------
    res : float
        The integral of func.
    abserr : float
        A crude estimate of the absolute error of the integration
     """
    pts1 = [list(y_data[:,i]) for i in range(0,n_pts/2)]
    pts2 = [list(y_data[:,i]) for i in range(n_pts/2, n_pts)]


    # Evaluate the function at these points
    A_1 = np.asarray(map(func, pts1))
    A_2 = np.asarray(map(func, pts2))

    # Number of arguments our function func takes
    n_vars = len(pts1[:][0])
    n_res = len(A_1[0])

    res1 = A_1.sum(axis=0) * volume
    res2 = A_2.sum(axis=0) * volume

    return [(res1+res2) / float(n_pts), np.abs(res1-res2) / float(n_pts)]



def integral_of_tensor_product_w_pts(func, y_data, n_pts, volume):
    """
    Compute a definite Monte Carlo integral.

    For a function func = (A_1 .. A_NUM_RES) and phase space point list
    rnd_pts = [y_1 ... y_d], return the Monte Carlo integrals
        I[i,j] = int A_i(y) * A_j(y) dy.
    The points y_1 .. y_d consist of a list of variables: y_i = [m_1, .. m_n].

    Parameters
    ----------
    func : function
        A Python function or method to integrate.
    pts_gen : Genetator that generates elements 
        of the type [x_1_i,x_2_i, ..., x_var_i].

    Returns
    -------
    res : float
        The integral of func.
    abserr : float
        A crude estimate of the absolute error of the integration
    """
    # A list containing points where we want to evaluate
    pts1 = [list(y_data[:,i]) for i in range(0,n_pts/2)]
    pts2 = [list(y_data[:,i]) for i in range(n_pts/2, n_pts)]

    # Evaluate the function at these points
    A_1 = map(func, pts1)
    A_2 = map(func, pts2)

    # Number of arguments our function func takes
    n_vars = len(pts1[:][0])
    n_res = len(A_1[0])

    # Define the function I[r,r'](y_d) = conj(A_r(y_d)) * A_r'(y_d)
    def I(A):
        return np.asarray([[np.conj(A[r]) * A[s] for r in range(n_res)] for s in range(n_res)])

    I_1 = np.asarray(map(I, A_1))
    I_2 = np.asarray(map(I, A_2))

    res1 = I_1.sum(axis=0) * volume
    res2 = I_2.sum(axis=0) * volume

    return [(res1+res2) / float(n_pts), np.abs(res1-res2)/float(n_pts)]





def integral_b_w_pts(func, amp, n_bins, y_data, n_pts, volume):
    """
    Compute a definite Monte Carlo integral for a 'binned function'.

    Consider the function func: double |-> int.
    With other words, for a data point the function returns the
    number of bin into which that point falls. Assume the total
    number of bins is num_bins. 'integral_b_w_pts' returns a vector
    res with num_bins entries. For each entry
              res[i] = #events_landed_in_bin_i / n_pts * volume 
                       product(amp(events_landed_in_bin_i)).

    amp(y_data_i) is the amplitude associated with each event;
    it is a function taking same arguments as func.
    Number of data points is n_pts; they are stored in the list y_data.
    An example entry of y_data is [y_1i, ... y_num_var_i]
    (num_var denotes the number of variables that 'func' takes).

    I am sorry for bad documentation, but its late and I have flu and
    deadline. Sorry.

    Returns
    -------
    res : list of floats
        The integral of func over each bin.
    abserr : list of floats
        A crude estimate of the absolute error of the integration
     """
    pts1 = [list(y_data[:,i]) for i in range(0,n_pts/2)]
    pts2 = [list(y_data[:,i]) for i in range(n_pts/2, n_pts)]

    # These vars are list of integers - e.g., A_1 = [1, 5, 8, 5, ..]
    # tells 1st event is in 1st bin, 2nd evt. is in 5th bin, etc.
    A_1_flat = np.asarray(map(func, pts1))
    A_2_flat = np.asarray(map(func, pts2))

    amp1 = np.asarray(map(amp, pts1))
    amp2 = np.asarray(map(amp, pts2))

    # Now repack these results: res1 = [52, 15, ...]
    # means that in the 1st bin, there were 52 events, in 2nd bin, 15 evts...
    res1 = np.zeros(n_bins, dtype=float)
    res2 = np.zeros(n_bins, dtype=float)

    for i in A_1_flat:
        res1[i-1] += amp1
    for i in A_2_flat:
        res1[i-1] += amp2

    # Rescale the results with the volume over which we are integrating
    res1 = res1 * volume
    res2 = res1 * volume

    return [(res1+res2) / float(n_pts), np.abs(res1-res2) / float(n_pts)]





def integral_b_tensor_func_w_pts(b, func_b, n_bins, 
                                 func, y_data, n_pts, volume):
    """
    Compute a definite Monte Carlo integral.

    The idea is to return the tensor product integral of two functions,
    just as before. Lets call these functions 'func' and 'func2'.
    The argument 'y_data' is a list of points where these functions
    are evaluated at.

    Now, assume that 'func2' has a weird structure: it is a binned function.
    What I mean by that is that this function is computed in the following 
    way: take an event y_i from the list y_data. The function 'b' takes
    'y_i' and tells, which bin it falls into. The function func_b also takes
    'y_i' and returns the value func2(y_i) - the actual value
    of the function. However, func_b does not know whether we are in the 
    right bin.

    Let the target space of 'func' be R-dimensional:
                  func(y_i) = [x1 ... xR],
    and the target space of 'func2' be n_bins-dimensional:
                  finc(y_i) = [x1 ... x_n_bins].
    This method returns an R x n_bins -dimensional matrix with following
    entries:

        res[r,n] = integral func(y)[r] * func2(y)[n] dy.

    And the integration is the Monte Carlo integration over the events
    stored in the list 'y_data'.

    Returns
    -------
    res : matrix of floats
        The integral of func cross func2.
    abserr : matrix of floats
        A crude estimate of the absolute error of the integration
    """
    # A list containing points where we want to evaluate
    pts1 = [list(y_data[:,i]) for i in range(0,n_pts/2)]
    pts2 = [list(y_data[:,i]) for i in range(n_pts/2, n_pts)]

    # Determine the bins corresponding to the events
    b_1 = map(b, pts1)
    b_2 = map(b, pts2)

    # Determine the amplitudes corresponding to the events for func_b
    func_b_1 = map(func_b, pts1)
    func_b_2 = map(func_b, pts2)

    # Determine the amplitudes corresponding to the events for func
    func_1 = map(func, pts1)
    func_2 = map(func, pts2)

    # Number of arguments our function func takes
    n_res = len(func_1[0])

    print('b_1')
    print(b_1[0])
    print('func_b_1')
    print(func_b_1[0])
    print('func_1')
    print(func_1[0])

    # Loop over all events, add the appropriate amplitude
    # to each bin (thus, taking the average for each bin)
    res1 = np.zeros([n_bins, n_res], dtype=complex)
    res2 = np.zeros([n_bins, n_res], dtype=complex)
    for i in range(len(b_1)):
        res1[b_1[i]-1, :] += np.multiply(np.conj(func_b_1[i]), func_1[i]) 
    for i in range(len(b_2)):
        res2[b_2[i]-1, :] += np.multiply(np.conj(func_b_2[i]), func_2[i]) 

    return [(res1+res2) / float(n_pts), np.abs(res1-res2)/float(n_pts)]




def integral_b_tensor_b_w_pts(b, func_b, n_bins_b,
                              c, func_c, n_bins_c, 
                              y_data, n_pts, volume):
    """
    Compute a definite Monte Carlo integral.

    As above, but with two binned functions.

    Returns
    -------
    res : matrix of floats
        The integral of func cross func2.
    abserr : matrix of floats
        A crude estimate of the absolute error of the integration
    """
    # A list containing points where we want to evaluate
    pts1 = [list(y_data[:,i]) for i in range(0,n_pts/2)]
    pts2 = [list(y_data[:,i]) for i in range(n_pts/2, n_pts)]

    # Determine the bins corresponding to the events
    b_1 = map(b, pts1)
    b_2 = map(b, pts2)

    # Determine the amplitudes corresponding to the events for func_b
    func_b_1 = map(func_b, pts1)
    func_b_2 = map(func_b, pts2)

    # The same for the functions b2, func_b2
    # Determine the bins corresponding to the events
    c_1 = map(c, pts1)
    c_2 = map(c, pts2)

    # Determine the amplitudes corresponding to the events for func_b
    func_c_1 = map(func_c, pts1)
    func_c_2 = map(func_c, pts2)

    # Loop over all events, add the appropriate amplitude
    # to each bin (thus, taking the average for each bin)
    res1 = np.zeros([n_bins_b, n_bins_c], dtype=float)
    res2 = np.zeros([n_bins_b, n_bins_c], dtype=float)

    for i in range(len(pts1)):
        res1[b_1[i]-1, c_1[i]-1] += np.conj(func_b_1[i]) * func_c_1[i] 
    for i in range(len(pts2)):
        res2[b_2[i]-1, c_2[i]-1] += np.conj(func_b_2[i]) * func_c_2[i] 

    return [(res1+res2) / float(n_pts), np.abs(res1-res2)/float(n_pts)]
