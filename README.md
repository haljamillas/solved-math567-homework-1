Download Link: https://assignmentchef.com/product/solved-math567-homework-1
<br>
1.   (Order of accuracy). The following three tables of errors were produced from different finite difference approximations to u0(x), for some function u(x). The first column of each table shows the size h used in the difference approximation, and the second column is the absolute error, e(h), in that approximation, as a function of h. Assume that the error can be modeled as

e(h) ≈ Chp,

where p is an integer. Approximate the order of accuracy p for each sequence. You can find the data for these sequences in the file errors.txt on the Blackboard website, under the homework 1 link.

————————— ————————–h err(h) h err(h)

—————————                                                               —————————7.8125000e-032.0844e-02                 7.8125000e-031.9059e-033.9062500e-031.1118e-02                 3.9062500e-034.3086e-041.9531250e-035.3455e-03                 1.9531250e-031.0318e-049.7656250e-042.7049e-03                 9.7656250e-042.6007e-054.8828125e-041.3469e-03                 4.8828125e-04

———————–h    err(h)

————————

1.00000e+00          1.3829e-02

5.00000e-01          1.8805e-03

2.50000e-01          1.3742e-04

1.25000e-01          8.9170e-06

6.25000e-02          5.6252e-07

3.12500e-02          3.5239e-08

1.56250e-02          2.2037e-09

7.81250e-03          1.3775e-10

3.90625e-03          8.6098e-126.5716e-062.   (Finite-volume scheme) An alternative way to approximate the second derivative than those discussed in the class and Chapter 1 of the book is to use the Divergence Theorem (or, in 1D, the Fundamental Theorem of Calculus). The technique is fundamental to finite-volume methods where cell-centered grids are used.

One convention for a cell-centered grid is to let integer indices indicate the cell boundaries and half-integer indices indicate the center of the cell, where the function (or it’s average value) is given. This is illustrated in the figure below, where the circle at xj+1/2 is the center of the cell

[xj,xj+1], i.e. :

The second derivative of some function u at xj+1/2 can be approximated as

.

Here we have approximated u00 at xj+1/2 by its average value over the cell [xj,xj+1] and then applied the Divergence Theorem. If we approximate the first derivative as

and

we obtain the discretization

.     (1)

(a)    Show that if we assume a uniform grid, i.e. xj = x0 +jh, h = 1/n, the stencil you obtain for the approximation to u00(xj+1/2) at the cell centers xj+1/2 is essentially the same second order approximation we obtained in class for u00(xj) on a vertex centered grid.

(b)   Show that on a non-equispaced grid, however, this approximation differs from the technique we discussed in class for approximating the second derivative (i.e. using Lagrange interpolation or method of undetermined coefficients).

You can show this by giving a specific example with non-equally spaced points. For example, let x0 = 0, x1 = 1, x2 = 7/4, and x3 = 3, and compute the stencils for u00(x3/2) using (1) and the technique discussed in class (or simply use Fornberg’s weights algorithm [2]; see course Blackboard page).

Alternatively, you can show that (1) is exact for u = 1 and u = x, but not u = x2, then argue that the technique we discussed in class has to give the exact value for these functions.

3.   (Compact finite difference formulas) It is often desirable in certain applications to keep the width (i.e. the number of points) in a finite difference (FD) formula relatively small. This unfortunately directly impacts the order of accuracy we can obtain with the formula. For example, the centered, 3-point FD formulas for the first and second derivative, given, respectively, as

(2)

,                                  (3)

are only second order accurate (i.e. O(h2)). The order of accuracy can often be increased without changing the width of the formula by using an implicit (or compact or Hermite) FD formula, albeit at the cost of solving a sparse and banded linear system [1]. For example, the following are the most popular centered, 3-point implicit FD formulas for the first and second derivative of u at ¯x:

(4)

.           (5)

(a)   Write a code that implements these two implicit formulas for a given function u sampled at m+2 equally spaced points xj over the interval [a,b] (i.e. xj = a+jh, j = 0,1,…,m+1, h = (b−a)/(m+1)). In your calculation, you may only use information about u0(x) and u00(x) at the boundary points x0 = a and xm+1 = b. This means to get the values of the derivatives at the interior you will need to solve two m-by-m tridiagonal systems of equations—one for the first derivative and one for the second. You can solve this system using any built-in linear solver to the software you are using (e.g. in Matlab I would use the spdiags and ‘backslash’ functions).

(b)   Apply your code to approximate the first and second derivatives of the function u(x) = x2e−x over [0,1]. Generate tables and plots showing how the relative two-norm of the errors in the approximations over the interval [0,1] decreases as the grid spacing h decreases (or N increases). Start with h = 1/8 (m = 7) and successively reduce h by a factor of 2 until h = 1/256 (m = 255). The plots should be made on a log-log scale (i.e., logarithmic on by the x and y axes). Determine the orders of accuracy of these two implicit formulas experimentally from the plots and table.

(c)    The formal orders of accuracy and exact truncation errors of these implicit FD formulas can be determined in a number of ways. For example, we can use the Taylor series approach discussed in class for the standard FD formulas. Computer Algebra Systems, such as Mathematica or Maple, make this approach considerably easy. For example, the Mathematica code (which can be used directly in WolframAlpha) for deriving the truncation error for the implicit FD formula (4) is given by: s1 = Series[(D[u[x – h] + 4 u[x] + u[x + h], {x, 1}]) – 3/h(-u[x – h] + u[x + h]), {h, 0, 6}] The corresponding Maple code is given by: convert(series(diff(u(x-h)+4*u(x)+u(x+h), x)-3/h*(-u(x-h)+u(x+h)), h = 0, 6), diff)

Use either of these codes to determine the formal truncation error and order of accuracy of (4) then modify the code and do the same for (5). Do your numerical results from part (a) match theoretical results for the order of accuracy? Explain.

If you have ever studied interpolation with splines then you may recognize the implicit FD formula for the first derivative. It is the same as the one that arises when working out the conditions on continuity of a cubic spline with equally spaced knots. Thus the order of accuracy for computing the derivative with a cubic spline interpolant at the nodes xj is four (assuming the boundaries are handled appropriately), this is one order higher than one might expect based on the accuracy of the interpolant and is sometimes called super-convergence.

4.   (Increasing the FD stencil width) As discussed in class, Fornberg’s algorithm [2] can be used to rapidly generate finite difference (FD) weights of any order for arbitrarily spaced points in one dimension. Matlab , Python, and Julia functions for this algorithm are posted on the course Blackboard page under the name weights.m, weights.py, and weights.jl, respectively (these are similar to the books Matlab function fdcoeffF). You can download any of these functions to perform the tasks listed below.

(a)   For the following three node sets (which are known as equispaced, Chebyshev, and Legendre, respectively), use Fornberg’s algorithm to compute the FD weights for approximating the first derivative at x = 0 and x = −1 + 3/14:

-0.937273392400706 -0.201194093997435  0.724417731360170

-0.848206583410427                   0                    0.848206583410427

-0.724417731360170 0.201194093997435  0.937273392400706

-0.570972172608539 0.394151347077563  0.987992518020485.

Ordering the Legendre points in (iii) from smallest (-0.987…) to largest (0.987…) will make things easier in what follows.

You should use all the nodes in the approximations, which means there should be 15 weights for each node set. Plot each set of weights for the x = 0 approximation on the same graph (use the index j as the horizontal or independent coordinate). Repeat this for the x = −1+3/14 approximation. Comment on how the weights differ between node sets and between approximation points.

(b) Use the weights from each of the three node sets from part (a) to approximate the derivative of u(x) = e−cos(2(x−1/5)) at x = 0 and x = −1 + 3/14. Report the error in these approximations in a nice table and comment on how the errors differ. Which set of points would you prefer to use to approximate derivatives?

Note that in the case of equally spaced and Chebyshev points, analytical expressions can be worked out for the FD weights (see the next problem for an example).

5. (Trigonometric interpolation) The trigonometric interpolant at N equally spaced nodes over [0,2π) (i.e. xj =  2Nπj, j = 0,1,…,N −1) to a function u can be written in Lagrange form as

N−1 pN(x) = X SN(x − xj)u(xj),

j=0

where

is even, is odd.

The function SN(t) is called the ‘periodic sinc function’.

(a)    Verify that pN(x) actually interpolates the data, i.e. show that pN(xk) = u(xk), for k = 0,1,…,N − 1. This can be done by simply showing SN(xj − xk) = 0 when k 6= j and SN(0) = 1.

(b)   For the case of N being even show that



, 0

.

These are the Fourier (or trigonometric) weights for approximating the first derivative of u at x = 0 from samples of u at xk. Since SN(t) is periodic, they can be shifted appropriately to approximate the derivative of u at any xk.