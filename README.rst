This repositary contains the solution to the coding homework for APMA
2810 Discontinuous Galerkin Methods by Prof. Chi-Wang Shu. Below is a
description of these homeworks(including 1 take-home exam). The written
homework is not included here.

**HW1**
~~~~~~~

Use DG to find the solution to the following equation:
:raw-latex:`\begin{equation}
\left\{
             \begin{array}{lr}
             u_x=\cos x, &  0\leq x \leq 1\\
             u(0)=0 &  
             \end{array}
\right.
\end{equation}` and plot :math:`L_1`, :math:`L_2`, :math:`L_{\infty}`
error tables.

**HW2**
~~~~~~~

| Use DG to find the solution to the following equation with two
  different initial conditions on :math:`[0,2\pi]\times[0,2\pi]`:
  :raw-latex:`\begin{equation}
  \left\{
               \begin{array}{lr}
               u_t+u_x=0, \\
               u(x,0)=\sin(x) &  
               \end{array}
  \right.
  \end{equation}`
| :raw-latex:`\begin{equation}
  \left\{
               \begin{array}{lr}
               u_t+u_x=0, \\
               u(x,0)=1, &  x\in(\frac{\pi}{2},\frac{3\pi}{2})\\
               u(x,0)=0, &  x\in(0,\frac{\pi}{2})\cup(\frac{3\pi}{2},2\pi)
               \end{array}
  \right.
  \end{equation}` and plot :math:`L_1`, :math:`L_2`, :math:`L_{\infty}`
  error tables.

**HW3**
~~~~~~~

Plot moment error tables of the problem last week.

**HW4**
~~~~~~~

Apply TVD limiter and TVB limiter with :math:`M = 0.1`, :math:`1`,
:math:`5`, :math:`10` to the problem in week 2. Plot error tables and
report the points changed by the limiters.

**HW5**
~~~~~~~

Apply Bound Preserving Limiter to the problem in week 2. Plot error
tables.

**HW6**
~~~~~~~

Use LDG with central flux and alternating flux to find the solution to
the following equation on :math:`[0,2\pi]\times[0,1]`:
:raw-latex:`\begin{equation}
\left\{
             \begin{array}{lr}
             u_t=u_{xx}, \\
             u(x,0)=\sin(x) &  
             \end{array}
\right.
\end{equation}` and plot :math:`L_1`, :math:`L_2`, :math:`L_{\infty}`
error tables.

**HW7**
~~~~~~~

Use Bauman-Oden, SIPG and Ultra-Weak scheme to solve the problem last
week and plot error tables.

**HW8**
~~~~~~~

Use LDG to find the solution to the following equation on
:math:`[0,2\pi]\times[0,1]`: :raw-latex:`\begin{equation}
\left\{
             \begin{array}{lr}
             u_t=u_{xxx}, \\
             u(x,0)=\sin(x) &  
             \end{array}
\right.
\end{equation}` and plot :math:`L_1`, :math:`L_2`, :math:`L_{\infty}`
error tables.

**HW9**
~~~~~~~

Use Ultra-weak scheme to solve the problem last week and plot error
tables.

**Final Exam**
~~~~~~~~~~~~~~

| Use LDG to find the solution to the following two equations on
  :math:`[0,2\pi]\times[0,1]`: :raw-latex:`\begin{equation}
  \left\{
               \begin{array}{lr}
               u_t+u_x=\epsilon u_{xx}, \\
               u(x,0)=\sin(x) &  
               \end{array}
  \right.
  \end{equation}`
| :raw-latex:`\begin{equation}
  \left\{
               \begin{array}{lr}
               u_t+u_x=\epsilon u_{xxx}, \\
               u(x,0)=\sin(x) &  
               \end{array}
  \right.
  \end{equation}` and plot :math:`L_1`, :math:`L_2`, :math:`L_{\infty}`
  error tables.
