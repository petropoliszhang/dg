This repositary contains the solution to the coding homework for APMA 2810 Discontinuous Galerkin Methods by Prof. Chi-Wang Shu. Below is a description of these homeworks(including 1 take-home exam). The written homework is not included here.

### **HW1**

Use DG to find the solution to the following equation:
<p align="center"><img alt="\begin{equation}&#10;\left\{&#10;             \begin{array}{lr}&#10;             u_x=\cos x, &amp;  0\leq x \leq 1\\&#10;             u(0)=0 &amp;  &#10;             \end{array}&#10;\right.&#10;\end{equation}" src="https://rawgit.com/zzhang222/dg/master/svgs/13344509d7232284f0fb074711415c7a.png" align="middle" width="447.71743994999997pt" height="39.452455349999994pt"/></p>
and plot <img alt="$L_1$" src="https://rawgit.com/zzhang222/dg/master/svgs/929ed909014029a206f344a28aa47d15.png" align="middle" width="17.73978854999999pt" height="22.465723500000017pt"/>, <img alt="$L_2$" src="https://rawgit.com/zzhang222/dg/master/svgs/4327ea69d9c5edcc8ddaf24f1d5b47e4.png" align="middle" width="17.73978854999999pt" height="22.465723500000017pt"/>, <img alt="$L_{\infty}$" src="https://rawgit.com/zzhang222/dg/master/svgs/986e40e0f11ab0c97cfd953d6e3fb747.png" align="middle" width="24.292324649999987pt" height="22.465723500000017pt"/> error tables.

### **HW2**

Use DG to find the solution to the following equation with two different initial conditions on <img alt="$[0,2\pi]\times[0,2\pi]$" src="https://rawgit.com/zzhang222/dg/master/svgs/0f1880c6e7bbc3ac10285066ac8ba019.png" align="middle" width="105.76486634999999pt" height="24.65753399999998pt"/>:
<p align="center"><img alt="\begin{equation}&#10;\left\{&#10;             \begin{array}{lr}&#10;             u_t+u_x=0, \\&#10;             u(x,0)=\sin(x) &amp;  &#10;             \end{array}&#10;\right.&#10;\end{equation}" src="https://rawgit.com/zzhang222/dg/master/svgs/825df8af772fed0c9d30e6b4442110fd.png" align="middle" width="429.42427605pt" height="39.452455349999994pt"/></p>
\
<p align="center"><img alt="\begin{equation}&#10;\left\{&#10;             \begin{array}{lr}&#10;             u_t+u_x=0, \\&#10;             u(x,0)=1, &amp;  x\in(\frac{\pi}{2},\frac{3\pi}{2})\\&#10;             u(x,0)=0, &amp;  x\in(0,\frac{\pi}{2})\cup(\frac{3\pi}{2},2\pi)&#10;             \end{array}&#10;\right.&#10;\end{equation}" src="https://rawgit.com/zzhang222/dg/master/svgs/f7d752fc0ab8308a0f68f8073ba3a8ac.png" align="middle" width="491.1943278pt" height="59.178683850000006pt"/></p>
and plot <img alt="$L_1$" src="https://rawgit.com/zzhang222/dg/master/svgs/929ed909014029a206f344a28aa47d15.png" align="middle" width="17.73978854999999pt" height="22.465723500000017pt"/>, <img alt="$L_2$" src="https://rawgit.com/zzhang222/dg/master/svgs/4327ea69d9c5edcc8ddaf24f1d5b47e4.png" align="middle" width="17.73978854999999pt" height="22.465723500000017pt"/>, <img alt="$L_{\infty}$" src="https://rawgit.com/zzhang222/dg/master/svgs/986e40e0f11ab0c97cfd953d6e3fb747.png" align="middle" width="24.292324649999987pt" height="22.465723500000017pt"/> error tables.

### **HW3**

Plot moment error tables of the problem last week.

### **HW4**

Apply TVD limiter and TVB limiter with <img alt="$M = 0.1$" src="https://rawgit.com/zzhang222/dg/master/svgs/5781b091b470317acd7f461f904dcf36.png" align="middle" width="60.66198929999999pt" height="22.465723500000017pt"/>, <img alt="$1$" src="https://rawgit.com/zzhang222/dg/master/svgs/034d0a6be0424bffe9a6e7ac9236c0f5.png" align="middle" width="8.219209349999991pt" height="21.18721440000001pt"/>, <img alt="$5$" src="https://rawgit.com/zzhang222/dg/master/svgs/9612eecfec9dadf1a81d296bd2473777.png" align="middle" width="8.219209349999991pt" height="21.18721440000001pt"/>, <img alt="$10$" src="https://rawgit.com/zzhang222/dg/master/svgs/b0c08f9b595a704efb907fc688034d80.png" align="middle" width="16.438418699999993pt" height="21.18721440000001pt"/> to the problem in week 2. Plot error tables and report the points changed by the limiters.

### **HW5**

Apply Bound Preserving Limiter to the problem in week 2. Plot error tables.

### **HW6**

Use LDG with central flux and alternating flux to find the solution to the following equation on <img alt="$[0,2\pi]\times[0,1]$" src="https://rawgit.com/zzhang222/dg/master/svgs/0bb8c0cc4fff220f3dc9ff3b13b2c985.png" align="middle" width="95.80477829999998pt" height="24.65753399999998pt"/>:
<p align="center"><img alt="\begin{equation}&#10;\left\{&#10;             \begin{array}{lr}&#10;             u_t=u_{xx}, \\&#10;             u(x,0)=\sin(x) &amp;  &#10;             \end{array}&#10;\right.&#10;\end{equation}" src="https://rawgit.com/zzhang222/dg/master/svgs/f73503b1bc5e3c3bf8cb90960bbdb18f.png" align="middle" width="429.42427605pt" height="39.452455349999994pt"/></p>
and plot <img alt="$L_1$" src="https://rawgit.com/zzhang222/dg/master/svgs/929ed909014029a206f344a28aa47d15.png" align="middle" width="17.73978854999999pt" height="22.465723500000017pt"/>, <img alt="$L_2$" src="https://rawgit.com/zzhang222/dg/master/svgs/4327ea69d9c5edcc8ddaf24f1d5b47e4.png" align="middle" width="17.73978854999999pt" height="22.465723500000017pt"/>, <img alt="$L_{\infty}$" src="https://rawgit.com/zzhang222/dg/master/svgs/986e40e0f11ab0c97cfd953d6e3fb747.png" align="middle" width="24.292324649999987pt" height="22.465723500000017pt"/> error tables.

### **HW7**

Use Bauman-Oden, SIPG and Ultra-Weak scheme to solve the problem last week and plot error tables.

### **HW8**

Use LDG to find the solution to the following equation on <img alt="$[0,2\pi]\times[0,1]$" src="https://rawgit.com/zzhang222/dg/master/svgs/0bb8c0cc4fff220f3dc9ff3b13b2c985.png" align="middle" width="95.80477829999998pt" height="24.65753399999998pt"/>:
<p align="center"><img alt="\begin{equation}&#10;\left\{&#10;             \begin{array}{lr}&#10;             u_t=u_{xxx}, \\&#10;             u(x,0)=\sin(x) &amp;  &#10;             \end{array}&#10;\right.&#10;\end{equation}" src="https://rawgit.com/zzhang222/dg/master/svgs/9554cd7d12947524b47a123734ad8cee.png" align="middle" width="429.42427605pt" height="39.452455349999994pt"/></p>
and plot <img alt="$L_1$" src="https://rawgit.com/zzhang222/dg/master/svgs/929ed909014029a206f344a28aa47d15.png" align="middle" width="17.73978854999999pt" height="22.465723500000017pt"/>, <img alt="$L_2$" src="https://rawgit.com/zzhang222/dg/master/svgs/4327ea69d9c5edcc8ddaf24f1d5b47e4.png" align="middle" width="17.73978854999999pt" height="22.465723500000017pt"/>, <img alt="$L_{\infty}$" src="https://rawgit.com/zzhang222/dg/master/svgs/986e40e0f11ab0c97cfd953d6e3fb747.png" align="middle" width="24.292324649999987pt" height="22.465723500000017pt"/> error tables.

### **HW9**

Use Ultra-weak scheme to solve the problem last week and plot error tables.

### **Final Exam**

Use LDG to find the solution to the following two equations on <img alt="$[0,2\pi]\times[0,1]$" src="https://rawgit.com/zzhang222/dg/master/svgs/0bb8c0cc4fff220f3dc9ff3b13b2c985.png" align="middle" width="95.80477829999998pt" height="24.65753399999998pt"/>:
<p align="center"><img alt="\begin{equation}&#10;\left\{&#10;             \begin{array}{lr}&#10;             u_t+u_x=\epsilon u_{xx}, \\&#10;             u(x,0)=\sin(x) &amp;  &#10;             \end{array}&#10;\right.&#10;\end{equation}" src="https://rawgit.com/zzhang222/dg/master/svgs/bd1e61bdc049b95a08d40047a42420a0.png" align="middle" width="429.42427605pt" height="39.452455349999994pt"/></p>
\
<p align="center"><img alt="\begin{equation}&#10;\left\{&#10;             \begin{array}{lr}&#10;             u_t+u_x=\epsilon u_{xxx}, \\&#10;             u(x,0)=\sin(x) &amp;  &#10;             \end{array}&#10;\right.&#10;\end{equation}" src="https://rawgit.com/zzhang222/dg/master/svgs/2c27fc6a7dbaac2b81ce8d7484c7e75a.png" align="middle" width="433.0896207pt" height="39.452455349999994pt"/></p>
and plot <img alt="$L_1$" src="https://rawgit.com/zzhang222/dg/master/svgs/929ed909014029a206f344a28aa47d15.png" align="middle" width="17.73978854999999pt" height="22.465723500000017pt"/>, <img alt="$L_2$" src="https://rawgit.com/zzhang222/dg/master/svgs/4327ea69d9c5edcc8ddaf24f1d5b47e4.png" align="middle" width="17.73978854999999pt" height="22.465723500000017pt"/>, <img alt="$L_{\infty}$" src="https://rawgit.com/zzhang222/dg/master/svgs/986e40e0f11ab0c97cfd953d6e3fb747.png" align="middle" width="24.292324649999987pt" height="22.465723500000017pt"/> error tables.
