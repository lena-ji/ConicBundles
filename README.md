# ConicBundles
This contains supplementary code for the paper "Rationality of real conic bundles with quartic discriminant curve" by Lena Ji and Mattie Ji.

## Files

The functions of the files are as follows:

- [Quadric-bundle-verifications.sage](Quadric-bundle-verifications.sage) can be used to verify the claims made about the examples constructed in Sections $4.1$ and $4.3$ of the paper. Given three quadrics $Q_1, Q_2, Q_3$ with rational coefficients, it checks smoothness of the curves $\Delta=(Q_1 Q_3-Q_2^2=0)$ and $\tilde{\Delta}$, computes the real branch points of the genus two curve $\Gamma$, and computes the signature of the quadric surface bundle $Y \dashrightarrow \mathbb{P}^1$ on each interval between real branch points of $\Gamma$. This code is a $\texttt{Sage}$ implementation of the $\texttt{Magma}$ code accompanying the paper "Curve classes on conic bundle threefolds and applications to rationality" by S. Frei, L. Ji, S. Sankar, B. Viray, and I. Vogt.

- [P1tilde.sage](P1tilde.sage) can be used to verify the claims about the intermediate Jacobian torsor $\tilde{P}^{(1)}$ for the examples in Section $4.2$. Given three quadrics $Q_1, Q_2, Q_3$, and a line in $\mathbb{P}^2$, this code first checks whether this line meets the curve $\Delta=(Q_1 Q_3 - Q_2^2=0)$ in four distinct complex points. If it does, the code then checks whether a $\operatorname{Gal}(\mathbb{C}/\mathbb{R})$-invariant set of four points over this intersection spans a $3$-plane in $\mathbb{P}^4$ by computing the rank of the $4 \times 4$ submatrices. If at least one of these has full rank, then this set of four points gives a real point on $\tilde{P}^{(1)}$.

- [P1tilde-bitangents.sage](P1tilde-bitangents.sage) also looks for real points on $\tilde{P}^{(1)}$. It takes $Q_1, Q_2, Q_3$ as inputs and, instead of checking over a given line, it checks whether $\tilde{P}^{(1)}$ has any real points mapping to a real bitangent of the quartic curve $\Delta=(Q_1 Q_3-Q_2^2=0)$. The portion of this code that computes the real bitangents of the quartic curve $\Delta$ was written by D. Plaumann, B. Sturmfels, and C. Vinzant, and is from the code "quartictype.sage" accompanying their paper "Quartic curves and their bitangents."

- [Singular-members.m2](Singular-members.m2) takes as inputs quadrics $Q_1, Q_2, Q_3$ where one of the coefficients of $Q_1=Q_1(s)$ varies in a one-parameter family parametrized by $s$. It outputs the values of $s$ at which the signature sequence of the fibers of the quadric surface bundle can change (by computing the discriminant of the polynomial defining the branch locus of $\Gamma$) and the values of $s$ for which the quartic curve $\Delta$ is singular. (If both of the curves $\Delta$ and $\tilde{\Delta}$ are smooth, then so is the threefold $Y$.) We do not use this code in our paper, but we include it in this GitHub repository for the interested reader who may want to construct some examples in families.

Note that throughout the paper we use the variables $u,v,w$ for the coordinates of $\mathbb{P}^2$, whereas in the code we use the variables $x,y,z$ instead.
