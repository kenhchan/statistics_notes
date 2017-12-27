
## Introduction

First we recall the following basic fact: given a standard normally distributed random variable, its square follows a chi-squared distribution with <em>1</em> degree of freedom. More generally, given <em>n</em> such independent random variables, the sum of their squares follow a chi-squared distribution with <em>n</em> degrees of freedom. The following theorem of Cochran shows that a decomposition of the canonical positive real quadratic form into components of minimal rank is actually a direct sum decomposition; a fact which yields the joint distribution of the components. 

### Notation

The <em>n</em> x <em>n</em> identity matrix will be denoted by <em>I<sub>n</sub></em>.
The <em>n</em> x <em>n</em> matrix where every entry is <em>1</em> will be denoted by <b>1</b><em><sub>n</sub></em>.

### Theorem (Cochran)
Let <em>Z</em> ~ <em>N</em>(<em>0</em>, <em>I<sub>n</sub></em>). Suppose we have a decomposition

<em>I<sub>n</sub></em> = <em>Q<sub>1</sub></em> + ... + <em>Q<sub>k</sub></em>

where <em>Q<sub>1</sub></em>, ..., <em>Q<sub>k</sub></em> are real symmetric matrices with ranks <em>r<sub>1</sub></em>, ..., <em>r<sub>k</sub></em>
satisfying <em>r<sub>1</sub></em> + ... + <em>r<sub>k</sub></em> = n. Then the random variable <em>U<sub>i</sub></em> = <em>Z</em><sup>T</sup> <em>Q<sub>i</sub></em> Z
follows a chi-squared distribution with <em>r<sub>i</sub></em> degrees of freedom for <em>i</em> = <em>1</em>, ..., <em>k</em>.
Moreover, <em>U<sub>1</sub></em>, ..., <em>U<sub>k</sub></em> are jointly independent.

#### Proof
First recall the following basic fact in linear algebra: given matrices <em>A</em> and <em>B</em> of the same 
shape, the rank of their sum is bounded above by the sum of their ranks, that is

rank(<em>A</em> + <em>B</em>) &le; rank(<em>A</em>) + rank(<em>B</em>).

Now let <em>Q</em> = <em>Q<sub>i</sub></em> for some <em>i</em> and <em>R</em> = <em>I<sub>n</sub></em> - <em>Q</em>. Then

rank(<em>R</em>)  &le; <em>n</em> - <em>r<sub>i</sub></em> 

By the spectral theory for symmetric matrices over ℝ, there exists an orthogonal matrix <em>P</em> and a 
diagonal matrix <em>D</em> such that <em>P</em><sup>T</sup><em>Q</em> <em>P</em> = <em>D</em>. 
The change of basis matrix <em>P</em> also diagonalises <em>R</em>, since <em>P</em><sup>T</sup><em>R</em> <em>P</em> = <em>I<sub>n</sub></em> - <em>D</em>.

Since <em>Q</em> has rank <em>r<sub>i</sub></em>, the matrix <em>D</em> contains an (<em>N</em> - <em>r<sub>i</sub></em>) x (<em>N</em> - <em>r<sub>i</sub></em>)
block of zeroes on the diagonal, so <em>I<sub>n</sub></em> - <em>D</em> has an <em>I<sub>N-r<sub>i</sub></sub></em> block on the
diagonal. Hence

rank(<em>T</em>) =  rank(<em>I<sub>n</sub></em> - <em>D</em>) &ge; <em>n - r<sub>i</sub></sub></em>

so in fact rank(<em>T</em>) = <em>n - r<sub>i</sub></sub></em>. 

This shows that the <em>I<sub>N-r<sub>i</sub></sub></em> diagonal block in <em>I<sub>n</sub></em> - <em>D</em> 
is the only nonzero block, and 
consequently <em>D</em> has an <em>I<sub>r<sub>i</sub></sub></em> block on the diagonal and is zero elsewhere. 
In other words 
<em>1</em> is an eigenvalue of <em>Q</em> with multiplicity <em>r<sub>i</sub></em>, 
and the eigenspaces of <em>Q</em> and <em>T</em> 
corresponding to nonzero eigenvalues are orthogonal.

Therefore <em>P</em> diagonalises all the <em>Q<sub>i</sub></em>s simultaneously, and that

<em>P</em><sup>T</sup> <em>Q<sub>i</sub></em> <em>P</em> 

is zero everywhere except for 1s along the following diagonal positions

(<em>x</em> + <em>1</em>, <em>x</em> + <em>1</em>) , ..., (<em>x</em> + <em>r<sub>i</sub></em>, <em>x</em> + <em>r<sub>i</sub></em>)

where x = <em>r<sub>1</sub></em> + ... + <em>r<sub>i-1</sub></em>.

Next we show that <em>P</em> <em>Z</em> ~ <em>N</em>(<em>0</em>, <em>I<sub>n</sub></em>). Since <em>P</em> is linear, 
we have 𝔼(<em>P</em> <em>Z</em>) = <em>P</em> 𝔼(<em>Z</em>) = <b>0</b>. Next we show that 
Cov(<em>P</em> <em>Z</em>) = <em>I<sub>n</sub></em>. 
The change of basis formula for the covariance matrix is 

Cov(<em>P</em> <em>Z</em>) = <em>P</em><sup>T</sup> Cov(<em>Z</em>) <em>P</em>.

Since Cov(<em>Z</em>) = <em>I<sub>n</sub></em> and <em>P</em> is orthogonal, the result follows. 
Since a linear combination of normally distributed random variables is again normal, 
we have <em>P</em> <em>Z</em> ~ <em>N</em>(<em>0</em>, <em>I<sub>n</sub></em>). This establishes
the joint distribution of <em>U<sub>1</sub></em>, ..., <em>U<sub>k</sub></em>.

### Application: One-way ANOVA

Suppose we have <em>i</em> = <em>1</em>, ..., <em>m</em> groups and each group has <em>n<sub>i</sub></em> observations.
Denote by <em>X<sub>ij</sub></em> the <em>j</em>-th observation in the <em>i</em>-th group, and let <em>N</em> = 
<em>n<sub>1</sub></em> + ... + <em>n<sub>m</sub></em> be the number of observations. We denote by <em>X<sub>i.</sub></em>
the mean of the observations in the <em>i</em>-th group, and <em>X<sub>..</sub></em> be the aggregate sample mean.

The usual null hypothesis in this situation is that the population means of each group are the same. Let <em>&mu;</em> be the common population mean of each group and <em>&sigma;<sup>2</sup></em> be the population variance. We will normalise the observed random variables by setting

<em>Z<sub>i</sub></em> = (<em>X<sub>i</sub></em> - <em>&mu;</em>) / <em>&sigma;</em>

for each <em>i</em>. The symbols <em>Z<sub>i.</sub></em> and <em>Z<sub>..</sub></em> are defined analogously. The three sums of squares usually computed in one-way ANOVA are 

SS<sub>total</sub>    = &Sigma; (<em>Z<sub>ij</sub></em> - <em>Z<sub>..</sub></em>)<sup><em>2</em></sup>

SS<sub>group</sub>    = &Sigma; (<em>Z<sub>i.</sub></em> - <em>Z<sub>..</sub></em>)<sup><em>2</em></sup>

SS<sub>residual</sub> = SS<sub>total</sub> - SS<sub>group</sub> 

where the first sum is over all possible <em>i</em>, <em>j</em> and ths second sum is over all possible <em>i</em>. We will show the following

### Proposition

SS<sub>group</sub> ~ <em>&chi;<sup>2</sup></em>(<em>m</em> - <em>1</em>), and

SS<sub>residual</sub> ~ <em>&chi;<sup>2</sup></em>(<em>N</em> - <em>m</em>).

#### Proof
Consider the following decomposition 

<em>I<sub>N</sub></em> = (<em>I<sub>N</sub></em> - <em>H</em>) + (<em>H</em> - <b>1</b><sub><em>N</em></sub>/N) + <b>1</b><sub><em>N</em></sub>/N

where <em>H</em> = diag(<b>1</b><sub><em>n<sub>1</sub></em></sub>/<em>n<sub>1</sub></em>, ..., <b>1</b><sub><em>n<sub>m</sub></em></sub>/<em>n<sub>m</sub></em>). The three summands in the decomposition of <em>I<sub>N</sub></em> are clearly symmetric matrices. It is elementary to verify the following calculations (<b>1</b><sub><em>N</em></sub>/N)<sup><em>2</em></sup> = (<b>1</b><sub><em>N</em></sub>/N), <em>H<sup>2</sup></em> = <em>H</em> and <em>H</em> <b>1</b><sub><em>N</em></sub>/N = <b>1</b><sub><em>N</em></sub>/N. These show that each of the three summands above are idempotent. Idempotent matrices have characteristic equation <em>x</em>(<em>x</em> - <em>1</em>) = 0, so they only have <em>0</em> or <em>1</em> as eigenvalues. Hence the rank of an idempotent matrix is equal to its trace. This allows us to determine the ranks of each summand above easily:

rank(<em>I<sub>N</sub></em> - <em>H</em>) = tr(<em>I<sub>N</sub></em>) - tr(<em>H</em>) = <em>N</em> - <em>m</em>

rank(<em>H</em> - <b>1</b><sub><em>N</em></sub>/N) = tr(<em>H</em>) - tr(<b>1</b><sub><em>N</em></sub>/N) = <em>m</em> - <em>1</em>

rank(<b>1</b><sub><em>N</em></sub>/N) = <em>1</em>
