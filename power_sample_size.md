Let <em>&theta;&#770;<sub>n</sub></em> be an estimator for a parameter <em>&theta;</em> (obtained from a sample of size <em>n</em>). 
We wish to determine whether there is significant evidence to reject a null hypothesis 
<em>H<sub>0</sub></em>: <em>&theta;</em> = <em>&theta;<sub>0</sub></em> in favour of a simple alternative 
<em>H<sub>1</sub></em>: <em>&theta;</em> = <em>&theta;<sub>1</sub></em> &gt; <em>&theta;<sub>0</sub></em>.
Suppose that <em>&theta;&#770;<sub>n</sub></em> ∼ <em>F<sub>0</sub></em> under <em>H<sub>0</sub></em> and 
<em>&theta;&#770;<sub>n</sub></em> ∼ <em>F<sub>1</sub></em> under <em>H<sub>1</sub></em>.

Given a significance level <em>&alpha;</em>, we define the rejection region <em>R</em> by 

<em>R</em> = [<em>F<sub>0</sub><sup>-1</sup></em>(<em>1</em> - <em>&alpha;</em>), &infin;)

The power of the hypothesis test is the measure of <em>R</em> under <em>H<sub>1</sub></em>, i.e. 

Power = <em>1</em> - <em>F<sub>1</sub></em>(<em>F<sub>0</sub><sup>-1</sup></em>(<em>1</em> - <em>&alpha;</em>))


## Example

Let <em>X&#772;<sub>n</sub></em> be the sample mean of a sample of size <em>n</em>. 
Suppose the population mean and variance are <em>&mu;</em> and <em>&sigma;<sup>2</sup></em>. 
By the central limit theorem <em>X&#772;<sub>n</sub></em> converges in distribution to
N(<em>&mu;</em>, <em>&sigma;<sup>2</sup></em>/<em>n</em>) as <em>n</em> &rarr; &infin;.
Let <em>&tau;<sub>n</sub></em> = <em>&sigma;</em>/<em>n</em><sup><em>1</em>/<em>2</em></sup> denote the standard
deviation of <em>X&#772;<sub>n</sub></em>.

We wish to test the null hypothesis <em>H<sub>0</sub></em>: <em>&mu;</em> = <em>&mu;<sub>0</sub></em> against 
<em>H<sub>1</sub></em>: <em>&mu;</em> = <em>&mu;<sub>1</sub></em> &gt; <em>&mu;<sub>0</sub></em> at a significance 
level of <em>&alpha;</em>.

Let &Phi; be the cumulative distribution function of the standard normal.
Since (<em>X&#772;<sub>n</sub></em> - <em>&mu;</em>)/<em>&tau;<sub>n</sub></em>
is approximately standard normal, we have

<em>R</em> = [<em>&mu;<sub>0</sub></em> + &Phi;<sup>-<em>1</em></sup>(<em>1</em> - <em>&alpha;</em>)/<em>&tau;<sub>n</sub></em>, &infin;)

The measure of this interval under <em>H<sub>1</sub></em> gives the power

Power = <em>1</em> - &Phi;((<em>&mu;<sub>0</sub></em> - <em>&mu;<sub>1</sub></em>)/<em>&tau;<sub>n</sub></em> + &Phi;<sup>-<em>1</em></sup>(<em>1</em> - <em>&alpha;</em>))
