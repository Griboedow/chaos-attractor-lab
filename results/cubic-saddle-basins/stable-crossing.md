# A stable arc crosses the passage box

A connected arc of the stable manifold of $`z=1`$ crosses both horizontal sides of the [phase-zero passage box](../cubic-saddle-passages/). The box half-width is 0.14. In its own coordinates the entire arc has $`|x|<0.04530`$, while its end heights satisfy $`y_-<-0.15471`$ and $`y_+>0.15471`$. Truncating it at the two sides gives a vertical crossing. [Interval bounds](stable-crossing.json)

## Building the arc in reverse time

Let $`f=P^{-1}`$ in the centered coordinate $`w=z-1`$. Its expanding eigenpair is $`(\lambda,v)`$, with $`|v|=1`$ and $`\lambda<-2^{63}`$. Use the [forward proof's local estimates](README.md#local-bounds) with the additional term $`\kappa(6\delta+3\delta^2)^2`$ in the reversed logarithmic-norm bound. The same constants $`H=2^{78}`$ and $`C=2^{164}`$ work on $`|w|\le2^{-84}`$.

Set $`s=1/(5\lambda)`$. For $`|a|\le |s|`$ consider $`\psi_n(a)=f^n(a\lambda^{-n}v)`$. The one-step quadratic bound and $`n`$ applications of the Lipschitz bound give

```math
|\psi_{n+1}(a)-\psi_n(a)|\le \frac{C|a|^2H^n}{|\lambda|^{2n+2}}.
```

The error induction from the [three-exit proof](README.md#three-exits-at-arbitrarily-small-scales) keeps every preterminal state in the local disk. In particular, $`|s/\lambda|+K|s/\lambda|^2<2^{-84}`$, with $`K=C/(2^{126}-H)`$. The geometric series gives a continuous limit $`\psi`$ within $`K|a|^2`$ of $`av`$. Since $`f(\psi(a/\lambda))=\psi(a)`$, each point converges to the saddle in forward time. The curve is injective: equality at two parameters implies equality at the rescaled parameters $`a/\lambda^n,b/\lambda^n`$ by uniqueness of the flow. Dividing their quadratic error bound by $`|\lambda|^{-n}`$ and letting $`n\to\infty`$ gives $`a=b`$.

## Crossing the sides

Apply $`f`$ once more to the whole arc $`\psi(ts)`$, $`-1\le t\le1`$. Use $`tV`$ as the nominal path, where $`V`$ follows the reversed linear equation from $`sv`$. In a time cell with $`|V|\le r`$, the uniform nonlinear contribution to the defect per unit time is at most

```math
\kappa\left[3r^2|b+2|+(6r+3r^2)(3+r)r^2\right].
```

This follows by expanding the reversed field $`\kappa\overline{3w(2+w)}[-b-2+(3+w)w^2]`$. The checker propagates the Taylor defect and this remainder with a logarithmic norm on the enclosing tube. The final Euclidean error is below 0.03369 for the entire parameter interval. The certified horizontal clearance and the two end heights above use the exact rational frame of the passage box.

```sh
python research/shape_search_20261121/frontiers/cubic/saddle_basins.py --stable-crossing --out stable-crossing.json
```
