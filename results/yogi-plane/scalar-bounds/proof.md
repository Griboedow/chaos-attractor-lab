# Boundedness of scalar adaptive momentum with a positive denominator floor

This note gives a direct boundedness argument for the scalar quadratic experiments in this repository. It covers every `0<=beta1<1`, rather than requiring a comparison with the second-memory parameter. The conclusion concerns finite forward orbits in exact real arithmetic. It does not establish convergence, chaos, a common global attractor, or historical priority of the argument.

## An abstract scalar statement

Let `e>0`, `sigma>0`, `0<=b<1`, and let finite real sequences satisfy

$$
u_{t+1}=a_tu_t+(1-a_t)x_t,\qquad
x_{t+1}=x_t-u_{t+1}/d_t,
$$

where, at every step,

$$
0\le a_t\le b,\qquad d_t\ge e,\qquad d_t\ge\sigma|x_t|.
$$

The denominator may depend on the entire preceding history. No monotonicity, continuity, or upper bound for it is needed. For any finite initial pair `(x0,u0)`, both sequences remain bounded.

Here is one explicit, deliberately loose bound. Put

$$
\kappa=\frac{\min(e,1)}8,
\quad L\in\mathbb N_{>0}\text{ with }b^L\le\kappa,
\quad K=\frac1{\sigma\kappa},
\quad R_0=\max(|x_0|,|u_0|),
$$

and choose

$$
B=\max\left\{2(1+1/e)^L R_0,
\frac{8LK}{e},\;4(L+1)K\right\}.
$$

Then `|x_t|<=B` and `|u_t|<=B` for all `t>=0`. For `b=0`, take `L=1`. For small `e` and `b` near one this bound is enormous; it is an existence bound, not an estimate of the pictured cloud's size.

## Proof by excluding the first exit

Write `S=kappa*B`. Suppose `N` is the first index for which `|x_N|>B`. Up to this proposed exit, the convex momentum recurrence gives `|u_t|<=B`, including `u_N`.

First, `R_t=max(|x_t|,|u_t|)` satisfies `R_(t+1)<=(1+1/e)R_t`. Consequently `|x_t|<=B/2` for `0<=t<=L`, and `N>L`.

Second, whenever `t<N` and `|x_t|>S`, the displacement obeys

$$
|x_{t+1}-x_t|\le\frac{B}{\sigma S}=K. \tag{1}
$$

Third, a state with `|x_t|<=S` cannot be followed by a state of magnitude greater than `B/2`. For `t<L`, the preceding initial bound already proves this. For `t>=L-1`, work backward through the last `L` gradient inputs. At each backward step, either the predecessor is inside `[-S,S]`, or (1) limits its distance from its successor by `K`. Thus every one of these inputs has magnitude at most `S+LK`. Unrolling the convex recurrence for `L` steps leaves an old-momentum weight

$$
q=\prod_{j=t-L+1}^{t}a_j\le b^L\le\kappa.
$$

The other nonnegative weights sum to `1-q`. It follows that

$$
|u_{t+1}|\le qB+S+LK\le2\kappa B+LK,
$$

and hence

$$
|x_{t+1}|\le\kappa B+\frac{2\kappa B+LK}{e}
\le\frac38B+\frac18B=\frac12B. \tag{2}
$$

Here `kappa*(1+2/e)<=3/8` follows separately for `e<=1` and `e>=1`, and `LK/e<=B/8` follows from the definition of `B`.

Consider a positive first exit, `x_N>B`. Equation (2) excludes `|x_(N-1)|<=S`; equation (1) therefore gives `x_(N-1)>B-K`. Its predecessor cannot lie inside `[-S,S]`, since that would force `|x_(N-1)|<=B/2` by (2). Repeating this argument backward gives

$$
x_{N-j}>B-jK\ge\frac34 B,
\qquad j=1,\ldots,L.
$$

All indices are nonnegative because `N>L`. Unrolling momentum over these `L` positive inputs yields

$$
u_N\ge(1-q)\frac34B-qB
\ge\left(\frac34-\frac74\kappa\right)B>0.
$$

The update then gives `x_N=x_(N-1)-u_N/d_(N-1)<x_(N-1)<=B`, a contradiction. Replacing `(x,u)` by `(-x,-u)` excludes a negative first exit by the same argument. Thus no first exit exists. Convexity gives the asserted bound for `u` as well.

## Application to the three scalar quadratic optimizers

Use the normalized coordinates of the linked experiment reports, so that the scalar gradient is `x`, and the first moment obeys `M=b*m+(1-b)*x`. All initial positions and first memories may be arbitrary finite real numbers; the second memory must be finite and nonnegative.

| Update | Denominator `d` | A valid `sigma` | Bound for the second memory after `|x|<=B` |
| --- | --- | --- | --- |
| [Yogi](yogi-structure.md) | `sqrt(V)+e` | `sqrt(min(c,1-c))`, `0<c<1` | `v_t<=max(v0,(2-c)*B^2)` |
| Adam | `sqrt(V)+e`, `V=c*v+(1-c)*x^2` | `sqrt(1-c)`, `0<=c<1` | `v_t<=max(v0,B^2)` |
| [Adamax](adamax-structure.md) | `U=max(c*u,abs(x)+e)` | `1`, `0<=c<1` | `u_t<=max(u0,B+e)` |

In the abstract theorem, `u` denotes the first moment; the conventional Adamax symbol `u` in the table denotes its maximum accumulator. They are different variables.

For Yogi, if `v>x^2`, then `V=v-(1-c)x^2>=c*x^2`; if `v<x^2`, then `V=v+(1-c)x^2>=(1-c)*x^2`; equality gives `V=x^2`. This proves the required current-gradient floor, including the switch itself. When `v>x^2` the accumulator decreases; otherwise its next value is at most `(2-c)*B^2`. The other two rows follow immediately from their recurrences.

Thus **every forward orbit of each stated autonomous scalar map is bounded** for positive `e` and the stated memory parameters. The result applies in particular to the high-first-memory Yogi cases and to Adamax case 4099, which were not covered by this repository's earlier, sharper restricted-parameter bounds. It also makes the finite-sheet omega-limit conclusion in the Adamax structure note unconditional on a separately assumed orbit bound.

### Continuously bias-corrected moments

For update indices `t>=1`, define `mhat_t=m_t/(1-b^t)`. Then

$$
\widehat m_{t+1}=a_t\widehat m_t+(1-a_t)x_t,
\qquad a_t=\frac{b(1-b^t)}{1-b^{t+1}}\in[0,b].
$$

For Adam and Yogi, correcting the nonnegative second memory replaces `sqrt(V)` by `sqrt(V/(1-c^(t+1)))`, which preserves both denominator lower bounds. Adamax corrects only the first moment. Starting the theorem after the first finite update therefore also bounds these continuously corrected scalar iterations, including nonzero initial memories. This statement refers to these exact denominator conventions; implementations that place epsilon elsewhere require a separate check.

## Relation to earlier results and limits of this statement

[Fong and Yang, arXiv:2608.20638v2](https://arxiv.org/abs/2608.20638v2), Lemma H.3, prove a momentum/second-moment estimate and trajectory bounds for diagonal-quadratic Adam under `beta1^2<beta2`. [Dereich, Do and Jentzen, arXiv:2603.18899v2](https://arxiv.org/abs/2603.18899v2), Corollary 2.7, give much broader objective and input-sequence bounds, with `beta2>beta1^2` and nonincreasing learning rates tending to zero. Neither stated result directly covers all parameter choices and constant-step Yogi recurrences considered here. This comparison is about hypotheses, not a claim that a scalar first-exit argument has never appeared elsewhere.

For the no-momentum RMSProp case, [Dimitrieski et al., arXiv:2603.15823v2](https://arxiv.org/abs/2603.15823v2) already establish global stability and boundedness for arbitrary bounded step-size rules under their objective assumptions. Boundedness of adaptive optimization is therefore an existing research subject, not a new setting introduced by this repository.

The proof above uses the scalar sign and convex averaging in an essential way. It is not a theorem for a general coupled quadratic objective or a stochastic gradient sequence. Positive `e` is essential to this proof. The bound depends on the initial state, and Yogi has an unbounded line of fixed points, so no common compact global attractor follows. Moreover, Yogi is discontinuous: compactness of an orbit's accumulation set alone does not establish its invariance under the map. None of these boundedness statements proves positive asymptotic Lyapunov exponents, an attracting chaotic invariant set, or a new topological type.
