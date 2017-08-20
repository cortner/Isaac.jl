
# Notes on Preconditioned Arnoldi

This is a brief note, mostly for myself as a record, on how the preconditioned Arnoldi
process is implemented.

## Standard Arnoldi

Let $A \in \mathbb{R}^{N \times N}$, for simplicity assume diagonalisable.
Then by Arnoldi iteration we generate an ONB $V = (v_1, \dots, v_M)$ as
well as
$$
  y_i = A \cdot v_i
$$
and by construction $y_i \in {\rm span} \{ v_1, \dots, v_{i-1} \}$ hence we
get a orthogonal projection of $A$ onto $V$ of the form
$$
  A \approx V \cdot H \cdot V^\top
$$
where $H$ is upper hessenberg. The eigenvalues of $H$ are the Ritz values
of $A$ and (are understood to) converge to the extremal eigenvalues of $A$.

## Preconditioned Arnoldi

If we precondition this process then we have
$$
  y_i = P\^{-1} A v_i
$$
The Krylov vectors $v_i$ may still be chosen $\ell^2$-orthogonal, $V^T V = I$;
but we may want to consider a generalisation to $V^T M V = I$. (TODO)

A linear system $Au = b$ can still be solved as before by just preconditioning
$b$ as well.

The question of eigenvalues becomes more subtle: we now have
$$
  P\^{-1} A \cdot V = Y = V \cdot H
$$
hence the eigenvalues of $H$ are now the Ritz-values of $P^{-1} A$. If these
are in fact the eigenvalues we are after (this depends on the context) then
we can simply diagonalise
$$
  H = Z \Lambda Z\^{-1}, \qquad Z = [z_1 \dots z_M],
$$
(assuming this can in fact be done) and setting $w_j = V z_j$ we obtain
$$
  A w_j = A V z_j= P V H z_j = \lambda_j P V z_j = \lambda_j P w_j;
$$
i.e. $w_j$ are the eigen-vectory approximations corresponding to the
Ritz values $\lambda_j$ for the generalised eigenvalue problem $A w = \lambda P w$.

If we insist that on standard eigenvalues, then we need to solve the generalised
eigenvalue problem
$$
  P_V H z_j = \lambda_j z_j, \qquad \text{where} \quad P_V = V^\top P V.
$$
Defining again $w_j = V z_j$ we obtain
$$
  V^\top A w_j = V^\top P V H z_j = \lambda_j z_j = \lambda_j V^\top w_j;
$$
that is, $\lambda_j$ are the Ritz values of $A$ (rather than those of $P^{-1} A$)
and $w_j$ the corresponding approximate eigenvectors.


## With metric

We can consider orthogonalising the Krylov vectors w.r.t. an IP $u^T M v$ where
$M$ is spd. It is important to stress, however, that no matter how we orthogonalise
we always get the same Krylov subspace ($P = I$ for simplicity)
$$
  K_M = {\rm span} \\{ v_0, A v_0, A\^2 v_0, \dots, A\^{M-1} v_0\\}
$$

A potential advantage of general $M$ for the case when $A$ is
self-adjoint (hermitian), is that choosing $M = P$ leads to $H$ being hermitian
as well, hence it is tri-diagonal instead of diagonal.

In general it is unclear however whether this has any advantage other than
possibly numerical stability, and it would come at the cost of a more
complex implementation?
