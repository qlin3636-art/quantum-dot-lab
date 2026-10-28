# DOT · Quantum Playground v2

Bilingual, browser-local quantum-dot laboratory. No AI inference service.

## Experiments
- Quantum-dot emission using a simplified spherical effective-mass model.
- 2D potential studio: presets, mathematical V(x,y), painted barriers, four low eigenstates and probability maps.
- QD electroluminescence: carrier injection balance, phenomenological losses, relative photon output and approximate emission colour. Replaces the original RGB drawing studio.

## Solver
Five-point finite differences on a 24x24 or 36x36 interior grid, x,y in (-8,8) nm; zero Dirichlet boundaries and constant effective mass. Block imaginary-time iteration with orthogonalization and Rayleigh-Ritz subspace rotation. Each state is normalized by dx^2. The maximum residual in eV and convergence status are shown. Finite grids, finite box size and incomplete convergence limit accuracy. Degenerate eigenvectors are not unique. This is a 2D single-particle teaching model, not a complete 3D exciton simulation.

## EL model
Dimensionless injection rates Je=D*se and Jh=D*sh. Balance gamma=2min(Je,Jh)/(Je+Jh); radiative fraction=1/(1+T+0.12D^2); EQE=gamma*radiative*0.20. All rate, loss and extraction coefficients are educational assumptions, not fitted experimental parameters. Real-device loss mechanisms vary. Source: https://www.nature.com/articles/s41467-024-52521-0

## Validation
Analytic discrete 2D box energies and degeneracy, normalization, orthogonality, eigen-equation residuals, a nonseparable xy potential and constant energy shifts; EL limits and no-injection darkness. Ten unit tests including the original solver tests pass.

## Deployment and rollback
Publish main branch root using GitHub Pages. Original v1 commit: 1b8772b865f602a8a229c6985007ffc35943b0fc. Revert the v2 commit through GitHub to restore v1. A local ZIP of the complete original source was also saved before editing.
