# DOT · Quantum Playground

A bilingual interactive quantum-dot laboratory for MECH 6045, built through vibe coding.

## Explore

- Quantum-dot emission: radius, material, approximate optical gap and spectrum.
- Potential studio: draw or enter V(x), solve the 1D stationary Schrodinger equation, inspect the five lowest states, and sample positions.
- Pixel studio: RGB mixing, 16 x 16 painting, subpixel inspection and PNG export.
- Field notes: selected original research and model assumptions.

All calculations run in the browser. No AI API or backend is required. Experiment state is saved only in browser local storage.

## Deploy

Publish the main branch root with GitHub Pages. index.html must remain at the repository root. The .nojekyll file disables Jekyll processing.

## Scientific scope

The potential solver uses 256 interior grid points, central differences, constant effective mass and zero Dirichlet boundaries. The separate light lab uses a simplified spherical Brus-type model, including the bulk band gap, electron/hole confinement and Coulomb attraction. The 1D E2-E1 spacing is not used as quantum-dot visible recombination energy. RGB artwork is an illustrative display model. Detailed limitations and references appear on the website.
