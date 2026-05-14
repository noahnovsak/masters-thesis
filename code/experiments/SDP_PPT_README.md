# SDP for Generating PPT Entangled States

This directory contains an implementation of the approach from the paper:

**"A simple class of bound entangled states based on the properties of the antisymmetric subspace"**  
by Enrico Sindici and Marco Piani

## Notebook: `sdp_ppt_entangled_states.ipynb`

### Overview

The notebook implements the semidefinite programming (SDP) approach for generating positive partial transpose (PPT) entangled states. These are bound entangled states that cannot be distilled into pure entanglement, but are still entangled.

### Key Concepts

1. **Symmetric and Antisymmetric Subspaces**: The Hilbert space of two d-dimensional systems can be decomposed as:
   - Symmetric subspace: $\mathcal{H}_S = \mathbb{C}^d \vee \mathbb{C}^d$ (dimension d(d+1)/2)
   - Antisymmetric subspace: $\mathcal{H}_A = \mathbb{C}^d \wedge \mathbb{C}^d$ (dimension d(d-1)/2)

2. **Swap Operator**: The operator that exchanges subsystems A and B, with eigenprojectors onto these subspaces:
   - $P_S = (I + V)/2$ projects onto symmetric subspace
   - $P_A = (I - V)/2$ projects onto antisymmetric subspace

3. **SDP Problem**: For a given antisymmetric state $\rho_A$, solve:
   ```
   max Tr(P_A σ)
   s.t. P_A σ P_A = Tr(P_A σ) ρ_A
        σ ≥ 0, Tr(σ) = 1, σ^Γ ≥ 0 (PPT)
   ```

4. **Main Result**: If $p^{PPT}(\rho_A) < 1/2$, the optimal PPT state σ* is entangled (by Lemma 1 in the paper).

### Theoretical Bounds

For any antisymmetric state $\rho_A$:
$$\frac{2}{d(d+1)+2} \leq p^{PPT}(\rho_A) \leq \frac{1}{2}$$

The lower bound is achieved by states of the form:
$$\sigma(p) = p \rho_A \oplus (1-p) \frac{P_S}{d_S}$$

### Notebook Sections

1. **Helper Functions**: Swap operators, projectors, partial transpose
2. **Random State Generators**: Generate random quantum states and antisymmetric states
3. **SDP Formulation**: JuMP model for the optimization problem
4. **Bounds**: Implementation of theoretical bounds from Theorem 1
5. **Main Procedure**: Algorithm for generating PPT entangled states
6. **Examples**: Generate PPT entangled states for d=3
7. **Batch Generation**: Generate and analyze multiple states
8. **Verification**: Check that theoretical bounds hold

### Usage

The notebook requires the Julia environment with the following packages installed:
- `LinearAlgebra`, `Random` (standard library)
- `JuMP`, `MosekTools` (optimization)
- `ProgressMeter` (progress bars)
- `ppt2` (custom quantum info utilities)

Run the notebook cells sequentially to:
1. Set up the environment
2. Define helper functions
3. Generate PPT entangled states
4. Verify they satisfy PPT and entanglement properties
5. Collect statistics on multiple states

### Main Functions

- `swap_operator(d)`: Constructs the swap operator for d-dimensional systems
- `symmetric_projector(d)`: Projector onto symmetric subspace
- `antisymmetric_projector(d)`: Projector onto antisymmetric subspace
- `partial_transpose_A(rho, d_A, d_B)`: Partial transpose on subsystem A
- `is_ppt(rho, d_A, d_B)`: Check if state is PPT
- `is_entangled(rho, d_A, d_B)`: Check if state is entangled
- `rand_antisymmetric_state(d)`: Generate random antisymmetric state
- `solve_ppt_sdp(rho_A, d)`: Solve the SDP for a given antisymmetric state
- `generate_ppt_entangled_state(d)`: Generate a PPT entangled state (with retries)
- `batch_generate_ppt_entangled_states(d, n_samples)`: Generate multiple states and collect statistics

### Example Output

For d=3:
- Dimension of symmetric subspace: d_S = 6
- Dimension of antisymmetric subspace: d_A = 3
- Lower bound on p^PPT: 2/(3·4+2) ≈ 0.1429
- Upper bound on p^PPT: 0.5

The notebook generates PPT entangled states where the optimal probability p^PPT(ρ_A) is between these bounds.

### References

- Sindici, E., & Piani, M. (2018). "A simple class of bound entangled states based on the properties of the antisymmetric subspace." arXiv preprint arXiv:0902.1834.
- Horodecki, R., Horodecki, P., Horodecki, M., & Horodecki, K. (2009). "Quantum entanglement." Reviews of Modern Physics, 81(2), 865.
