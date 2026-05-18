import numpy as np
from numpy.linalg import norm
from scipy.stats import unitary_group

def random_entangled_choi(dim_in=2, dim_out=2):
    """
    Generate a random PSD Choi matrix that is entangled.
    
    Args:
        dim_in (int): Input dimension of the channel.
        dim_out (int): Output dimension of the channel.
    
    Returns:
        np.ndarray: PSD Choi matrix (dim_in*dim_out x dim_in*dim_out).
    """
    d_total = dim_in * dim_out

    # Step 1: Generate a random unitary
    U = unitary_group.rvs(d_total)

    # Step 2: Apply it to a maximally entangled state
    # |Φ+> = sum_i |i>_in ⊗ |i>_out / sqrt(dim_in)
    phi_plus = np.zeros((dim_in, dim_out), dtype=complex)
    for i in range(min(dim_in, dim_out)):
        phi_plus[i, i] = 1
    phi_plus = phi_plus.flatten() / np.sqrt(dim_in)

    # Apply random unitary to create an entangled state
    psi = U @ phi_plus

    # Step 3: Create density matrix (pure state projector)
    rho = np.outer(psi, np.conjugate(psi))

    # Step 4: Ensure PSD (it is, since it's a projector)
    # Interpret rho as Choi matrix
    choi = rho

    # Step 5: Normalize trace to dim_in (Choi normalization)
    choi *= dim_in / np.trace(choi)

    return choi

# Example usage
choi_matrix = random_entangled_choi(2, 2)

print(choi_matrix)

# Check PSD property
eigvals = np.linalg.eigvalsh(choi_matrix)
print("Eigenvalues:", eigvals)
print("Is PSD?", np.all(eigvals >= -1e-12))

# Check entanglement via partial transpose (Peres-Horodecki criterion)
def partial_transpose(rho, dim_in, dim_out):
    rho_reshaped = rho.reshape(dim_in, dim_out, dim_in, dim_out)
    rho_pt = rho_reshaped.swapaxes(1, 2).reshape(dim_in*dim_out, dim_in*dim_out)
    return rho_pt

rho_pt = partial_transpose(choi_matrix, 2, 2)
eigvals_pt = np.linalg.eigvalsh(rho_pt)
print("Partial transpose eigenvalues:", eigvals_pt)
print("Is entangled?", np.any(eigvals_pt < -1e-12))