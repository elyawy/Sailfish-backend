"""
Tests for auto-discrete-gamma correlation model implementation.
"""

import numpy as np
from typing import List
from msasim.correlation import (
    build_auto_gamma_transition_matrix,
    calculate_discrete_gamma_correlation,
)


def test_transition_matrix_properties():
    """Test A: Verify mathematical properties of the transition matrix."""
    
    alpha = 1.0
    K = 4
    rho = 0.5
    
    M = build_auto_gamma_transition_matrix(alpha, K, rho)
    
    # Property 1: Rows sum to 1 (valid probability distribution)
    row_sums = M.sum(axis=1)
    assert np.allclose(row_sums, 1.0, atol=1e-6), \
        f"Rows don't sum to 1: {row_sums}"
    
    # Property 2: All entries are non-negative
    assert np.all(M >= 0), "Transition matrix has negative entries"
    
    # Property 3: Matrix is symmetric (special property of this model)
    # M[i,j] should equal M[j,i] due to symmetry of bivariate normal
    assert np.allclose(M, M.T, atol=1e-6), \
        "Transition matrix is not symmetric"
    
    # Property 4: Stationary distribution is uniform
    # π @ M = π where π = [1/K, 1/K, ..., 1/K]
    stationary = np.ones(K) / K
    result = stationary @ M
    assert np.allclose(result, stationary, atol=1e-6), \
        f"Stationary distribution not preserved: {result} vs {stationary}"
    
    print("✓ All matrix properties verified")


def test_boundary_case_rho_zero():
    """Test that rho=0 gives independent rates (uniform transition matrix)."""
    
    alpha = 1.0
    K = 4
    rho = 0.0
    
    M = build_auto_gamma_transition_matrix(alpha, K, rho)
    
    # When rho=0, all transitions should be equally likely
    # M[i,j] = 1/K for all i,j (independence)
    expected = np.ones((K, K)) / K
    
    assert np.allclose(M, expected, atol=0.01), \
        f"rho=0 should give uniform transition matrix"
    
    # ρ_dG should also be 0
    rho_dG = calculate_discrete_gamma_correlation(M, alpha, K)
    assert abs(rho_dG) < 0.01, \
        f"rho=0 should give rho_dG≈0, got {rho_dG}"
    
    print("✓ Boundary case rho=0 verified")


def test_correlation_monotonicity():
    """Test that ρ_dG increases monotonically with ρ."""
    
    alpha = 1.0
    K = 4
    rho_values = np.linspace(-0.9, 0.9, 10)
    
    rho_dG_values = []
    for rho in rho_values:
        M = build_auto_gamma_transition_matrix(alpha, K, rho)
        rho_dG = calculate_discrete_gamma_correlation(M, alpha, K)
        rho_dG_values.append(rho_dG)
    
    # Check monotonicity
    diffs = np.diff(rho_dG_values)
    assert np.all(diffs > -1e-6), \
        "ρ_dG is not monotonically increasing with ρ"
    
    print("✓ Monotonicity verified")


def test_correlation_bounds():
    """Test that |ρ_dG| < |ρ| for most cases (transformation attenuates)."""
    
    alpha = 0.5
    K = 8
    
    test_cases = [0.3, 0.5, 0.7, -0.3, -0.5, -0.7]
    
    for rho in test_cases:
        M = build_auto_gamma_transition_matrix(alpha, K, rho)
        rho_dG = calculate_discrete_gamma_correlation(M, alpha, K)
        
        # For small alpha, transformation typically reduces correlation magnitude
        # This is a qualitative test based on Figure 3
        assert abs(rho_dG) <= abs(rho) + 0.1, \
            f"For α={alpha}, rho={rho}: |ρ_dG|={abs(rho_dG)} > |ρ|={abs(rho)}"
    
    print("✓ Correlation bounds verified")


def test_large_alpha_approaches_diagonal():
    """Test that as α→∞, ρ_dG → ρ (less transformation)."""
    
    K = 8
    rho = 0.5
    
    alpha_values = [1.0, 2.0, 5.0, 10.0]
    rho_dG_values = []
    
    for alpha in alpha_values:
        M = build_auto_gamma_transition_matrix(alpha, K, rho)
        rho_dG = calculate_discrete_gamma_correlation(M, alpha, K)
        rho_dG_values.append(rho_dG)
    
    # ρ_dG should approach ρ as alpha increases
    assert rho_dG_values[-1] > rho_dG_values[0], \
        "ρ_dG should increase with α"
    
    # For large alpha, should be close to rho
    assert abs(rho_dG_values[-1] - rho) < 0.1, \
        f"For large α, ρ_dG should approach ρ: {rho_dG_values[-1]} vs {rho}"
    
    print("✓ Large alpha behavior verified")


def reproduce_yang_figure3(
    alpha_values: List[float] = [0.1, 0.3, 0.5, 1.0, 5.0],
    rho_range: tuple = (-1, 1),
    n_points: int = 30,
    K: int = 8,
    save_path: str = "yang_figure3_reproduction.png"
):
    """
    Reproduce Yang (1995) Figure 3: ρ_dG vs ρ for different α values.
    
    This validates the implementation by comparing to published results.
    
    Args:
        alpha_values: List of gamma shape parameters to plot
        rho_range: Range of ρ values to plot
        n_points: Number of points to compute
        K: Number of rate categories
        save_path: Where to save the figure
    """
    import matplotlib.pyplot as plt
    
    # Avoid exact boundaries where covariance matrix becomes singular
    rho_min = max(rho_range[0], -0.99)
    rho_max = min(rho_range[1], 0.99)
    rho_values = np.linspace(rho_min, rho_max, n_points)
    
    plt.figure(figsize=(10, 7))
    
    for alpha in alpha_values:
        rho_dG_values = []
        
        print(f"Computing for α = {alpha}...")
        
        for rho in rho_values:
            # Build transition matrix
            M = build_auto_gamma_transition_matrix(alpha, K, rho)
            
            # Calculate discrete gamma correlation
            rho_dG = calculate_discrete_gamma_correlation(M, alpha, K)
            rho_dG_values.append(rho_dG)
        
        plt.plot(rho_values, rho_dG_values, label=f'α = {alpha}', linewidth=2)
    
    plt.xlabel('ρ', fontsize=14)
    plt.ylabel('ρ_dG', fontsize=14)
    plt.legend(fontsize=12, loc='best')
    plt.grid(True, alpha=0.3)
    plt.title('Yang (1995) Figure 3: Relationship between ρ and ρ_dG', fontsize=14)
    
    # Add reference lines
    plt.axhline(y=0, color='k', linestyle='-', linewidth=0.5, alpha=0.5)
    plt.axvline(x=0, color='k', linestyle='-', linewidth=0.5, alpha=0.5)
    
    # Add diagonal reference (ρ_dG = ρ)
    plt.plot(rho_values, rho_values, 'k--', alpha=0.3, linewidth=1, label='ρ_dG = ρ')
    
    plt.xlim(-1, 1)
    plt.ylim(-1, 1)
    plt.legend(fontsize=11, loc='best')
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    print(f"\nFigure saved to: {save_path}")
    plt.show()


def test_reproduce_figure3_visual():
    """
    Generate Yang (1995) Figure 3 for visual comparison.
    This is not an automated test - requires manual inspection.
    """
    print("\nGenerating Yang (1995) Figure 3 reproduction...")
    print("This will take a few minutes...")
    
    reproduce_yang_figure3(
        alpha_values=[0.1, 0.3, 0.5, 1.0, 5.0],
        rho_range=(-1.0, 1.0),
        n_points=20,  # Fewer points for faster testing
        K=8,
        save_path="test_outputs/test_figure3.png"
    )
    
    print("\n✓ Figure generated successfully")
    print("  Visual inspection required: Compare with Yang (1995) Figure 3")
    print("  Expected properties:")
    print("    1. All curves pass through origin")
    print("    2. Larger α → closer to diagonal")
    print("    3. Smaller α → more nonlinearity")
    print("    4. All curves monotonically increasing")


if __name__ == "__main__":
    print("Running Auto-Discrete-Gamma Correlation Tests\n")
    print("=" * 60)
    
    test_transition_matrix_properties()
    test_boundary_case_rho_zero()
    test_correlation_monotonicity()
    test_correlation_bounds()
    test_large_alpha_approaches_diagonal()
    test_reproduce_figure3_visual()
    
    print("\n" + "=" * 60)
    print("All tests passed! ✓")
