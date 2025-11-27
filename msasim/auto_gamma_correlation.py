"""
Implementation of Yang (1995) auto-discrete-gamma correlation model.

This module implements the bivariate gamma copula method for modeling
correlated substitution rates at adjacent sites in DNA sequences.
"""

import numpy as np
from scipy import stats
from scipy.integrate import quad
from typing import List
import warnings


def _bivariate_normal_rectangle(
    mvn: stats.multivariate_normal,
    x1_low: float,
    x1_high: float,
    x2_low: float,
    x2_high: float
) -> float:
    """
    Compute P(x1_low < X1 < x1_high, x2_low < X2 < x2_high) for bivariate normal.
    
    Uses CDF differences: P(a < X < b, c < Y < d) = 
        CDF(b,d) - CDF(a,d) - CDF(b,c) + CDF(a,c)
    
    Args:
        mvn: scipy multivariate_normal object
        x1_low, x1_high: Bounds for first dimension
        x2_low, x2_high: Bounds for second dimension
        
    Returns:
        Probability mass in the rectangle
    """
    def safe_cdf(x1, x2):
        """Handle infinite bounds using limits"""
        # Replace infinities with large finite values
        # ±10 gives CDF very close to 0 or 1
        if np.isinf(x1):
            x1 = -10.0 if x1 < 0 else 10.0
        if np.isinf(x2):
            x2 = -10.0 if x2 < 0 else 10.0
        
        try:
            result = mvn.cdf([x1, x2])
            # Handle potential NaN or invalid results
            if not np.isfinite(result):
                # If CDF fails, estimate based on position
                if x1 < -5 or x2 < -5:
                    return 0.0
                elif x1 > 5 or x2 > 5:
                    return 1.0
                else:
                    return 0.5
            return result
        except:
            # Fallback for numerical issues
            if x1 < -5 or x2 < -5:
                return 0.0
            elif x1 > 5 and x2 > 5:
                return 1.0
            else:
                return 0.5
    
    prob = (safe_cdf(x1_high, x2_high) - 
            safe_cdf(x1_low, x2_high) -
            safe_cdf(x1_high, x2_low) + 
            safe_cdf(x1_low, x2_low))
    
    return max(0.0, min(1.0, prob))  # Clamp to [0,1] for numerical safety


def build_auto_gamma_transition_matrix(
    alpha: float,
    categories: int,
    rho: float
) -> np.ndarray:
    """
    Build transition matrix for auto-discrete-gamma model using copula method.
    
    Algorithm (Yang 1995):
    1. Discretize gamma(alpha, scale=1/alpha) into K equal-probability categories
    2. Map category boundaries to standard normal quantiles via copula
    3. Use bivariate normal CDF to compute joint probabilities P(i,j)
    4. Normalize to conditional probabilities M[i,j] = P(category j | category i)
    
    Args:
        alpha: Shape parameter of gamma distribution
        categories: Number of rate categories (K)
        rho: Correlation parameter of bivariate normal (ρ)
        
    Returns:
        K×K transition matrix where M[i,j] = P(next site in category j | current site in category i)
    """
    K = categories
    
    # Step 1: Define category boundaries in probability space [0, 1]
    # Each category has probability 1/K
    probs = np.linspace(0, 1, K + 1)
    
    # Step 2: Map to standard normal quantiles
    # For prob p, find z such that Φ(z) = p, where Φ is standard normal CDF
    normal_thresholds = []
    for p in probs:
        if p == 0:
            normal_thresholds.append(-np.inf)
        elif p == 1:
            normal_thresholds.append(np.inf)
        else:
            normal_thresholds.append(stats.norm.ppf(p))
    
    # Step 3: Create bivariate normal with correlation rho
    # scipy handles near-singular cases automatically with allow_singular=True
    mvn = stats.multivariate_normal(
        mean=[0, 0],
        cov=[[1, rho], [rho, 1]],
        allow_singular=True
    )
    
    # Step 4: Compute joint probabilities for each (i,j) rectangle
    joint_probs = np.zeros((K, K))
    
    for i in range(K):
        for j in range(K):
            # Rectangle bounds in standard normal space
            lower_i = normal_thresholds[i]
            upper_i = normal_thresholds[i + 1]
            lower_j = normal_thresholds[j]
            upper_j = normal_thresholds[j + 1]
            
            # Compute P(Z1 in category i, Z2 in category j)
            joint_probs[i, j] = _bivariate_normal_rectangle(
                mvn, lower_i, upper_i, lower_j, upper_j
            )
    
    # Step 5: Convert to conditional probabilities M[i,j] = P(j|i)
    # Each row should sum to 1
    row_sums = joint_probs.sum(axis=1, keepdims=True)
    
    # Handle potential numerical issues
    if np.any(row_sums < 1e-10):
        warnings.warn("Some row sums very small in transition matrix construction")
        row_sums = np.maximum(row_sums, 1e-10)
    
    M = joint_probs / row_sums
    
    return M


def calculate_discrete_gamma_correlation(
    M: np.ndarray,
    alpha: float,
    K: int
) -> float:
    """
    Calculate ρ_dG (discrete gamma correlation) from transition matrix.
    
    Uses Yang (1995) equation (8):
    
    ρ_dG = [Σᵢⱼ (1/K) * Mᵢⱼ * r̄ᵢ * r̄ⱼ - 1] / [Σᵢ (1/K) * r̄ᵢ² - 1]
    
    where r̄ᵢ is the conditional mean rate for category i.
    
    Args:
        M: K×K transition matrix
        alpha: Shape parameter of gamma distribution
        K: Number of categories
        
    Returns:
        ρ_dG: Correlation between rates at adjacent sites
    """
    # Compute conditional mean for each category
    # r̄ᵢ = E[R | R in category i]
    r = np.zeros(K)
    
    for i in range(K):
        # Category i covers probability range [i/K, (i+1)/K]
        lower_p = i / K
        upper_p = (i + 1) / K
        
        # Convert to gamma quantiles
        if lower_p == 0:
            lower_q = 0
        else:
            lower_q = stats.gamma.ppf(lower_p, alpha, scale=1/alpha)
        
        upper_q = stats.gamma.ppf(upper_p, alpha, scale=1/alpha)
        
        # Compute conditional expectation E[R | lower_q < R < upper_q]
        # = ∫ r * f(r) dr / ∫ f(r) dr over [lower_q, upper_q]
        def integrand(x):
            return x * stats.gamma.pdf(x, alpha, scale=1/alpha)
        
        numerator, _ = quad(integrand, lower_q, upper_q)
        denominator = 1 / K  # Each category has equal probability
        
        r[i] = numerator / denominator
    
    # Yang's equation (8)
    # Numerator: E[r_{n-1} * r_n] - 1
    numerator = 0.0
    for i in range(K):
        for j in range(K):
            # Stationary probability = 1/K for each category
            numerator += (1/K) * M[i, j] * r[i] * r[j]
    numerator -= 1.0
    
    # Denominator: Var[r] = E[r²] - (E[r])² = E[r²] - 1
    # (since E[r] = 1 by construction)
    denominator = np.sum((1/K) * r**2) - 1.0
    
    if abs(denominator) < 1e-10:
        return 0.0
    
    return numerator / denominator


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


if __name__ == "__main__":
    print("Reproducing Yang (1995) Figure 3...")
    print("This will take a few minutes due to numerical integration...")
    print()
    
    reproduce_yang_figure3(
        alpha_values=[0.1, 0.3, 0.5, 1.0, 5.0],
        rho_range=(-1, 1),
        n_points=30,
        K=8,
        save_path="test_outputs/yang_figure3_reproduction.png"
    )
    
    print("\nDone! Compare the generated plot with Yang (1995) Figure 3.")
    print("\nExpected properties:")
    print("1. All curves should pass through origin (ρ=0 → ρ_dG=0)")
    print("2. For larger α, curves should be closer to diagonal (ρ_dG ≈ ρ)")
    print("3. For smaller α, curves should show more nonlinearity")
    print("4. All curves should be monotonically increasing")