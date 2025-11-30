"""
msasim/correlation.py

Auto-discrete-gamma correlation model (Yang 1995) for correlated substitution rates.
"""

import numpy as np
import warnings

try:
    from scipy import stats
    from scipy.integrate import quad
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False


def build_auto_gamma_transition_matrix(
    alpha: float,
    categories: int,
    rho: float
) -> np.ndarray:
    """
    Build transition matrix for auto-discrete-gamma model using copula method.
    
    Implements Yang (1995) equations for correlated rate categories.
    Uses bivariate normal copula to generate correlated gamma variates.
    
    Args:
        alpha: Shape parameter of gamma distribution
        categories: Number of discrete rate categories (K)
        rho: Correlation parameter for bivariate normal (-1 to 1)
            This is the input parameter ρ.
            The realized discrete correlation ρ_dG will differ somewhat.
        
    Returns:
        K×K transition matrix M where M[i,j] = P(category j at next site | category i at current site)
        
    Raises:
        ImportError: If scipy is not installed
        ValueError: If parameters are invalid
        
    Reference:
        Yang, Z. (1995). A space-time process model for the evolution of DNA sequences.
        Genetics, 139(2), 993-1005.
    """
    if not SCIPY_AVAILABLE:
        raise ImportError(
            "scipy is required for rate correlation. "
            "Install with: pip install scipy"
        )
    
    if categories < 2:
        raise ValueError(f"Need at least 2 categories, got {categories}")
    if not -1 <= rho <= 1:
        raise ValueError(f"Correlation rho must be in [-1, 1], got {rho}")
    if alpha <= 0:
        raise ValueError(f"Alpha must be positive, got {alpha}")
    
    K = categories
    
    # Handle perfect correlation as special cases (avoid numerical issues)
    if abs(abs(rho) - 1.0) < 1e-6:
        if rho > 0:
            # Perfect positive correlation: X2 = X1 exactly
            # Sites stay in same category: diagonal matrix
            M = np.eye(K)
        else:
            # Perfect negative correlation: X2 = -X1 exactly  
            # Sites flip to opposite category: anti-diagonal matrix
            M = np.eye(K)[:, ::-1]
        return M
    
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
            
            # Compute P(lower_i < X1 < upper_i, lower_j < X2 < upper_j)
            joint_probs[i, j] = _bivariate_normal_rectangle(
                mvn, lower_i, upper_i, lower_j, upper_j
            )
    
    # Step 5: Normalize to get conditional probabilities M[i,j] = P(j|i)
    # M[i,j] = P(i,j) / P(i) = P(i,j) / (1/K)
    row_sums = joint_probs.sum(axis=1, keepdims=True)
    
    # Handle potential numerical issues
    if np.any(row_sums < 1e-10):
        warnings.warn("Some row sums very small in transition matrix construction")
        row_sums = np.maximum(row_sums, 1e-10)
    
    M = joint_probs / row_sums
    
    return M


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
        
    Raises:
        ImportError: If scipy is not installed
    """
    if not SCIPY_AVAILABLE:
        raise ImportError("scipy is required for correlation calculation")
    
    # Compute conditional mean for each category
    r = np.zeros(K)
    
    for i in range(K):
        lower_p = i / K
        upper_p = (i + 1) / K
        
        if lower_p == 0:
            lower_q = 0
        else:
            lower_q = stats.gamma.ppf(lower_p, alpha, scale=1/alpha)
        
        upper_q = stats.gamma.ppf(upper_p, alpha, scale=1/alpha)
        
        def integrand(x):
            return x * stats.gamma.pdf(x, alpha, scale=1/alpha)
        
        numerator, _ = quad(integrand, lower_q, upper_q)
        r[i] = numerator / (1/K)
    
    # Yang's equation (8)
    numerator = 0.0
    for i in range(K):
        for j in range(K):
            numerator += (1/K) * M[i, j] * r[i] * r[j]
    numerator -= 1.0
    
    denominator = np.sum((1/K) * r**2) - 1.0
    
    if abs(denominator) < 1e-10:
        return 0.0  # No variance, no correlation
    
    return numerator / denominator