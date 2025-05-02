import numpy as np
import pandas as pd

# Configuration
num_points = 1000  # Number of data points
mu, sigma = 0, 0.5  # Lognormal parameters (for ln(X))
alpha = 0.2         # Power law exponent
log_param = 1.5     # Logarithmic parameter

# Generate data
ranks = np.arange(1, num_points + 1)

# Lognormal distribution (probability density)
x = np.linspace(0.01, 10, num_points)
lognormal = (1 / (x * sigma * np.sqrt(2 * np.pi))) * \
            np.exp(-(np.log(x) - mu)**2 / (2 * sigma**2))

# Power law distribution (Grinder's implementation)
powerlaw = ranks ** -alpha

# Logarithmic distribution (Grinder's implementation)
logarithmic = 1 / np.log(ranks + log_param)

# Create DataFrame and save to CSV
df = pd.DataFrame({
    'Rank': ranks,
    'Lognormal': lognormal,
    'Powerlaw_0.2': powerlaw,
    'Logarithmic_1.5': logarithmic
})

# Normalize for better visual comparison
df['Lognormal'] = df['Lognormal'] / df['Lognormal'].max()
df['Powerlaw_0.2'] = df['Powerlaw_0.2'] / df['Powerlaw_0.2'].max()
df['Logarithmic_1.5'] = df['Logarithmic_1.5'] / df['Logarithmic_1.5'].max()

df.to_csv('distribution_comparison.csv', index=False)
