import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import convolve
import warnings
warnings.filterwarnings('ignore')

# Model parameters
c = 5.0      # wave velocity
rho = 1.0    # density
h = 5        # source depth
r = 10       # receiver distance
z = h

# Time window
dt = 0.01
n = 700
nbefore = 100

# Constants
scale = -1 / (4 * np.pi * rho * c**2)
eps = 1e-8  # small imaginary part for causality

# Calculate geometry
R = np.sqrt(r**2 + h**2)
t0 = R / c

print(f"R = {R}, t0 = {t0}")

# Time vector: t = t0 + dt*(-nbefore:n-nbefore-1)' + 0.5*dt
t = t0 + dt * np.arange(-nbefore, n - nbefore) + 0.5 * dt

# Compute p (slowness) - this is complex
# p = (r/R^2)*t - (z/R^2)*sign(t0-t).*sqrt(t0^2-t.^2) + eps*i
p_real = (r / R**2) * t
p_imag = -(z / R**2) * np.sign(t0 - t) * np.sqrt(np.abs(t0**2 - t**2))
# Note: In MATLAB, sqrt of negative gives imaginary, here we handle it differently
# Actually we should use complex sqrt: sqrt(t0^2 - t^2 + 0i)
p_complex = (r / R**2) * t - (z / R**2) * np.sign(t0 - t) * np.sqrt(t0**2 - t**2 + 0j) + eps * 1j
p = p_complex

# Compute eta = sqrt(1/c^2 - p.^2)
eta = np.sqrt(1 / c**2 - p**2)

# Compute dpdt = sign(t0-t).*eta./sqrt(t0^2-t.^2)
# Handle division by zero near t = t0
sqrt_term = np.sqrt(t0**2 - t**2 + 0j)
# Avoid division by zero
#sqrt_term_safe = np.where(np.abs(sqrt_term) < 1e-12, 1e-12 + 0j, sqrt_term)

dpdt = np.sign(t0 - t) * eta / sqrt_term

# Compute jt = imag((sqrt(p)./eta).*dpdt)
jt = np.imag((np.sqrt(p) / eta) * dpdt)

# Replace NaN/Inf
#jt = np.nan_to_num(jt, nan=0.0, posinf=0.0, neginf=0.0)

# Compute a = 2*sqrt(t) for convolution
a = 2 * np.sqrt(np.arange(0, n) * dt)

# Convolution: conv(a, jt)
conv_result = np.convolve(a, jt, mode='full')

# Derivative: diff(conv_result)
diff_conv = np.diff(conv_result)

# Apply scaling
phi = (scale * np.sqrt(2 / r) / np.pi) * diff_conv

# Take first n-2 points (as in MATLAB)
phi_trimmed = phi[:n-2]
t_phi = t[1:n-1]

# Also compute the original p for plotting (real and imag)
p_for_plot = p

# Plotting
fig, axes = plt.subplots(3, 1, figsize=(12, 10))

# Subplot 1: plot(p,'+') - plot real and imaginary parts
axes[0].plot(np.real(p_for_plot),np.imag(p_for_plot), 'b+', label='Real(p)', markersize=3)
#axes[0].plot(np.imag(p_for_plot), 'r+', label='Imag(p)', markersize=3)
axes[0].set_xlabel('Index')
axes[0].set_ylabel('p')
axes[0].set_title('p (slowness)')
axes[0].legend()
axes[0].grid(True, alpha=0.3)

# Subplot 2: plot(t, jt)
axes[1].plot(t, jt, 'b-', linewidth=1.5)
axes[1].set_xlabel('Time (s)')
axes[1].set_ylabel('jt')
axes[1].set_title('j(t)')
axes[1].grid(True, alpha=0.3)

# Subplot 3: plot(t(2:n-1), phi(1:n-2))
axes[2].plot(t_phi, phi_trimmed, 'b-', linewidth=1.5)
axes[2].set_xlabel('Time (s)')
axes[2].set_ylabel('phi')
axes[2].set_title('phi (pressure)')
axes[2].grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

# Print diagnostics
print(f"\nDiagnostics:")
print(f"t shape: {t.shape}")
print(f"p shape: {p.shape}")
print(f"jt shape: {jt.shape}")
print(f"a shape: {a.shape}")
print(f"conv_result shape: {conv_result.shape}")
print(f"phi_trimmed shape: {phi_trimmed.shape}")
print(f"t_phi shape: {t_phi.shape}")
print(f"\njt statistics:")
print(f"  min: {np.min(jt):.6f}")
print(f"  max: {np.max(jt):.6f}")
print(f"  mean: {np.mean(jt):.6f}")
print(f"  std: {np.std(jt):.6f}")
print(f"\nphi statistics:")
print(f"  min: {np.min(phi_trimmed):.6f}")
print(f"  max: {np.max(phi_trimmed):.6f}")
print(f"  mean: {np.mean(phi_trimmed):.6f}")
print(f"  std: {np.std(phi_trimmed):.6f}")