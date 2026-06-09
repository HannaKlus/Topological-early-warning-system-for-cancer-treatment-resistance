from Model import model
import numpy as np
from gtda.time_series import SingleTakensEmbedding
from gtda.plotting import plot_diagram
import matplotlib.pyplot as plt

ess_vals = np.array([6060.60606060606,7575.75757575758,1.e-09,27272.7272727273])
c_m_plus = np.array([1.0, 0.7, 0.8])
c_m_P = np.array([0.4, 1.0, 0.5])
c_m_minus = np.array([0.6, 0.9, 1.0])

t_regular, psa_regular = model(ess_vals, c_m_plus, c_m_P, c_m_minus)

#Takens Embedding
embedder = SingleTakensEmbedding(parameters_type='search', time_delay=200, dimension = 5, n_jobs=-1)
point_cloud = embedder.fit_transform(psa_regular.reshape(-1, 1))
tau = embedder.time_delay_
d = embedder.dimension_
print(tau, d)

X_pc = point_cloud[np.newaxis, :, :]
X_pc = X_pc[:, ::4, :]

fig, ax = plt.subplots(figsize=(8, 8))
ax.plot(point_cloud[:, 0], point_cloud[:, 1],
        color='#102f51', linewidth=2, alpha=0.7)
ax.set_xlabel(f'PSA(t)', fontsize=20)
ax.set_ylabel(f'PSA(t + {tau})', fontsize=20)
ax.set_title('Phase Portrait', fontsize=20, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()

plt.savefig('phase_portrait.png')
plt.show()
