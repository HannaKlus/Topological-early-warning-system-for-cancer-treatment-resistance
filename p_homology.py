import numpy as np
from gtda.homology import VietorisRipsPersistence
from gtda.time_series import SlidingWindow
from gtda.time_series import SingleTakensEmbedding
from Model import model
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

#Sliding Window
step_size = 100
window_size = 1800
SW = SlidingWindow(size=window_size, stride=step_size)
X_windows = SW.fit_transform(point_cloud)

X_windows = X_windows[:, ::4, :]

#Persitent homology
VR = VietorisRipsPersistence(homology_dimensions=[1], n_jobs=-1)
diagrams = VR.fit_transform(X_windows)

h1_persistence = []
for diag in diagrams:
    h1_pts = diag[diag[:, 2] == 1]
    h1_pts = h1_pts[h1_pts[:, 1] != np.inf]
    if len(h1_pts) == 0:
        h1_persistence.append(0.0)
    else:
        persistence = h1_pts[:, 1] - h1_pts[:, 0]
        h1_persistence.append(np.max(persistence))

h1_persistence = np.array(h1_persistence)

n_windows = X_windows.shape[0]
time_axis = [
    t_regular[start + window_size // 2]
    for start in range(0, n_windows * step_size, step_size)
]

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(18, 13), gridspec_kw={'hspace': 0.4})

ax1.plot(time_axis, h1_persistence, color='#cf2128', linewidth=3)
ax1.set_title('Topological Mutation Detector (Early Warning)', fontsize=26, fontweight='bold', pad=20)
ax1.set_xlabel('Time [days]', fontsize=25)
ax1.set_ylabel('Max $H_1$ Persistence', fontsize=25)


ax2.plot(t_regular, psa_regular, color='#102f51', linewidth=3, label='PSA Level')
ax2.set_title('Tumor Dynamics (PSA Concentration)', fontsize=26, fontweight='bold', pad=20)
ax2.set_xlabel('Time [days]', fontsize=25)
ax2.set_ylabel('PSA Concentration', fontsize=25)


for ax in [ax1, ax2]:
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.tick_params(axis='both', which='major', labelsize=18)
    ax.axvline( x = 9200, linestyle = '--', color = '#2babe2', linewidth=3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


plt.tight_layout()
plt.savefig('mutation_detect.png')
plt.show()