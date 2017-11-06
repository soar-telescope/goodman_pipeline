import matplotlib.pyplot as plt
import matplotlib.colors as mcol
import matplotlib.cm as cm
import numpy as np

fig, ax = plt.subplots()
ax.axis('equal')
width = 0.3





core_pie, _ = ax.pie([1000], radius=0.2, labels=['Core'], colors=['w'], labeldistance=-0.001)
plt.setp(core_pie, width=.2, edgecolor='white')

cm2 = plt.get_cmap('Blues')
cmid = cm2(np.arange(3, 6) * 50)
print(cmid)

mid_labels = ['data_classifier',
              'image_processor',
              'wavmode_translator',
              'night_organizer',
              'redspec',
              'wavelength']
mid_pie, _ = ax.pie([100, 100, 100, 100, 100, 100], radius=0.7, labels=mid_labels, colors=cmid, labeldistance=0.4, explode=[.1,.1,.1,.1,.1,.1])
plt.setp(mid_pie, width=.5, edgecolor='white')


#
cm3 = plt.get_cmap("Reds")
cin = cm3(np.array([1,2,5,6,9,10])*10)
labels = map("".join, zip(list("aabbcc"), map(str, [1, 2] * 3)))
pie2, _ = ax.pie([60,60,37,40,29,10], radius=1.3, labels=labels,
                                      labeldistance=.9, colors=cin)
plt.setp(pie2, width=width, edgecolor='white')
plt.show()

# Accent
# Accent_r
# Blues
# Blues_r
# BrBG
# BrBG_r
# BuGn
# BuGn_r
# BuPu
# BuPu_r
# CMRmap
# CMRmap_r
# Dark2
# Dark2_r
# GnBu
# GnBu_r
# Greens
# Greens_r
# Greys
# Greys_r
# OrRd
# OrRd_r
# Oranges
# Oranges_r
# PRGn
# PRGn_r
# Paired
# Paired_r
# Pastel1
# Pastel1_r
# Pastel2
# Pastel2_r
# PiYG
# PiYG_r
# PuBu
# PuBuGn
# PuBuGn_r
# PuBu_r
# PuOr
# PuOr_r
# PuRd
# PuRd_r
# Purples
# Purples_r
# RdBu
# RdBu_r
# RdGy
# RdGy_r
# RdPu
# RdPu_r
# RdYlBu
# RdYlBu_r
# RdYlGn
# RdYlGn_r
# Reds
# Reds_r
# Set1
# Set1_r
# Set2
# Set2_r
# Set3
# Set3_r
# Spectral
# Spectral_r
# Wistia
# Wistia_r
# YlGn
# YlGnBu
# YlGnBu_r
# YlGn_r
# YlOrBr
# YlOrBr_r
# YlOrRd
# YlOrRd_r
# afmhot
# afmhot_r
# autumn
# autumn_r
# binary
# binary_r
# bone
# bone_r
# brg
# brg_r
# bwr
# bwr_r
# cool
# cool_r
# coolwarm
# coolwarm_r
# copper
# copper_r
# cubehelix
# cubehelix_r
# flag
# flag_r
# gist_earth
# gist_earth_r
# gist_gray
# gist_gray_r
# gist_heat
# gist_heat_r
# gist_ncar
# gist_ncar_r
# gist_rainbow
# gist_rainbow_r
# gist_stern
# gist_stern_r
# gist_yarg
# gist_yarg_r
# gnuplot
# gnuplot2
# gnuplot2_r
# gnuplot_r
# gray
# gray_r
# hot
# hot_r
# hsv
# hsv_r
# inferno
# inferno_r
# jet
# jet_r
# magma
# magma_r
# nipy_spectral
# nipy_spectral_r
# ocean
# ocean_r
# pink
# pink_r
# plasma
# plasma_r
# prism
# prism_r
# rainbow
# rainbow_r
# seismic
# seismic_r
# spectral
# spectral_r
# spring
# spring_r
# summer
# summer_r
# terrain
# terrain_r
# viridis
# viridis_r
# winter
# winter_r