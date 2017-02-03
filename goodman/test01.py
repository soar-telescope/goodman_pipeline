from astropy.io import fits
import matplotlib.pyplot as plt
import wsbuilder


def test_func(message='Hola', **kwargs):
    print message
    print len(kwargs)


test_func()

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, gridspec_kw={'width_ratios': [5, 1]}, figsize=(25, 10))
ax1.set_title('Ax 1')
ax2.set_title('Ax 2')
ax3.set_title('Ax 3')
ax4.set_title('Ax 4')
plt.tight_layout()
plt.show()