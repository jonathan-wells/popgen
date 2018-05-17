from matplotlib import pyplot as plt
import numpy as np

def plot_dominance(het=0.5, sel=0.1):
    """Plot the effect of selection and heterozygosity effect on the change in
    allele frequencies.
    """
    fig = plt.figure()
    ax = fig.add_axes([0, 0, 1, 1])

    # prepare the axes limits
    ax.set_xlim((0.0, 1.0))
    ax.set_ylim((-0.025, 0.05))
    
    def deltap(p, h=het, s=sel):
        """Compute deltap"""
        q = 1 - p
        return (p*q*s*(p*h + q*(1-h)))/(1 - 2*p*q*h*s - s*(q**2)) 
    
    p = np.linspace(0.0, 1.0)
    dp = deltap(p, het, sel)
    lines = ax.plot(p, dp, '-')
    plt.xlabel('p')
    plt.ylabel(r'$\Delta_s p$')
    plt.show()
    return p, dp
