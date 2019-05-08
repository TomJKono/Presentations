#!/usr/bin/env python
"""Demo of some of scikit-allel functionality. Adapted from
http://alimanfoo.github.io/2016/06/10/scikit-allel-tour.html"""

import allel
import numpy as np
import scipy
import pandas
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('white')
sns.set_style('ticks')


def plot_variant_hist(f, bins=30):
    """Plot a distribution of a given variant annotation score"""
    k = 'variants/' + f
    x = calls[k][:]
    # Remove NA
    x = x[~np.isnan(x)]
    fig, ax = plt.subplots(figsize=(7, 5))
    sns.despine(ax=ax, offset=10)
    ax.hist(x, bins=bins)
    ax.set_xlabel(f)
    ax.set_ylabel('Count')
    ax.set_title('Variant %s Distribution' % f)


def locate_ti(x):
    """Return a boolean array that is True for transitions and False for
    transversions."""
    x = np.asarray(x)
    return (x == 'AG') | (x == 'GA') | (x == 'CT') | (x == 'TC')


def ti_tv(x):
    """Return the Ti/Tv ratio."""
    if len(x) == 0:
        return np.nan
    is_ti = locate_ti(x)
    n_ti = np.count_nonzero(is_ti)
    n_tv = np.count_nonzero(~is_ti)
    if n_tv > 0:
        return n_ti / n_tv
    else:
        return np.nan


def plot_ti_tv(f, downsample, bins):
    """Plot Ti/Tv against some other specified variant annotation"""
    fig, ax = plt.subplots(figsize=(7, 5))
    sns.despine(ax=ax, offset=10)
    k = 'variants/' + f
    x = calls[k][:]
    keep = ~np.isnan(x)
    x = x[keep]
    x = x[::downsample]
    nona_mut = mut[keep]
    # plot a histogram
    ax.hist(x, bins=bins)
    ax.set_xlabel(f)
    ax.set_ylabel('Count')
    # plot Ti/Tv
    ax = ax.twinx()
    sns.despine(ax=ax, bottom=True, left=True, right=False, offset=10)
    values = nona_mut[::downsample]
    with np.errstate(over='ignore'):
        # binned_statistic generates an annoying overflow warning which we can ignore
        y1, _, _ = scipy.stats.binned_statistic(x, values, statistic=ti_tv, bins=bins)
    bx = (bins[1:] + bins[:-1]) / 2
    ax.plot(bx, y1, color='k')
    ax.set_ylabel('Ti/Tv')
    ax.set_ylim(0.5, 3)
    ax.set_title('Variant %s and Ti/Tv' % f)


def plot_joint_ti_tv(f1, f2, downsample, gridsize=20, mincnt=20, vmin=0.6, vmax=1.4, extent=None):
    fig, ax = plt.subplots()
    sns.despine(ax=ax, offset=10)
    k1 = 'variants/' + f1
    k2 = 'variants/' + f2
    x = calls[k1][:][::downsample]
    y = calls[k2][:][::downsample]
    m = mut[::downsample]
    keepx = ~np.isnan(x)
    keepy = ~np.isnan(y)
    x = x[keepx & keepy]
    y = y[keepx & keepy]
    C = m[keepx & keepy]
    im = ax.hexbin(x, y, C=C, reduce_C_function=ti_tv, mincnt=mincnt, extent=extent,
                   gridsize=gridsize, cmap='jet', vmin=vmin, vmax=vmax)
    fig.colorbar(im)
    ax.set_xlabel(f1)
    ax.set_ylabel(f2)
    ax.set_title('Variant %s versus %s and Ti/Tv' % (f1, f2))


# Read the sample VCF
calls = allel.read_vcf('sample.vcf.bgz', fields='*')
# View the keys in the calls dict
print(list(calls))

# Plot a histogram of DP
plot_variant_hist('DP', bins=100)
plt.savefig('DP_Hist.png', dpi=300)

# And MQ
plot_variant_hist('MQ', bins=100)
plt.savefig('MQ_Hist.png', dpi=300)

# And QD
plot_variant_hist('QD', bins=100)
plt.savefig('QD_Hist.png', dpi=300)

# Transition-transversion ratio
# Make a vector of "mutations"
mut = np.asarray(
    [''.join(x)
     for x
     in zip(calls['variants/REF'], calls['variants/ALT'][:, 0])
     ])
is_ti = locate_ti(mut)
print(is_ti)
# And calculate the proportion of ti/tv
titv = ti_tv(mut)
print(titv)

# Plot Ti:Tv against some other quality metrics
plot_ti_tv('DP', downsample=1, bins=np.arange(0, 2000, 50))
plt.savefig('Ti-Tv_DP.png', dpi=300)
plot_ti_tv('MQ', downsample=1, bins=np.arange(0, 100, 2))
plt.savefig('Ti-Tv_MQ.png', dpi=300)
plot_ti_tv('QD', downsample=1, bins=np.arange(0, 40, 1))
plt.savefig('Ti-Tv_QD.png', dpi=300)

# Make a joint distribution
plot_joint_ti_tv('QD', 'DP', downsample=1, mincnt=25, extent=(0, 40, 0, 2000))
plt.savefig('QD_DP_Ti-Tv.png', dpi=300)

# Convert the calls to a GenotypeArray
geno_arr = allel.GenotypeArray(calls['calldata/GT'])
# Then to allele counts
ac = geno_arr.count_alleles()
# Estimate pi in windows
pi_est, win, nbases, segsites = allel.windowed_diversity(
    calls['variants/POS'],
    ac,
    size=1000,
    start=1,
    stop=558416963)
print(pi_est)

# Plot it
# Mean of windows gives the X-coordinate for the plot
coord = np.mean(win, axis=1)
fig, ax = plt.subplots()
ax.plot(coord/1000000, pi_est)
ax.set(
    xlabel='Position (Mb)',
    ylabel='Pi per bp',
    title='Diversity in 1kb windows')
plt.savefig('Pi_est.png', dpi=300)

# Make a distance heatmap!
nalt = geno_arr.to_n_alt()
dist = allel.pairwise_distance(nalt, metric='cityblock')
distplot = allel.plot_pairwise_distance(dist, labels=list(calls['samples']))
plt.savefig('Pairwise_Distance_Heatmap.png', dpi=300)

# What about inbreeding and FST?
# the first six samples are BA, then 7 MN, then 8 ND
subpops = [
    [0, 1, 2, 3, 4, 5],
    [6, 7, 8, 9, 10, 11, 12],
    [13, 14, 15, 16, 17, 18, 19, 20]]
a, b, c = allel.weir_cockerham_fst(geno_arr, subpops)
fst = (np.sum(a, axis=1) / (np.sum(a, axis=1) + np.sum(b, axis=1) + np.sum(c, axis=1)))

print(fst)
