"""Directional intersection-union and conservative partial-conjunction tests."""
import numpy as np
from scipy.stats import false_discovery_control


def checked_pvalues(values):
    p=np.asarray(values,dtype=float)
    if not np.isfinite(p).all() or (p<0).any() or (p>1).any():raise ValueError('Invalid component P-value')
    return p


def intersection_union(values,axis=-1):
    p=checked_pvalues(values)
    if p.ndim==0 or p.shape[axis]==0:raise ValueError('Empty intersection-union test')
    return p.max(axis=axis)


def partial_conjunction(values,r,axis=-1):
    p=checked_pvalues(values)
    if p.ndim==0:raise ValueError('Need component hypotheses')
    n=p.shape[axis]
    if not isinstance(r,int) or not 1<=r<=n:raise ValueError('Invalid replicability level')
    return np.minimum(1,(n-r+1)*np.take(np.sort(p,axis=axis),r-1,axis=axis))


def adjust_joint_family(values,method='by'):
    p=checked_pvalues(values)
    if not p.size:raise ValueError('Empty hypothesis family')
    if method not in ['by','bh']:raise ValueError('Unsupported multiplicity adjustment')
    return false_discovery_control(p,axis=None,method=method).reshape(p.shape)


def directional_tail_from_two_sided(two_sided,t,direction):
    p=checked_pvalues(two_sided);t=np.asarray(t,dtype=float)
    if p.shape!=t.shape or not np.isfinite(t).all() or direction not in ['up','down']:raise ValueError('Invalid directional statistics')
    if ((t==0)&(p!=1)).any():raise ValueError('Zero t requires two-sided P=1')
    aligned=t>=0 if direction=='up' else t<=0
    return np.where(aligned,p/2,1-p/2)
