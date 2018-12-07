import pandas as pd
import numpy as np
import cobra

from pyefm.ElementaryFluxModes import EFMToolWrapper
from tqdm import tqdm

class EFVWrapper(EFMToolWrapper):

    def create_matrices(self, extra_g=None, extra_h=None):
        """ Initialize the augmented stoichiometric matrix.

        extra_g: (n x nr) array
            Extra entries in the constraint matrix. postive values for lower
            bounds, negative values for upper bounds
        extra_h: (n) array
            Corresponding bounds for the extra entries matrix

        """

        # Create stoichiometric matrix, get key dimensions
        N = cobra.util.create_stoichiometric_matrix(self.model)
        nm, nr = N.shape
        self.nm = nm
        self.nr = nr

        # Construct full G and h matrices, then drop homogeneous (or near
        # homogeneous) entries
        g_full = np.vstack([np.eye(nr), -np.eye(nr)])
        h_full = np.array([(r.lower_bound, -r.upper_bound)
                           for r in self.model.reactions]).T.flatten()

        inhomogeneous = ~((h_full <= -1000) | np.isclose(h_full, 0))
        h_full = h_full[inhomogeneous]
        g_full = g_full[inhomogeneous]

        if extra_g is not None:
            assert extra_g.shape[1] == nr
            assert extra_g.shape[0] == len(extra_h)

            g_full = np.vstack([g_full, extra_g])
            h_full = np.hstack([h_full, extra_h])

        G = g_full
        h = h_full

        self.nt = nt = len(h)

        self.D = np.vstack([
            np.hstack([N, np.zeros((nm, nt)), np.zeros((nm, 1))]),
            np.hstack([G, -np.eye(nt), np.atleast_2d(-h).T])
        ])

    def create_model_files(self, temp_dir):

        # Stoichiometric Matrix
        np.savetxt(temp_dir + '/stoich.txt', self.D, delimiter='\t')

        # Reaction reversibilities
        np.savetxt(
            temp_dir + '/revs.txt', np.hstack([
                np.array([r.lower_bound < 0 for r in self.model.reactions]),
                np.zeros((self.nt + 1))]),
            delimiter='\t', fmt='%d', newline='\t')

        # Reaction Names
        r_names = np.hstack([
                np.array([r.id for r in self.model.reactions]),
                np.array(['s{}'.format(i) for i in range(self.nt)]),
                np.array(['lambda'])
            ])
        with open(temp_dir + '/rnames.txt', 'w') as f:
            f.write('\t'.join(('"{}"'.format(name) for name in r_names)))

        # Metabolite Names
        m_names = np.hstack([
                np.array([m.id for m in self.model.metabolites]),
                np.array(['s{}'.format(i) for i in range(self.nt)]),
            ])
        with open(temp_dir + '/mnames.txt', 'w') as f:
            f.write('\t'.join(('"{}"'.format(name) for name in m_names)))

        pass

    def read_double_out(self, out_file):
        with open(out_file, 'rb') as f:
            out_arr = np.fromstring(f.read()[13:], dtype='>d').reshape(
                (-1, self.nt + self.nr + 1)).T
            out_arr = np.asarray(out_arr, dtype=np.float64).T

        # Sort by the absolute value of the stoichiometry
        sort_inds= np.abs(out_arr[:, :self.nr]).sum(1).argsort()
        out_arr = out_arr[sort_inds]

        unbounded = out_arr[np.isclose(out_arr[:,-1], 0.)]
        bounded = out_arr[~np.isclose(out_arr[:,-1], 0.)]

        if bounded.size:  # Test if its empty
            bounded /= np.atleast_2d(bounded[:,-1]).T

        unbounded_df = pd.DataFrame(
            unbounded[:, :self.nr], 
            columns=[r.id for r in self.model.reactions],
            index=['UEV{}'.format(i) 
                   for i in range(1, unbounded.shape[0] + 1)])

        bounded_df = pd.DataFrame(
            bounded[:, :self.nr], 
            columns=[r.id for r in self.model.reactions],
            index=('BEV{}'.format(i) 
                   for i in range(1, bounded.shape[0] + 1)))

        return unbounded_df.append(bounded_df)
 

def calculate_elementary_vectors(cobra_model, opts=None, verbose=True,
                                 java_args=None, extra_g=None, extra_h=None):
    """Calculate elementary flux vectors, which capture arbitrary linear
    constraints. Approach as detailed in S. Klamt et al., PLoS Comput Biol. 13,
    e1005409–22 (2017).

    Augmented constraints as a hacky workaround for implementing more
    complicated constraints without using optlang.

    java_args: string
        Extra command-line options to pass to the java virtual machine.
        Eg. '-Xmx1g' will set the heap space to 1 GB.
    
    extra_g: (n x nr) array
        Extra entries in the constraint matrix. postive values for lower
        bounds, negative values for upper bounds
    extra_h: (n) array
        Corresponding bounds for the extra entries matrix

    """
    efv_wrap =  EFVWrapper(cobra_model, opts, verbose, java_args=java_args)
    efv_wrap.create_matrices(extra_g=extra_g, extra_h=extra_h)
    return efv_wrap()
    

def get_support_minimal(efvs):
    """Return only those elementary flux vectors whose support is not a proper
    superset of another EFV"""
    
    bool_df = pd.DataFrame(np.isclose(efvs, 0),
                           columns=efvs.columns, index=efvs.index)
    set_df = bool_df.apply(lambda x: set(x.index[~x]), 1)
    set_df = set_df[set_df != set()]  # Drop the empty set EFV
    set_dict = set_df.to_dict()    
    
    is_support_minimal = _get_support_minimal_list(set_dict)

    return efvs.loc[is_support_minimal]


def _get_support_minimal_list(set_dict):

    all_keys = set(set_dict.keys())
    is_support_minimal = []

    for this_key, val in tqdm(set_dict.items()):

        for key in all_keys.difference(set([this_key])):
            if val.issuperset(set_dict[key]):
                break
        else:
            is_support_minimal.append(this_key)

    return is_support_minimal
