"""HNN utils to build the L5 cell"""

import neuron
from neuron import h
from neuron.units import ms, mV

def _get_nseg(L):
    nseg = 1
    if L > 100.:  # 100 um
        nseg = int(L / 50.)
        # make dend.nseg odd for all sections
        if not nseg % 2:
            nseg += 1
    return nseg


def _get_dends(params, cell_type, section_names):
    """Convert a flat dictionary to a nested dictionary.

    Returns
    -------
    sections : dict
        Dictionary of sections. Keys are section names
    """
    prop_names = ['L', 'diam', 'Ra', 'cm']
    sections = dict()
    for section_name in section_names:
        dend_prop = dict()
        for key in prop_names:
            if key in ['Ra', 'cm']:
                middle = 'dend'
            else:
                # map apicaltrunk -> apical_trunk etc.
                middle = section_name.replace('_', '')
            dend_prop[key] = params[f'{cell_type}_{middle}_{key}']
        sections[section_name] = Section(L=dend_prop['L'],
                                        diam=dend_prop['diam'],
                                        Ra=dend_prop['Ra'],
                                        cm=dend_prop['cm'])
    return sections

def _get_pyr_soma(p_all, cell_type):
    """Get somatic properties."""
    return Section(
        L=p_all[f'{cell_type}_soma_L'],
        diam=p_all[f'{cell_type}_soma_diam'],
        cm=p_all[f'{cell_type}_soma_cm'],
        Ra=p_all[f'{cell_type}_soma_Ra']
    )

class Section:

    def __init__(self, L, diam, Ra, cm, end_pts=None):

        self._L = L
        self._diam = diam
        self._Ra = Ra
        self._cm = cm
        if end_pts is None:
            end_pts = list()
        self._end_pts = end_pts

        # For distance functionality
        self.nseg = _get_nseg(self.L)

    def __repr__(self):
        return f'L={self.L}, diam={self.diam}, cm={self.cm}, Ra={self.Ra}'

    @property
    def L(self):
        return self._L

    @property
    def diam(self):
        return self._diam

    @property
    def cm(self):
        return self._cm

    @property
    def Ra(self):
        return self._Ra

    @property
    def end_pts(self):
        return self._end_pts
    
class Cell: 

    def __init__(self, sections, cell_tree, gid):

        self._gid = gid
        self.sections = sections

        self._nrn_sections = dict()

        self.cell_tree = cell_tree
        
    def __repr__(self):
        return 'L5_pyramidal[{}]'.format(self._gid)
    
    def _create_sections(self): 
        
        for sec_name in self.sections:
            sec = h.Section(name=f'L5_pyramidal_{sec_name}', cell=self)
            self._nrn_sections[sec_name] = sec

            h.pt3dclear(sec=sec)
            h.pt3dconst(0, sec=sec)  
            for pt in self.sections[sec_name].end_pts:
                h.pt3dadd(pt[0], pt[1], pt[2], 1, sec=sec)
            # with pt3dconst==0, these will alter the 3d points defined above!
            sec.L = self.sections[sec_name].L
            sec.diam = self.sections[sec_name].diam
            sec.Ra = self.sections[sec_name].Ra
            sec.cm = self.sections[sec_name].cm
            sec.nseg = self.sections[sec_name].nseg

        if self.cell_tree is None:
            self.cell_tree = dict()

        # Connects sections of THIS cell together.
        for parent_node in self.cell_tree:
            for child_node in self.cell_tree[parent_node]:
                parent_sec = self._nrn_sections[parent_node[0]]
                child_sec = self._nrn_sections[child_node[0]]
                if parent_sec == child_sec:
                    continue
                parent_loc = parent_node[1]
                child_loc = child_node[1]
                child_sec.connect(parent_sec, parent_loc, child_loc)

        h.define_shape()

    def build(self):

        self._create_sections()



def cell_params():

    section_params = {
        # Soma
        'L5Pyr_soma_L': 39.,
        'L5Pyr_soma_diam': 28.9,
        'L5Pyr_soma_cm': 0.85,
        'L5Pyr_soma_Ra': 495.73,

        # Dendrite
        'L5Pyr_dend_cm': 0.85,
        'L5Pyr_dend_Ra': 495.73,

        'L5Pyr_apicaltrunk_L': 102.,
        'L5Pyr_apicaltrunk_diam': 10.2,

        'L5Pyr_apical1_L': 680.,
        'L5Pyr_apical1_diam': 7.48,

        'L5Pyr_apical2_L': 680.,
        'L5Pyr_apical2_diam': 4.93,

        'L5Pyr_apicaltuft_L': 425.,
        'L5Pyr_apicaltuft_diam': 3.4,

        'L5Pyr_apicaloblique_L': 255.,
        'L5Pyr_apicaloblique_diam': 5.1,

        'L5Pyr_basal1_L': 85.,
        'L5Pyr_basal1_diam': 6.8,

        'L5Pyr_basal2_L': 255.,
        'L5Pyr_basal2_diam': 8.5,

        'L5Pyr_basal3_L': 255.,
        'L5Pyr_basal3_diam': 8.5,
    }

    end_pts = {
            'soma': [[0, 0, 0], [0, 0, 23]],
            'apical_trunk': [[0, 0, 23], [0, 0, 83]],
            'apical_oblique': [[0, 0, 83], [-150, 0, 83]],
            'apical_1': [[0, 0, 83], [0, 0, 483]],
            'apical_2': [[0, 0, 483], [0, 0, 883]],
            'apical_tuft': [[0, 0, 883], [0, 0, 1133]],
            'basal_1': [[0, 0, 0], [0, 0, -50]],
            'basal_2': [[0, 0, -50], [-106, 0, -156]],
            'basal_3': [[0, 0, -50], [106, 0, -156]]
        }

    cell_tree = {
            ('apical_trunk', 0): [('apical_trunk', 1)],
            ('apical_1', 0): [('apical_1', 1)],
            ('apical_2', 0): [('apical_2', 1)],
            ('apical_tuft', 0): [('apical_tuft', 1)],
            ('apical_oblique', 0): [('apical_oblique', 1)],
            ('basal_1', 0): [('basal_1', 1)],
            ('basal_2', 0): [('basal_2', 1)],
            ('basal_3', 0): [('basal_3', 1)],
            # Different sections connected
            ('soma', 0): [('soma', 1), ('basal_1', 0)],
            ('soma', 1): [('apical_trunk', 0)],
            ('apical_trunk', 1): [('apical_1', 0), ('apical_oblique', 0)],
            ('apical_1', 1): [('apical_2', 0)],
            ('apical_2', 1): [('apical_tuft', 0)],
            ('basal_1', 1): [('basal_2', 0), ('basal_3', 0)]
        }

    return section_params, end_pts, cell_tree