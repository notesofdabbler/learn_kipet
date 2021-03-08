
import re
import numpy as np
import pandas as pd
import pyomo.environ as pyo
from collections import defaultdict

def stoicoeff(x):
    if '*' in x:
        stoi = (x.split('*')[1], int(x.split('*')[0]))
    else:
        stoi = (x, 1)
    return stoi

def get_rxnstoi(rxnlist):
    rxnstoi = dict()
    cmplist = []

    for k, rxn in rxnlist.items():
        rxn2 = re.split('->|=', rxn)
        rxntype = 'eqb' if '=' in rxn else 'fwd'
        reactants = [x.strip() for x in rxn2[0].split('+')]
        products = [x.strip() for x in rxn2[1].split('+')]
        stoi = {stoicoeff(x)[0]: -stoicoeff(x)[1] for x in reactants}
        stoi_products = {stoicoeff(x)[0]: stoicoeff(x)[1] for x in products}
        stoi.update(stoi_products)
        #stoi.update({'rxntype': rxntype})
        cmplist = cmplist + reactants
        cmplist = cmplist + products
        rxnstoi[k] = stoi
        
    cmplist = list(set(cmplist))
    cmplist.sort()
    
    return cmplist, rxnstoi

def get_rxnrates(r1, rxndict, T = 325):

    cmplist, rxnstoi = get_rxnstoi(rxndict['rxn'])
    cmp = dict()
    for c in cmplist:
        cmp[c] = r1.component(c, value = 0.0)
    
    kf = dict()
    Keq = dict()
    Eaf = dict()
    Ear = dict()
    for r in rxndict['rxn'].keys():
        kf[r] = r1.parameter('k_'+ r, value = rxndict['kf'][r])
        if not np.isnan(rxndict['Keq'][r]):
            Keq[r] = r1.parameter('Keq_' + r, value = rxndict['Keq'][r])
        if not np.isnan(rxndict['Eaf'][r]):
            Eaf[r] = r1.parameter('Eaf_' + r, value = rxndict['Eaf'][r])
        if not np.isnan(rxndict['Ear'][r]):
            Ear[r] = r1.parameter('Ear_' + r, value = rxndict['Ear'][r])

    Tref = 298 # K
    R = 8.314 # J/mol.K

    rates = dict()
    for r in rxndict['rxn'].keys():
        ratesf = kf[r]
        if r in Eaf.keys():
            ratesf = ratesf * pyo.exp(-Eaf[r]/R * (1/T - 1/Tref))
        for c, s in rxnstoi[r].items():
            if s < 0:
                ratesf = ratesf * cmp[c] ** (-s)
        rates[r] = ratesf        
        
        if r in Keq.keys():
            ratesr = kf[r] / Keq[r]
            if r in Ear.keys():
                ratesr = ratesr * pyo.exp(-Ear[r]/R * (1/T - 1/Tref))
            for c, s in rxnstoi[r].items():
                if s > 0:
                    ratesr = ratesr * cmp[c] ** s
            rates[r] = rates[r] - ratesr   

    r_dict = defaultdict(lambda: 0)
    for r in rxndict['rxn'].keys():
        for c, s in rxnstoi[r].items():
            r_dict[c] = r_dict[c] + s * rates[r]

    r1.add_odes(r_dict)

    return r1