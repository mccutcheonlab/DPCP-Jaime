# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 11:45:57 2019

@author: admin
"""

df = pd.DataFrame([sessions[s].rat for s in sessions], columns=['rat'])
df['session'] = [sessions[s].session for s in sessions]

df['condition'] = [sessions[s].condition for s in sessions]

df['sexp_alpha'] = [sessions[s].sexp_alpha for s in sessions]
df['sexp_rsq'] = [sessions[s].sexp_rsq for s in sessions]

df['dexp_alpha'] = [sessions[s].dexp_alpha for s in sessions]
df['dexp_beta'] = [sessions[s].dexp_beta for s in sessions]
df['dexp_tau'] = [sessions[s].dexp_tau for s in sessions]
df['dexp_rsq'] = [sessions[s].dexp_rsq for s in sessions]





