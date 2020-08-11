from EJR_junction_class import *
from EJR_junction_auxfns import *

def get_EJR_junction_pool(original_seq, leftnick, rightnick, max_deln):
    """
    """
    deln_pool = all_valid_delns(original_seq, leftnick, rightnick, max_deln)
    EJR_junctions = []
    for deln in deln_pool: #deln with be a tuple of deletion indices
        jnxn = EJR_junction(original_seq, deln, leftnick, rightnick)
        EJR_junctions.append(jnxn)
    return EJR_junctions


##def get_ABJs(EJR_junction_pool):
##    ABJs = []
##    for ejrj in EJR_junction_pool:
##        if ejrj.jnxn_type == "ABJ": ABJs.append(ejrj)
##    return ABJs
##
##
##def num_ABJs(EJR_junction_pool):
##    ABJs = get_ABJs(EJR_junction_pool)
##    return len(ABJs)
##
##
##def get_MHJs(EJR_junction_pool):
##    MHJs = []
##    for ejrj in EJR_junction_pool:
##        if ejrj.jnxn_type == "MHJ": MHJs.append(ejrj)
##    return MHJs
##
##
##def num_MHJs(EJR_junction_pool):
##    MHJs = get_MHJs(EJR_junction_pool)
##    return len(MHJs)
##



    
        
    
    
