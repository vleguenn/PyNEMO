# 
# Library of some functions used by multiple classes
# Written by John Kazimierz Farey, Sep 2012
#


def sub2ind(shap, subx, suby):
    ind = (subx * shap[0]) + suby
    return ind
