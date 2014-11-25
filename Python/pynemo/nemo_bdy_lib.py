# 
# Library of some functions used by multiple classes
# Written by John Kazimierz Farey, Sep 2012
#


def sub2ind(shap, subx, suby):
    ind = (subx * shap[0]) + suby
    return ind
    
    # THIS FUNCTION MAY BE BROKEN
def rot_rep(pxin, pyin, cd_type, cd_todo, gcos, gsin):
    if cd_todo.lower() in ['en to i', 'ij to e']:
        x,y = pxin, pyin
    elif cd_todo.lower() in ['en to j', 'ij to n']:
        x,y = pyin, pxin*-1
    else:
        raise SyntaxError('rot_rep cd_todo %s is invalid' %cd_todo)
    # cd_type = superfluous?
    return x * gcos + y * gsin
