"""
 Library of some functions used by multiple classes
 Written by John Kazimierz Farey, Sep 2012
"""


def sub2ind(shap, subx, suby):
    """subscript to index of a 1d array"""
    ind = (subx * shap[0]) + suby
    return ind

    # THIS FUNCTION MAY BE BROKEN
def rot_rep(pxin, pyin, dummy, cd_todo, gcos, gsin):
    """rotate function"""
    if cd_todo.lower() in ['en to i', 'ij to e']:
        x_var, y_var = pxin, pyin
    elif cd_todo.lower() in ['en to j', 'ij to n']:
        x_var, y_var = pyin, pxin*-1
    else:
        raise SyntaxError('rot_rep cd_todo %s is invalid' %cd_todo)
    return x_var * gcos + y_var * gsin

def get_output_filename(setup_var, year, month, var_type):
    """This returns a output filename constructed for a given var_type, year and month"""
    if var_type == 'ice':
        return setup_var.settings['dst_dir']+setup_var.settings['fn']+'_bdyT_y'+str(year)+ \
               'm'+str(month)+'.nc'
    elif var_type == 'bt':
        return setup_var.settings['dst_dir']+setup_var.settings['fn']+'_bt_bdyT_y'+str(year)+ \
               'm'+str(month)+'.nc'
    elif var_type == 'u':
        return setup_var.settings['dst_dir'] + setup_var.settings['fn'] + '_bdyU_y' + \
               str(year) + 'm' + str(month) + '.nc'
    elif var_type == 'v':
        return setup_var.settings['dst_dir'] + setup_var.settings['fn'] + '_bdyV_y' + \
               str(year) + 'm' + str(month) + '.nc'

def get_output_tidal_filename(setup_var, const_name, grid_type):
    """This method returns a output filename constructed for a given tidal constituent and
    grid type"""
    return setup_var.settings['dst_dir']+setup_var.settings['fn']+"_bdytide_rotT_"+const_name+ \
           "_grid_"+grid_type.upper()+".nc"
