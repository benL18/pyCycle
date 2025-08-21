import sys

import numpy as np

# protection incase env doesn't have matplotlib installed, since its not strictly required
try:
    import matplotlib
    import matplotlib.pyplot as plt
except ImportError:
  plt = None


def get_val(prob, point, element, var_name, units=None):
    """
    Get the value of the requested element from the OpenMDAO model.

    Parameters
    ----------
    prob: <Problem>
        OpenMDAO problem that contains the pycycle model.
    point: str
        OpenMDAO pathname of the cycle point.
    element: str
        Name of the pcycle element system.
    var_name: str
        Name of the variable to get.
    units: str or None
        Units for the return value. Default is None, which will return using the delared units.

    Returns
    float
        Value of requested element.
    """
    try:
        val = prob.get_val(f"{point}.{element}.{var_name}", units=units)
    except KeyError:
        # This is a cycle parameter.
        val = prob.get_val(f"{element}.{var_name}", units=units)

    return val[0]

def print_flow_station(prob, fs_names, file=sys.stdout, units='Default'):
    # user can specify the variable and units to be used below 
    # leaving the cell empty ('') means to use the default 
    names       = ['tot:P', 'tot:T', 'tot:h', 'tot:S', 'stat:P', 'stat:W', 'stat:MN', 'stat:V', 'stat:area']
    if units == 'Default':
        names_units = [''   ,''     ,   ''  ,    ''  ,    ''   ,   ''    ,    ''    ,    ''   ,     ''     ]
    # define the SI units for each variables
    elif units=='SI': 
        names_units = ['Pa'   ,'K'     ,'kJ/kg','kJ/kg/K','Pa'     ,   'kg/s',        '',    'm/s',  'm**2']
    else:
        print('UNITS not recognized')
        return None
    
    n_names = len(names)
    line_tmpl = '{:<23}|  '+'{:>13}'*n_names
    len_header = 27+13*n_names

    print("-"*len_header, file=file, flush=True)
    print("                            FLOW STATIONS", file=file, flush=True)
    print("-"*len_header, file=file, flush=True)

    # header_line
    vals = ['Flow Station'] + names
    # print('-'*len_header, file=file, flush=True)
    print(line_tmpl.format(*vals), file=file, flush=True)
    # print(line_tmpl.format(*vals_units), file=file, flush=True)

    # Processing Units: extract the units of the variables from the full variables data structure
    # step1: create a dictionnary to link the full variable names and promoted names
    prom_name_to_abs = {
        meta['prom_name']: abs_name
        for abs_name, meta in prob.model.get_io_metadata(iotypes=('input', 'output')).items()
    }
    # step2: extract the units for the first flow station
    line_tmpl_units = '{:<23.23}|  ' + '{:>13}'*n_names
    fs_name=fs_names[0]
    data_units = []
    ind = 0
    for name in names:
        full_name = '{}:{}'.format(fs_name, name)
        # within the list of units: only use the default units when specified by ''
        if names_units[ind] != '':
            data_units.append(names_units[ind])
        else:
            data_units.append(prob.model.get_io_metadata(iotypes=('input','output'))[prom_name_to_abs[full_name]]['units'])
        # index ind is used to pick the corresponding value in the names_unit list
        ind = ind + 1 
        
    vals_units = ['         ' + units + ' units'] + data_units
    vals_units_clean = [str(v) if v is not None else '-' for v in vals_units]
    # step3: print the units as part of the output table
    print(line_tmpl_units.format(*vals_units_clean), file=file, flush=True)

# print separator line
    print('-'*len_header, file=file, flush=True)

# loop through all elements and variable names to collect the data
    line_tmpl = '{:<23.23}|  ' + '{:13.3f}'*n_names
    for fs_name in fs_names:
        data = []
        ind = 0
        for name in names:
            full_name = '{}:{}'.format(fs_name, name)

            if names_units[ind] != '':
                data.append(prob.get_val(full_name,units=names_units[ind])[0])
            else: 
                data.append(prob[full_name][0])
            ind = ind + 1

        vals = [fs_name] + data
        print(line_tmpl.format(*vals), file=file, flush=True)
        
        
    print('-'*len_header, file=file, flush=True)



def print_compressor(prob, element_names, file=sys.stdout, units='Default'):

    # user can specify the variable and units to be used below 
    # leaving the cell empty ('') means to use the default 
    display_names = ['Wc'  , 'PR', 'eta_a', 'eta_p'   , 'Nc' , 'pwr'  , 'RlineMap'    , 'NcMap'    , 'PRmap'    ,  'WcMap'   , 'alphaMap'    , 'SMN', 'SMW', 'effMap'    ]
    names         = ['Wc'  , 'PR',  'eff' , 'eff_poly', 'Nc' , 'power', 'map.RlineMap', 'map.NcMap', 'map.PRmap', 'map.WcMap', 'map.alphaMap', 'SMN', 'SMW', 'map.effMap']
    
    if units == 'Default':
        names_units   = [''    ,  '' ,   ''   ,    ''     , ''   ,  ''    ,      ''       ,  ''        ,    ''      ,   ''       ,   ''          , ''   ,  ''  ,      ''     ]    
    # define the SI units for each variables
    elif units=='SI': 
        names_units   = ['kg/s',  '' ,   ''   ,    ''     , 'rpm',  'kW'  ,      ''       ,  'rpm'     ,    ''      ,   'kg/s'   ,   ''          , ''   ,  ''  ,      ''     ]    
    else:
        print('UNITS not recognized')
        return None    


    # configure & print header
    len_header = 23+14*13
    # print("-"*len_header)
    print("-"*len_header, file=file, flush=True)
    print("                            COMPRESSOR PROPERTIES", file=file, flush=True)
    print("-"*len_header, file=file, flush=True)

    # configure & print the variables and units
    line_tmpl = '{:<23}|  '+'{:>11}'*14
    print(line_tmpl.format(*(['Compressor'] + display_names)),
          file=file, flush=True)

    # NEW: extract the units of the variables from the full variables data structure
    # step1: create a dictionnary to link the full variable names and promoted names
    prom_name_to_abs = {
        meta['prom_name']: abs_name
        for abs_name, meta in prob.model.get_io_metadata(iotypes=('input', 'output')).items()
    }
    # step2: extract the units for the first flow station
    line_tmpl_units = '{:<23}|  ' + '{:>11}'*14
    e_name=element_names[0]
    data_units = []
    
    # loop through names of variables to get units
    ind = 0
    for name in names:
        full_name = '{}.{}'.format(e_name, name)
        if names_units[ind] != '':
            data_units.append(names_units[ind])
        else:
            data_units.append(prob.model.get_io_metadata(iotypes=('input','output'))[prom_name_to_abs[full_name]]['units'])

        ind = ind + 1 
        
    vals_units = ['         ' + units + ' units'] + data_units
    vals_units_clean = [str(v) if v is not None else '-' for v in vals_units]
    # step3: print the units as part of the output table
    print(line_tmpl_units.format(*vals_units_clean), file=file, flush=True)

    # process data
    print("-"*len_header, file=file, flush=True)

    line_tmpl = '{:<23}|  '+'{:11.3f}'*14
    for e_name in element_names:
        sys = prob.model._get_subsystem(e_name)
        if sys.options['design']:
          PR_temp = prob[e_name+'.map.scalars.PR'][0]
          eff_temp = prob[e_name+'.map.scalars.eff'][0]
        else:
          PR_temp = prob[e_name+'.PR'][0]
          eff_temp = prob[e_name+'.eff'][0]

        data = []
        ind = 0
        for name in names:
            full_name = '{}.{}'.format(e_name, name)
            
            if names_units[ind] != '':
                data.append(prob.get_val(full_name,units=names_units[ind])[0])
            else:
                data.append(prob[full_name][0])
            ind = ind +1 

        vals = [e_name] + data
        print(line_tmpl.format(*vals), file=file, flush=True)
        
    print("-"*len_header, file=file, flush=True)


def print_burner(prob, element_names, file=sys.stdout, units='Default'):
    
    # user can specify the variable and units to be used below 
    # leaving the cell empty ('') means to use the default 
    # note dPqP is not listed as it is not found via the method used for units, also it hs no units
    # FAR is calculated based on massflow through burner and fuel flow (Wfuel)
    display_names = ['dPqP',    'TtOut'  , 'Wfuel', 'FAR']
    names         = [        'Fl_O:tot:T', 'Wfuel']
    
    
    if units == 'Default':
        print('Default')
        names_units   = [     ''   ,   ''   ]    
    # define the SI units for each variables
    elif units=='SI': 
        names_units   = [  'K'     , 'kg/s']
    else:
        print('UNITS not recognized')
        return None  
    
    # configure & print header
    len_header = 27+4*13
    print("-"*len_header, file=file, flush=True)
    print("                            BURNER PROPERTIES", file=file, flush=True)
    print("-"*len_header, file=file, flush=True)

    line_tmpl = '{:<23}|  '+'{:>13}'*4
    print(line_tmpl.format(*(['Burner'] + display_names)), file=file, flush=True)
    
    # NEW: extract the units of the variables from the full variables data structure
    # step1: create a dictionnary to link the full variable names and promoted names
    prom_name_to_abs = {
        meta['prom_name']: abs_name
        for abs_name, meta in prob.model.get_io_metadata(iotypes=('input', 'output')).items()
    }
    # step2: extract the units for the first flow station
    line_tmpl_units = '{:<23}|  ' + '{:>13}'*4
    e_name=element_names[0]
    data_units = []
    
    # loop through names of variables to get units
    ind = 0
    for name in names:
        full_name = '{}.{}'.format(e_name, name)
        if names_units[ind] != '':
            data_units.append(names_units[ind])
        else:
            data_units.append(prob.model.get_io_metadata(iotypes=('input','output'))[prom_name_to_abs[full_name]]['units'])
            # temp_debug
            # print(data_units)
        ind = ind + 1 
        
    vals_units = ['         ' + units + ' units'] + ["-"] + data_units  + [' - '] # extra for dPqP and FAR that are defined as input (dPqP) and calculated (FAR)
    vals_units_clean = [str(v) if v is not None else '-' for v in vals_units]
    # step3: print the units as part of the output table
    print(line_tmpl_units.format(*vals_units_clean), file=file, flush=True)

    
    # extract data for each element
    print("-"*len_header, file=file, flush=True)

    # line_tmpl = '{:<20}|  '+'{:13.3f}'*4
    line_tmpl = '{:<23}|  {:13.4f}{:13.2f}{:13.4f}{:13.5f}'

    for e_name in element_names:

        point, _, element = e_name.rpartition(".")

        dPqP = get_val(prob, point, element, 'dPqP')
        T_tot = get_val(prob, point, element, 'Fl_O:tot:T',units = data_units[0])
        W_fuel = get_val(prob, point, element, 'Wfuel',units = data_units[1])
        W_tot = get_val(prob, point, element, 'Fl_O:stat:W',units = data_units[1])
        W_air = W_tot - W_fuel
        FAR = W_fuel/W_air

        print(line_tmpl.format(e_name, dPqP, T_tot, W_fuel, FAR),
              file=file, flush=True)


def print_turbine(prob, element_names, file=sys.stdout, units='Default'):
    # user can specify the variable and units to be used below 
    # leaving the cell empty ('') means to use the default 
    display_names = ['Wp'  , 'PR', 'eff_a', 'eff_p'   , 'Np', 'pwr'  , 'NpMap'    , 'PRmap'    , 'alphaMap'    ]
    names         = ['Wp'  , 'PR', 'eff'  , 'eff_poly', 'Np', 'power', 'map.NpMap', 'map.PRmap', 'map.alphaMap']

    
    if units == 'Default':
        print('Default')
        names_units   = [  ''  , ''  ,   ''   ,    ''     , ''  ,   ''   ,    ''      ,      ''    ,         ''    ]      
    # define the SI units for each variables
    elif units=='SI': 
        names_units   = ['kg/s', ''  ,   ''   ,    ''     ,'rpm',   'kW' ,   'rpm'    ,      ''    ,         ''    ]      
    else:
        print('UNITS not recognized')
        return None    

    # configure & print header
    len_header = 27+9*13
    print("-"*len_header, file=file, flush=True)
    print("                            TURBINE PROPERTIES", file=file, flush=True)
    print("-"*len_header, file=file, flush=True)

    line_tmpl = '{:<23}|  '+'{:>13}'*9
    print(line_tmpl.format(*(['Turbine'] + display_names)),file=file, flush=True)

    # UNITS --------------------------------
    # NEW: extract the units of the variables from the full variables data structure 
    # step1: create a dictionnary to link the full variable names and promoted names
    prom_name_to_abs = {
        meta['prom_name']: abs_name
        for abs_name, meta in prob.model.get_io_metadata(iotypes=('input', 'output')).items()
    }
    # step2: extract the units for the first flow station
    line_tmpl_units = '{:<23}|  ' + '{:>13}'*9
    e_name=element_names[0]
    data_units = []
       
    # loop through names of variables to get units
    ind = 0
    for name in names:
        full_name = '{}.{}'.format(e_name, name)
        if names_units[ind] != '':
            data_units.append(names_units[ind])
        else:
            data_units.append(prob.model.get_io_metadata(iotypes=('input','output'))[prom_name_to_abs[full_name]]['units'])
        ind = ind + 1 
        
    vals_units = ['         ' + units + ' units'] + data_units
    vals_units_clean = [str(v) if v is not None else '-' for v in vals_units]
    # step3: print the units as part of the output table
    print(line_tmpl_units.format(*vals_units_clean), file=file, flush=True)
    # end UNITS --------------------------------
    print("-"*len_header, file=file, flush=True)

    # extract data and print output
    line_tmpl = '{:<23}|  '+'{:13.3f}'*9
    for e_name in element_names:

        point, _, element = e_name.rpartition(".")

        sys = prob.model._get_subsystem(e_name)
        if sys.options['design']:
            PR_temp = get_val(prob, point, element, 'map.scalars.PR')
            eff_temp = get_val(prob, point, element, 'map.scalars.eff')
        else:
            PR_temp = get_val(prob, point, element, 'PR')
            eff_temp = get_val(prob, point, element, 'eff')


        data = []
        ind = 0
        for name in names:
            full_name = '{}.{}'.format(e_name, name)
            
            if names_units[ind] != '':
                data.append(prob.get_val(full_name,units=names_units[ind])[0])
            else:
                data.append(prob[full_name][0])
            ind = ind +1 

        vals = [e_name] + data
        print(line_tmpl.format(*vals), file=file, flush=True)


def print_nozzle(prob, element_names, file=sys.stdout, units='Default'):
    # user can specify the variable and units to be used below 
    # leaving the cell empty ('') means to use the default 
    display_names = [ 'PR', 'Cv', 'Cfg', 'Ath'             , 'MNth'          , 'MNout'       , 'V'          , 'Fg']
    names         = [ 'PR', 'Cv', 'Cfg', 'Throat:stat:area', 'Throat:stat:MN', 'Fl_O:stat:MN', 'Fl_O:stat:V', 'Fg']
     
    if units == 'Default':
         print('Default')
         names_units   = [  '' ,  '' , ''   ,   ''              ,     ''          ,       ''      ,    ''        ,  '' ]
    # define the SI units for each variables
    elif units=='SI': 
         names_units   = [  '' ,  '' , ''   , 'm**2'            ,     ''          ,       ''      ,   'm/s'      , 'kN']
    else:
         print('UNITS not recognized')
         return None    
     
     # configure & print header
    len_header = 27+8*13
    print("-"*len_header, file=file, flush=True)
    print("                            NOZZLE PROPERTIES", file=file, flush=True)
    print("-"*len_header, file=file, flush=True)

    line_tmpl = '{:<23}|  '+'{:>13}'*8
    print(line_tmpl.format(*(['Nozzle'] + display_names)), file=file, flush=True)

    # UNITS --------------------------------
    # NEW: extract the units of the variables from the full variables data structure       
    # step1: create a dictionnary to link the full variable names and promoted names
    prom_name_to_abs = {
        meta['prom_name']: abs_name
        for abs_name, meta in prob.model.get_io_metadata(iotypes=('input', 'output')).items()
    }
    # step2: extract the units for the first flow station
    line_tmpl_units = '{:<23}|  ' + '{:>13}'*8
    e_name=element_names[0]
    data_units = []
       
    # loop through names of variables to get units
    ind = 0
    for name in names:
        
        if name == "Cv" or name == "Cfg":
            data_units.append('')
        else:
            full_name = '{}.{}'.format(e_name, name)
            if names_units[ind] != '':
                data_units.append(names_units[ind])
            else:
                data_units.append(prob.model.get_io_metadata(iotypes=('input','output'))[prom_name_to_abs[full_name]]['units'])    
        ind = ind +1
            
    vals_units = ['         ' + units + ' units'] + data_units
    vals_units_clean = [str(v) if v is not None else '-' for v in vals_units]
    
    # step3: print the units as part of the output table
    print(line_tmpl_units.format(*vals_units_clean), file=file, flush=True)
    # end UNITS --------------------------------
    print("-"*len_header, file=file, flush=True)

    for e_name in element_names:

        point, _, element = e_name.rpartition(".")

        sys = prob.model._get_subsystem(e_name)
        if sys.options['lossCoef'] == 'Cv':

            Cv_val = get_val(prob, point, element, 'Cv')
            Cfg_val = '        N/A  '
            line_tmpl = '{:<23}|  ' + '{:13.3f}'*2 + '{}' + '{:13.3f}'*5

        else:
            Cv_val = '        N/A  '
            Cfg_val = get_val(prob, point, element, 'Cfg')
            line_tmpl = '{:<23}|  ' + '{:13.3f}'*1 + '{}' + '{:13.3f}'*6

        data = []
        ind = 0
        for name in names:
            full_name = '{}.{}'.format(e_name, name)
            
            if name == "Cv":
                data.append(Cv_val)
            else:
                if name == "Cfg":
                    data.append(Cfg_val)
                else:
                    if names_units[ind] != '':
                        data.append(prob.get_val(full_name,units=names_units[ind])[0])
                    else:
                        data.append(prob[full_name][0])
            ind = ind +1 

        vals = [e_name] + data
        print(line_tmpl.format(*vals), file=file, flush=True)


def print_bleed(prob, element_names, file=sys.stdout, units='Default'):

    # user can specify the variable and units to be used below 
    # leaving the cell empty ('') means to use the default 
    display_names = ['Wb/Win', 'Pfrac', 'Workfrac',  'W'   , 'Tt'   , 'ht'   , 'Pt'   ]
    names         = [                              'stat:W', 'tot:T', 'tot:h', 'tot:P']
    
    if units == 'Default':
        print('Default')
        names_units   = [                           ''     ,  ''    , ''     ,  ''    ]   
    # define the SI units for each variables
    elif units=='SI': 
        names_units   = [                           'kg/s' ,  'K'   , 'kJ/kg', 'Pa'   ]   
    else:
        print('UNITS not recognized')
        return None        

    # configure & print header
    # get max name length:
    max_name_len = 0
    for e_name in element_names:
        # print('foo', e_name)
        bleed = prob.model._get_subsystem(e_name)
        for bn in bleed.options['bleed_names']:
            max_name_len = max(max_name_len, len(e_name+bn))

    max_name_len += 2

    len_header = max_name_len+3+7*13
    print("-"*len_header, file=file, flush=True)
    print("                            BLEED PROPERTIES", file=file, flush=True)
    print("-"*len_header, file=file, flush=True)

    # configure & print the variables and units
    max_name_len = str(max_name_len)
    line_tmpl = '{:<'+max_name_len+'}|  '+'{:>13}'*7
    print(line_tmpl.format(*(['Bleed'] + display_names)), file=file, flush=True)

    # NEW: extract the units of the variables from the full variables data structure
    # step1: create a dictionnary to link the full variable names and promoted names
    prom_name_to_abs = {
        meta['prom_name']: abs_name
        for abs_name, meta in prob.model.get_io_metadata(iotypes=('input', 'output')).items()
    }

    # step2: extract the units for the first flow station
    line_tmpl_units = '{:<'+max_name_len+'}|  ' + '{:>13}'*7
    e_name=element_names[0]
    data_units = []
    
    # loop through names of variables to get units
    ind = 0
    for name in names:
        full_name = '{}.{}'.format(e_name, name)
        
        # need to find the name of the bleed
        bleed = prob.model._get_subsystem(e_name)
        bleed_names = bleed.options['bleed_names']
        full_name = bleed.pathname + '.' + bleed_names[0] + ':' + name
        
        if names_units[ind] != '':
            data_units.append(names_units[ind])
        else:
            try: 
                data_units.append(prob.model.get_io_metadata(iotypes=('input','output'))[prom_name_to_abs[full_name]]['units'])
            except KeyError: # for data not described in the metadata - try the direct way
                data_units.append(prob.model.get_io_metadata(iotypes=('input','output'))[full_name]['units'])
            
        ind = ind + 1 
        
    vals_units = ['         ' + units + ' units'] + ["-"]*3 + data_units
    vals_units_clean = [str(v) if v is not None else '-' for v in vals_units]
    # step3: print the units as part of the output table
    print(line_tmpl_units.format(*vals_units_clean), file=file, flush=True)

    # process data
    print("-"*len_header, file=file, flush=True)
    # extract data for each element
    line_tmpl = '{:<'+max_name_len+'}|  '+'{:13.3f}'*7
    for e_name in element_names:
        bleed = prob.model._get_subsystem(e_name)

        bleed_names = bleed.options['bleed_names']

        for bn in bleed_names:

            full_bleed_name = bleed.pathname + '.' + bn

            try: # this one will work for interstage bleeds on compressor
                bld_pwr_path = bleed.pathname + '.blds_pwr.' + bn
                frac_p = prob[bld_pwr_path+':frac_P'][0]
                frac_work = prob[bld_pwr_path+':frac_work'][0]
                frac_W = prob[bld_pwr_path +':frac_W'][0]

            except KeyError: # for stand alone bleeds
                bld_clacs_path = bleed.pathname + '.bld_calcs.' + bn
                frac_W = prob[bld_clacs_path +':frac_W'][0]
                frac_p = 1.0
                frac_work = 1.0


            data = [frac_W, frac_p,frac_work]
            ind = 0
            for name in names:
                full_name = '{}:{}'.format(full_bleed_name, name)
                
                if names_units[ind] != '':
                    data.append(prob.get_val(full_name,units=names_units[ind])[0])
                else:
                    data.append(prob[full_name][0])
                ind = ind +1 


            vals = [full_bleed_name] + data
            print(line_tmpl.format(*vals), file=file, flush=True)
            
            # print(line_tmpl.format(full_bleed_name, frac_W, frac_p,
            #                        frac_work, prob[full_bleed_name+':stat:W'][0], prob[full_bleed_name+':tot:T'][0],
            #                        prob[full_bleed_name+':tot:h'][0], prob[full_bleed_name+':tot:P'][0]),
            #       file=file, flush=True)
        print("-"*len_header, file=file, flush=True)


def print_shaft(prob, element_names, file=sys.stdout, units='Default'):
    # user can specify the variable and units to be used below 
    # leaving the cell empty ('') means to use the default 
    display_names = ['Nmech', 'trq_in', 'trq_out', 'pwr_in', 'pwr_out']
    names         = display_names
    
    if units == 'Default':
        # print('Default')
        #  !!!! automated method to find default units does not work for shaft properties - define manually !!!!
        # names_units   = [ 'rpm' , 'ft*lbf', 'ft*lbf' ,   'hp'  ,    'hp'  ]    
        names_units   = [ ''    , ''      , ''       ,   ''    ,    ''    ]    
    # define the SI units for each variables
    elif units=='SI': 
        names_units   = [ 'rpm' , 'N*m'   ,  'N*m'   ,   'kW'  ,    'kW'  ]
    else:
        print('UNITS not recognized')
        return None    
    
    # configure & print header
    len_header = len_header = 27+20*5

    print("-"*len_header, file=file, flush=True)
    print("                            SHAFT PROPERTIES", file=file, flush=True)
    print("-"*len_header, file=file, flush=True)

    line_tmpl = '{:<23}|  '+'{:>20}'*5
    print(line_tmpl.format(*(['Shaft'] + names)), file=file)
    
    # NEW: extract the units of the variables from the full variables data structure
    # step1: create a dictionnary to link the full variable names and promoted names
    prom_name_to_abs = {
         meta['prom_name']: abs_name
         for abs_name, meta in prob.model.get_io_metadata(iotypes=('input', 'output')).items()
     }
    
    # step2: extract the units for the first flow station
    line_tmpl_units = '{:<23}|  '+'{:>20}'*5
    e_name=element_names[0]
    data_units = []
    
    # loop through names of variables to get units
    ind = 0
    for name in names:
        full_name = '{}.{}'.format(e_name, name)
        if names_units[ind] != '':
            data_units.append(names_units[ind])
        else:
            try: 
                data_units.append(prob.model.get_io_metadata(iotypes=('input','output'))[prom_name_to_abs[full_name]]['units'])
            except KeyError: 
                data_units.append(prob.model.get_io_metadata(iotypes=('input','output'))[full_name]['units'])
        ind = ind + 1 
        
    vals_units = ['         ' + units + ' units'] + data_units
    vals_units_clean = [str(v) if v is not None else '-' for v in vals_units]
    # step3: print the units as part of the output table
    print(line_tmpl_units.format(*vals_units_clean), file=file, flush=True)

    # process data
    print("-"*len_header, file=file, flush=True)

    line_tmpl = '{:<23}|  '+'{:20.3f}'*5
    for e_name in element_names:
        data = []
        ind = 0
        for name in names:
            full_name = '{}.{}'.format(e_name, name)
            
            if names_units[ind] != '':
                data.append(prob.get_val(full_name,units=names_units[ind])[0])
            else:
                data.append(prob[full_name][0])
            ind = ind +1 
        
        vals = [e_name] + data
        print(line_tmpl.format(*vals), file=file, flush=True)

    print("-"*len_header, file=file, flush=True)
    


def print_mixer(prob, element_names, file=sys.stdout, units='Default'):
    # units not implemented
    len_header = len_header = 27+20*6

    print("-"*len_header, file=file, flush=True)
    print("                            MIXER PROPERTIES", file=file, flush=True)
    print("-"*len_header, file=file, flush=True)

    line_tmpl = '{:<23}|  '+'{:>20}'*6
    print(line_tmpl.format('Mixer', 'balance.P_tot', 'designed_stream', 'Fl_calc:stat:P', 'Fl_calc:stat:area', 'Fl_calc:stat:MN', 'ER'),
          file=file, flush=True)

    line_tmpl = '{:<23}|  {:20.3f}{:^20}'+'{:20.3f}'*3
    for e_name in element_names:
        mixer = prob.model._get_subsystem(e_name)
        ds = mixer.options['designed_stream']
        if ds == 1:
            print(line_tmpl.format(e_name, prob[e_name+'.balance.P_tot'][0], 1,
                                   prob[e_name+'.Fl_I1_calc:stat:P'][0],
                                   prob[e_name+'.Fl_I1_calc:stat:area'][0],
                                   prob[e_name+'.Fl_I1_calc:stat:MN'][0]),
                                   prob[e_name+'.ER'][0],
                  file=file, flush=True)
        else:
            print(line_tmpl.format(e_name, prob[e_name+'.balance.P_tot'][0], 2,
                                   prob[e_name+'.Fl_I2_calc:stat:P'][0],
                                   prob[e_name+'.Fl_I2_calc:stat:area'][0],
                                   prob[e_name+'.Fl_I2_calc:stat:MN'][0]),
                                   prob[e_name+'.ER',][0],
                  file=file, flush=True)


def plot_compressor_maps(prob, element_names, eff_vals=np.array([0,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0]),alphas=[0]):

    for e_name in element_names:
        comp = prob.model._get_subsystem(e_name)
        map_data = comp.options['map_data']

        s_Wc = prob[e_name+'.s_Wc']
        s_PR = prob[e_name+'.s_PR']
        s_eff = prob[e_name+'.s_eff']
        s_Nc = prob[e_name+'.s_Nc']

        RlineMap, NcMap = np.meshgrid(map_data.RlineMap, map_data.NcMap, sparse=False)

        for alpha in alphas:
          scaled_PR = (map_data.PRmap[alpha,:,:] - 1.) * s_PR + 1.

          plt.figure(figsize=(11,8))
          Nc = plt.contour(map_data.WcMap[alpha,:,:]*s_Wc,scaled_PR,NcMap*s_Nc,colors='k',levels=map_data.NcMap*s_Nc)
          R = plt.contour(map_data.WcMap[alpha,:,:]*s_Wc,scaled_PR,RlineMap,colors='k',levels=map_data.RlineMap)
          eff = plt.contourf(map_data.WcMap[alpha,:,:]*s_Wc,scaled_PR,map_data.effMap[alpha,:,:]*s_eff,levels=eff_vals)

          plt.colorbar(eff)

          plt.plot(prob[e_name+'.Wc'], prob[e_name+'.map.scalars.PR'][0], 'ko')

          plt.clabel(Nc, fontsize=9, inline=False)
          plt.clabel(R, fontsize=9, inline=False)
          # plt.clabel(eff, fontsize=9, inline=True)

          plt.xlabel('Wc, lbm/s')
          plt.ylabel('PR')
          plt.title(e_name)
          # plt.show()
          plt.savefig(e_name+'.pdf')


def plot_turbine_maps(prob, element_names, eff_vals=np.array([0,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0]),alphas=[0]):

    for e_name in element_names:
        comp = prob.model._get_subsystem(e_name)
        map_data = comp.options['map_data']

        s_Wp = prob[e_name+'.s_Wp']
        s_PR = prob[e_name+'.s_PR']
        s_eff = prob[e_name+'.s_eff']
        s_Np = prob[e_name+'.s_Np']

        PRmap, NpMap = np.meshgrid(map_data.PRmap, map_data.NpMap, sparse=False)

        for alpha in alphas:
          scaled_PR = (PRmap - 1.) * s_PR + 1.

          plt.figure(figsize=(11,8))
          Nc = plt.contour(map_data.WpMap[alpha,:,:]*s_Wp,scaled_PR,NpMap*s_Np,colors='k',levels=map_data.NpMap*s_Np)
          eff = plt.contourf(map_data.WpMap[alpha,:,:]*s_Wp,scaled_PR,map_data.effMap[alpha,:,:]*s_eff,levels=eff_vals)

          plt.colorbar(eff)

          plt.plot(prob[e_name+'.Wp'], prob[e_name+'.map.scalars.PR'][0], 'ko')

          plt.clabel(Nc, fontsize=9, inline=False)
          # plt.clabel(eff, fontsize=9, inline=True)

          plt.xlabel('Wc, lbm/s')
          plt.ylabel('PR')
          plt.title(e_name)
          # plt.show()
          plt.savefig(e_name+'.pdf')


