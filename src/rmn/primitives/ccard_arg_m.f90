    character *8192 function ccard_arg_m(nom,ord)
    character *(*) nom
    integer ord
    
    character * 40 keyname
    character *8192 value
    integer lng, c_get_appl_var,L
    external c_get_appl_var
    
    L = len(trim(nom))
    if ((nom(L:L) == '.') .or. (nom(L:L) == ':') .or. (nom(L:L) == '_')) L = L-1
    write(keyname,77) '%%'//nom(1:L),ord,'%%'
 77 format(a,i4.4,a)
    lng = c_get_appl_var(keyname,value)
    ccard_arg_m = value
    return
    end
    
    character *8192 function ccard_arg(nom)
    character *(*) nom
    integer ord
    
    character * 40 keyname
    character *8192 value
    integer lng, c_get_appl_var
    external c_get_appl_var
    
    L = len(trim(nom))
    if ((nom(L:L) == '.') .or. (nom(L:L) == ':') .or. (nom(L:L) == '_')) L = L-1
    write(keyname,77) '%%'//nom(1:L),0,'%%'
 77 format(a,i4.4,a)
    lng = c_get_appl_var(keyname,value)
    ccard_arg = value
    return
    end
    