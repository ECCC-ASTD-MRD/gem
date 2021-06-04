function ccard_arg_m(nom, ord) result(value)
    character(len = *), intent(in) :: nom
    integer, intent(in) :: ord

    character(len = 40) :: keyname
    character(len = 8192) :: value
    integer :: lng
    integer :: c_get_appl_var
    integer :: lnom
    external c_get_appl_var

    lnom = len(trim(nom))
    if ((nom(lnom:lnom) == '.') .or. (nom(lnom:lnom) == ':') .or. (nom(lnom:lnom) == '_')) lnom = lnom-1

    write(keyname, '(a,i4.4,a)') '%%'//nom(1:lnom), ord, '%%'
    lng = c_get_appl_var(keyname, value)
end function ccard_arg_m


function ccard_arg(nom) result(value)
    character(len=*), intent(in) :: nom

    character(len = 40) :: keyname
    character(len = 8192) :: value
    integer :: lng

    integer :: c_get_appl_var
    external c_get_appl_var

    lnom = len(trim(nom))
    if ((nom(lnom:lnom) == '.') .or. (nom(lnom:lnom) == ':') .or. (nom(lnom:lnom) == '_')) lnom = lnom-1

    write(keyname, '(a,i4.4,a)') '%%'//nom(1:lnom), 0, '%%'
    lng = c_get_appl_var(keyname, value)
end function ccard_arg
