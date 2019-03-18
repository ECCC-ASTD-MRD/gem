subroutine exp_cal_vor ( F_QR,F_QQ, F_uu,F_vv          , &
                         F_filtqq, F_coefqq, F_absvor_L, &
                         Minx,Maxx,Miny,Maxy,Nk )
logical  F_absvor_L
integer  F_filtqq, Minx,Maxx,Miny,Maxy,Nk
real     F_QR (Minx:Maxx,Miny:Maxy,Nk), &
         F_QQ (Minx:Maxx,Miny:Maxy,Nk), &
         F_uu (Minx:Maxx,Miny:Maxy,Nk), &
         F_vv (Minx:Maxx,Miny:Maxy,Nk), F_coefqq
return
end


subroutine exp_dynstep()
return
end

integer function exp_nml(F_namelistf_S)
character(len=*) F_namelistf_S
exp_nml = 0
return
end

subroutine exp_set_vt()
return
end
