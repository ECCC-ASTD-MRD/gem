module levels
   use dimout, only: MAXELEM, MAXSET
   implicit none
   public
   save
!______________________________________________________________________
!                                                                      |
!  VARIABLES FOR DEFINITION OF THE OUTPUT LEVELS (set_level)           |
!______________________________________________________________________|
!                    |                                                 |
! NAME               | DESCRIPTION                                     |
!--------------------|-------------------------------------------------|
! MAXLEV             | maximum number of output levels per set         |
! Level_sets         | total number of sets of defined output levels   |
! Level_id           | Level_id(i) are the id of each defined level set|
! Level              | contains the values of the defined output levels|
!                    | in mb and/or eta, ie: Level(i,j)=1000           |
! Level_max          | Level_max(i) contains the number of output      |
!                    | levels defined in each set                      |
! Level_typ_S        | level type: "P" for pressure in mb, "M" for eta |
!                    | which corresponds to each defined Level set     |
! Level_kind_ip1     | value of KIND to be used in CONVIP function     |
! Level_kind_diag    | value of KIND for the diagnostic level          |
! Level_momentum     | # of levels in momentum including the surface   |
! Level_thermo       | # of levels in momentum including the surface   |
! Level_allpres      | all pressure levels requested across all sets   |
! Level_npres        | number of Level_allpres levels                  |
! Level_version      | Level version 2 is Staggered CP Girard V4       |
!----------------------------------------------------------------------
!
!

      integer, parameter :: MAXLEV = 201
      real Level(MAXLEV,MAXSET),Level_allpres(MAXSET*MAXLEV)
      integer Level_id(MAXSET),Level_sets,Level_thermo,Level_momentum
      integer Level_max(MAXSET),Level_kind_ip1,Level_kind_diag
      integer :: Level_version,Level_npres
      character(len=1) Level_typ_S(MAXSET)

end module levels
