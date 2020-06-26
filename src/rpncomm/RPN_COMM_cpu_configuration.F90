! RPN_COMM - Library of useful routines for C and FORTRAN programming
! Copyright (C) 2019  Division de Recherche en Prevision Numerique
!                     Environnement Canada
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation,
! version 2.1 of the License.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the
! Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! Boston, MA 02111-1307, USA.
!
module node_lcpu_configuration
  save
  integer :: logical_cpus     = 0          ! hyperthreads * number of physical cores
  integer :: sockets_per_node = 2          ! number of sockets in node
  integer :: numa_spaces      = 1          ! number of numa spaces per socket
  integer :: hyperthreads     = 1          ! number of logical threads per physical core
  integer :: thread_mapping   = 0          ! 0  : all ht 0, all ht 1 , ...
                                           ! 1  : hyperthreads of same core have consecutive numbers
end module node_lcpu_configuration

subroutine print_logical_cpu_configuration
  use node_lcpu_configuration
  implicit none
  if(logical_cpus == 0) call init_logical_cpu_configuration
  print 1,logical_cpus,' logical cpus,',sockets_per_node,' socket(s),',numa_spaces,' numa space(s),',hyperthreads,' thread(s) per core, mapping type ',thread_mapping
1 format(4(I4,A),I2)
end subroutine print_logical_cpu_configuration

subroutine set_logical_cpu_configuration(lcpus, nsockets, nnuma, ht, map)
  use node_lcpu_configuration
  implicit none
  integer, intent(IN) :: lcpus, nsockets, nnuma, ht, map
  logical_cpus = lcpus
  sockets_per_node = nsockets
  numa_spaces = nnuma
  hyperthreads = ht
  thread_mapping = map
end subroutine set_logical_cpu_configuration

subroutine get_logical_cpu_configuration(lcpus, nsockets, nnuma, ht, map)
  use node_lcpu_configuration
  implicit none
  integer, intent(OUT) :: lcpus, nsockets, nnuma, ht, map
  if(logical_cpus == 0) call init_logical_cpu_configuration
  lcpus = logical_cpus
  nsockets = sockets_per_node
  nnuma = numa_spaces
  ht = hyperthreads
  map = thread_mapping
end subroutine get_logical_cpu_configuration

subroutine init_logical_cpu_configuration
  use node_lcpu_configuration
  implicit none
  character(len=1024) :: node_cpu_config
  integer :: length, status
  call get_environment_variable("NODE_CPU_CONFIG",node_cpu_config,length, status)
  if(status .ne. 0) return
  read(node_cpu_config,*,err=1)logical_cpus, sockets_per_node, numa_spaces, hyperthreads, thread_mapping
1 continue
end subroutine init_logical_cpu_configuration
