ifneq (,$(DEBUGMAKE))
$(info ## ====================================================================)
$(info ## File:/tmp/viv001/31061/dependencies+.mk)
$(info ## )
endif
TMPL90DECKS=
INCLUDEDECKS=
PTN90DECKS=
FDECKS=
HFDECKS= \
	RPN_COMM_adj_halo.hf  
F03DECKS=
PTNDECKS=
INCDECKS= \
	RPN_COMM.inc  RPN_COMM_constants.inc  RPN_COMM_ftoc.inc  RPN_COMM_int.inc  RPN_COMM_is_null.inc   \
	RPN_COMM_legacy.inc  RPN_COMM_spread.inc  RPN_COMM_times.inc  RPN_COMM_times_post.inc  RPN_COMM_types.inc   \
	RPN_COMM_types_int.inc  RPN_COMM_version.inc  
HDECKS= \
	RPN_COMM_mpi.h  RPN_COMM_mpif.h  
FTN90DECKS=
CDECKS= \
	RPN_COMM_c_io.c  RPN_COMM_chdir.c  RPN_COMM_f2c.c  RPN_COMM_getenv.c  RPN_COMM_manage_halo.c   \
	RPN_COMM_shmget.c  RPN_COMM_sock.c  RPN_COMM_unbind_process.c  RPN_COMM_wtime.c  
F95DECKS=
F90DECKS= \
	RPN_COMM_C_interface.F90  RPN_COMM_WTIME.F90  RPN_COMM_Xpose1.F90  RPN_COMM_Xpose148.F90  RPN_COMM_Xpose2.F90   \
	RPN_COMM_Xpose248.F90  RPN_COMM_adj_halo.F90  RPN_COMM_adj_halo48.F90  RPN_COMM_allgather.F90  RPN_COMM_allreduce.F90   \
	RPN_COMM_alltoall.F90  RPN_COMM_barrier.F90  RPN_COMM_bcast.F90  RPN_COMM_bcst_world.F90  RPN_COMM_bcst_x.F90   \
	RPN_COMM_bloc.F90  RPN_COMM_bloctopo.F90  RPN_COMM_carac.F90  RPN_COMM_coll.F90  RPN_COMM_colors.F90   \
	RPN_COMM_comm.F90  RPN_COMM_compress.F90  RPN_COMM_const.F90  RPN_COMM_datyp.F90  RPN_COMM_debug.F90   \
	RPN_COMM_defo.F90  RPN_COMM_diag.F90  RPN_COMM_dist.F90  RPN_COMM_ezwin.F90  RPN_COMM_finalize.F90   \
	RPN_COMM_gather.F90  RPN_COMM_gather_e1y3.F90  RPN_COMM_gatherv.F90  RPN_COMM_globalsum.F90  RPN_COMM_grank.F90   \
	RPN_COMM_grid_2d.F90  RPN_COMM_grid_redist.F90  RPN_COMM_group.F90  RPN_COMM_halo.F90  RPN_COMM_haloflip.F90   \
	RPN_COMM_init.F90  RPN_COMM_io.F90  RPN_COMM_io_dist_coll.F90  RPN_COMM_io_pe_list.F90  RPN_COMM_io_pes.F90   \
	RPN_COMM_irecv.F90  RPN_COMM_is_null.F90  RPN_COMM_isend.F90  RPN_COMM_limit.F90  RPN_COMM_local_dist.F90   \
	RPN_COMM_low2up.F90  RPN_COMM_mod.F90  RPN_COMM_move.F90  RPN_COMM_mype.F90  RPN_COMM_openmp.F90   \
	RPN_COMM_oper.F90  RPN_COMM_options.F90  RPN_COMM_optn.F90  RPN_COMM_petopo.F90  RPN_COMM_progatate_boundary.F90   \
	RPN_COMM_ptr.F90  RPN_COMM_qadl.F90  RPN_COMM_rank.F90  RPN_COMM_recv.F90  RPN_COMM_reduce.F90   \
	RPN_COMM_scatterv.F90  RPN_COMM_send.F90  RPN_COMM_sendrecv.F90  RPN_COMM_size.F90  RPN_COMM_split.F90   \
	RPN_COMM_split_by_node.F90  RPN_COMM_spread.F90  RPN_COMM_status_size.F90  RPN_COMM_swapns.F90  RPN_COMM_topo.F90   \
	RPN_COMM_topo_xy.F90  RPN_COMM_transpose.F90  RPN_COMM_transpose48.F90  RPN_COMM_transpose_utils.F90  RPN_COMM_types.F90   \
	RPN_COMM_unit.F90  RPN_COMM_version.F90  RPN_COMM_waitall_nostat.F90  RPN_COMM_window.F90  RPN_COMM_xch_halo.F90   \
	RPN_COMM_xch_halo_8.F90  RPN_COMM_xch_haloew.F90  RPN_COMM_xch_halon.F90  RPN_COMM_xch_halons.F90  RPN_COMM_xch_halosl.F90   \
	RPN_COMM_xch_halox.F90  
FHDECKS=
CDK90DECKS=
FTNDECKS=
CDKDECKS=
OBJECTS= \
	RPN_COMM_C_interface.o  RPN_COMM_WTIME.o  RPN_COMM_Xpose1.o  RPN_COMM_Xpose148.o  RPN_COMM_Xpose2.o   \
	RPN_COMM_Xpose248.o  RPN_COMM_adj_halo.o  RPN_COMM_adj_halo48.o  RPN_COMM_allgather.o  RPN_COMM_allreduce.o   \
	RPN_COMM_alltoall.o  RPN_COMM_barrier.o  RPN_COMM_bcast.o  RPN_COMM_bcst_world.o  RPN_COMM_bcst_x.o   \
	RPN_COMM_bloc.o  RPN_COMM_bloctopo.o  RPN_COMM_c_io.o  RPN_COMM_carac.o  RPN_COMM_chdir.o   \
	RPN_COMM_coll.o  RPN_COMM_colors.o  RPN_COMM_comm.o  RPN_COMM_compress.o  RPN_COMM_const.o   \
	RPN_COMM_datyp.o  RPN_COMM_debug.o  RPN_COMM_defo.o  RPN_COMM_diag.o  RPN_COMM_dist.o   \
	RPN_COMM_ezwin.o  RPN_COMM_f2c.o  RPN_COMM_finalize.o  RPN_COMM_gather.o  RPN_COMM_gather_e1y3.o   \
	RPN_COMM_gatherv.o  RPN_COMM_getenv.o  RPN_COMM_globalsum.o  RPN_COMM_grank.o  RPN_COMM_grid_2d.o   \
	RPN_COMM_grid_redist.o  RPN_COMM_group.o  RPN_COMM_halo.o  RPN_COMM_haloflip.o  RPN_COMM_init.o   \
	RPN_COMM_io.o  RPN_COMM_io_dist_coll.o  RPN_COMM_io_pe_list.o  RPN_COMM_io_pes.o  RPN_COMM_irecv.o   \
	RPN_COMM_is_null.o  RPN_COMM_isend.o  RPN_COMM_limit.o  RPN_COMM_local_dist.o  RPN_COMM_low2up.o   \
	RPN_COMM_manage_halo.o  RPN_COMM_mod.o  RPN_COMM_move.o  RPN_COMM_mype.o  RPN_COMM_openmp.o   \
	RPN_COMM_oper.o  RPN_COMM_options.o  RPN_COMM_optn.o  RPN_COMM_petopo.o  RPN_COMM_progatate_boundary.o   \
	RPN_COMM_ptr.o  RPN_COMM_qadl.o  RPN_COMM_rank.o  RPN_COMM_recv.o  RPN_COMM_reduce.o   \
	RPN_COMM_scatterv.o  RPN_COMM_send.o  RPN_COMM_sendrecv.o  RPN_COMM_shmget.o  RPN_COMM_size.o   \
	RPN_COMM_sock.o  RPN_COMM_split.o  RPN_COMM_split_by_node.o  RPN_COMM_spread.o  RPN_COMM_status_size.o   \
	RPN_COMM_swapns.o  RPN_COMM_topo.o  RPN_COMM_topo_xy.o  RPN_COMM_transpose.o  RPN_COMM_transpose48.o   \
	RPN_COMM_transpose_utils.o  RPN_COMM_types.o  RPN_COMM_unbind_process.o  RPN_COMM_unit.o  RPN_COMM_version.o   \
	RPN_COMM_waitall_nostat.o  RPN_COMM_window.o  RPN_COMM_wtime.o  RPN_COMM_xch_halo.o  RPN_COMM_xch_halo_8.o   \
	RPN_COMM_xch_haloew.o  RPN_COMM_xch_halon.o  RPN_COMM_xch_halons.o  RPN_COMM_xch_halosl.o  RPN_COMM_xch_halox.o  

TOPDIRLIST=
TOPDIRLIST_NAMES=
OBJECTS_all= \
	RPN_COMM_C_interface.o  RPN_COMM_WTIME.o  RPN_COMM_Xpose1.o  RPN_COMM_Xpose148.o  RPN_COMM_Xpose2.o   \
	RPN_COMM_Xpose248.o  RPN_COMM_adj_halo.o  RPN_COMM_adj_halo48.o  RPN_COMM_allgather.o  RPN_COMM_allreduce.o   \
	RPN_COMM_alltoall.o  RPN_COMM_barrier.o  RPN_COMM_bcast.o  RPN_COMM_bcst_world.o  RPN_COMM_bcst_x.o   \
	RPN_COMM_bloc.o  RPN_COMM_bloctopo.o  RPN_COMM_c_io.o  RPN_COMM_carac.o  RPN_COMM_chdir.o   \
	RPN_COMM_coll.o  RPN_COMM_colors.o  RPN_COMM_comm.o  RPN_COMM_compress.o  RPN_COMM_const.o   \
	RPN_COMM_datyp.o  RPN_COMM_debug.o  RPN_COMM_defo.o  RPN_COMM_diag.o  RPN_COMM_dist.o   \
	RPN_COMM_ezwin.o  RPN_COMM_f2c.o  RPN_COMM_finalize.o  RPN_COMM_gather.o  RPN_COMM_gather_e1y3.o   \
	RPN_COMM_gatherv.o  RPN_COMM_getenv.o  RPN_COMM_globalsum.o  RPN_COMM_grank.o  RPN_COMM_grid_2d.o   \
	RPN_COMM_grid_redist.o  RPN_COMM_group.o  RPN_COMM_halo.o  RPN_COMM_haloflip.o  RPN_COMM_init.o   \
	RPN_COMM_io.o  RPN_COMM_io_dist_coll.o  RPN_COMM_io_pe_list.o  RPN_COMM_io_pes.o  RPN_COMM_irecv.o   \
	RPN_COMM_is_null.o  RPN_COMM_isend.o  RPN_COMM_limit.o  RPN_COMM_local_dist.o  RPN_COMM_low2up.o   \
	RPN_COMM_manage_halo.o  RPN_COMM_mod.o  RPN_COMM_move.o  RPN_COMM_mype.o  RPN_COMM_openmp.o   \
	RPN_COMM_oper.o  RPN_COMM_options.o  RPN_COMM_optn.o  RPN_COMM_petopo.o  RPN_COMM_progatate_boundary.o   \
	RPN_COMM_ptr.o  RPN_COMM_qadl.o  RPN_COMM_rank.o  RPN_COMM_recv.o  RPN_COMM_reduce.o   \
	RPN_COMM_scatterv.o  RPN_COMM_send.o  RPN_COMM_sendrecv.o  RPN_COMM_shmget.o  RPN_COMM_size.o   \
	RPN_COMM_sock.o  RPN_COMM_split.o  RPN_COMM_split_by_node.o  RPN_COMM_spread.o  RPN_COMM_status_size.o   \
	RPN_COMM_swapns.o  RPN_COMM_topo.o  RPN_COMM_topo_xy.o  RPN_COMM_transpose.o  RPN_COMM_transpose48.o   \
	RPN_COMM_transpose_utils.o  RPN_COMM_types.o  RPN_COMM_unbind_process.o  RPN_COMM_unit.o  RPN_COMM_version.o   \
	RPN_COMM_waitall_nostat.o  RPN_COMM_window.o  RPN_COMM_wtime.o  RPN_COMM_xch_halo.o  RPN_COMM_xch_halo_8.o   \
	RPN_COMM_xch_haloew.o  RPN_COMM_xch_halon.o  RPN_COMM_xch_halons.o  RPN_COMM_xch_halosl.o  RPN_COMM_xch_halox.o  
$(LIBDIR)/liball.a: $(OBJECTS_all) $(LIBDEP_all) $(LIBDEP_ALL)
	rm -f $@; ar r $@_$$$$ $(OBJECTS_all); mv $@_$$$$ $@
liball.a: $(LIBDIR)/liball.a

ALL_LIBS= \
	$(LIBDIR)/liball.a  

FORTRAN_MODULES= \
	rpn_comm  rpn_comm_barrier_priv  rpn_comm_bloc_mgt  rpn_comm_ezwin_mod  rpn_comm_grids   \
	rpn_comm_halos  rpn_comm_io  rpn_comm_io_pe_tables  rpn_comm_localdist  rpn_comm_test   \
	rpn_comm_transpose_utils  rpn_comm_types  rpn_comm_windows  rpncomm_com  split_by_node  

FMOD_FILE_rpn_comm = RPN_COMM_mod.F90
FMOD_FILE_rpn_comm_barrier_priv = RPN_COMM_barrier.F90
FMOD_FILE_rpn_comm_bloc_mgt = RPN_COMM_bloc.F90
FMOD_FILE_rpn_comm_ezwin_mod = RPN_COMM_ezwin.F90
FMOD_FILE_rpn_comm_grids = RPN_COMM_grid_2d.F90
FMOD_FILE_rpn_comm_halos = RPN_COMM_halo.F90
FMOD_FILE_rpn_comm_io = RPN_COMM_io.F90
FMOD_FILE_rpn_comm_io_pe_tables = RPN_COMM_io_pes.F90
FMOD_FILE_rpn_comm_localdist = RPN_COMM_local_dist.F90
FMOD_FILE_rpn_comm_test = RPN_COMM_progatate_boundary.F90
FMOD_FILE_rpn_comm_transpose_utils = RPN_COMM_transpose_utils.F90
FMOD_FILE_rpn_comm_types = RPN_COMM_types.F90
FMOD_FILE_rpn_comm_windows = RPN_COMM_window.F90
FMOD_FILE_rpncomm_com = RPN_COMM_comm.F90
FMOD_FILE_split_by_node = RPN_COMM_split_by_node.F90

FMOD_LIST_RPN_COMM_barrier.F90= \
	rpn_comm_barrier_priv  
FMOD_LIST_RPN_COMM_bloc.F90= \
	rpn_comm_bloc_mgt  
FMOD_LIST_RPN_COMM_comm.F90= \
	rpncomm_com  
FMOD_LIST_RPN_COMM_ezwin.F90= \
	rpn_comm_ezwin_mod  
FMOD_LIST_RPN_COMM_grid_2d.F90= \
	rpn_comm_grids  
FMOD_LIST_RPN_COMM_halo.F90= \
	rpn_comm_halos  
FMOD_LIST_RPN_COMM_io.F90= \
	rpn_comm_io  
FMOD_LIST_RPN_COMM_io_pes.F90= \
	rpn_comm_io_pe_tables  
FMOD_LIST_RPN_COMM_local_dist.F90= \
	rpn_comm_localdist  
FMOD_LIST_RPN_COMM_mod.F90= \
	rpn_comm  
FMOD_LIST_RPN_COMM_progatate_boundary.F90= \
	rpn_comm_test  
FMOD_LIST_RPN_COMM_split_by_node.F90= \
	split_by_node  
FMOD_LIST_RPN_COMM_transpose_utils.F90= \
	rpn_comm_transpose_utils  
FMOD_LIST_RPN_COMM_types.F90= \
	rpn_comm_types  
FMOD_LIST_RPN_COMM_window.F90= \
	rpn_comm_windows  


$(eval MYVAR2 = $$($(MYVAR)))
echo_mydepvar:
	echo $(MYVAR2)


RPN_COMM_C_interface.o:	RPN_COMM_C_interface.F90 \
	RPN_COMM_mod.o  
RPN_COMM_WTIME.o:	RPN_COMM_WTIME.F90
RPN_COMM_Xpose1.o:	RPN_COMM_Xpose1.F90 \
	RPN_COMM_mod.o  
RPN_COMM_Xpose148.o:	RPN_COMM_Xpose148.F90 \
	RPN_COMM_mod.o  
RPN_COMM_Xpose2.o:	RPN_COMM_Xpose2.F90 \
	RPN_COMM_mod.o  
RPN_COMM_Xpose248.o:	RPN_COMM_Xpose248.F90 \
	RPN_COMM_mod.o  
RPN_COMM_adj_halo.o:	RPN_COMM_adj_halo.F90 \
	RPN_COMM_mod.o  
RPN_COMM_adj_halo48.o:	RPN_COMM_adj_halo48.F90 \
	RPN_COMM_adj_halo.hf  RPN_COMM_mod.o  
RPN_COMM_allgather.o:	RPN_COMM_allgather.F90
RPN_COMM_allreduce.o:	RPN_COMM_allreduce.F90
RPN_COMM_alltoall.o:	RPN_COMM_alltoall.F90
RPN_COMM_barrier.o:	RPN_COMM_barrier.F90
RPN_COMM_bcast.o:	RPN_COMM_bcast.F90
RPN_COMM_bcst_world.o:	RPN_COMM_bcst_world.F90 \
	RPN_COMM_mod.o  
RPN_COMM_bcst_x.o:	RPN_COMM_bcst_x.F90
RPN_COMM_bloc.o:	RPN_COMM_bloc.F90 \
	RPN_COMM_comm.o  RPN_COMM_mod.o  
RPN_COMM_bloctopo.o:	RPN_COMM_bloctopo.F90 \
	RPN_COMM_mod.o  
RPN_COMM_c_io.o:	RPN_COMM_c_io.c
RPN_COMM_carac.o:	RPN_COMM_carac.F90 \
	RPN_COMM_mod.o  
RPN_COMM_chdir.o:	RPN_COMM_chdir.c
RPN_COMM_coll.o:	RPN_COMM_coll.F90 \
	RPN_COMM_mod.o  
RPN_COMM_colors.o:	RPN_COMM_colors.F90 \
	RPN_COMM_mod.o  
RPN_COMM_comm.o:	RPN_COMM_comm.F90 \
	RPN_COMM_mod.o  
RPN_COMM_compress.o:	RPN_COMM_compress.F90
RPN_COMM_const.o:	RPN_COMM_const.F90 \
	RPN_COMM_mod.o  
RPN_COMM_datyp.o:	RPN_COMM_datyp.F90 \
	RPN_COMM_mod.o  
RPN_COMM_debug.o:	RPN_COMM_debug.F90
RPN_COMM_defo.o:	RPN_COMM_defo.F90 \
	RPN_COMM_comm.o  
RPN_COMM_diag.o:	RPN_COMM_diag.F90 \
	RPN_COMM_mod.o  
RPN_COMM_dist.o:	RPN_COMM_dist.F90 \
	RPN_COMM_mod.o  
RPN_COMM_ezwin.o:	RPN_COMM_ezwin.F90 \
	RPN_COMM.inc  RPN_COMM_constants.inc  RPN_COMM_is_null.inc  RPN_COMM_types.inc  
RPN_COMM_f2c.o:	RPN_COMM_f2c.c
RPN_COMM_finalize.o:	RPN_COMM_finalize.F90 \
	RPN_COMM_mod.o  
RPN_COMM_gather.o:	RPN_COMM_gather.F90
RPN_COMM_gather_e1y3.o:	RPN_COMM_gather_e1y3.F90 \
	RPN_COMM_mod.o  
RPN_COMM_gatherv.o:	RPN_COMM_gatherv.F90
RPN_COMM_getenv.o:	RPN_COMM_getenv.c
RPN_COMM_globalsum.o:	RPN_COMM_globalsum.F90 \
	RPN_COMM_mod.o  
RPN_COMM_grank.o:	RPN_COMM_grank.F90
RPN_COMM_grid_2d.o:	RPN_COMM_grid_2d.F90 \
	RPN_COMM_constants.inc  RPN_COMM_mod.o  
RPN_COMM_grid_redist.o:	RPN_COMM_grid_redist.F90 \
	RPN_COMM_mod.o  
RPN_COMM_group.o:	RPN_COMM_group.F90 \
	RPN_COMM_mod.o  
RPN_COMM_halo.o:	RPN_COMM_halo.F90
RPN_COMM_haloflip.o:	RPN_COMM_haloflip.F90 \
	RPN_COMM_mod.o  
RPN_COMM_init.o:	RPN_COMM_init.F90 \
	RPN_COMM_mod.o  
RPN_COMM_io.o:	RPN_COMM_io.F90
RPN_COMM_io_dist_coll.o:	RPN_COMM_io_dist_coll.F90 \
	RPN_COMM_io_pes.o  RPN_COMM_mod.o  
RPN_COMM_io_pe_list.o:	RPN_COMM_io_pe_list.F90
RPN_COMM_io_pes.o:	RPN_COMM_io_pes.F90 \
	RPN_COMM.inc  RPN_COMM_constants.inc  RPN_COMM_is_null.inc  RPN_COMM_types.inc  RPN_COMM_mod.o  
RPN_COMM_irecv.o:	RPN_COMM_irecv.F90
RPN_COMM_is_null.o:	RPN_COMM_is_null.F90 \
	RPN_COMM_types.inc  
RPN_COMM_isend.o:	RPN_COMM_isend.F90
RPN_COMM_limit.o:	RPN_COMM_limit.F90
RPN_COMM_local_dist.o:	RPN_COMM_local_dist.F90
RPN_COMM_low2up.o:	RPN_COMM_low2up.F90
RPN_COMM_manage_halo.o:	RPN_COMM_manage_halo.c
RPN_COMM_mod.o:	RPN_COMM_mod.F90 \
	RPN_COMM_constants.inc  RPN_COMM_types.inc  RPN_COMM_types_int.inc  
RPN_COMM_move.o:	RPN_COMM_move.F90 \
	RPN_COMM_mod.o  
RPN_COMM_mype.o:	RPN_COMM_mype.F90 \
	RPN_COMM_mod.o  
RPN_COMM_openmp.o:	RPN_COMM_openmp.F90
RPN_COMM_oper.o:	RPN_COMM_oper.F90 \
	RPN_COMM_mod.o  
RPN_COMM_options.o:	RPN_COMM_options.F90 \
	RPN_COMM_mod.o  
RPN_COMM_optn.o:	RPN_COMM_optn.F90 \
	RPN_COMM_mod.o  
RPN_COMM_petopo.o:	RPN_COMM_petopo.F90 \
	RPN_COMM_mod.o  
RPN_COMM_progatate_boundary.o:	RPN_COMM_progatate_boundary.F90 \
	RPN_COMM_mod.o  
RPN_COMM_ptr.o:	RPN_COMM_ptr.F90 \
	RPN_COMM_types.inc  
RPN_COMM_qadl.o:	RPN_COMM_qadl.F90 \
	RPN_COMM_mod.o  
RPN_COMM_rank.o:	RPN_COMM_rank.F90
RPN_COMM_recv.o:	RPN_COMM_recv.F90
RPN_COMM_reduce.o:	RPN_COMM_reduce.F90
RPN_COMM_scatterv.o:	RPN_COMM_scatterv.F90
RPN_COMM_send.o:	RPN_COMM_send.F90
RPN_COMM_sendrecv.o:	RPN_COMM_sendrecv.F90
RPN_COMM_shmget.o:	RPN_COMM_shmget.c
RPN_COMM_size.o:	RPN_COMM_size.F90
RPN_COMM_sock.o:	RPN_COMM_sock.c
RPN_COMM_split.o:	RPN_COMM_split.F90 \
	RPN_COMM_mod.o  
RPN_COMM_split_by_node.o:	RPN_COMM_split_by_node.F90 \
	RPN_COMM_constants.inc  
RPN_COMM_spread.o:	RPN_COMM_spread.F90 \
	RPN_COMM_int.inc  RPN_COMM_constants.inc  RPN_COMM_types.inc  RPN_COMM_types_int.inc  RPN_COMM_spread.inc  
RPN_COMM_status_size.o:	RPN_COMM_status_size.F90
RPN_COMM_swapns.o:	RPN_COMM_swapns.F90 \
	RPN_COMM_mod.o  
RPN_COMM_topo.o:	RPN_COMM_topo.F90 \
	RPN_COMM_mod.o  
RPN_COMM_topo_xy.o:	RPN_COMM_topo_xy.F90 \
	RPN_COMM_mod.o  
RPN_COMM_transpose.o:	RPN_COMM_transpose.F90 \
	RPN_COMM_mod.o  
RPN_COMM_transpose48.o:	RPN_COMM_transpose48.F90 \
	RPN_COMM_mod.o  
RPN_COMM_transpose_utils.o:	RPN_COMM_transpose_utils.F90 \
	RPN_COMM_mod.o  
RPN_COMM_types.o:	RPN_COMM_types.F90 \
	RPN_COMM_types.inc  
RPN_COMM_unbind_process.o:	RPN_COMM_unbind_process.c
RPN_COMM_unit.o:	RPN_COMM_unit.F90 \
	RPN_COMM_mod.o  
RPN_COMM_version.o:	RPN_COMM_version.F90 \
	RPN_COMM_version.inc  
RPN_COMM_waitall_nostat.o:	RPN_COMM_waitall_nostat.F90
RPN_COMM_window.o:	RPN_COMM_window.F90 \
	RPN_COMM.inc  RPN_COMM_constants.inc  RPN_COMM_is_null.inc  RPN_COMM_types.inc  RPN_COMM_types_int.inc  
RPN_COMM_wtime.o:	RPN_COMM_wtime.c
RPN_COMM_xch_halo.o:	RPN_COMM_xch_halo.F90 \
	RPN_COMM.inc  RPN_COMM_constants.inc  RPN_COMM_is_null.inc  RPN_COMM_types.inc  RPN_COMM_mod.o  
RPN_COMM_xch_halo_8.o:	RPN_COMM_xch_halo_8.F90
RPN_COMM_xch_haloew.o:	RPN_COMM_xch_haloew.F90 \
	RPN_COMM_mod.o  RPN_COMM_times.inc  
RPN_COMM_xch_halon.o:	RPN_COMM_xch_halon.F90 \
	RPN_COMM_mod.o  
RPN_COMM_xch_halons.o:	RPN_COMM_xch_halons.F90 \
	RPN_COMM_mod.o  
RPN_COMM_xch_halosl.o:	RPN_COMM_xch_halosl.F90 \
	RPN_COMM_mod.o  
RPN_COMM_xch_halox.o:	RPN_COMM_xch_halox.F90 \
	RPN_COMM_mod.o  
.PHONY: _invdep_.RPN_COMM.inc
_invdep_.RPN_COMM.inc: \
	RPN_COMM_io_pes.o  RPN_COMM_xch_halo.o  RPN_COMM_ezwin.o  RPN_COMM_window.o  RPN_COMM_io_dist_coll.o  
.PHONY: _invdep_.RPN_COMM_adj_halo.hf
_invdep_.RPN_COMM_adj_halo.hf: \
	RPN_COMM_adj_halo48.o  
.PHONY: _invdep_.RPN_COMM_comm.F90
_invdep_.RPN_COMM_comm.F90: \
	RPN_COMM_bloc.o  RPN_COMM_defo.o  
.PHONY: _invdep_.RPN_COMM_constants.inc
_invdep_.RPN_COMM_constants.inc: \
	RPN_COMM_io_pes.o  RPN_COMM_gather_e1y3.o  RPN_COMM_split.o  RPN_COMM_grid_2d.o  RPN_COMM_adj_halo.o   \
	RPN_COMM_datyp.o  RPN_COMM_colors.o  RPN_COMM_Xpose148.o  RPN_COMM_bloc.o  RPN_COMM_oper.o   \
	RPN_COMM_topo.o  RPN_COMM_unit.o  RPN_COMM_xch_halon.o  RPN_COMM_C_interface.o  RPN_COMM_xch_halo.o   \
	RPN_COMM_coll.o  RPN_COMM_Xpose1.o  RPN_COMM_diag.o  RPN_COMM_xch_halons.o  RPN_COMM_defo.o   \
	RPN_COMM_globalsum.o  RPN_COMM_move.o  RPN_COMM_adj_halo48.o  RPN_COMM_transpose48.o  RPN_COMM_bloctopo.o   \
	RPN_COMM_const.o  RPN_COMM_comm.o  RPN_COMM_options.o  RPN_COMM_mod.o  RPN_COMM_xch_halosl.o   \
	RPN_COMM_bcst_world.o  RPN_COMM_group.o  RPN_COMM_transpose_utils.o  RPN_COMM_spread.o  RPN_COMM_ezwin.o   \
	RPN_COMM_topo_xy.o  RPN_COMM_mype.o  RPN_COMM_qadl.o  RPN_COMM_Xpose2.o  RPN_COMM_transpose.o   \
	RPN_COMM_xch_halox.o  RPN_COMM_window.o  RPN_COMM_optn.o  RPN_COMM_carac.o  RPN_COMM_grid_redist.o   \
	RPN_COMM_finalize.o  RPN_COMM_haloflip.o  RPN_COMM_swapns.o  RPN_COMM_progatate_boundary.o  RPN_COMM_petopo.o   \
	RPN_COMM_Xpose248.o  RPN_COMM_xch_haloew.o  RPN_COMM_io_dist_coll.o  RPN_COMM_init.o  RPN_COMM_split_by_node.o   \
	RPN_COMM_dist.o  
.PHONY: _invdep_.RPN_COMM_int.inc
_invdep_.RPN_COMM_int.inc: \
	RPN_COMM_spread.o  
.PHONY: _invdep_.RPN_COMM_io_pes.F90
_invdep_.RPN_COMM_io_pes.F90: \
	RPN_COMM_io_dist_coll.o  
.PHONY: _invdep_.RPN_COMM_is_null.inc
_invdep_.RPN_COMM_is_null.inc: \
	RPN_COMM_io_pes.o  RPN_COMM_xch_halo.o  RPN_COMM_ezwin.o  RPN_COMM_window.o  RPN_COMM_io_dist_coll.o  
.PHONY: _invdep_.RPN_COMM_mod.F90
_invdep_.RPN_COMM_mod.F90: \
	RPN_COMM_io_pes.o  RPN_COMM_gather_e1y3.o  RPN_COMM_split.o  RPN_COMM_grid_2d.o  RPN_COMM_adj_halo.o   \
	RPN_COMM_datyp.o  RPN_COMM_colors.o  RPN_COMM_Xpose148.o  RPN_COMM_bloc.o  RPN_COMM_oper.o   \
	RPN_COMM_topo.o  RPN_COMM_unit.o  RPN_COMM_xch_halon.o  RPN_COMM_C_interface.o  RPN_COMM_xch_halo.o   \
	RPN_COMM_coll.o  RPN_COMM_Xpose1.o  RPN_COMM_diag.o  RPN_COMM_xch_halons.o  RPN_COMM_defo.o   \
	RPN_COMM_globalsum.o  RPN_COMM_move.o  RPN_COMM_adj_halo48.o  RPN_COMM_transpose48.o  RPN_COMM_bloctopo.o   \
	RPN_COMM_const.o  RPN_COMM_comm.o  RPN_COMM_options.o  RPN_COMM_xch_halosl.o  RPN_COMM_bcst_world.o   \
	RPN_COMM_group.o  RPN_COMM_transpose_utils.o  RPN_COMM_topo_xy.o  RPN_COMM_mype.o  RPN_COMM_qadl.o   \
	RPN_COMM_Xpose2.o  RPN_COMM_transpose.o  RPN_COMM_xch_halox.o  RPN_COMM_optn.o  RPN_COMM_carac.o   \
	RPN_COMM_grid_redist.o  RPN_COMM_finalize.o  RPN_COMM_haloflip.o  RPN_COMM_swapns.o  RPN_COMM_progatate_boundary.o   \
	RPN_COMM_petopo.o  RPN_COMM_Xpose248.o  RPN_COMM_xch_haloew.o  RPN_COMM_io_dist_coll.o  RPN_COMM_init.o   \
	RPN_COMM_dist.o  
.PHONY: _invdep_.RPN_COMM_spread.inc
_invdep_.RPN_COMM_spread.inc: \
	RPN_COMM_spread.o  
.PHONY: _invdep_.RPN_COMM_times.inc
_invdep_.RPN_COMM_times.inc: \
	RPN_COMM_xch_haloew.o  
.PHONY: _invdep_.RPN_COMM_types.inc
_invdep_.RPN_COMM_types.inc: \
	RPN_COMM_io_pes.o  RPN_COMM_gather_e1y3.o  RPN_COMM_split.o  RPN_COMM_grid_2d.o  RPN_COMM_adj_halo.o   \
	RPN_COMM_datyp.o  RPN_COMM_colors.o  RPN_COMM_Xpose148.o  RPN_COMM_bloc.o  RPN_COMM_ptr.o   \
	RPN_COMM_oper.o  RPN_COMM_topo.o  RPN_COMM_unit.o  RPN_COMM_xch_halon.o  RPN_COMM_C_interface.o   \
	RPN_COMM_types.o  RPN_COMM_xch_halo.o  RPN_COMM_coll.o  RPN_COMM_Xpose1.o  RPN_COMM_diag.o   \
	RPN_COMM_xch_halons.o  RPN_COMM_defo.o  RPN_COMM_globalsum.o  RPN_COMM_move.o  RPN_COMM_adj_halo48.o   \
	RPN_COMM_transpose48.o  RPN_COMM_bloctopo.o  RPN_COMM_const.o  RPN_COMM_comm.o  RPN_COMM_options.o   \
	RPN_COMM_mod.o  RPN_COMM_xch_halosl.o  RPN_COMM_bcst_world.o  RPN_COMM_group.o  RPN_COMM_transpose_utils.o   \
	RPN_COMM_spread.o  RPN_COMM_ezwin.o  RPN_COMM_topo_xy.o  RPN_COMM_mype.o  RPN_COMM_qadl.o   \
	RPN_COMM_Xpose2.o  RPN_COMM_transpose.o  RPN_COMM_xch_halox.o  RPN_COMM_window.o  RPN_COMM_optn.o   \
	RPN_COMM_carac.o  RPN_COMM_grid_redist.o  RPN_COMM_finalize.o  RPN_COMM_haloflip.o  RPN_COMM_swapns.o   \
	RPN_COMM_is_null.o  RPN_COMM_progatate_boundary.o  RPN_COMM_petopo.o  RPN_COMM_Xpose248.o  RPN_COMM_xch_haloew.o   \
	RPN_COMM_io_dist_coll.o  RPN_COMM_init.o  RPN_COMM_dist.o  
.PHONY: _invdep_.RPN_COMM_types_int.inc
_invdep_.RPN_COMM_types_int.inc: \
	RPN_COMM_io_pes.o  RPN_COMM_gather_e1y3.o  RPN_COMM_split.o  RPN_COMM_grid_2d.o  RPN_COMM_adj_halo.o   \
	RPN_COMM_datyp.o  RPN_COMM_colors.o  RPN_COMM_Xpose148.o  RPN_COMM_bloc.o  RPN_COMM_oper.o   \
	RPN_COMM_topo.o  RPN_COMM_unit.o  RPN_COMM_xch_halon.o  RPN_COMM_C_interface.o  RPN_COMM_xch_halo.o   \
	RPN_COMM_coll.o  RPN_COMM_Xpose1.o  RPN_COMM_diag.o  RPN_COMM_xch_halons.o  RPN_COMM_defo.o   \
	RPN_COMM_globalsum.o  RPN_COMM_move.o  RPN_COMM_adj_halo48.o  RPN_COMM_transpose48.o  RPN_COMM_bloctopo.o   \
	RPN_COMM_const.o  RPN_COMM_comm.o  RPN_COMM_options.o  RPN_COMM_mod.o  RPN_COMM_xch_halosl.o   \
	RPN_COMM_bcst_world.o  RPN_COMM_group.o  RPN_COMM_transpose_utils.o  RPN_COMM_spread.o  RPN_COMM_topo_xy.o   \
	RPN_COMM_mype.o  RPN_COMM_qadl.o  RPN_COMM_Xpose2.o  RPN_COMM_transpose.o  RPN_COMM_xch_halox.o   \
	RPN_COMM_window.o  RPN_COMM_optn.o  RPN_COMM_carac.o  RPN_COMM_grid_redist.o  RPN_COMM_finalize.o   \
	RPN_COMM_haloflip.o  RPN_COMM_swapns.o  RPN_COMM_progatate_boundary.o  RPN_COMM_petopo.o  RPN_COMM_Xpose248.o   \
	RPN_COMM_xch_haloew.o  RPN_COMM_io_dist_coll.o  RPN_COMM_init.o  RPN_COMM_dist.o  
.PHONY: _invdep_.RPN_COMM_version.inc
_invdep_.RPN_COMM_version.inc: \
	RPN_COMM_version.o  
INVDEP_LIST_RPN_COMM.inc= \
	RPN_COMM_io_pes.o  RPN_COMM_xch_halo.o  RPN_COMM_ezwin.o  RPN_COMM_window.o  RPN_COMM_io_dist_coll.o  
INVDEP_LIST_RPN_COMM_adj_halo.hf= \
	RPN_COMM_adj_halo48.o  
INVDEP_LIST_RPN_COMM_comm.F90= \
	RPN_COMM_bloc.o  RPN_COMM_defo.o  
INVDEP_LIST_RPN_COMM_constants.inc= \
	RPN_COMM_io_pes.o  RPN_COMM_gather_e1y3.o  RPN_COMM_split.o  RPN_COMM_grid_2d.o  RPN_COMM_adj_halo.o   \
	RPN_COMM_datyp.o  RPN_COMM_colors.o  RPN_COMM_Xpose148.o  RPN_COMM_bloc.o  RPN_COMM_oper.o   \
	RPN_COMM_topo.o  RPN_COMM_unit.o  RPN_COMM_xch_halon.o  RPN_COMM_C_interface.o  RPN_COMM_xch_halo.o   \
	RPN_COMM_coll.o  RPN_COMM_Xpose1.o  RPN_COMM_diag.o  RPN_COMM_xch_halons.o  RPN_COMM_defo.o   \
	RPN_COMM_globalsum.o  RPN_COMM_move.o  RPN_COMM_adj_halo48.o  RPN_COMM_transpose48.o  RPN_COMM_bloctopo.o   \
	RPN_COMM_const.o  RPN_COMM_comm.o  RPN_COMM_options.o  RPN_COMM_mod.o  RPN_COMM_xch_halosl.o   \
	RPN_COMM_bcst_world.o  RPN_COMM_group.o  RPN_COMM_transpose_utils.o  RPN_COMM_spread.o  RPN_COMM_ezwin.o   \
	RPN_COMM_topo_xy.o  RPN_COMM_mype.o  RPN_COMM_qadl.o  RPN_COMM_Xpose2.o  RPN_COMM_transpose.o   \
	RPN_COMM_xch_halox.o  RPN_COMM_window.o  RPN_COMM_optn.o  RPN_COMM_carac.o  RPN_COMM_grid_redist.o   \
	RPN_COMM_finalize.o  RPN_COMM_haloflip.o  RPN_COMM_swapns.o  RPN_COMM_progatate_boundary.o  RPN_COMM_petopo.o   \
	RPN_COMM_Xpose248.o  RPN_COMM_xch_haloew.o  RPN_COMM_io_dist_coll.o  RPN_COMM_init.o  RPN_COMM_split_by_node.o   \
	RPN_COMM_dist.o  
INVDEP_LIST_RPN_COMM_int.inc= \
	RPN_COMM_spread.o  
INVDEP_LIST_RPN_COMM_io_pes.F90= \
	RPN_COMM_io_dist_coll.o  
INVDEP_LIST_RPN_COMM_is_null.inc= \
	RPN_COMM_io_pes.o  RPN_COMM_xch_halo.o  RPN_COMM_ezwin.o  RPN_COMM_window.o  RPN_COMM_io_dist_coll.o  
INVDEP_LIST_RPN_COMM_mod.F90= \
	RPN_COMM_io_pes.o  RPN_COMM_gather_e1y3.o  RPN_COMM_split.o  RPN_COMM_grid_2d.o  RPN_COMM_adj_halo.o   \
	RPN_COMM_datyp.o  RPN_COMM_colors.o  RPN_COMM_Xpose148.o  RPN_COMM_bloc.o  RPN_COMM_oper.o   \
	RPN_COMM_topo.o  RPN_COMM_unit.o  RPN_COMM_xch_halon.o  RPN_COMM_C_interface.o  RPN_COMM_xch_halo.o   \
	RPN_COMM_coll.o  RPN_COMM_Xpose1.o  RPN_COMM_diag.o  RPN_COMM_xch_halons.o  RPN_COMM_defo.o   \
	RPN_COMM_globalsum.o  RPN_COMM_move.o  RPN_COMM_adj_halo48.o  RPN_COMM_transpose48.o  RPN_COMM_bloctopo.o   \
	RPN_COMM_const.o  RPN_COMM_comm.o  RPN_COMM_options.o  RPN_COMM_xch_halosl.o  RPN_COMM_bcst_world.o   \
	RPN_COMM_group.o  RPN_COMM_transpose_utils.o  RPN_COMM_topo_xy.o  RPN_COMM_mype.o  RPN_COMM_qadl.o   \
	RPN_COMM_Xpose2.o  RPN_COMM_transpose.o  RPN_COMM_xch_halox.o  RPN_COMM_optn.o  RPN_COMM_carac.o   \
	RPN_COMM_grid_redist.o  RPN_COMM_finalize.o  RPN_COMM_haloflip.o  RPN_COMM_swapns.o  RPN_COMM_progatate_boundary.o   \
	RPN_COMM_petopo.o  RPN_COMM_Xpose248.o  RPN_COMM_xch_haloew.o  RPN_COMM_io_dist_coll.o  RPN_COMM_init.o   \
	RPN_COMM_dist.o  
INVDEP_LIST_RPN_COMM_spread.inc= \
	RPN_COMM_spread.o  
INVDEP_LIST_RPN_COMM_times.inc= \
	RPN_COMM_xch_haloew.o  
INVDEP_LIST_RPN_COMM_types.inc= \
	RPN_COMM_io_pes.o  RPN_COMM_gather_e1y3.o  RPN_COMM_split.o  RPN_COMM_grid_2d.o  RPN_COMM_adj_halo.o   \
	RPN_COMM_datyp.o  RPN_COMM_colors.o  RPN_COMM_Xpose148.o  RPN_COMM_bloc.o  RPN_COMM_ptr.o   \
	RPN_COMM_oper.o  RPN_COMM_topo.o  RPN_COMM_unit.o  RPN_COMM_xch_halon.o  RPN_COMM_C_interface.o   \
	RPN_COMM_types.o  RPN_COMM_xch_halo.o  RPN_COMM_coll.o  RPN_COMM_Xpose1.o  RPN_COMM_diag.o   \
	RPN_COMM_xch_halons.o  RPN_COMM_defo.o  RPN_COMM_globalsum.o  RPN_COMM_move.o  RPN_COMM_adj_halo48.o   \
	RPN_COMM_transpose48.o  RPN_COMM_bloctopo.o  RPN_COMM_const.o  RPN_COMM_comm.o  RPN_COMM_options.o   \
	RPN_COMM_mod.o  RPN_COMM_xch_halosl.o  RPN_COMM_bcst_world.o  RPN_COMM_group.o  RPN_COMM_transpose_utils.o   \
	RPN_COMM_spread.o  RPN_COMM_ezwin.o  RPN_COMM_topo_xy.o  RPN_COMM_mype.o  RPN_COMM_qadl.o   \
	RPN_COMM_Xpose2.o  RPN_COMM_transpose.o  RPN_COMM_xch_halox.o  RPN_COMM_window.o  RPN_COMM_optn.o   \
	RPN_COMM_carac.o  RPN_COMM_grid_redist.o  RPN_COMM_finalize.o  RPN_COMM_haloflip.o  RPN_COMM_swapns.o   \
	RPN_COMM_is_null.o  RPN_COMM_progatate_boundary.o  RPN_COMM_petopo.o  RPN_COMM_Xpose248.o  RPN_COMM_xch_haloew.o   \
	RPN_COMM_io_dist_coll.o  RPN_COMM_init.o  RPN_COMM_dist.o  
INVDEP_LIST_RPN_COMM_types_int.inc= \
	RPN_COMM_io_pes.o  RPN_COMM_gather_e1y3.o  RPN_COMM_split.o  RPN_COMM_grid_2d.o  RPN_COMM_adj_halo.o   \
	RPN_COMM_datyp.o  RPN_COMM_colors.o  RPN_COMM_Xpose148.o  RPN_COMM_bloc.o  RPN_COMM_oper.o   \
	RPN_COMM_topo.o  RPN_COMM_unit.o  RPN_COMM_xch_halon.o  RPN_COMM_C_interface.o  RPN_COMM_xch_halo.o   \
	RPN_COMM_coll.o  RPN_COMM_Xpose1.o  RPN_COMM_diag.o  RPN_COMM_xch_halons.o  RPN_COMM_defo.o   \
	RPN_COMM_globalsum.o  RPN_COMM_move.o  RPN_COMM_adj_halo48.o  RPN_COMM_transpose48.o  RPN_COMM_bloctopo.o   \
	RPN_COMM_const.o  RPN_COMM_comm.o  RPN_COMM_options.o  RPN_COMM_mod.o  RPN_COMM_xch_halosl.o   \
	RPN_COMM_bcst_world.o  RPN_COMM_group.o  RPN_COMM_transpose_utils.o  RPN_COMM_spread.o  RPN_COMM_topo_xy.o   \
	RPN_COMM_mype.o  RPN_COMM_qadl.o  RPN_COMM_Xpose2.o  RPN_COMM_transpose.o  RPN_COMM_xch_halox.o   \
	RPN_COMM_window.o  RPN_COMM_optn.o  RPN_COMM_carac.o  RPN_COMM_grid_redist.o  RPN_COMM_finalize.o   \
	RPN_COMM_haloflip.o  RPN_COMM_swapns.o  RPN_COMM_progatate_boundary.o  RPN_COMM_petopo.o  RPN_COMM_Xpose248.o   \
	RPN_COMM_xch_haloew.o  RPN_COMM_io_dist_coll.o  RPN_COMM_init.o  RPN_COMM_dist.o  
INVDEP_LIST_RPN_COMM_version.inc= \
	RPN_COMM_version.o  
ifneq (,$(DEBUGMAKE))
$(info ## ==== /tmp/viv001/31061/dependencies+.mk [END] =========================================)
endif
