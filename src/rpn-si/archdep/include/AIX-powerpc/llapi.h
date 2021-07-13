#ifndef _llapi_h__
#define _llapi_h__
/*========================================================================*/
/*                                                                        */
/* Module: llapi.a                                                        */
/*                                                                        */
/*========================================================================*/
/* Do not remove the following line! Required for SCCS. */
/* static char sccsid[] = "@(#) src/ll/h/llapi.h, includes, comm_rsnep, rsneps006a 1.47.2.6 5/11/12 16:17:55"; */

#include <stdio.h>
#include <netdb.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <sys/resource.h>

#if defined(__linux__)
typedef int  crid_t;
#endif /* _linux__ */

#if defined(__64BIT__) && !defined(__linux__)
typedef signed long tid_t;   /* thread ID */
#else
typedef signed int tid_t;
#endif /* __64BIT__ */

#ifdef	__cplusplus
extern "C" {
#endif

	typedef struct {
		void     *rst_buffer;     /* Checkpoint-restart interface buffer */
		int       rst_len;        /* Length of interface buffer */
		int       rsvd[8];
	} rst_state_t;

	/* Max data length for rst_state_t.rst_len */
#define MAX_RESTART_DATA	(128*1024)

	typedef struct {
		int      *chk_fdp;         /* Pointer to list of fds to be inherited*/
		int       chk_nfd;         /* No of fds to be inherited */
		pid_t     chk_pid;	   /* Process Id of process     */
		int       rsvd[7];
	} chk_state_t;

	struct chkpt_handler {
		void    (*ch_handler)(int);	/* Ptr to the chkpt handler fn */
		void     *ch_data;		/* Addr of checkpnt_pending variable */
		sigset_t  ch_mask;		/* Signals to be blocked  */
		uint      ch_flags;		/* Style of handler. */
		void	 *ch_attr;		/* Attributes for dedicated handler */
		void     (*ch_advh)(int);	/* Reserved          */
		int	  rsvd[2];		/* Reserved	     */
	};

	/* chkpt_handler.ch_flags = Checkpoint handler style flags */
#define CHK_HANDLER_BASIC      	0x00000001
#define CHK_HANDLER_ADVANCED    0x00000002
#define CHK_HANDLER_POLLING     0x00000004
#define CHK_HANDLER_MASK	0x0000000f

#define CHK_HANDLER_THREAD      0x00000010

#define CHK_MAX_ERRORSZ		2048	/* Maximum length of error data */


	/* User error structure */
	typedef struct {
		char      *error_data;    /* Null-terminated user error data */
		int       Py_error;       /* From errno.h */
		int       Sy_error;       /* Secondary error */
		int       Xtnd_error;     /* Extended error number,  if any */
		int       epid;           /* Pid of process that caused error */
		int       error_len;      /* Length of user error data */
		int       rsvd[4];
	} cr_error_t;

	int checkpnt(char *path, id_t id, unsigned int cflag,
	             chk_state_t *cstate, char *epath);
	int checkpnt_commit(unsigned int ccflag, chk_state_t *cstate,
	                    rst_state_t *rst_common, rst_state_t *rst_primary);
	int checkpnt_register(struct chkpt_handler *newh, struct chkpt_handler *old);
	int checkpnt_fail(cr_error_t *errbuf);
	int checkpnt_block(unsigned int option);
	int checkpnt_deliver(void);
	int checkpnt_wait(id_t id, unsigned int flags, int *status);
	int restart(char *path, id_t id, unsigned int rflag,
	            rst_state_t *rstate, char *epath);
	int restart_wait(unsigned int rwflag, rst_state_t *rstate, chk_state_t *cstate);
	int restart_data(rst_state_t *rstate, unsigned int which, int *len);
	crid_t mkcrid(unsigned int token);
	crid_t getcrid(pid_t pid);
	int setcrid(crid_t crid, pid_t pid, int flags, unsigned int token);
	pid_t getrpid(pid_t vpid, crid_t rcrid);
	tid_t getrtid(void);

	/*
	 * checkpnt() - values for cflag
	 *
	 * CHKPNT_AND_TERMINATE and CHKPNT_AND_STOP should not both be specified.
	 * CHKPNT_IGNORE_SHMEM and CHKPNT_ALWAYS_SHMEM should not both be specified.
	 *
	 * NOTE : these flags should NOT be equal to any of checkpnt_commit flags
	 */
#define CHKPNT_CRID		0x00000001
#define CHKPNT_NODELAY		0x00000002
#define CHKPNT_FORCE		0x00000004
#define CHKPNT_AND_TERMINATE	0x00000010
#define CHKPNT_AND_STOP		0x00000020
#define CHKPNT_AND_STOPTRC      0x00000040
#define CHKPNT_IGNORE_SHMEM	0x00000100
#define CHKPNT_ALWAYS_SHMEM	0x00000200
#define CHKPNT_KEEP_SHMEM	0x00001000
#define CHKPNT_KEEP_SEMAPHORES	0x00002000
#define CHKPNT_KEEP_MSGQS	0x00004000

	/*
	 * checkpnt_commit() - values for ccflag
	 *
	 * NOTE : these flags should NOT be equal to any of checkpnt cflag values
	 */
#define CHKPNT_IGNORE_SOCKETS	0x00010000
#define CHKPNT_IGNORE_TIMERS	0x00020000
#define CHKPNT_IGNORE_PIPEDATA	0x00040000
#define CHKPNT_AND_CONTINUE	0x00100000

	/* checkpnt_block() - values for option */
#define CHK_BLOCK              	1
#define CHK_UNBLOCK             2

	/* restart() - values for rflag */
#define RESTART_OVER_CRID       0x00000001
#define RESTART_IMMED           0x00000002
#define RESTART_NONBLOCK	0x00000004
#define RESTART_AND_STOP	0x00000008
#define RESTART_AND_STOPTRC	0x00000010
#define RESTART_AND_STOPTRC_ALL	0x00000020
#define RESTART_IGNORE_BADSC  0x00000040

	/* source of extended filesystem credentials -- mutually exclusive */
#define	RESTART_CALLER_PAG	0x00000100
#define	RESTART_WAITER_PAG	0x00000200

	/* propogation rules for extended filesystem credentials */
#define	RESTART_PAG_PRIMARY	0x00000400
#define	RESTART_PAG_NONPRIMARY	0x00000800
#define	RESTART_PAG_ALL		(RESTART_PAG_PRIMARY|RESTART_PAG_NONPRIMARY)

	/*
	 * restart_wait() - values for rwflag
	 *
	 * NOTE : these flags should NOT be equal to any of restart rflag values.
	 */
#define RESTART_INHERIT_PGL     0x00010000

	/* restart_wait() - values for source of data */
#define RST_WAITER_DATA         1
#define RST_CALLER_DATA         2

	/* Return values from checkpnt() and checkpnt_commit() */
#define CHECKPOINT_OK           0
#define CHECKPOINT_FAILED      -1
#define RESTART_OK              1

	/* setcrid() flags */
#define SETCRID_FLAGS_NONE 	0
#define SETCRID_FLAGS_PRIMARY	1

#ifdef	__cplusplus
}
#endif

#ifdef   __cplusplus
extern "C" {
#endif

	struct cpr_header {
		int      flags;      /* Checkpoint status flags */
		int      version; /* Version of checkpoint file */
		crid_t      crid;    /* CRID or 0 if single-process chkpt */
		uid_t           owner;     /* Owner of c/r group      */
		int      n_primary;  /* No. of primaries     */
		int             n_secondary;  /* No of secondary processes  */
		int      n_zombie;   /* No. of zombie processes */
		int      reserved[4];
	};

	/*
	 * cpr_header.flags
	 */

#define CHKPT_FILE_OK      0x00000000       /* Checkpoint completed ok */
#define CHKPT_FILE_TRUNC   0x00000001       /* Middle of checkpoint */
#define CHKPT_FILE_ERROR   0x00000100       /* Error file header */
#define RESTART_FILE_ERROR 0x00000200       /* Error file header */


	extern   void  *cp_open( const char * );
	extern  void   cp_header(void *, struct cpr_header *);
	/* cp_error( void *, cr_error_t *, int ); */  /* vukhac : must return int ... */
	extern int cp_error( void *, cr_error_t *, int );
	extern   int   cp_close( void * );


#ifdef   __cplusplus
}
#endif

/***********************************************************************
 * Support for ll_ckpt API.
 **********************************************************************/

#ifdef	__cplusplus
extern "C" {
#endif

#ifndef CKPT_DATA
#define CKPT_DATA
	typedef enum ckpt_type {CKPT_AND_CONTINUE, CKPT_AND_TERMINATE, CKPT_AND_HOLD} CkptType_t;
	typedef enum wait_option {CKPT_NO_WAIT, CKPT_WAIT} WaitOption_t;
#endif

	enum CkptStart {CKPT_YES, CKPT_NO, CKPT_FAIL};

	/* Structure for invoking checkpoint on a specific job step                 */
	/* This structure is also used by ll_init_ckpt to return error information  */
	/*     When used with ll_init_ckpt, the version should be filled in by the  */
	/*     caller, an address to the cp_error_data structure should be passed,  */
	/*     error data information will be filled in when rc from ll_init_ckpt   */
	/*     is -7, all other values should be set to NULL;			    */
	typedef struct LL_ckpt_info {
		int             version;	/* version of the api compiled with */
		char            *step_id;	/* id of step being checkpointed    */
		enum ckpt_type  ckptType;	/* action after success     */
		enum wait_option waitType;	/* identify if blocking enabled */
		int             abort_sig;	/* signal to abort ckpt, default is SIGINT*/
		cr_error_t 	*cp_error_data;  /* structure containing error info from ckpt */
		int		ckpt_rc; 	/* return code from checkpnt()	    */
		int             soft_limit;	/* ckpt soft time limit, in seconds */
		int             hard_limit;	/* ckpt hard time limit, in seconds */
	} LL_ckpt_info;

	/***********************************************************************
	 * Support for ll_(un)set_ckpt_callbacks APIs.
	 **********************************************************************/
	typedef struct {
		void (*checkpoint_callback) (void);
		void (*restart_callback) (void);
		void (*resume_callback) (void);
	} callbacks_t;

#ifdef	__cplusplus
}
#endif

#if defined(__linux__)

#if !defined(LL_RUSAGE64)
#define LL_RUSAGE64
struct   rusage64 {
	struct timeval ru_utime;   /* user time used */
	struct timeval ru_stime;   /* system time used */
	long long   ru_maxrss;
	long long   ru_ixrss;   /* integral shared memory size */
	long long   ru_idrss;   /* integral unshared data */
	long long   ru_isrss;   /* integral unshared stack */
	long long   ru_minflt;  /* page reclaims */
	long long   ru_majflt;  /* page faults */
	long long   ru_nswap;   /* swaps */
	long long   ru_inblock; /* block input operations */
	long long   ru_oublock; /* block output operations */
	long long   ru_msgsnd;  /* messages sent */
	long long   ru_msgrcv;  /* messages received */
	long long   ru_nsignals;   /* signals received */
	long long   ru_nvcsw;   /* voluntary context switches */
	long long   ru_nivcsw;  /* involuntary */
};
#endif /* !LL_RUSAGE64 */

#endif  /* __linux__ */

#if !defined(TRUE)
#define TRUE 1
#endif /* TRUE */

#if !defined(FALSE)
#define FALSE 0
#endif /* FALSE */

#ifndef MAXLEN_HOST
#define MAXLEN_HOST 256
#endif /* MAXLEN_HOST */

#define LL_API_VERSION 410	/* use to keep track of the version */
/* of the api code is compiled with */

typedef void LL_element;

enum LL_Daemon { LL_STARTD, LL_SCHEDD, LL_CM, LL_MASTER,
                 LL_STARTER, LL_HISTORY_FILE, LL_RESOURCE_MANAGER
               };
enum QueryType    { JOBS, MACHINES, PERF,
                    CLUSTERS, WLMSTAT, UNUSED_MATRIX,
                    CLASSES, RESERVATIONS, MCLUSTERS,
                    BLUE_GENE, FAIRSHARE, REGIONS,
                    REGISTERED_HOST_NAMES, MACHINE_GROUP,
                    JOBQ_SUMMARY, ACROSS_SUPERNODE, SUPERNODE
                  };
enum QueryFlags { QUERY_ALL = (1<<0), QUERY_JOBID = (1<<1),
                  QUERY_STEPID = (1<<2), QUERY_USER = (1<<3),
                  QUERY_GROUP = (1<<4), QUERY_CLASS = (1<<5),
                  QUERY_HOST = (1<<6), QUERY_PERF = (1<<7),
                  QUERY_STARTDATE = (1<<8), QUERY_ENDDATE = (1<<9),
                  QUERY_PROCID = (1<<10),
                  QUERY_RESERVATION_ID = (1<<11),
                  QUERY_LOCAL = (1<<12),
                  QUERY_BG_JOB = (1<<13),
                  QUERY_BG_BASE_PARTITION = (1<<14),
                  QUERY_BG_PARTITION = (1<<15),
                  QUERY_TOP_DOG = (1<<16),
                  QUERY_CLUSTER = (1<<17),
                  QUERY_MEDIUM = (1<<18),
                  QUERY_TOPOLOGY = (1<<19),
                  QUERY_STEPID_WITH_VERIFY = (1<<23)
                };

enum DataFilter { ALL_DATA, STATUS_LINE, Q_LINE };

enum JobType { SET_BATCH, SET_INTERACTIVE };

enum JobStepType { BATCH_JOB, INTERACTIVE_JOB };

enum StepParallelMode { SERIAL_TYPE, PARALLEL_TYPE, BLUE_GENE_TYPE, MPICH_TYPE };

enum Usage { SHARED, NOT_SHARED, SLICE_NOT_SHARED };

enum CommLevel { LOW, AVERAGE, HIGH, COMMLVL_UNSPECIFIED };

enum EventType { ERROR_EVENT=-1, STATUS_EVENT, TIMER_EVENT };

enum HoldType  { NO_HOLD, HOLDTYPE_USER, HOLDTYPE_SYSTEM, HOLDTYPE_USERSYS };

enum StepState { STATE_IDLE, STATE_PENDING, STATE_STARTING, STATE_RUNNING,
                 STATE_COMPLETE_PENDING, STATE_REJECT_PENDING, STATE_REMOVE_PENDING,
                 STATE_VACATE_PENDING, STATE_COMPLETED, STATE_REJECTED, STATE_REMOVED,
                 STATE_VACATED, STATE_CANCELED, STATE_NOTRUN, STATE_TERMINATED,
                 STATE_UNEXPANDED, STATE_SUBMISSION_ERR, STATE_HOLD, STATE_DEFERRED,
                 STATE_NOTQUEUED, STATE_PREEMPTED, STATE_PREEMPT_PENDING,
                 STATE_RESUME_PENDING
               };

enum SessionType { BATCH_SESSION, INTERACTIVE_SESSION,
                   INTERACTIVE_HOSTLIST_SESSION
                 };

enum SpawnFlags { MARK_ALL_TASKS_RUNNING = (1<<0) };

enum QueueTypes { QUEUE_SYS_PREEMPTED = (1<<0),
                  QUEUE_GLOBAL_WAIT   = (1<<1),
                  QUEUE_FASTPATH      = (1<<2)
                };


/* Ranges are hardcoded for each object's enumerations. This will enable
   future updates to be grouped with the object and not change
   compatibility. Please add any new enums to the end of the range for each
   object */
enum LLAPI_Specification {
	LL_JobManagementInteractiveClass=0, /* char **: LL interactive class */
	LL_JobManagementListenSocket,  /* int * : socket Schedd sends info on */
	LL_JobManagementAccountNo,     /* char ** : returns LOADL_ACCOUNT_NO */
	LL_JobManagementSessionType,   /* int * : session type info */
	LL_JobManagementPrinterFILE,   /* LL_element * : Default printer */
	LL_JobManagementRestorePrinter,/* LL_element * : restore previous printer */

	/* Job object data */
	LL_JobGetFirstStep=200, /* LL_element * (Step) : first step       */
	LL_JobGetNextStep,      /* LL_element * (Step) : next step        */
	LL_JobCredential,       /* LL_element * (Credential): credentials */
	LL_JobName,             /* char ** : job name                     */
	LL_JobStepCount,        /* int *   : number of steps in job       */
	LL_JobStepType,         /* int *   : INTERACTIVE_JOB or BATCH_JOB */
	LL_JobSubmitHost,       /* char ** : host job was submitted from  */
	LL_JobSubmitTime,       /* time_t  : time job was queued          */
	LL_JobVersionNum,       /* int *   : job version number           */
	LL_JobSchedd,           /* char ** : Schedd managing the job      */
	LL_JobJobQueueKey,      /* int *   : key used to write job to spool */
	LL_JobIsRemote,         /* int *   : is job remote - 1 true       */
	/* The following data is available for remote jobs. Local jobs - NULL */
	LL_JobSchedulingCluster,/* char ** : cluster scheduling remote job*/
	LL_JobSubmittingCluster,/* char ** : original cluster job was submitted on*/
	LL_JobSubmittingUser,   /* char ** : user name job was submitted under */
	LL_JobSendingCluster,   /* char ** : cluster which sent the remote job */
	LL_JobRequestedCluster, /* char ** : cluster_list JCF keyword value */
	LL_JobLocalOutboundSchedds, /* char *** : list of local outbound schedds, last schedd in list is current outbound */
	LL_JobScheddHistory,        /* char *** : list of scheduling schedds, last schedd in list is current */
	/*  End of remote job values */

	LL_JobGetFirstClusterInputFile,  /* LL_element *(ClusterFile): first input ClusterFile object */
	LL_JobGetNextClusterInputFile,   /* LL_element *(ClusterFile): next input ClusterFile object */
	LL_JobGetFirstClusterOutputFile, /* LL_element *(ClusterFile): first output ClusterFile object */
	LL_JobGetNextClusterOutputFile,  /* LL_element *(ClusterFile): next output ClusterFile file object */
	LL_JobUsersJCF, 		     /* char **: users JCF statements */

	/* Step object data */
	LL_StepNodeCount=400,          /* int * : number of nodes in step	*/
	LL_StepGetFirstNode,           /* LL_element * (Node) : first node	*/
	LL_StepGetNextNode,            /* LL_element * (Node) : next node	*/
	LL_StepMachineCount,           /* int * : number of machines assigned 	*/
	LL_StepGetFirstMachine,        /* LL_element * (LlMachine): first machine*/
	LL_StepGetNextMachine,         /* LL_element * (LlMachine): next machine */
	LL_StepGetFirstSwitchTable,    /* LL_element * (SwitchTable): 	*/
	LL_StepGetNextSwitchTable,     /* LL_element * (SwitchTable): 	*/
	LL_StepGetMasterTask,          /* LL_element * (Task): master task on step */
	LL_StepTaskInstanceCount,      /* int * : number of task instances	*/
	LL_StepAccountNumber,          /* char ** : associated account number   */
	LL_StepAdapterUsage,           /* int * :  requested adapter usaqe 	*/
	LL_StepComment,                /* char ** : user defined comment	*/
	LL_StepCompletionCode,         /* int * : exit status 			*/
	LL_StepCompletionDate,         /* time_t * : UNIX time step was completed */
	LL_StepEnvironment,            /* char ** : environment			*/
	LL_StepErrorFile,              /* char ** : name of error file 		*/
	LL_StepExecSize,               /* int * : size of the master executable */
	LL_StepHostName,               /* char * : name of host to be scheduled*/
	LL_StepID,                     /* char ** : scheddhostname.jobid.stepid */
	LL_StepInputFile,              /* char ** : name of input file 		*/
	LL_StepImageSize,              /* int * : size of the virtual image in K */
	LL_StepImmediate,              /* int * : immediate scheduling of job 	*/
	LL_StepIwd,                    /* char ** : initial directory 		*/
	LL_StepJobClass,               /* char ** : defined class for step 	*/
	LL_StepMessages,               /* char ** : string of messages from LL  */
	LL_StepName,                   /* char ** : name assigned to step 	*/
	LL_StepNodeUsage,              /* int * : requested node usage 	*/
	LL_StepOutputFile,             /* char ** : name of output file 	*/
	LL_StepParallelMode,           /* int * : mode of the step  		*/
	LL_StepPriority,               /* int * : User priority for step 	*/
	LL_StepShell,                  /* char ** : shell to be used 		*/
	LL_StepStartDate,              /* time_t * : date step was started 	*/
	LL_StepDispatchTime,           /* time_t * : time step was started 	*/
	LL_StepState,                  /* int * : current state of the step	*/
	LL_StepStartCount,             /* int * : times step has been started 	*/
	LL_StepCpuLimitHard,           /* int * : cpu hard limit		*/
	LL_StepCpuLimitSoft,           /* int * : cpu soft limit		*/
	LL_StepCpuStepLimitHard,       /* int * : cpu hard limit for entire job */
	LL_StepCpuStepLimitSoft,       /* int * : cpu soft limit for entire job */
	LL_StepCoreLimitHard,          /* int * : core file hard limit 	*/
	LL_StepCoreLimitSoft,          /* int * : core file soft limit		*/
	LL_StepDataLimitHard,          /* int * : data hard limit 		*/
	LL_StepDataLimitSoft,          /* int * : data soft limit 		*/
	LL_StepFileLimitHard,          /* int * : file size hard limit		*/
	LL_StepFileLimitSoft,          /* int * : file size soft limit		*/
	LL_StepRssLimitHard,           /* int * : resident set size hard limit */
	LL_StepRssLimitSoft,           /* int * : resident set size soft limit */
	LL_StepStackLimitHard,         /* int * : stack size hard limit	*/
	LL_StepStackLimitSoft,         /* int * : stack size soft limit	*/
	LL_StepWallClockLimitHard,     /* int * : wall clock hard limit 	*/
	LL_StepWallClockLimitSoft,     /* int * : wall clock soft limit 	*/
	LL_StepHostList,               /* char *** : hosts in host.list file */
	LL_StepHoldType,               /* int * : hold state of step           */
	LL_StepLoadLevelerGroup,       /* char ** : LL group specified by step  */
	LL_StepGetFirstAdapterReq,     /* LL_element*(AdapterReq):first AdapterReq*/
	LL_StepGetNextAdapterReq,      /* LL_element*(AdapterReq):next AdapterReq */
	LL_StepRestart,                /* int * : restartable job? */
	LL_StepBlocking,               /* int * : blocking factor requested */
	LL_StepTaskGeometry,           /* char ** : task geometry requested */
	LL_StepTotalTasksRequested,    /* int * : total tasks requested */
	LL_StepTasksPerNodeRequested,  /* int * : tasks per node requested */
	LL_StepTotalNodesRequested,    /* char ** : total nodes requested */
	LL_StepSystemPriority,         /* int * : system priority for step 	*/
	LL_StepClassSystemPriority,    /* int * : class system priority for step 	*/
	LL_StepGroupSystemPriority,    /* int * : group system priority for step 	*/
	LL_StepUserSystemPriority,     /* int * : user system priority for step 	*/
	LL_StepQueueSystemPriority,    /* int * : queue system priority for step 	*/
	Unused_LL_StepExecutionFactor, /* reserved for future use  */
	LL_StepImageSize64,            /* int64_t * : size of the virtual image in K */
	LL_StepCpuLimitHard64,         /* int64_t * : cpu hard limit		*/
	LL_StepCpuLimitSoft64,         /* int64_t * : cpu soft limit		*/
	LL_StepCpuStepLimitHard64,     /* int64_t*: cpu hard limit for entire job */
	LL_StepCpuStepLimitSoft64,     /* int64_t*: cpu soft limit for entire job */
	LL_StepCoreLimitHard64,        /* int64_t * : core file hard limit 	*/
	LL_StepCoreLimitSoft64,        /* int64_t * : core file soft limit	*/
	LL_StepDataLimitHard64,        /* int64_t * : data hard limit 	*/
	LL_StepDataLimitSoft64,        /* int64_t * : data soft limit 	*/
	LL_StepFileLimitHard64,        /* int64_t * : file size hard limit	*/
	LL_StepFileLimitSoft64,        /* int64_t * : file size soft limit	*/
	LL_StepRssLimitHard64,         /* int64_t*: resident set size hard limit */
	LL_StepRssLimitSoft64,         /* int64_t*: resident set size soft limit */
	LL_StepStackLimitHard64,       /* int64_t * : stack size hard limit	*/
	LL_StepStackLimitSoft64,       /* int64_t * : stack size soft limit	*/
	LL_StepWallClockLimitHard64,   /* int64_t * : wall clock hard limit*/
	LL_StepWallClockLimitSoft64,   /* int64_t * : wall clock soft limit*/
	LL_StepStepUserTime64,         /* int64_t * : step user time - microseconds */
	LL_StepStepSystemTime64,       /* int64_t*: step system time - microseconds */
	LL_StepStepMaxrss64,           /* int64_t * : step ru_maxrss */
	LL_StepStepIxrss64,            /* int64_t * : step ru_ixrss */
	LL_StepStepIdrss64,            /* int64_t * : step ru_idrss */
	LL_StepStepIsrss64,            /* int64_t * : step ru_isrss */
	LL_StepStepMinflt64,           /* int64_t * : step ru_minflt */
	LL_StepStepMajflt64,           /* int64_t * : step ru_majflt */
	LL_StepStepNswap64,            /* int64_t * : step ru_nswap */
	LL_StepStepInblock64,          /* int64_t * : step ru_inblock */
	LL_StepStepOublock64,          /* int64_t * : step ru_oublock */
	LL_StepStepMsgsnd64,           /* int64_t * : step ru_msgsnd */
	LL_StepStepMsgrcv64,           /* int64_t * : step ru_msgrcv */
	LL_StepStepNsignals64,         /* int64_t * : step ru_nsignals */
	LL_StepStepNvcsw64,            /* int64_t * : step ru_nvcsw */
	LL_StepStepNivcsw64,           /* int64_t * : step ru_nivcsw */
	LL_StepStarterUserTime64,/* int64_t*: starter user time-microseconds */
	LL_StepStarterSystemTime64,/* int64_t*: starter system time-microseconds */
	LL_StepStarterMaxrss64,        /* int64_t * : starter ru_maxrss */
	LL_StepStarterIxrss64,         /* int64_t * : starter ru_ixrss */
	LL_StepStarterIdrss64,         /* int64_t * : starter ru_idrss */
	LL_StepStarterIsrss64,         /* int64_t * : starter ru_isrss */
	LL_StepStarterMinflt64,        /* int64_t * : starter ru_minflt */
	LL_StepStarterMajflt64,        /* int64_t * : starter ru_majflt */
	LL_StepStarterNswap64,         /* int64_t * : starter ru_nswap */
	LL_StepStarterInblock64,       /* int64_t * : starter ru_inblock */
	LL_StepStarterOublock64,       /* int64_t * : starter ru_oublock */
	LL_StepStarterMsgsnd64,        /* int64_t * : starter ru_msgsnd */
	LL_StepStarterMsgrcv64,        /* int64_t * : starter ru_msgrcv */
	LL_StepStarterNsignals64,      /* int64_t * : starter ru_nsignals */
	LL_StepStarterNvcsw64,         /* int64_t * : starter ru_nvcsw */
	LL_StepStarterNivcsw64,        /* int64_t * : starter ru_nivcsw */
	LL_StepMachUsageCount,         /* int * : count of machine usages */
	LL_StepGetFirstMachUsage,      /* LL_element * (MachineUsage): first machine usage */
	LL_StepGetNextMachUsage,       /* LL_element * (MachineUsage): next machine usage */
	LL_StepCheckpointable,         /* int * : 0-job not checkpointable */
	/*         1-job is checkpointable */
	LL_StepCheckpointing,          /* Boolean * : step is being checkpointed */
	LL_StepCkptAccumTime,          /* int * : accumulated ckpt time */
	LL_StepCkptFailStartTime,      /* time_t * : time of last failed ckpt*/
	LL_StepCkptFile,               /* char ** : name of ckpt file */
	LL_StepCkptGoodElapseTime,     /* time_t * : time taken to complete last good ckpt */
	LL_StepCkptGoodStartTime,      /* time_t * : time of last good ckpt */
	LL_StepCkptTimeHardLimit,      /* int * : ckpt time hard limit */
	LL_StepCkptTimeHardLimit64,    /* int64_t * : ckpt time hard limit */
	LL_StepCkptTimeSoftLimit,      /* int * : ckpt time soft limit */
	LL_StepCkptTimeSoftLimit64,    /* int64_t * : ckpt time soft limit */
	LL_StepCkptRestart,            /* int * : 0|1 job restarted from ckpt*/
	LL_StepCkptRestartSameNodes,   /* int * : 0|1 restart job on same nodes */
	LL_StepWallClockUsed,          /* int * : wallclock already used by the step */
	LL_StepLargePage,              /* char ** : Large Page requirement (= "M", "Y", or "N") */
	LL_StepMaxProtocolInstances,   /* int * : largest number of instances allowed on network stmt */
	LL_StepBulkXfer,		/* int * : 0|1 Step Requests Bulk Transfer */
	LL_StepTotalRcxtBlocks,	/* int * : total number of RCXT blocks application needs */
	LL_StepStartTime,       /* time_t *: time the starter process for the job step started  */
	LL_StepUserRcxtBlocks,	/* int * : number of user RCXT blocks application requests */
	LL_StepRequestedReservationID, /* char ** : The step's requested reservation ID */
	LL_StepReservationID,          /* char ** : The step's reservation ID */
	LL_StepPreemptable,            /* int * : 0|1 step is preemptable */
	LL_StepPreemptWaitList,      /* char *** list of steps to preempted */
	LL_StepRsetName,            /* char ** rset name */
	LL_StepCkptExecuteDirectory,   /* char ** : Executable directory for a checkpointable step */
	LL_StepAcctKey,                /* int64_t * : accounting key for step */
	LL_StepDependency,             /* char ** : Step Dependency */
	LL_StepFavoredJob,             /* int * : whether the step is favored using llfavorjob */
	LL_StepBgJobId,                /* char ** : Blue Gene ID for the step */
	LL_StepBgJobState,             /* int * : Blue Gene state for the step */
	LL_StepBgSizeRequested,        /* int * : size requested for Blue Gene step */
	LL_StepBgSizeAllocated,        /* int * : size allocated for Blue Gene step */
	LL_StepBgShapeRequested,       /* int ** : shape requested for Blue Gene step */
	LL_StepBgShapeAllocated,       /* int ** : shape allocated for Blue Gene step */
	LL_StepBgConnectionRequested,  /* int * : type of wiring requested for Blue Gene step */
	LL_StepBgConnectionAllocated,  /* int * : type of wiring allocated for Blue Gene step */
	LL_StepBgPartitionRequested,   /* char ** : Blue Gene partition requested for step */
	LL_StepBgPartitionAllocated,   /* char ** : Blue Gene partition allocated for step */
	LL_StepBgPartitionState,       /* int * : state of Blue Gene partition allocated for step */
	LL_StepBgErrorText,            /* char ** : error text from Blue Gene system for step */
	LL_StepMcmAffinityOptions,     /* char ** mcm affinity options */
	LL_StepCoschedule,             /* int * : 0|1 coschedule option */
	LL_StepSMTRequired,            /* int * : if SMT is required to be turned on , off or as_is for step */
	LL_StepMetaClusterJobID,       /* int * : MetaCluster job ID LoadLeveler has assigned to job step */
	LL_StepMetaClusterJob,         /* int * : 1|0 = 1 if metacluster_job keyword is set to "yes". */
	LL_StepMasterVirtualIP,        /* char ** : MetaCluster virtual IP associated with POE process */
	LL_StepMasterRealIP,           /* char ** : MetaCluster: real IP of adapter associated with POE process */
	LL_StepMasterNetmask,          /* char ** : MetaCluster: real netmask of adapter associated with POE process */
	LL_StepVipNetmask,             /* char ** : MetaCluster vipserver virtual netmask */
	LL_StepMetaClusterPoeHostname,     /* char ** : MetaCluster:  POE host name */
	LL_StepMetaClusterPoePmdPhysnet,   /* char ** : MetaCluster:  Physnet of POE/PMD adapters */
	LL_StepCkptSubDir,             /* char ** : name of ckpt_subdir file */
	LL_StepAsLimitHard64,          /* int64_t * : as hard limit         */
	LL_StepAsLimitSoft64,          /* int64_t * : as soft limit         */
	LL_StepNprocLimitHard64,       /* int64_t * : nproc hard limit      */
	LL_StepNprocLimitSoft64,       /* int64_t * : nproc soft limit      */
	LL_StepMemlockLimitHard64,     /* int64_t * : memlock hard limit    */
	LL_StepMemlockLimitSoft64,     /* int64_t * : memlock soft limit    */
	LL_StepLocksLimitHard64,       /* int64_t * : locks hard limit      */
	LL_StepLocksLimitSoft64,       /* int64_t * : locks soft limit      */
	LL_StepNofileLimitHard64,      /* int64_t * : nofile hard limit     */
	LL_StepNofileLimitSoft64,      /* int64_t * : nofile soft limit     */
	LL_StepTaskAffinity,           /* char** : task affinity request    */
	LL_StepCpusPerCore,            /* int * : cpus per core             */
	LL_StepIsTopDog,               /* Boolean * : step is a top dog     */
	LL_StepConsideredAt,           /* time_t *  : time the top dog is considered by CM */
	LL_StepEstimatedStartTime,     /* time_t *  : estimated time for the top dog */
	LL_StepUserHoldTime,           /* time_t *  : time of in UserHold state	*/
	LL_StepQueueId,                /* int * : id of queue which the step is hold */
	LL_StepQueueIndex,             /* int * : index in the queue */
	LL_StepClusterOption,          /* char ** : cluster option */
	LL_StepScaleAcrossClusterCount,    /* int * : number of scale across clusters where this step is being scheduled or is dispatched to */
	LL_StepGetFirstScaleAcrossCluster, /* LL_element * (LlMCluster) : first cluster */
	LL_StepGetNextScaleAcrossCluster,  /* LL_element * (LlMCluster) : next cluster */

	LL_StepBgPartitionType,		/* int * : Blue Gene partition type */
	LL_StepJobKey,           	/* int * : job key */
	/*topology scheduling*/
	LL_StepTopologyName,            /*char * : topology group name */
	LL_StepTopologyRequire,         /*char * : Require topology type */
	LL_StepFlexibleReservationID,   /* char ** : The step's flexible reservation ID */
	LL_StepEligibilityTime,         /* time_t * : time step became eligible for dispatch */
	LL_StepGetFirstStepResourceRequirement,	/* LL_element * (ResourceReq) */
	LL_StepGetNextStepResourceRequirement,  /* LL_element * (ResourceReq) */

	/* Machine object data */
	LL_MachineAdapterList=800,   /* char *** : adapters defined	*/
	LL_MachineArchitecture,      /* char ** : machine architecture	*/
	LL_MachineAvailableClassList,/* char *** : classes defined */
	LL_MachineCPUs,              /* int * : number of cpus		*/
	LL_MachineDisk,              /* int * : avail space (KB) in execute dir */
	LL_MachineFeatureList,       /* char *** : features defined	*/
	LL_MachineConfiguredClassList,/* char ***: initiators defined */
	LL_MachineKbddIdle,          /* int * : seconds kbdd is idle  	*/
	LL_MachineLoadAverage,       /* double * : load average 		*/
	LL_MachineMachineMode,       /* char ** : configured machine mode 	*/
	LL_MachineMaxTasks,          /* int * : max number of tasks allowed 	*/
	LL_MachineName,              /* char * : official hostname 		*/
	LL_MachineOperatingSystem,   /* char ** : machine operating system */
	LL_MachinePoolList,          /* int ** : list of configured pools */
	LL_MachineRealMemory,        /* int * : physical memory 		*/
	LL_MachineScheddRunningJobs, /* int * : running jobs assigned schedd */
	LL_MachineScheddState,       /* int * : schedd state			*/
	LL_MachineScheddTotalJobs,   /* int * : total jobs assigned schedd */
	LL_MachineSpeed,             /* double * : speed associated with machine */
	LL_MachineStartdState,       /* char ** : startd state		*/
	LL_MachineStartdRunningJobs, /* int * : running jobs assigned startd */
	LL_MachineStepList,          /* char *** : stepids scheduled to run*/
	LL_MachineTimeStamp,         /* time_t *: time when machine data received */
	LL_MachineVirtualMemory,     /* int * : available swap space in kilobytes */
	LL_MachinePoolListSize,      /* int * : number of configured pools   */
	LL_MachineFreeRealMemory,    /* int * : free real memory in Mbytes */
	LL_MachinePagesScanned,      /* int * : pages scanned/sec by page replacement algorithm */
	LL_MachinePagesFreed,        /* int * : pages freed/sec by page replacement algorithm */
	LL_MachinePagesPagedIn,      /* int * : pages paged in from paging space */
	LL_MachinePagesPagedOut,     /* int * : pages paged out to paging space  */
	LL_MachineGetFirstResource,  /* LL_element * (Resource):first resource */
	LL_MachineGetNextResource,   /* LL_element * (Resource):next  resource */
	LL_MachineGetFirstAdapter,   /* LL_element * (Adapter): first adapter */
	LL_MachineGetNextAdapter,    /* LL_element * (Adapter): next adapter  */
	LL_MachineDrainingClassList, /* char *** : draining class list */
	LL_MachineDrainClassList,    /* char *** : drain class list */
	LL_MachineStartExpr,         /* char ** : START expression */
	LL_MachineSuspendExpr,       /* char ** : SUSPEND expression */
	LL_MachineContinueExpr,      /* char ** : CONTINUE expression */
	LL_MachineVacateExpr,        /* char ** : VACATE expression */
	LL_MachineKillExpr,          /* char ** : KILL expression */
	LL_MachineDisk64,            /* int64_t * : avail space (KB) in execute dir */
	LL_MachineRealMemory64,      /* int64_t * : physical memory 		*/
	LL_MachineVirtualMemory64,   /* int64_t * : available swap space in Kb */
	LL_MachineFreeRealMemory64,  /* int64_t * : free real memory in Mbytes */
	LL_MachinePagesScanned64,    /* int64_t * : pages scanned/sec by page replacement algorithm */
	LL_MachinePagesFreed64,      /* int64_t * : pages freed/sec by page replacement algorithm */
	LL_MachinePagesPagedIn64,    /* int64_t*:pages paged in from paging space*/
	LL_MachinePagesPagedOut64,   /* int64_t*:pages paged out to paging space*/
	LL_MachineLargePageSize64,   /* int64_t * : Size of Large Page memory block (bytes) */
	LL_MachineLargePageCount64,  /* int64_t * : Total number of Large Page memory blocks */
	LL_MachineLargePageFree64,   /* int64_t * : Number of Large Page memory blocks on freelist */
	LL_MachineReservationPermitted, /* int * : boolean, can this machine be reserved */
	LL_MachineReservationList,   /* char *** : list of IDs of reservations using this machine */
	LL_MachinePrestartedStarters,/* int * : Prestarted Starters to be started */
	LL_MachineCPUList,           /* int ** : list of cpus */
	LL_MachineUsedCPUs,          /* int *  : used cpus on this machine */
	LL_MachineUsedCPUList,       /* int ** : list of used cpus on this machine */
	LL_MachineGetFirstMCM,       /* LL_element * (MCM): first mcm */
	LL_MachineGetNextMCM,        /* LL_element * (MCM): next mcm  */
	LL_MachineConfigTimeStamp,   /* int * : Time of last reconfig */
	LL_MachineRSetSupport,       /* int * : The RSet Support configured */
	LL_MachineSMTState,          /* int * : The SMT state of a node */
	LL_MachineMaxDstgStarters,   /* int * : max number of data staging tasks allowed */
	LL_MachineMachineGroupName,  /* char ** : configured machine group name */
	LL_MachineScheddRunsHere,    /* int * : Whether schedd daemon is configured on this machine */
	LL_MachineStartdRunsHere,    /* int * : Whether startd daemon is configured on this machine */
	LL_MachineScheddStepPending,     /* int * : Number of steps of this schedd machine which has state pending */
	LL_MachineScheddStepRemovePending,     /* int * : Number of steps of this schedd machine which has state remove pending */
	LL_MachineScheddStepStarting,    /* int * : Number of steps of this schedd machine which has state Starting */
	LL_MachineScheddStepRunning,     /* int * : Number of steps of this schedd machine which has state running */

	/*topology scheduling*/
	LL_MachineSuperNode,             /* int * : Supernode index */
	LL_MachineSuperSegment,          /* int * : Supernode segment index */
	LL_MachineShuffleExchangeSegment, /* int * : Shuffle exchange segment index*/
	LL_MachineSector,                /* int * : Sector index */

	/* Node object data */
	LL_NodeTaskCount=1000, /* int * : number of task instances	*/
	LL_NodeGetFirstTask,   /* LL_element * (Task) : first task	*/
	LL_NodeGetNextTask,    /* LL_element * (Task) : next task 	*/
	LL_NodeMaxInstances,   /* int * : maximum # requested nodes    */
	LL_NodeMinInstances,   /* int * : minimum # requested nodes    */
	LL_NodeRequirements,   /* char ** : defined requirements 	*/
	LL_NodeInitiatorCount, /* int * : initiator count      	*/
	LL_NodeGetFirstResourceRequirement, /* LL_element * (ResourceReq) */
	LL_NodeGetNextResourceRequirement,  /* LL_element * (ResourceReq) */

	LL_SwitchTableJobKey=1200,/* int * : job key			*/

	/* Task object data */
	LL_TaskTaskInstanceCount=1400,      /* int * : number of task instances*/
	LL_TaskGetFirstTaskInstance,        /* LL_element * (TaskInstance)	*/
	LL_TaskGetNextTaskInstance,         /* LL_element * (TaskInstance)	*/
	LL_TaskExecutable,                  /* char ** : executable 		*/
	LL_TaskExecutableArguments,         /* char ** : executable arguments	*/
	LL_TaskIsMaster,                    /* int * : boolean, is this the master task */
	LL_TaskGetFirstResourceRequirement, /* LL_element * (ResourceReq) */
	LL_TaskGetNextResourceRequirement,  /* LL_element * (ResourceReq) */

	/* Task instance object data */
	LL_TaskInstanceAdapterCount=1600,    /* int * : number of adapters	*/
	LL_TaskInstanceGetFirstAdapter,      /* LL_element * (Adapter)	*/
	LL_TaskInstanceGetNextAdapter,       /* LL_element * (Adapter)	*/
	LL_TaskInstanceGetFirstAdapterUsage, /* LL_element * (AdapterUsage) */
	LL_TaskInstanceGetNextAdapterUsage,  /* LL_element * (AdapterUsage) */
	LL_TaskInstanceMachineName,          /*  char ** : machine assigned   */
	LL_TaskInstanceTaskID,               /*  int * : task id	 	*/
	LL_TaskInstanceMachineAddress,       /*  char * : machine IP address */
	LL_TaskInstanceMachine,              /*  LL_element * (LlMachine): machine Object */
	LL_TaskInstanceCpuList,              /* cpus used by task */
	LL_TaskInstanceMachineVirtualIP,     /* char ** : MetaCluster virtual IP associated PMD process */

	/* Adapter object data */
	LL_AdapterInterfaceAddress=1800,/* char ** : interface IP address*/
	LL_AdapterMode,              /* char ** : Use LL_AdapterUsageMode		*/
	LL_AdapterName,              /* char ** : Adapter name		*/
	Unused_LL_AdapterUsageWindow,       /* int * : window assigned to task 	*/
	Unused_LL_AdapterUsageProtocol,     /* char ** : Protocol used by task	*/
	LL_AdapterCommInterface=1806,     /* int * : communication interface	*/
	Unused_LL_AdapterUsageMode,         /* char ** : Used for css IP or US	*/
	LL_AdapterTotalWindowCount=1811,  /* int * : # of windows on adapter  */
	LL_AdapterAvailWindowCount,  /* int * : # of windows not in use  */
	Unused_LL_AdapterUsageAddress,      /* char ** : IP Address to use adapter  */
	Unused_LL_AdapterUsageCommunicationInterface, /* int * : comm interface  */
	Unused_LL_AdapterUsageDevice,       /* char ** : Name of adapter device being used */
	Unused_LL_AdapterUsageInstanceNumber,/* int * :Unique ID for multiple instances */
	Unused_LL_AdapterUsageNetworkId,	/* int *: Network ID of adapter being used */
	Unused_LL_AdapterWindowList,        /* int ** : Array of window numbers */
	Unused_LL_AdapterUsageWindowMemory64,/* uint64_t * : bytes used by window    */
	LL_AdapterMinWindowSize64,   /* uint64_t * : min allocatable window memory */
	LL_AdapterMaxWindowSize64,   /* uint64_t * : max allocatable window memory */
	LL_AdapterMemory64,          /* uint64_t * : Total adapter memory		*/
	Unused_LL_AdapterUsageTag,		 /* char** : Tag that identifies switch table for usage */
	LL_AdapterMCMId,						 /* int * : mcm this adpter connected to */
	Unused_LL_AdapterUsageRcxtBlocks,   /* int * : number of rCxt blocks used by window */
	LL_AdapterRcxtBlocks,         /* int * : number of rCxt blocks available on adapter */
	Unused_LL_AdapterUsageExclusive, /* int * : 1 if the usage is exclusive, 0 if it is not */
	Unused_LL_AdapterUsagePortNumber,    /* int *: Port number of adapter being used */
	Unused_LL_AdapterUsageLmc,           /* int *: lmc of adapter being used */
	LL_AdapterPortNumber,         /* int *: Port number of adapter */
	LL_AdapterLmc,                /* int *: lmc of adapter */
	Unused_LL_AdapterUsageNetworkId64,	  /* uint64_t *: 64 bit Network ID of adapter being used */
	Unused_LL_AdapterUsageDeviceDriver,  /* char ** : Name of adapter device driver being used */
	Unused_LL_AdapterUsageDeviceType,    /* int *: Device type of adapter being used */
	LL_AdapterInterfaceNetmask,   /* char ** : Netmask of adapter */
	Unused_LL_AdapterUsageVirtualIP,     /* char ** : MetaCluster virtual IP associated with an application task */
	Unused_LL_AdapterUsageNetmask,       /* char ** : Netmask of adapter in AdapterUsage object */

	/* Credential object data */
	LL_CredentialGid=2000,  /* int * : unix group id of submitter 	*/
	LL_CredentialGroupName, /* char ** : User group for job		*/
	LL_CredentialUid,       /* int * : unix userid of submitter 	*/
	LL_CredentialUserName,  /* char ** : login of person submitting job */

	LL_StartdPerfJobsRunning=2200,  /* All of the StartdPerf are of */
	LL_StartdPerfJobsPending,       /* int * data type		*/
	LL_StartdPerfJobsSuspended,
	LL_StartdPerfCurrentJobs,
	LL_StartdPerfTotalJobsReceived,
	LL_StartdPerfTotalJobsCompleted,
	LL_StartdPerfTotalJobsRemoved,
	LL_StartdPerfTotalJobsVacated,
	LL_StartdPerfTotalJobsRejected,
	LL_StartdPerfTotalJobsSuspended,
	LL_StartdPerfTotalConnections,
	LL_StartdPerfFailedConnections,
	LL_StartdPerfTotalOutTransactions,
	LL_StartdPerfFailedOutTransactions,
	LL_StartdPerfTotalInTransactions,
	LL_StartdPerfFailedInTransactions,

	LL_ScheddPerfJobsIdle=2400,	/* All of the ScheddPerf are of */
	LL_ScheddPerfJobsPending,	/* int * data type 		*/
	LL_ScheddPerfJobsStarting,
	LL_ScheddPerfJobsRunning,
	LL_ScheddPerfCurrentJobs,
	LL_ScheddPerfTotalJobsSubmitted,
	LL_ScheddPerfTotalJobsCompleted,
	LL_ScheddPerfTotalJobsRemoved,
	LL_ScheddPerfTotalJobsVacated,
	LL_ScheddPerfTotalJobsRejected,
	LL_ScheddPerfTotalConnections,
	LL_ScheddPerfFailedConnections,
	LL_ScheddPerfTotalOutTransactions,
	LL_ScheddPerfFailedOutTransactions,
	LL_ScheddPerfTotalInTransactions,
	LL_ScheddPerfFailedInTransactions,

	LL_VersionCheck=2600,         /* used by POE for release checking */

	/* AdapterReq object data */
	LL_AdapterReqCommLevel=2700,  /* int * : communication level      */
	LL_AdapterReqUsage,           /* int * : requested adapter usage */
	LL_AdapterReqInstances,       /* int * : requested number of instances for protocol */
	LL_AdapterReqRcxtBlks,        /* int * : requested number of user rCxt Blocks for protocol */
	LL_AdapterReqProtocol,        /* char ** : requested adapter protocol */
	LL_AdapterReqMode,            /* char ** : requested adapter mode */
	LL_AdapterReqTypeName,        /* char ** : requested adapter type */
	LL_AdapterReqCollectiveGroups,    /* int * : requested collective group indexes of HFI Adapter for protocol */
	LL_AdapterReqImmSendBuffers,    /* int * : requested immediate send buffers of HFI Adapter for protocol */

	/* Cluster object data */
	LL_ClusterGetFirstResource=2800,  /* LL_element * (Resource): first */
	LL_ClusterGetNextResource,        /* LL_element * (Resource): next  */
	LL_ClusterSchedulingResources,    /* char ***:scheduling resources*/
	LL_ClusterDefinedResources,       /* char *** : resources defined   */
	LL_ClusterSchedulingResourceCount,/* int *: # of scheduling resources*/
	LL_ClusterDefinedResourceCount,   /* int *: number of defined resources*/
	LL_ClusterEnforcedResources,      /* char ***: enforced resources */
	LL_ClusterEnforcedResourceCount,  /* int * : # of enforced resources   */
	LL_ClusterEnforceSubmission,      /* int * : Boolean, are resources required at submission time    */
	LL_ClusterSchedulerType,          /* char **: scheduler type     */
	LL_ClusterPreemptionEnabled,      /* int * : Boolean, is the preemption function enabled */
	LL_ClusterSysPrioThreshold,       /* int *: value of SYSPRIO_THRESHOLD_TO_IGNORE_STEP */
	LL_ClusterMusterEnvironment,      /* int *: 1-Muster Environment Enabled */
	LL_ClusterClusterMetric,          /* char**: CLUSTER_METRIC string */
	LL_ClusterClusterUserMapper,      /* char**: CLUSTER_USER_MAPPER string */
	LL_ClusterClusterRemoteJobFilter, /* char**: CLUSTER_REMOTE_JOB_FILTER string */
	LL_ClusterEnforceMemory,          /* int * : Boolean, is the Absolute Memory Limit enabled */
	LL_ClusterScaleAcrossEnv,	      /* int * : Boolean, is scale-across environment enabled */
	LL_ClusterStartdTotal,		  /* int * : number of machines known to have startd configured */
	LL_ClusterStartdAvailable,	  /* int * : number of machines known to have startd running */
	LL_ClusterScheddTotal,		  /* int * : number of machines known to have schedd configured */
	LL_ClusterScheddAvailable,	  /* int * : number of machines known to have schedd running */
	LL_ClusterJobStepsInQueue,	  /* int * : The total number of job steps queued. */
	LL_ClusterTasksRunning,		  /* int * : The total number of tasks of all job steps running in the cluster. */
	LL_ClusterCentralManager,	  /* char ** : The name of the machine currently running central manager */
	LL_ClusterMachineAbsentList,	  /* char *** : machines which are present in the LoadLeveler configuration but on which no daemons are running. */
	LL_ClusterGetFirstUnavailableMachine,/* LL_element* (Machine) : A pointer to the first machine object which has either an unavailable startd or an unavailable schedd. */
	LL_ClusterGetNextUnavailableMachine, /* LL_element* (Machine) : A pointer to the next machine object which has either an unavailable startd or an unavailable schedd. */

	/* Resource object data */
	LL_ResourceName=2900,            /* char ** : Resource Name   */
	LL_ResourceInitialValue,         /* int * : # of initial resources  */
	LL_ResourceAvailableValue,       /* int * : # of available resources*/
	LL_ResourceInitialValue64,       /* int64_t * # of initial resources   */
	LL_ResourceAvailableValue64,     /* int64_t * # of available resources */

	/* ResourceReq object data */
	LL_ResourceRequirementName=3000, /* char **:job's resource requirement*/
	LL_ResourceRequirementValue,     /* int  *:# of resources requested  */
	LL_ResourceRequirementValue64,   /* int64_t*: # of resources requested */

	/* WlmStat object data */
	LL_WlmStatCpuTotalUsage=3100, /* int64_t* :WLM reported total cpu usage*/
	LL_WlmStatCpuSnapshotUsage,   /* int*     :WLM snapshot cpu usage      */
	LL_WlmStatMemoryHighWater,    /* int64_t* :WLM real memory high water  */
	LL_WlmStatMemorySnapshotUsage,/* int*:WLM real memory snapshot usage   */
	LL_WlmStatVMemoryHighWater,   /* int64_t* :WLM virtual memory high water     */
	LL_WlmStatVMemorySnapshotUsage,/* int64_t* :WLM virtual memory usage         */
	LL_WlmStatLargePageMemorySnapshotUsage,/* int* :WLM large page memory usage  */

	/* MachineUsage object data */
	LL_MachUsageMachineName = 3400,   /* char **  */
	LL_MachUsageMachineSpeed,         /* double *  */
	LL_MachUsageDispUsageCount,       /* int    *  */
	LL_MachUsageGetFirstDispUsage,    /* LL_element * (DispUsage): first dispatch usage */
	LL_MachUsageGetNextDispUsage,     /* LL_element * (DispUsage): next  dispatch usage */

	/* DispatchUsage object data */
	LL_DispUsageEventUsageCount=3500, /* int  *  */
	LL_DispUsageGetFirstEventUsage,   /* LL_element * (EventUsage): first event usage */
	LL_DispUsageGetNextEventUsage,    /* LL_element * (EventUsage): next  event usage */
	LL_DispUsageStepUserTime64,       /* int64_t * : dispatch usage step user time - microseconds */
	LL_DispUsageStepSystemTime64,     /* int64_t * : dispatch usage step system time - microseconds */
	LL_DispUsageStepMaxrss64,         /* int64_t * : dispatch usage step ru_maxrss */
	LL_DispUsageStepIxrss64,          /* int64_t * : dispatch usage step ru_ixrss */
	LL_DispUsageStepIdrss64,          /* int64_t * : dispatch usage step ru_idrss */
	LL_DispUsageStepIsrss64,          /* int64_t * : dispatch usage step ru_isrss */
	LL_DispUsageStepMinflt64,         /* int64_t * : dispatch usage step ru_minflt */
	LL_DispUsageStepMajflt64,         /* int64_t * : dispatch usage step ru_majflt */
	LL_DispUsageStepNswap64,          /* int64_t * : dispatch usage step ru_nswap */
	LL_DispUsageStepInblock64,        /* int64_t * : dispatch usage step ru_inblock */
	LL_DispUsageStepOublock64,        /* int64_t * : dispatch usage step ru_oublock */
	LL_DispUsageStepMsgsnd64,         /* int64_t * : dispatch usage step ru_msgsnd */
	LL_DispUsageStepMsgrcv64,         /* int64_t * : dispatch usage step ru_msgrcv */
	LL_DispUsageStepNsignals64,       /* int64_t * : dispatch usage step ru_nsignals */
	LL_DispUsageStepNvcsw64,          /* int64_t * : dispatch usage step ru_nvcsw */
	LL_DispUsageStepNivcsw64,         /* int64_t * : dispatch usage step ru_nivcsw */
	LL_DispUsageStarterUserTime64,    /* int64_t * : dispatch usage starter user time - microseconds */
	LL_DispUsageStarterSystemTime64,  /* int64_t * : dispatch usage starter system time - microseconds */
	LL_DispUsageStarterMaxrss64,      /* int64_t * : dispatch usage starter ru_maxrss */
	LL_DispUsageStarterIxrss64,       /* int64_t * : dispatch usage starter ru_ixrss */
	LL_DispUsageStarterIdrss64,       /* int64_t * : dispatch usage starter ru_idrss */
	LL_DispUsageStarterIsrss64,       /* int64_t * : dispatch usage starter ru_isrss */
	LL_DispUsageStarterMinflt64,      /* int64_t * : dispatch usage starter ru_minflt */
	LL_DispUsageStarterMajflt64,      /* int64_t * : dispatch usage starter ru_majflt */
	LL_DispUsageStarterNswap64,       /* int64_t * : dispatch usage starter ru_nswap */
	LL_DispUsageStarterInblock64,     /* int64_t * : dispatch usage starter ru_inblock */
	LL_DispUsageStarterOublock64,     /* int64_t * : dispatch usage starter ru_oublock */
	LL_DispUsageStarterMsgsnd64,      /* int64_t * : dispatch usage starter ru_msgsnd */
	LL_DispUsageStarterMsgrcv64,      /* int64_t * : dispatch usage starter ru_msgrcv */
	LL_DispUsageStarterNsignals64,    /* int64_t * : dispatch usage starter ru_nsignals */
	LL_DispUsageStarterNvcsw64,       /* int64_t * : dispatch usage starter ru_nvcsw */
	LL_DispUsageStarterNivcsw64,      /* int64_t * : dispatch usage starter ru_nivcsw */

	/* EventUsage object data */
	LL_EventUsageEventID = 3600,      /* int   *  */
	LL_EventUsageEventName,           /* char  **  */
	LL_EventUsageEventTimestamp,      /* int   *  */
	LL_EventUsageStepUserTime64,      /* int64_t * : event usage step user time - microseconds */
	LL_EventUsageStepSystemTime64,    /* int64_t * : event usage step system time - microseconds */
	LL_EventUsageStepMaxrss64,        /* int64_t * : event usage step ru_maxrss */
	LL_EventUsageStepIxrss64,         /* int64_t * : event usage step ru_ixrss */
	LL_EventUsageStepIdrss64,         /* int64_t * : event usage step ru_idrss */
	LL_EventUsageStepIsrss64,         /* int64_t * : event usage step ru_isrss */
	LL_EventUsageStepMinflt64,        /* int64_t * : event usage step ru_minflt */
	LL_EventUsageStepMajflt64,        /* int64_t * : event usage step ru_majflt */
	LL_EventUsageStepNswap64,         /* int64_t * : event usage step ru_nswap */
	LL_EventUsageStepInblock64,       /* int64_t * : event usage step ru_inblock */
	LL_EventUsageStepOublock64,       /* int64_t * : event usage step ru_oublock */
	LL_EventUsageStepMsgsnd64,        /* int64_t * : event usage step ru_msgsnd */
	LL_EventUsageStepMsgrcv64,        /* int64_t * : event usage step ru_msgrcv */
	LL_EventUsageStepNsignals64,      /* int64_t * : event usage step ru_nsignals */
	LL_EventUsageStepNvcsw64,         /* int64_t * : event usage step ru_nvcsw */
	LL_EventUsageStepNivcsw64,        /* int64_t * : event usage step ru_nivcsw */
	LL_EventUsageStarterUserTime64,   /* int64_t * : event usage starter user time - microseconds */
	LL_EventUsageStarterSystemTime64, /* int64_t * : event usage starter system time - microseconds */
	LL_EventUsageStarterMaxrss64,     /* int64_t * : event usage starter ru_maxrss */
	LL_EventUsageStarterIxrss64,      /* int64_t * : event usage starter ru_ixrss */
	LL_EventUsageStarterIdrss64,      /* int64_t * : event usage starter ru_idrss */
	LL_EventUsageStarterIsrss64,      /* int64_t * : event usage starter ru_isrss */
	LL_EventUsageStarterMinflt64,     /* int64_t * : event usage starter ru_minflt */
	LL_EventUsageStarterMajflt64,     /* int64_t * : event usage starter ru_majflt */
	LL_EventUsageStarterNswap64,      /* int64_t * : event usage starter ru_nswap */
	LL_EventUsageStarterInblock64,    /* int64_t * : event usage starter ru_inblock */
	LL_EventUsageStarterOublock64,    /* int64_t * : event usage starter ru_oublock */
	LL_EventUsageStarterMsgsnd64,     /* int64_t * : event usage starter ru_msgsnd */
	LL_EventUsageStarterMsgrcv64,     /* int64_t * : event usage starter ru_msgrcv */
	LL_EventUsageStarterNsignals64,   /* int64_t * : event usage starter ru_nsignals */
	LL_EventUsageStarterNvcsw64,      /* int64_t * : event usage starter ru_nvcsw */
	LL_EventUsageStarterNivcsw64,     /* int64_t * : event usage starter ru_nivcsw */
	/* Class object data */
	LL_ClassName = 3700,         /* char ** : The name of the class */
	LL_ClassPriority,            /* int *  : The class system priority */
	LL_ClassExcludeUsers,        /* char *** : users not permitted to use the class */
	LL_ClassIncludeUsers,        /* char *** : users permitted to use the class */
	LL_ClassExcludeGroups,       /* char *** : groups not permitted to use the class */
	LL_ClassIncludeGroups,       /* char *** : groups permitted to use the class */
	LL_ClassAdmin,               /* char *** : administrators for the class */
	Unused_LL_ClassNqsClass,     /* Reserved for future use  */
	Unused_LL_ClassNqsSubmit,    /* Reserved for future use  */
	Unused_LL_ClassNqsQuery,     /* Reserved for future use  */
	LL_ClassMaxProcessors,       /* Obsolete  - replaced by LL_ClassMaxNode */
	LL_ClassMaxJobs,             /* int *    : The maximum number of job steps that can run at any time */
	LL_ClassGetFirstResourceRequirement, /* LL_element * (ResourceReq) */
	LL_ClassGetNextResourceRequirement,  /* LL_element * (ResourceReq) */
	LL_ClassComment,             /* char **  : Class comment */
	LL_ClassCkptDir,             /* char **  : The directory for checkpoint files */
	LL_ClassCkptTimeHardLimit,   /* int64_t *: ckpt time hard limit */
	LL_ClassCkptTimeSoftLimit,   /* int64_t *: ckpt time soft limit */
	LL_ClassWallClockLimitHard,  /* int64_t *: wall clock hard limit */
	LL_ClassWallClockLimitSoft,  /* int64_t *: wall clock soft limit */
	LL_ClassCpuStepLimitHard,    /* int64_t *: Hard Job_cpu_limit */
	LL_ClassCpuStepLimitSoft,    /* int64_t *: Soft Job_cpu_limit */
	LL_ClassCpuLimitHard,        /* int64_t *: cpu hard limit */
	LL_ClassCpuLimitSoft,        /* int64_t *: cpu soft limit */
	LL_ClassDataLimitHard,       /* int64_t *: data hard limit */
	LL_ClassDataLimitSoft,       /* int64_t *: data soft limit */
	LL_ClassCoreLimitHard,       /* int64_t *: core file hard limit */
	LL_ClassCoreLimitSoft,       /* int64_t *: core file soft limit */
	LL_ClassFileLimitHard,       /* int64_t *: file size hard limit */
	LL_ClassFileLimitSoft,       /* int64_t *: file size soft limit */
	LL_ClassStackLimitHard,      /* int64_t *: stack size hard limit */
	LL_ClassStackLimitSoft,      /* int64_t *: stack size soft limit */
	LL_ClassRssLimitHard,        /* int64_t *: resident set size hard limit */
	LL_ClassRssLimitSoft,        /* int64_t *: resident set size soft limit */
	LL_ClassNice,                /* int *    : The nice value */
	LL_ClassFreeSlots,           /* int *    : The number of available initiators */
	LL_ClassMaximumSlots,        /* int *    : The total number of configured initiators */
	LL_ClassConstraints,         /* int *    : Whether values of Maximum and Free Slots are constrained by MAX_STARTERS and MAXJOBS */
	Unused_LL_ClassExecutionFactor,/* reserved for future use */
	LL_ClassMaxTotalTasks,       /* int *    : The value for Max_total_tasks */
	LL_ClassPreemptClass,        /* char **  : The PREEMPT_CLASS rule */
	LL_ClassStartClass,          /* char **  : The START_CLASS rule */
	LL_ClassMaxProtocolInstances,/* int *    : The maximum number of windows per protocol per task */
	LL_ClassGetFirstUser,        /* LL_element *: Get the first user within the class */
	LL_ClassGetNextUser,         /* LL_element *: Get the next user within the class */
	LL_ClassDefWallClockLimitHard,  /* int64_t *: default wall clock hard limit */
	LL_ClassDefWallClockLimitSoft,  /* int64_t *: default wall clock soft limit */
	LL_ClassGetFirstNodeResourceRequirement, /* LL_element * (ResourceReq) */
	LL_ClassGetNextNodeResourceRequirement,  /* LL_element * (ResourceReq) */
	LL_ClassAsLimitHard,         /* int64_t *: as hard limit */
	LL_ClassAsLimitSoft,         /* int64_t *: as soft limit */
	LL_ClassNprocLimitHard,      /* int64_t *: nproc hard limit */
	LL_ClassNprocLimitSoft,      /* int64_t *: nproc soft limit */
	LL_ClassMemlockLimitHard,    /* int64_t *: memlock hard limit */
	LL_ClassMemlockLimitSoft,    /* int64_t *: memlock soft limit */
	LL_ClassLocksLimitHard,      /* int64_t *: locks hard limit */
	LL_ClassLocksLimitSoft,      /* int64_t *: locks soft limit */
	LL_ClassNofileLimitHard,     /* int64_t *: nofile hard limit */
	LL_ClassNofileLimitSoft,     /* int64_t *: nofile soft limit */
	LL_ClassExcludeBg,          /* char *** : list of bg not permitted to use in class */
	LL_ClassIncludeBg,          /* char *** : list of bg permitted to use in class */
	LL_ClassAllowScaleAcrossJobs, /* int * : The boolean allow scale across jobs status for this class */
	LL_ClassGetFirstMaxResourceRequirement, /* LL_element* (ResourceReq) max consumable resource per task */
	LL_ClassGetNextMaxResourceRequirement, /* LL_element* (ResourceReq) max consumable resource per task */
	LL_ClassGetFirstMaxNodeResourceRequirement, /* LL_element* (ResourceReq) max consumable resource per node */
	LL_ClassGetNextMaxNodeResourceRequirement, /* LL_element* (ResourceReq) max consumable resource per node */
	LL_ClassStripingWithMinimumNetworks, /* int *:The boolean allow striping on nodes with more than half networks in READY state */
	LL_ClassMaxNode,             /* int *: Maximum number of nodes that can be requested */
	LL_ClassCollectiveGroups, /* int * : collective groups */
	LL_ClassImmSendBuffers, /* int * : Immediate Send Buffers of HFI */
	LL_ClassRestart, /* char **  : Restart */
	/* end of Class object data */

	/* Reservation object data */
	LL_ReservationID = 3800,     /* char ** : The ID of the reservation */
	LL_ReservationStartTime,     /* time_t *: The beginning time of the reservation */
	LL_ReservationDuration,      /* int *   : Reservation duration in the unit of minutes */
	LL_ReservationMachines,      /* char ***: Machines reserved by the reservation */
	LL_ReservationJobs,          /* char ***: Job steps bound to the reservation */
	LL_ReservationModeShared,    /* int *   : RESERVATION_SHARED mode is on if 1; off if 0 */
	LL_ReservationModeRemoveOnIdle, /* int *: RESERVATION_REMOVE_ON_IDLE mode is on if 1; off if 0 */
	LL_ReservationStatus,        /* int *   : The state of the reservation */
	LL_ReservationOwner,         /* char ** : The owner of the reservation */
	LL_ReservationGroup,         /* char ** : The LoadLeveler group which owns the reservation */
	LL_ReservationCreateTime,    /* time_t *: The creation time of the reservation */
	LL_ReservationModifiedBy,    /* char ** : The userid who last modified the reservation */
	LL_ReservationModifyTime,    /* time_t *: The last modification time */
	LL_ReservationUsers,         /* char ***: The users who may run jobs in the reservation */
	LL_ReservationGroups,        /* char ***: The LoadLeveler groups whose users may run jobs in the reservation */
	LL_ReservationBgCNodes,      /* int *   : The number of Blue Gene c-nodes reserved */
	LL_ReservationBgConnection,  /* int *   : The BG connection (if any) of the reserved Blue Gene resources */
	LL_ReservationBgShape,       /* int **  : The BG shape (if any) of the reserved Blue Gene resources */
	LL_ReservationBgBPs,         /* char ***: The Blue Gene BPs reserved by the reservation; if a BP is partially reserved, the reserved node cards will also be listed */
	LL_ReservationBgPartition,   /* char ** : The predefined partition whose resources are reserved */
	LL_ReservationExpiration,	/*time_t*: The expiration date and time of the recurring reservation */
	LL_ReservationCanceledOccurrences,	/*int**: An NULL-terminated array of occurrence IDs indicating which occurrences of a recurring reservation have been canceled. */
	LL_ReservationCanceledOccurrencesCount,	/* int*: The number of canceled occurrences corresponds to the length of the arrary retured by LL_ReservationCanceledOccurrences.*/
	LL_ReservationRecurrenceString,	/*char**, The specification of the reservation's recurrence as a string.*/
	LL_ReservationRecurrenceStructure,		/*LL_crontab_time*: The specification of the reservation's recurrence as a structure. */
	LL_ReservationBindingMethod,	/*int *: The binding method ofr the reservation. Bits are set using Reservation_mode_t values.*/
	LL_ReservationGetNextOccurrence,		/*LL_Element*. A pointer to the next occurrence of the Reservation. "This" Reservation can be thought of as the "First" occurrence in the GetFirst/GetNext sequence of travesing a list of LoadLeveler date objects. */
	LL_ReservationOccurrenceID,	/* int*: The occurrence ID as an integer. Always 0 for one-time reservation.*/
	LL_StepReservationBindingMethod,	/*int*: The method used to bind the step to its reservation. The value is 0 if the step is not bound to any reservation, or RESERVATION_BIND_SOFT or RESERVATION_BIND_FIRM from Reservation_mode_t.*/
	LL_StepReservationFirstOidStepBoundTo,	/* int*: for one-time job, this is the oid to which the job step is bound to, for a recurring job step, this is the first oid the job step is bound to*/
	LL_ReservationJobOids, /* int * : occurrence ids of job steps bound to the reservation  */
	LL_ReservationResType,	/* int * : The reservation type */
	LL_ReservationFlexibleJobId, /*  char **: The flexible job id of flexible reservation */
	LL_ReservationNotificationPgm,  /* char **: The notification program of the reservation */
	LL_ReservationNotificationArgs,  /* char **: The arguments specified by user */
	LL_ReservationFlexibleUserSelectionMethod, /* int *: For Flexible Reservation, user specified selection method for create */
	LL_ReservationFlexibleUserNumNodes,        /* int *: For Flexible Reservation, user specified number of nodes */
	LL_ReservationFlexibleUserHostList,        /* char ***: For Flexible Reservation, user specified list of hosts */
	LL_ReservationFlexibleUserJCF,             /* char **: For Flexible Reservation, user specified job command file */
	LL_ReservationFlexibleUserFloatingResList, /* char **: For Flexible Reservation, user specified floating consumable resource list */
	LL_ReservationGetFirstFloatingResourceRequested,  /* LL_element *: A pointer to the element associated with the first requested floating resource */
	LL_ReservationGetNextFloatingResourceRequested, /* LL_element *: A pointer to the element associated with the next requested floating resource */
	LL_ReservationGetFirstFloatingResource, /* LL_element *: A pointer to the element associated with the first floating resource */
	LL_ReservationGetNextFloatingResource, /* LL_element *: A pointer to the element associated with the next floating resource */

	LL_ReservationFlexibleUserHostFile,             /* char **: For Flexible Reservation, user specified host file */
	/* end of Reservation object data */

	/* Multicluster object data */
	LL_MClusterName = 3900,    /* char ** : The name of the multi cluster */
	LL_MClusterInboundScheddPort, /* int * : The inbound Schedd port for the multi cluster */
	LL_MClusterLocal,          /* int * : The boolean local status for this multi cluster */
	LL_MClusterInboundHosts,   /* char ** : The inbound schedd hosts(clusters) */
	LL_MClusterOutboundHosts,  /* char ** : The outbound schedd hosts(clusters) */
	LL_MClusterIncludeUsers,   /* char ** : The include users(clusters) */
	LL_MClusterExcludeUsers,   /* char ** : The exclude users(clusters) */
	LL_MClusterIncludeGroups,  /* char ** : The include groups(clusters) */
	LL_MClusterExcludeGroups,  /* char ** : The exclude groups(clusters) */
	LL_MClusterIncludeClasses, /* char ** : The include classes(clusters) */
	LL_MClusterExcludeClasses, /* char ** : The exclude classes(clusters) */
	LL_MClusterSecureScheddPort, /* int * : The secure schedd port for the cluster */
	LL_MClusterMulticlusterSecurity, /* char ** : The security method for mutlicluster */
	LL_MClusterSslCipherList, /* char ** : The list of cipher for SSL */
	LL_MClusterAllowScaleAcrossJobs, /* int * : The boolean allow scale across jobs status for this multi cluster */
	LL_MClusterMainScaleAcrossCluster, /* int * : The boolean main scale across cluster status for this multi cluster */
	/* end of Multicluster object data */

	/* MCM object data */
	LL_MCMID = 4000,					  	/* int *  : mcm id */
	LL_MCMCPUs,					          /* int *  : total available cpus on a mcm */
	LL_MCMCPUList,					      /* int ** : list of available cpus on a mcm */
	/* end of MCM object data */

	/* Blue Gene machine data */
	LL_BgMachineBPSize = 4100,     /* int ** : Size of base partitions in c-nodes in each dimension */
	LL_BgMachineSize,              /* int ** : Size of system in base partitions in each dimension */
	LL_BgMachineSwitchCount,       /* int * : Number of switches in the system */
	LL_BgMachineWireCount,         /* int * : Number of wires in the system */
	LL_BgMachinePartitionCount,    /* int * : Number of partitions defined */
	LL_BgMachineGetFirstBP,        /* LL_element * : First element in the base partition list */
	LL_BgMachineGetNextBP,         /* LL_element * : Next element in the base partition list */
	LL_BgMachineGetFirstSwitch,    /* LL_element * : First element in the switch list */
	LL_BgMachineGetNextSwitch,     /* LL_element * : Next element in the switch list */
	LL_BgMachineGetFirstWire,      /* LL_element * : First element in the wire list */
	LL_BgMachineGetNextWire,       /* LL_element * : Next element in the wire list */
	LL_BgMachineGetFirstPartition, /* LL_element * : First element in the partition list */
	LL_BgMachineGetNextPartition,  /* LL_element * : Next element in the partition list */

	/* Blue Gene base partition data */
	LL_BgBPId = 4200,                     /* char ** : Id of base partition */
	LL_BgBPState,                  /* int * : State of the base partition */
	LL_BgBPLocation,               /* int ** : Location of base partition in system in each dimension */
	LL_BgBPSubDividedBusy,         /* int * : Flag indicates small partiton active in BP */
	LL_BgBPCurrentPartition,       /* char * : Id of assigned partition */
	LL_BgBPCurrentPartitionState,  /* int * : State of assigned partition */
	LL_BgBPNodeCardCount,          /* int * : Number of node cards defined */
	LL_BgBPGetFirstNodeCard,       /* LL_element * : First element in the node card list */
	LL_BgBPGetNextNodeCard,        /* LL_element * : Next element in the node card list */
	LL_BgBPCnodeMemory,            /* int * : C-node memory of the base partition */
	LL_BgBPIONodeCount,            /* int * : The total number of I/O nodes in the base partition */

	/* Blue Gene switch data */
	LL_BgSwitchId = 4300,          /* char ** : Id of switch */
	LL_BgSwitchBasePartitionId,    /* char ** : Id of base partition connected to switch */
	LL_BgSwitchState,              /* int * : State of the switch */
	LL_BgSwitchDimension,          /* int * : Dimension the switch is associated with */
	LL_BgSwitchConnCount,          /* int * : Number of connections in the switch */
	LL_BgSwitchGetFirstConn,       /* LL_element * : First element in the switch connection list */
	LL_BgSwitchGetNextConn,        /* LL_element * : Next element in the switch connection list */

	/* Blue Gene switch connection data */
	LL_BgPortConnToSwitchPort = 4400, /* int * : Id of to switch port */
	LL_BgPortConnFromSwitchPort,   /* int * : Id of from switch port */
	LL_BgPortConnCurrentPartition, /* char ** : Id of partition the connection is assigned to */
	LL_BgPortConnCurrentPartitionState, /* int * : State of partition the connection is assigned to */

	/* Blue Gene wire data */
	LL_BgWireId = 4500,             /* char ** : Id of the wire */
	LL_BgWireState,                 /* int * : State of the wire */
	LL_BgWireFromPortCompId,        /* char ** : Id of the source component */
	LL_BgWireFromPortId,            /* int * : Id of the source port  */
	LL_BgWireToPortCompId,          /* char ** : Id of the destination component */
	LL_BgWireToPortId,              /* int * : Id of the destination port  */
	LL_BgWireCurrentPartition,      /* char ** : Id of partiton which wire is assigned to */
	LL_BgWireCurrentPartitionState, /* int * : State of partiton which wire is assigned to */

	/* Blue Gene partition data */
	LL_BgPartitionId = 4600,        /* char ** : Partition Id */
	LL_BgPartitionState,            /* int * : Partition State */
	LL_BgPartitionBPCount,          /* int * : Number of base partitions in partition */
	LL_BgPartitionSwitchCount,      /* int * : Number of switches in partition */
	LL_BgPartitionBPList,           /* char *** : List of base partition ids assigned to partition */
	LL_BgPartitionGetFirstSwitch,   /* LL_element * : First element in the switch list */
	LL_BgPartitionGetNextSwitch,    /* LL_element * : Next element in the switch list */
	LL_BgPartitionNodeCardList,     /* char *** : List of node card ids assigned to partition */
	LL_BgPartitionConnection,       /* int * : Connection type */
	LL_BgPartitionOwner,            /* char ** : User who owns the partition */
	LL_BgPartitionMode,             /* int * : Node mode of the partition (For Blue Gene/L only) */
	LL_BgPartitionSmall,            /* int * : Partition is smaller than base partition */
	LL_BgPartitionMLoaderImage,     /* char ** : File name of machine loader image */
	LL_BgPartitionBLRTSImage,       /* char ** : File name of cnode's kernel image (For Blue Gene/L only) */
	LL_BgPartitionLinuxImage,       /* char ** : File name of I/O nodes Linux image (For Blue Gene/L only) */
	LL_BgPartitionRamDiskImage,     /* char ** : File name of ramdisk image (For Blue Gene/L only) */
	LL_BgPartitionDescription,      /* char ** : Partition description */
	LL_BgPartitionSize,             /* int *   : Size of the partition */
	LL_BgPartitionShape,            /* int **  : Shape of the partition */
	LL_BgPartitionUserList,         /* char ***: List of users for the partition */
	LL_BgPartitionIONodeCount,      /* int *   : The total number of I/O nodes in the Partition */
	LL_BgPartitionCnLoadImage,      /* char ** : Comma-separated list of images to load on the compute nodes (Not for Blue Gene/L) */
	LL_BgPartitionIoLoadImage,      /* char ** : Comma-separated list of images to load on the I/O nodes (Not for Blue Gene/L) */
	LL_BgPartitionIONodeList,       /* char ***: List of I/O nodes used by a BP in the partition (Not for Blue Gene/L) */
	LL_BgPartitionType,			/* int * : Type of the partition (Not for Blue Gene/L) */

	/* Blue Gene node card data */
	LL_BgNodeCardId = 4700,            /* char ** : Id of the wire */
	LL_BgNodeCardState,                /* int * : State of the wire */
	LL_BgNodeCardQuarter,              /* int * : quarter of BP which node card is in */
	LL_BgNodeCardCurrentPartition,     /* char ** : Id of partiton which node card is assigned to */
	LL_BgNodeCardCurrentPartitionState,/* int * : State of partiton which node card is assigned to */
	LL_BgNodeCardSubDividedBusy,       /* int * : Flag indicates at least one partition is using part of the node card (Not for Blue Gene/L) */
	LL_BgNodeCardIONodeCount,          /* int * : Total number of I/O nodes in the node card */
	LL_BgNodeCardGetFirstIONode,       /* LL_element * : First I/O node in the node card */
	LL_BgNodeCardGetNextIONode,        /* LL_element * : Next I/O node in the node card */

	/* Blue Gene IONode data (Not for Blue Gene/L) */
	LL_BgIONodeId = 4750,              /* char ** : Id of the I/O Node */
	LL_BgIONodeIPAddr,                 /* char ** : IP address of the I/O Node */
	LL_BgIONodeCurrentPartition,       /* char ** : Id of the partiton the I/O node belongs to */
	LL_BgIONodeCurrentPartitionState,  /* int *   : State of the partiton the I/O node belongs to */

	/* ClusterFile data */
	LL_ClusterFileLocalPath = 4800, /* char **: expanded local file pathname */
	LL_ClusterFileRemotePath, /* char **: expanded remote file pathname */
	/* end of ClusterFile data */

	/* Fair Share object data */
	LL_FairShareCurrentTime = 4900, /* time_t * : The time of query          */
	LL_FairShareTotalShares,        /* int * : FAIR_SHARE_TOTAL_SHARES value */
	LL_FairShareInterval,           /* int * : FAIR_SHARE_INTERVAL value     */
	LL_FairShareNumberOfEntries,    /* int * : Number of users and groups    */
	LL_FairShareEntryNames,         /* char *** : Names of users and groups  */
	LL_FairShareEntryTypes,         /* int ** : Types:0 for user,1 for group */
	LL_FairShareAllocatedShares,    /* int ** : Number of allocated shares   */
	LL_FairShareUsedShares,         /* int ** : Number of used shares        */
	LL_FairShareUsedBgShares,       /* int ** : Number of used BG shares     */
	/* end of Fair Share object data */

	/* ClassUser data */
	LL_ClassUserName = 5000,         /* char ** : Name of the user in the class */
	LL_ClassUserMaxIdle,             /* int * : The maximum number of idle job steps */
	LL_ClassUserMaxQueued,           /* int * : The maximum number of total job steps in the queue */
	LL_ClassUserMaxJobs,             /* int * : The maximum number of running job steps */
	LL_ClassUserMaxTotalTasks,       /* int * : The maximum number of running tasks */
	/* end of ClassUser data */

	/* Class NetworkUsage */
	LL_NetworkUsageNetworkId = 5100, /* uint64_t * : network ID */
	LL_NetworkUsageInstances,	 /* int * : number of Instances */
	LL_NetworkUsageProtocol,	 /* char ** : Names of requested protocol */
	LL_NetworkUsageMode,		 /* char ** : Names of modeIP or US */
	LL_NetworkUsageWindows,		 /* int * : window ids count */
	LL_NetworkUsageRcxtBlocks,	 /* int * : number of rCxt blocks */
	LL_StepGetFirstNetworkUsage,	 /* LL_element* : First element in the network list */
	LL_StepGetNextNetworkUsage,	 /* LL_element* : Next element in the network list */
	LL_CollectiveGroups,		 /* int * : Number of collective groups. */
	LL_ImmSendBuffers,		 /* int * : Number of immediate send buffers. */
	/* end of Class NetworkUsage */

	/* LlMachineGroup objects */
	LL_MachineGroupName = 5200,		/* char ** : Name of the machine group */
	LL_MachineGroupScheddAvail,		/* int * : Number of available schedds within the machine group */
	LL_MachineGroupScheddDown,		/* int * : Number of down schedds within the machine group */
	LL_MachineGroupScheddDrained,		/* int * : Number of drained schedds within the machine group */
	LL_MachineGroupScheddDraining,		/* int * : Number of draining schedds within the machine group */
	LL_MachineGroupScheddRunning,		/* int * : Number of running schedds within the machine group */
	LL_MachineGroupScheddTotal,		/* int * : Number of machines within the machine groups which are configed to have schedd daemon */
	LL_MachineGroupScheddTotalJobSteps,	/* int * : Number of job steps in the job queues of Schedd machines within the machine group */
	LL_MachineGroupScheddStepStarting,	/* int * : Number of starting job steps in the job queues of Schedd machines within the machine group */
	LL_MachineGroupScheddStepIdle,		/* int * : Number of idle job steps in the job queues of Schedd machines within the machine group */
	LL_MachineGroupScheddStepHeld,		/* int * : Number of held job steps in the job queues of Schedd machines within the machine group */
	LL_MachineGroupScheddStepPending,	/* int * : Number of pending job steps in the job queues of Schedd machines within the machine group */
	LL_MachineGroupScheddStepRemovePending,	/* int * : Number of remove_pending job steps in the job queues of Schedd machines within the machine group */
	LL_MachineGroupScheddStepUnexpanded,	/* int * : Number of unexpanded job steps in the job queues of Schedd machines within the machine group */
	LL_MachineGroupScheddStepRemoved,	/* int * : Number of removed job steps in the job queues of Schedd machines within the machine group */
	LL_MachineGroupScheddStepCompleted,	/* int * : Number of completed job steps in the job queues of Schedd machines within the machine group */
	LL_MachineGroupScheddRunningJobs,	/* int * : Number of running job steps in the job queues of Schedd machines within the machine group */
	LL_MachineGroupStartdAvail,		/* int * : Number of available startds within the machine group */
	LL_MachineGroupStartdDown,		/* int * : Number of down startds within the machine group */
	LL_MachineGroupStartdDrained,		/* int * : Number of drained startds within the machine group */
	LL_MachineGroupStartdDraining,		/* int * : Number of draining startds within the machine group */
	LL_MachineGroupStartdFlush,		/* int * : Number of flush startds within the machine group */
	LL_MachineGroupStartdSuspend,		/* int * : Number of suspend startds within the machine group */
	LL_MachineGroupStartdIdle,		/* int * : Number of idle startds within the machine group */
	LL_MachineGroupStartdRunning,		/* int * : Number of running startds within the machine group */
	LL_MachineGroupStartdBusy,		/* int * : Number of busy startds within the machine group */
	LL_MachineGroupStartdTotal,		/* int * : Number of machines within the machine groups which are configed to have startd daemon */
	LL_MachineGroupStartdTotalRunningTasks,	/* int * : Number of tasks running on startd machines within the machine group */
	LL_MachineGroupNumMachs,		/* int * : Number of machines within the machine group */
	LL_MachineGroupRange,			/* char ** : The expression which defines the range of machines included in the group */
	LL_MachineGroupFeatureList,		/* char *** : A NULL-terminated array of the features defined in the machine group */
	LL_MachineGroupConfiguredClassList,	/* char *** : A NULL-terminated array of the class initiators configured in the machine group */
	LL_MachineGroupMachineMode,		/* char ** : The mode configured for the machine group */
	LL_MachineGroupMaxStarters,		/* int * : The value of the max_starters keyword */
	LL_MachineGroupPoolList,		/* int ** : An array indicating the pool numbers to which machines in this machine group belong */
	LL_MachineGroupPoolListSize,		/* int * : The number of pools configured for the machine group */
	LL_MachineGroupGetFirstResource,	/* LL_element* : A pointer to the first element associated with the machine groups configured list of resources */
	LL_MachineGroupGetNextResource,		/* LL_element* : A pointer to the next element associated with the machine groups configured list of resources */
	LL_MachineGroupExplicitlyDefinedMachines,	/* char *** : A NULL-terminated list of hostnames of machines which have explicit definitions (substanzas) within the machine group */
	LL_MachineGroupReservationPermitted,	/* int * : A boolean which indicates whether machines in this group can be reserved */
	LL_MachineGroupPrestartedStarters,	/* int * : The number of prestarted starters to start on each machine in the machine group */
	LL_MachineGroupRegion,			/* char ** : The region of the machine group */
	LL_MachineGroupMaxDstgStarters,		/* int * : The value of the dstg_max_starters keyword */
	LL_MachineGroupMachineSpeed,		/* double * : The value of machine_speed */
	/* end of LlMachineGroup objects */

	/* JobSummary objects */
	LL_JobSummaryName = 5300,		/* char ** : Name of the job summary */
	LL_JobSummaryRunningCount,		/* int * : Number of job steps in Running state. */
	LL_JobSummaryPendingCount,		/* int * : Number of job steps in Pending state. */
	LL_JobSummaryWaitingCount,		/* int * : Number of job steps in Waiting state. */
	LL_JobSummaryHeldCount,			/* int * : Number of job steps in Held state. */
	LL_JobSummaryPreemptedCount,		/* int * : Number of job steps in Preempted state. */
	/* end of JobSummary objects */

	LL_LastGetDataSpecification
};

#define FREE_SLOTS_LIMITED_BY_MAX_STARTERS    1
#define MAXIMUM_SLOTS_LIMITED_BY_MAX_STARTERS 2
#define FREE_SLOTS_LIMITED_BY_MAX_JOBS        4

enum SummaryReportType { NUMERIC = (1<<0),
                         RESOURCE = (1<<1),
                         AVGTHROUGHPUT = (1<<2),
                         MAXTHROUGHPUT= (1<<3),
                         MINTHROUGHPUT=(1<<4),
                         THROUGHPUT=(AVGTHROUGHPUT|MAXTHROUGHPUT|MINTHROUGHPUT),
                         REPORT_ALL = (NUMERIC|RESOURCE|AVGTHROUGHPUT|MAXTHROUGHPUT|MINTHROUGHPUT),
                         REPORT_DEFAULT = (RESOURCE)
                       };

enum SummarySectionType { USER = (1<<0),
                          SECTION_GROUP = (1<<1),
                          CLASS = (1<<2),
                          ACCOUNT = (1<<3),
                          UNIXGROUP = (1<<4),
                          DAY = (1<<5),
                          WEEK = (1<<6),
                          MONTH = (1<<7),
                          JOBID = (1<<8),
                          JOBNAME = (1<<9),
                          ALLOCATED = (1<<10),
                          SECTION_ALL = (USER|SECTION_GROUP|CLASS|ACCOUNT|UNIXGROUP|DAY|WEEK|MONTH|JOBID|JOBNAME|ALLOCATED),
                          SECTION_DEFAULT=(USER|SECTION_GROUP|CLASS|ACCOUNT),
                          TIME_MASK=(DAY|WEEK|MONTH)
                        };


enum SummaryDisplayFormat {
	EXTENDED_FORMAT = (1<<0),
	SUMMARY_FORMAT  = (1<<1),
	QUERY_FORMAT    = (1<<2),
	GUI_FORMAT      = (1<<3)
};

enum LL_control_op {
	LL_CONTROL_RECYCLE, LL_CONTROL_RECONFIG, LL_CONTROL_START,
	LL_CONTROL_STOP, LL_CONTROL_DRAIN, LL_CONTROL_DRAIN_STARTD,
	LL_CONTROL_DRAIN_SCHEDD, LL_UNUSED_1, LL_CONTROL_FLUSH,
	LL_CONTROL_SUSPEND, LL_CONTROL_RESUME, LL_CONTROL_RESUME_STARTD,
	LL_CONTROL_RESUME_SCHEDD, LL_CONTROL_FAVOR_JOB, LL_CONTROL_UNFAVOR_JOB,
	LL_CONTROL_FAVOR_USER, LL_CONTROL_UNFAVOR_USER,
	LL_CONTROL_HOLD_USER, LL_CONTROL_HOLD_SYSTEM, LL_CONTROL_HOLD_RELEASE,
	LL_CONTROL_PRIO_ABS, LL_CONTROL_PRIO_ADJ, LL_CONTROL_START_DRAINED,
	LL_CONTROL_DUMP_LOGS, LL_CONTROL_DUMP_LOCKS
};

#define LL_CONTROL_VERSION_22  22
#define LL_CONTROL_VERSION_310 310
#define LL_CONTROL_VERSION     310

/*   Structures to support API interfaces */
typedef struct {
	int		cluster;
	int		proc;
	char		*from_host;	/* name of the schedd host */
}
LL_STEP_ID;

typedef struct {
	int		nqs_flags;	/* flags for controlling NQS step submission */
	char		*nqs_submit;	/* NQS submit queue */
	char		*nqs_query;	/* NQS query queues */
	char		*umask;		/* value of umask on submitting machine */
}
LL_NQS;

typedef struct {
	int  cpu_hard_limit;            /* cpu time cannot exceed this */
	int  cpu_soft_limit;            /* value set by user that is LE hard limit */
	int  data_hard_limit;           /* data size cannot exceed this */
	int  data_soft_limit;           /* value set by user that is LE hard limit */
	int  core_hard_limit;           /* size of core file cannot exceed this */
	int  core_soft_limit;           /* value set by user that is LE hard limit */
	int  file_hard_limit;           /* file size cannot exceed this */
	int  file_soft_limit;           /* value set by user that is LE hard limit */
	int  rss_hard_limit;            /* resident set size upper limit */
	int  rss_soft_limit;            /* value set by user that is LE hard limit */
	int  stack_hard_limit;          /* stack size cannot exceed this */
	int  stack_soft_limit;          /* value set by user that is LE hard limit */
	int  hard_cpu_step_limit;       /* hard CPU limit for the whole job step */
	int  soft_cpu_step_limit;       /* soft CPU limit for the whole job step */
	int  hard_wall_clock_limit;     /* hard limit for elapsed time */
	int  soft_wall_clock_limit;     /* soft limit for elapsed time */
	int  ckpt_time_hard_limit;	  /* ckpt time cannot exceed this */
	int  ckpt_time_soft_limit;	  /* value set by user that is LE hard limit */
}
LL_LIMITS;

typedef struct {
	int64_t cpu_hard_limit;        /* cpu time cannot exceed this */
	int64_t cpu_soft_limit;        /* set by user that is LE hard limit */
	int64_t data_hard_limit;       /* data size cannot exceed this */
	int64_t data_soft_limit;       /* set by user that is LE hard limit */
	int64_t core_hard_limit;       /* size of core file cannot exceed this */
	int64_t core_soft_limit;       /* set by user that is LE hard limit */
	int64_t file_hard_limit;       /* file size cannot exceed this */
	int64_t file_soft_limit;       /* set by user that is LE hard limit */
	int64_t rss_hard_limit;        /* resident set size upper limit */
	int64_t rss_soft_limit;        /* set by user that is LE hard limit */
	int64_t stack_hard_limit;      /* stack size cannot exceed this */
	int64_t stack_soft_limit;      /* set by user that is LE hard limit */
	int64_t hard_cpu_step_limit;   /* hard CPU limit for the whole job step */
	int64_t soft_cpu_step_limit;   /* soft CPU limit for the whole job step */
	int64_t hard_wall_clock_limit; /* hard limit for elapsed time */
	int64_t soft_wall_clock_limit; /* soft limit for elapsed time */
	int64_t ckpt_time_hard_limit;  /* ckpt time cannot exceed this */
	int64_t ckpt_time_soft_limit;  /* set by user that is LE hard limit */
	int64_t as_hard_limit;	       /* address space size cannot exceed this */
	int64_t as_soft_limit;	       /* set by user that is LE hard limit */
	int64_t nproc_hard_limit;      /* number of process cannot exceed this */
	int64_t nproc_soft_limit;      /* set by user that is LE hard limit */
	int64_t memlock_hard_limit;    /* size of memory that may be locked cannot exceed this */
	int64_t memlock_soft_limit;    /* set by user that is LE hard limit */
	int64_t locks_hard_limit;      /* number of file locks cannot exceed this */
	int64_t locks_soft_limit;      /* set by user that is LE hard limit */
	int64_t nofile_hard_limit;     /* number of open file descriptors cannot exceed this */
	int64_t nofile_soft_limit;     /* set by user that is LE hard limit */
}
LL_LIMITS64;

typedef struct ll_event_usage {
	int  event;                       /* event identifier             */
	char *name;                       /* event name			*/
	int  time;                        /* timestamp of this event      */
	struct rusage starter_rusage;     /* usage by starter at this event */
	struct rusage step_rusage;        /* usage by user's job step at this event */
	struct ll_event_usage *next;      /* next event			*/
}
LL_EVENT_USAGE;

typedef struct ll_event_usage64 {
	int   event;                      /* event identifier */
	char  *name;                      /* event name */
	int   time;                       /* timestamp of this event */
	struct rusage64 starter_rusage64; /* usage by starter at this event */
	struct rusage64 step_rusage64;    /* usage by user's job step at this event */
	struct ll_event_usage64 *next;    /* next event */
}
LL_EVENT_USAGE64;

typedef struct ll_dispatch_usage {
	int dispatch_num;                 /* # of event usages for this dispatch */
	struct rusage starter_rusage;     /* accumulated usage by starter */
	struct rusage step_rusage;        /* accumulated usage by user's job step */
	LL_EVENT_USAGE *event_usage;      /* per event usage detail    */
	struct ll_dispatch_usage *next;   /* next dispatch	*/
}
LL_DISPATCH_USAGE;

typedef struct ll_dispatch_usage64 {
	int dispatch_num;                 /* # of event usages for this dispatch */
	struct rusage64 starter_rusage64; /* accumulated usage by starter */
	struct rusage64 step_rusage64;    /* accumulated usage by user's job step */
	LL_EVENT_USAGE64 *event_usage64;  /* per event usage detail */
	struct ll_dispatch_usage64 *next; /* next dispatch  */
}
LL_DISPATCH_USAGE64;

typedef struct ll_mach_usage {
	char *name;                        /* machine name			*/
	float machine_speed;               /* machine speed		*/
	int   dispatch_num;                /* # of dispatches for this machine */
	LL_DISPATCH_USAGE *dispatch_usage; /* per dispatch usage detail*/
	struct ll_mach_usage *next;        /* next machine			*/
}
LL_MACH_USAGE;

typedef struct ll_mach_usage64 {
	char     *name;                        /* machine name */
	float    machine_speed;                /* machine speed */
	int      dispatch_num;                 /* # of dispatches for this machine */
	LL_DISPATCH_USAGE64 *dispatch_usage64; /* per dispatch usage detail */
	struct ll_mach_usage64 *next;          /* next machine */
}
LL_MACH_USAGE64;

typedef struct {
	struct rusage starter_rusage; /* accumulated usage by starters */
	struct rusage step_rusage;	/* accumulated usage by user's step */
	LL_MACH_USAGE	*mach_usage;    /* detail usage                   */
}
LL_USAGE;

typedef struct {
	struct rusage64 starter_rusage64; /* accumulated usage by starters */
	struct rusage64 step_rusage64;    /* accumulated usage by user's step */
	LL_MACH_USAGE64 *mach_usage64;    /* detail usage */
}
LL_USAGE64;

typedef struct {
	int		x;          /* Number of compute nodes in x-direction. */
	int		y;          /* Number of compute nodes in y-direction. */
	int		z;          /* Number of compute nodes in z-direction. */
}
LL_BG_SHAPE;

typedef struct {
	/* The following are inputs needed before scheduling is performed */

	char		*step_name;	/* step name */
	char		*requirements;	/* step requirements */
	char		*preferences;	/* step preferences */
	int		prio;		/* user step priority */
	char		*dependency;	/* step step dependancy */
	char		*group_name;	/* step group name */
	char		*stepclass;	/* step class */
	int		start_date;	/* Don't start before this date */
	int		flags;		/* Special step characteristics */
	int		min_processors; /* minimum # of requested processors */
	int		max_processors; /* maximum # of requested processors */
	char		*account_no;	/* Account number associated with step */
	char		*comment;	/* users comment about the step */

	/* The following are valid after queuing has occured */

	LL_STEP_ID	id;		/* step id */
	int		q_date;		/* UNIX time step was submitted */
	int		status;		/* Running, unexpanded, completed,.. */
	/* The following are valid after scheduling has occured */

	int		num_processors; /* actual number of assigned processors */
	char		**processor_list;/* list of processors on which to run step */
	char		**smt_status_list;/* list for saving smt status of machines
					     running steps*/

	/* The following are inputs needed to actually start an executable */

	char		*cmd;		/* a.out file */
	char		*args;		/* command line args */
	char		*env;		/* environment */
	char		*in;		/* file for stdin */
	char		*out;		/* file for stdout */
	char		*err;		/* file for stderr */
	char		*iwd;		/* Initial working directory */
	char		*notify_user;	/* target to be used when sending mail */
	char		*shell;		/* shell to be used */
	char		*tracker;	/* user's step tracking exit */
	char		*tracker_arg;	/* argument to tracking exit */
	int		notification;	/* Notification options */
	int		image_size;	/* Size of the virtual image in K */
	int		exec_size;	/* size of the executable */
	LL_LIMITS	limits;		/* step resource limits */
	LL_NQS	nqs_info;	/* info for a NQS step */
	int			smt_required;	/* if the SMT function is required to be turned on, off or just as_is before we start the executable */

	/* The following are valid after the executable has started */
	int     eligible_time; /* UNIX time Job eligibility time */
	int		dispatch_time; /* UNIX time negotiator dispatched job */
	int		start_time;	/* UNIX time the starter started */
	int		unused1;	/* reserved for future use */

	/* The following are valid after the executable has completed/terminated */

	int		completion_code;/* step exit status */
	int		completion_date;/* UNIX time step was completed */
	int		start_count;	/* times step has been started */
	LL_USAGE	usage_info;	/* step usage information */

	/* Priorities set from the admin file stanzas */
	int		user_sysprio;	/* user component of system  job priority */
	int		group_sysprio;	/* group component of system job priority */
	int		class_sysprio;	/* class component of system job priority */
	int		number;   	/* user number. Substitutable in LoadL config macros */

	/* Consumable resources requested and adapter pinned memory used */
	int		cpus_requested;           /* ConsumableCpus requested */
	int		virtual_memory_requested; /* VirtualMemory requested */
	int		memory_requested;         /* Memory (real) requested */
	int		adapter_used_memory;      /* pinned memory used by adapters */

	int   adapter_req_count;   /* number of adapter_req records */
	void  **adapter_req; /* adapter requirements - step->getFirstAdapterReq() ...  */

	/* 64-bit elements and structures */
	int64_t  image_size64;                /* Size of the virtual image in K */
	int64_t  exec_size64;                 /* size of the executable */
	int64_t  virtual_memory_requested64;  /* VirtualMemory requested */
	int64_t  memory_requested64;          /* Memory (real) requested */
	LL_LIMITS64 limits64;                 /* 64-bit step resource limits */
	LL_USAGE64  usage_info64;             /* 64-bit step usage information */

	/* Checkpoint statistics  */
	int	good_ckpt_start_time;		/* Time stamp of last successful ckpt */
	int	accum_ckpt_time;		/* accumulatd time job has spent checkpointing */
	char	*ckpt_dir;			/* Checkpoint directory */
	char	*ckpt_file;			/* Checkpoint file	*/

	/* Large Page Data/Heap support */
	char   *large_page;      /* Large Page requirement */

	/* RSet Support */
	char *rset;              /* RSet name requested by the step*/
	char *mcm_affinity_options; /* mcm affinity options requested by the step*/

	/* RDMA Support */
	int           bulkxfer;         /* Did step request BLOCKXFER? */
	int           rcxtblocks;       /* User rCxt blocks */
	int           userxfer;         /* Did step request USERXFER? */

	/* Advance Reservation */
	char   *reservation_id;              /* Reservation id */
	char   *requested_reservation_id;    /* Requested reservaiton id */

	/* AIX Advanced Accounting   */
	int64_t	acct_key;		     /* Job accounting key */

	/* Blue Gene Support   */
	int     bg_req_size;		     /* Requested size of Blue Gene partition in units of compute nodes.*/
	int     bg_alloc_size;		     /* Allocated size of Blue Gene partition in units of compute nodes.*/
	LL_BG_SHAPE 	bg_req_shape;	     /* Requested shape of Blue Gene Partition in units of cubical compute nodes.*/
	LL_BG_SHAPE 	bg_alloc_shape;	     /* Allocated shape of Blue Gene Partition in units of cubical compute nodes.*/
	char 	*bg_req_connection;	     /* Requested type of wiring for the partition.*/
	char 	*bg_alloc_connection;	     /* Allocated type of wiring for the partition.*/
	char 	*bg_mode;		     /* Requested mode of the partition. */
	char 	*bg_rotate;		     /* Whether the scheduler is free to rotate the requested Blue Gene shape.*/
	char 	*bg_job_id;		     /* Id of the Blue Gene job. */
	char 	*bg_partition_id;	     /* Id of the partition allocated for the job. */
	char 	*bg_alloc_partition;	     /* Id of the partition allocated for the job. */
	char 	*bg_req_partition;	     /* The partition Id requested for the job. */
	char 	*bg_error_text;		     /* Error text in the Blue Gene  job record. */
	char	*bg_requirements;            /* Requirements for Blue Gene machine  */
	char	*bg_partition_type;          /* Type of Blue Gene partition */

	/* Task affinity support */
	char *task_affinity;        /* task_affinity keyword value */
	int cpus_per_core;          /* cpus_per_core keyword value */
	int parallel_threads;          /* parallel_threads keyword value */

	/* Reserved fields */
	LL_element   *reserved001;      	/* Reserved */

	/* Metacluster data */
	int     metacluster_job_id;          /* The Metacluster ID of the job      */

	/* Time of job in User Hold state*/
	time_t		user_hold_time;

	/* LOADL_BG_BPS and LOADL_BG_IONODES values for a Blue Gene job */
	char *loadl_bg_bps;         /* BPs used by a Blue Gene job */
	char *loadl_bg_ionodes;     /* I/O nodes used by a Blue Gene BP */

	/* recurring job */
	int recurring;				/* Flag for recurring job */
	/* cluster option */
	char    *cluster_option;    /* cluster_option keyword value */

	char *soft;
	int trace;
	/*topology scheduling*/
	char *topology_group;
	char *topology_require;
	char   *flexible_reservation_id;    /* Flexible reservaiton id */
}
LL_job_step;

typedef struct LL_job {
	int		version_num;	/* LL_JOB_VERSION */
	char		*job_name;	/* job name */
	char		*owner;		/* login of person submitting job */
	char		*groupname;	/* group name of owner's login group */
	uid_t		uid;		/* unix userid of submitter */
	gid_t		gid;		/* unix group id of submitter */
	char		*submit_host;	/* Host of job submission */
	int		steps;		/* number of steps in job */
	LL_job_step	**step_list;	/* ptr to array of ptrs to job steps */
}
LL_job;

typedef struct LL_node {
	char          *nodename;              /* Name of this node */
	int           version_num;            /* PROC_VERSION */
	int           configtimestamp;        /* Date and time of last reconfig */
	int		time_stamp;		/* Time stamp of data */
	int           virtual_memory;         /* Available swap space in kilobytes */
	int           memory;                 /* Physical memory */
	int           disk;                   /* Avail space (KB) in execute directory */
	float         loadavg;                /* Berkeley one minute load average */
	float         speed;                  /* Speed associated with the node */
	int           max_starters;           /* Max number of jobs allowed */
	int           dstg_max_starters;      /* Max number of data staging steps allowed */
	int           pool;           	/* Pool number associated with node */
	int           cpus;           	/* number of CPUs */
	char          *state;         	/* Startd state */
	int           keywordidle;   		/* seconds since keyboard activity */
	int           totaljobs;      	/* total number of submitted jobs */
	char          *arch;          	/* Hardware Architecture */
	char          *opsys;         	/* Operating system */
	char          **adapter;      	/* Names of available adapters */
	char          **feature;      	/* set of all features */
	char          **job_class;    	/* Job classes allowed to run */
	char          **initiators;   	/* Initiators available */
	LL_STEP_ID    *steplist;    		/* steps allocated to run */
	int64_t       virtual_memory64;/* Available swap space in kilobytes */
	int64_t       memory64;        /* Physical memory */
	int64_t       disk64;          /* Avail space (KB) in execute directory */
}
LL_node;

/* This macro allows the mem field to be accessed by the name api_rcxtblocks */
/* so that when rCxt blocks are supported by the adapters the name of the    */
/* field is consistent with its meaing                                       */
#define api_rcxtblocks mem
typedef struct ll_adapter_usage {
	char * dev_name;		/* Adapter device name */
	char * protocol;		/* MPI, LAPI or MPI_LAPI */
	char * subsystem;		/* IP or US */
	int    wid;			/* window id */
	uint64_t mem;			/* adapter window memory */
}
LL_ADAPTER_USAGE;

typedef struct ll_network_usage {
	uint64_t network_id;
	char     *network_type;		/* switch, multilink, ethernet */

	char     *protocol;		/* MPI, LAPI, MPI_LAPI */
	char     *subsystem;		/* IP or US */

	int      windows_per_instance;
	int      instances_per_task;

	int      exclusive;		/* 0: false | 1: true */
}
LL_NETWORK_USAGE;


typedef struct LL_start_job_info_ext {
	int             version_num;
	LL_STEP_ID      StepId;
	char            **nodeList;
	int adapterUsageCount;
	LL_ADAPTER_USAGE * adapterUsage;
	int networkUsageCount;
	LL_NETWORK_USAGE * networkUsage;
}
LL_start_job_info_ext;

typedef struct LL_terminate_job_info {
	int             version_num;
	LL_STEP_ID      StepId;
	char            *msg;
}
LL_terminate_job_info;

/*
 *  Notification options.
 */
#define LL_NOTIFY_ALWAYS	0
#define LL_NOTIFY_COMPLETE	1
#define LL_NOTIFY_ERROR		2
#define LL_NOTIFY_NEVER		3
#define LL_NOTIFY_START		4

/*
 *  Status values.
 */
#define LL_IDLE			0
#define LL_STARTING		1
#define LL_RUNNING		2
#define LL_REMOVED		3
#define LL_COMPLETED		4
#define LL_HOLD			5
#define LL_DEFERRED		6
#define LL_SUBMISSION_ERR	7
#define LL_VACATE		8
#define LL_NOTRUN		9
#define LL_NOTQUEUED            10
#define LL_MAX_STATUS		10

/*
 *  Step flags.
 */
#define LL_CHECKPOINT		(1<<0)
#define LL_SYSTEM_HOLD		(1<<1)
#define LL_USER_HOLD		(1<<2)
#define LL_RESTART		(1<<3)
#define LL_CPU_LIMIT_USER	(1<<4)
#define LL_CORE_LIMIT_USER	(1<<5)
#define LL_DATA_LIMIT_USER	(1<<6)
#define LL_FILE_LIMIT_USER	(1<<7)
#define LL_RSS_LIMIT_USER	(1<<8)
#define LL_STACK_LIMIT_USER	(1<<9)
#define LL_NQS_STEP		(1<<10)
#define LL_STEP_PARALLEL	(1<<11)
#define LL_STEP_PVM3		(1<<12)
#define LL_IMMEDIATE		(1<<13)
#define LL_NO_ALLOCATE		(1<<14)
#define LL_INTERACTIVE		(1<<15)
#define LL_API_ACTIVE		(1<<16)
#define LL_API_SYNC_START	(1<<17)
#define LL_NODE_USAGE_NOT_SHARED (1<<18)
#define LL_RESTART_FROM_CKPT	(1<<19)
#define LL_CHECKPOINT_INTERVAL	(1<<20)
#define LL_RESTART_SAME_NODES	(1<<21)
#define LL_STEP_BLUEGENE	(1<<22)
#define LL_STEP_MPICH           (1<<23)
#define LL_METACLUSTER_JOB (1<<24)
#define LL_AS_LIMIT_USER	(1<<25)		/* add 11/2006 */
#define LL_NPROC_LIMIT_USER	(1<<26)		/* add 11/2006 */
#define LL_MEMLOCK_LIMIT_USER	(1<<27)		/* add 11/2006 */
#define LL_LOCKS_LIMIT_USER	(1<<28)		/* add 11/2006 */
#define LL_NOFILE_LIMIT_USER	(1<<29)		/* add 11/2006 */


#define LL_JOB_VERSION		(210)
#define LL_JOB_PROC_VERSION	9

/*
 * The following completion codes are
 * used when status is LL_SUBMISSION_ERR
 */

#define LL_NO_STORAGE		(1)
#define LL_BAD_STATUS		(2)
#define LL_BAD_NOTIFY		(3)
#define LL_BAD_CMD		(4)
#define LL_BAD_EXEC		(5)
#define LL_BAD_REQUIREMENTS	(6)
#define LL_BAD_PREFERENCES	(7)
#define LL_BAD_DEPENDENCY	(8)
#define LL_BAD_ACCOUNT_NO	(9)
#define LL_BAD_PRIO		(10)
#define LL_BAD_GROUP_CONFIG	(11)
#define LL_BAD_GROUP_NAME	(12)
#define LL_BAD_CLASS_CONFIG	(13)
#define LL_BAD_CLASS		(14)
#define LL_BAD_TRANSMIT		(15)

/*
 *      Values for accounting events
 */
#define	LL_LOADL_EVENT		1
#define	LL_INSTALLATION_EVENT	2

/*
 *      Values for scheduling API
 */
#define	LL_PROC_VERSION		9

/***********************************************************************
 * Status codes to support external scheduler.
 **********************************************************************/
#define     API_OK                   0  /* API call runs to complete */
#define     API_INVALID_INPUT        -1 /* Invalid input */
#define     API_INPUT_NOT_VALID      -1 /* input not valid */
#define     API_CANT_CONNECT         -2 /* can't connect to daemon */
#define     API_CANT_MALLOC          -3	/* out of memory */
#define     API_CONFIG_ERR           -4	/* Error from init_params() */
#define     API_CANT_FIND_PROC       -5	/* can't find proc */
#define     API_CANT_TRANSMIT        -6	/* xdr error */
#define     API_CANT_AUTH            -7 /* can't authorize */
#define     API_WRNG_PROC_VERSION    -8 /* Wrong proc version */
#define     API_WRNG_PROC_STATE      -9 /* Wrong proc state */
#define     API_MACH_NOT_AVAIL      -10	/* Machine not available */
#define     API_CANT_FIND_RUNCLASS  -11 /* Can't find run class */
#define     API_REQ_NOT_MET         -12 /* Can't meet requirements */
#define     API_WRNG_MACH_NO   	    -13 /* Wrong machine number */
#define     API_LL_SCH_ON     	    -14 /* Internal scheduler on */
#define     API_MACH_DUP      	    -15 /* Duplicate machine found */
#define     API_NO_DCE_ID           -16 /* User does not have a valid dce id */
#define     API_NO_DCE_CRED         -17 /* No dce credentials */
#define     API_INSUFFICIENT_DCE_CRED -18 /* Insufficient dce credential lifetime */
#define     API_64BIT_DCE_ERR       -19   /* 64-bit API is not supported when DCE is enabled. */

#define     API_BAD_ADAPTER_USAGE   -20	/* Job Start Adapter usage info is inconsistent */
#define     API_BAD_ADAPTER_DEVICE  -21 /* JobStart adapter usage info specified an adapter not on machine */
#define     API_BAD_ADAPTER_USAGE_COUNT -22 /* Wrong number of entries in adapter usage information to JobStart */
#define     API_BAD_ADAPTER_USAGE_PATTERN -23 /* The same protocol pattern is not specified for each task in the JobStart information */
#define     API_BAD_PROTOCOL -24        /* The JobStart information includes an unrecognized protocol string */
#define     API_INCOMPATIBLE_PROTOCOL -25 /* The JobStart information adapter usage information includes protocols that cannot be specified together (eg. PVM and MPI) */
#define     API_BAD_COMMUNICATION_SUBSYSTEM -26 /* The JobStart information adapter usage information includes a communication subsystem that is not IP or US */

#define     API_NO_DCE_SUPPORT_ERR     -27  /* This version of LL does not support DCE security. */
#define     API_NO_CTSEC_SUPPORT_ERR   -28  /* This version of LL does not support CTSEC security. */
#define     API_NO_GANG_SUPPORT_ERR    -29  /* This version of LL does not support GANG scheduling. */
#define     API_NO_PVM_SUPPORT_ERR     -30  /* This version of LL does not support PVM. */
#define     API_NO_NQS_SUPPORT_ERR     -31  /* This version of LL does not support NQS. */
#define     API_STEP_NOT_IDLE     -32 	/* A specified step is not in an idle state */
#define     API_JOB_NOT_FOUND     -33 	/* A specified was not found */
#define     API_JOBQ_ERR     -34 	/* Error occured writing to job queue */
#define     API_CANT_LISTEN     -35 /* Error occured creating listen socket */
#define     API_TIMEOUT    -36 	/* Timed out waiting for response */
#define     API_SSL_ERR    -37 	/* SSL communication error */
#define     API_SCHEDD_DOWN     -38 /* One or more Schedd is down     */
#define     API_NOT_SUPPORTED   -39 /* The function is not supported  */
#define     API_NOT_ENABLED     -40 /* The function is not enabled    */
#define     API_NO_PERMISSION   -41 /* The operation is not permitted */
#define     API_CANT_READ_FILE  -42 /* Can't read the file            */
#define     API_CANT_WRITE_FILE -43 /* Can't write the file           */
#define     API_SCHEDD_NOT_FENCED -44 /* Schedd not currently fenced  */

#define     API_BAD_NETWORK_USAGE_COUNT -45	/* input network usage must equal to the number of adapter requests. */
#define     API_BAD_NETWORK_USAGE -46
#define     API_BAD_NETWORK_ID    -47	/* does not exist in the cluster */
#define     API_BAD_NETWORK_TYPE  -48	/* switch, multilink, ethernet */


/***********************************************************************
 * Support for Performance Monitor APIs
 **********************************************************************/

#define     LL_INVALID_PTR           -1
#define     LL_INVALID_DAEMON_ID     -2
#define     LL_DAEMON_NOT_CONFIG     -3
#define     LL_HOST_NOT_CONFIG       -4
#define     LL_CANNOT_CONTACT_DAEMON -5
#define     LL_DATA_NOT_RECEIVED     -6
#define     LL_INVALID_FIELD_ID      -7
#define     LL_CONFIG_NOT_FOUND      -8

/***********************************************************************
 * Support for ll_control API.
 **********************************************************************/
#define  LL_CONTROL_OK                 0 /* Command successfully sent to appropriate LL daemon. */
#define  LL_CONTROL_CM_ERR            -2 /* Cannot send command to central manager. */
#define  LL_CONTROL_MASTER_ERR        -3 /* Cannot send command to one of LoadL_master daemons. */
#define  LL_CONTROL_CONFIG_ERR        -4 /* Errors encountered while processing the LL admin/config files. */
#define  LL_CONTROL_XMIT_ERR          -6 /* A data transmission failure occurred. */
#define  LL_CONTROL_AUTH_ERR          -7 /* Calling program does not have LL administrator authority. */
#define  LL_CONTROL_VERSION_ERR      -19 /* An incorrect ll_control version has been specified. */
#define  LL_CONTROL_SYSTEM_ERR       -20 /* A system error occurred. */
#define  LL_CONTROL_MALLOC_ERR       -21 /* Unable to allocate memory. */
#define  LL_CONTROL_INVALID_OP_ERR   -22 /* An invalid control_op operation has been specified. */
#define  LL_CONTROL_JOB_LIST_ERR     -23 /* job_list argument contains one or more errors. */
#define  LL_CONTROL_HOST_LIST_ERR    -24 /* host_list argument contains one or more errors. */
#define  LL_CONTROL_USER_LIST_ERR    -25 /* user_list argument contains one or more errors. */
#define  LL_CONTROL_HOLD_ERR         -26 /* Incompatible arguments specified for HOLD operation. */
#define  LL_CONTROL_PRIO_ERR         -27 /* Incompatible arguments specified for PRIORITY operation. */
#define  LL_CONTROL_FAVORJOB_ERR     -28 /* Incompatible arguments specified for FAVORJOB operation. */
#define  LL_CONTROL_FAVORUSER_ERR    -29 /* Incompatible arguments specified for FAVORUSER operation. */
#define  LL_CONTROL_SYS_ERR          -30 /* An error occurred while trying to start a child process. */
#define  LL_CONTROL_START_ERR        -31 /* An error occurred while trying to start the LoadL_master daemon. */
#define  LL_CONTROL_UNUSED_1_ERR     -32
#define  LL_CONTROL_CLASS_ERR        -33 /* class_list argument contains incompatible information. */
#define  LL_CONTROL_TMP_ERR          -34 /* Unable to create a file in /tmp directory. */
#define  LL_CONTROL_ERR              -35 /* Miscellaneous incompatible input specifications. */
#define  LL_CONTROL_NO_DCE_ID        -36 /* DCE identity can not be determined. */
#define  LL_CONTROL_NO_DCE_CRED      -37 /* No DCE credentials. */
#define  LL_CONTROL_INSUFFICIENT_DCE_CRED -38 /* DCE credentials within 300 secs of expiration. */
#define  LL_CONTROL_64BIT_DCE_ERR    -39 /* 64-bit API not supported when DCE is enabled */
#define  LL_CONTROL_NO_DCE_SUPPORT_ERR   -40  /* This version of LL does not support DCE security. */
#define  LL_CONTROL_NO_CTSEC_SUPPORT_ERR -41  /* This version of LL does not support CTSEC security. */
#define  LL_CONTROL_NO_GANG_SUPPORT_ERR  -42  /* This version of LL does not support GANG scheduling. */
#define  LL_CONTROL_NO_PVM_SUPPORT_ERR   -43  /* This version of LL does not support PVM. */
#define  LL_CONTROL_NO_NQS_SUPPORT_ERR   -44  /* This version of LL does not support NQS. */

static const int flush_ckpt_failure = 0xfcbad;

/***********************************************************************
 * Support for ll_modify API.
 **********************************************************************/
enum LL_modify_op {
	UNUSED_EXECUTION_FACTOR,/* reserved for future use */
	CONSUMABLE_CPUS,        /* use int * for data      */
	CONSUMABLE_MEMORY,      /* use int64 * for data    */
	WCLIMIT_ADD_MIN,        /* use int * for data      */
	JOB_CLASS,              /* use char * for data     */
	ACCOUNT_NO,             /* use char * for data     */
	STEP_PREEMPTABLE,       /* use int * for data      */
	SYSPRIO,                /* use int * for data      */
	BG_SIZE,                /* use int * for data      */
	BG_SHAPE,               /* use char * for data     */
	BG_CONNECTION,          /* use int * for data      */
	BG_PARTITION,           /* use char * for data     */
	BG_ROTATE,              /* use int * for data      */
	BG_REQUIREMENTS,        /* use char * for data     */
	RESOURCES,              /* use char * for data     */
	NODE_RESOURCES,         /* use char * for data     */
	CLUSTER_OPTION,         /* use char * for data     */
	DSTG_RESOURCES,         /* use char * for data     */
	BG_PARTITION_TYPE,      /* use int * for data	   */
	BG_USER_LIST,           /* use char * for data     */
	STARTDATE,              /* use char * for data     */
	WALL_CLOCK_LIMIT,       /* use char * for data     */
	JOB_TRACE,              /* use int * for data     */
	DYNAMIC_TRACE,          /* use array int * for data     */
	MAX_MODIFY_OP
} ;

typedef struct {
	enum LL_modify_op type;
	void *data;
}
LL_modify_param;

#define MODIFY_SUCCESS	             0
#define MODIFY_INVALID_PARAM        -1 	/* Invalid param specified */
#define MODIFY_CONFIG_ERROR         -2 	/* Configuration error     */
#define MODIFY_NOT_IDLE             -3 	/* Joblist has non-idle step */
#define MODIFY_WRONG_STATE          -4 	/* Joblist has step in wrong state */
#define MODIFY_NOT_AUTH	            -5 	/* Caller not authorized     */
#define MODIFY_SYSTEM_ERROR         -6 	/* Internal system error     */
#define MODIFY_CANT_TRANSMIT        -7 	/* Communication error       */
#define MODIFY_CANT_CONNECT         -8 	/* Connection error          */
#define MODIFY_NO_DCE_SUPPORT_ERR   -9  /* No support for DCE security */
#define MODIFY_NO_CTSEC_SUPPORT_ERR -10 /* No support for CTSEC security */
#define MODIFY_NO_GANG_SUPPORT_ERR  -11 /* No support for GANG scheduling */
#define MODIFY_NO_PVM_SUPPORT_ERR   -12 /* No support for PVM */
#define MODIFY_NO_NQS_SUPPORT_ERR   -13 /* No support for NQS */
#define MODIFY_OVERLAP_RESERVATION  -14 /* Would overlap with reservation */
#define MODIFY_BAD_BG_SHAPE         -15 /* BlueGene partition shape bad */
#define MODIFY_WRONG_JOB_TYPE       -16 /* Modify operation not valid for job type*/
#define MODIFY_BAD_BG_SIZE          -17 /* BlueGene partition size request not positive */
#define MODIFY_BAD_BG_CONNECTION    -18 /* BlueGene connection request not recognized */
#define MODIFY_EMPTY_BG_PARTITION   -19 /* BlueGene requested partition name blank */
#define MODIFY_CANT_MODIFY          -20 /* BlueGene partition exists, cant modify this value; or cannot modify this attribute for a scale across step */
#define MODIFY_BAD_RESOURCES        -28 /* Syntax error in resources or node_resources keyword */
#define MODIFY_INFO                 -29 /* Informational messages */
#define MODIFY_CLUSTER_OPTION_ERROR -30 /* Keyword cluster_option modification error */


/***********************************************************************
 * Support for ll_run_scheduler API.
 **********************************************************************/
#define RUN_SCHEDULER_SUCCESS	0
#define RUN_SCHEDULER_INVALID_PARAM -1		/* Invalid param specified */
#define RUN_SCHEDULER_CONFIG_ERROR  -2		/* Configuration error     */
#define RUN_SCHEDULER_NOT_AUTH	    -3		/* Caller not authorized     */
#define RUN_SCHEDULER_SYSTEM_ERROR  -4		/* Internal system error     */
#define RUN_SCHEDULER_CANT_TRANSMIT -5		/* Communication error       */
#define RUN_SCHEDULER_CANT_CONNECT  -6		/* Connection error          */
#define RUN_SCHEDULER_NEGOTIATOR_INTERVAL_NON_ZERO  -7	/* Non-zero negotiator interval */

/***********************************************************************
 * Support for ll_cluster API.
 **********************************************************************/
typedef enum LL_cluster_op {
	CLUSTER_SET,		/* Set the multicluster environment to cluster_list */
	CLUSTER_UNSET		/* Unset the multicluster environment */
} ClusterOp_t;

typedef struct {
	ClusterOp_t action;	/* CLUSTER_SET or CLUSTER_UNSET */
	char **cluster_list;	/* NULL terminated list of cluster names */
}
LL_cluster_param;

#define CLUSTER_SUCCESS 	0  /* Success */
#define CLUSTER_SYSTEM_ERROR 	-1  /* System error */
#define CLUSTER_INVALID_CLUSTER_PARAM	-2 /* cluster_list param is not valid */
#define CLUSTER_INVALID_ACTION_PARAM	-3 /* action param is not valid */

/***********************************************************************
 * Support for Preemption
 **********************************************************************/
typedef enum LL_preempt_op {
	PREEMPT_STEP,
	RESUME_STEP,
	SYSTEM_PREEMPT_STEP
} PreemptOp_t;

/* Values for the preempt method */
typedef enum preempt_method {
	LL_PREEMPT_SUSPEND,
	LL_PREEMPT_VACATE,
	LL_PREEMPT_REMOVE,
	LL_PREEMPT_SYS_HOLD,
	LL_PREEMPT_USER_HOLD
} PreemptMethod_t;

typedef struct {
	PreemptOp_t type;
	PreemptMethod_t method;
	char **user_list;
	char **host_list;
	char **job_list;
}
LL_preempt_param;

typedef struct {
	char *cluster_name;
	char *job_id;
}
LL_move_job_param;


/***********************************************************************
 * Support for poe API's(ll_spawn_connect,ll_spawn_connect_ext,ll_spawn_read and ll_spawn_write)
 **********************************************************************/
enum LL_JobManagement_RC {
	JOBMGMT_IO_COMPLETE=1, 		/* LoadLeveler I/O is  completed. */
	JOBMGMT_IO_PENDING=0, 		/* LoadLeveler I/O is pending. */
	JOBMGMT_BAD_JOBMGMT_OBJECT=-1, 	/* JobMgmtObject specified is not valid. */
	JOBMGMT_FAILED_CONNECT=-3, 	/* Can not connect to LoadLeveler Daemon. */
	JOBMGMT_SYSTEM=-5,			/* System Error. */
	JOBMGMT_NULL_EXECUTABLE=-6,	/* NULL executable specified. */
	JOBMGMT_TASKMGR_RUNNING=-7,	/* Parallel Task Manager for this job Step is alreay running on the targeted node. */
	JOBMGMT_INCOMPATABLE_NODES=-8,	/* All the nodes targeted to run the parallel job cannot support this interface. */
	JOBMGMT_BAD_MACHINE_OBJECT=-9,	/* Invalid Machine Object. */
	JOBMGMT_BAD_STEP_OBJECT=-10,	/* Invalid Step Object. */
	JOBMGMT_BAD_SEQUENCE=-11,		/* Function is called out of order. */
	JOBMGMT_BAD_FD=-12,			/* Invalid Socket descriptor. */
	JOBMGMT_JOB_PREEMPTED=-13           /* Starter is still in PREEMPTED state, job resume will be deferred. */
};

typedef struct {
	char *stepid;                      /*Const string of step id representing the step that the task belongs to. This step has to have been previously submitted via the ll_request function.*/
	char *machine_name;                /*Const string of machine name representing the machine where task Manager to be started. This argument is used to identify which OSI interface to connect to.*/
	char *executable;          /*character string for the name of executable to be started, which will serve as the parallel task manager.*/
} LL_spawn_connect_param;

/***********************************************************************
 * Support for Advance Reservation
 **********************************************************************/
typedef enum reservation_state_code {
	RESERVATION_WAITING,
	RESERVATION_SETUP,
	RESERVATION_ACTIVE,
	RESERVATION_ACTIVE_SHARED,
	RESERVATION_CANCEL,
	RESERVATION_COMPLETE
} Reservation_State_t;

typedef enum LL_reservation_mode {
	RESERVATION_DEFAULT_MODE = 0,
	RESERVATION_SHARED = (1<<0),
	RESERVATION_REMOVE_ON_IDLE = (1<<1),
	RESERVATION_BIND_FIRM = (1<<2),
	RESERVATION_BIND_SOFT = (1<<3)
} Reservation_Mode_t;

/* enums for Advance Reservation */

/* LL_reservation data identifies the type of data */
enum LL_reservation_data {
	RESERVATION_START_TIME,      /* char *         */
	RESERVATION_ADD_START_TIME,  /* int *          */
	RESERVATION_DURATION,        /* int *          */
	RESERVATION_ADD_DURATION,    /* int *          */
	RESERVATION_BY_NODE,         /* int *          */
	RESERVATION_ADD_NUM_NODE,    /* int *          */
	RESERVATION_BY_HOSTLIST,     /* char **, NULL terminated  */
	RESERVATION_ADD_HOSTS,       /* char **, NULL terminated  */
	RESERVATION_DEL_HOSTS,       /* char **, NULL terminated  */
	RESERVATION_BY_JOBSTEP,      /* char *      */
	RESERVATION_BY_JCF,          /* char *      */
	RESERVATION_USERLIST,        /* char **, NULL terminated  */
	RESERVATION_ADD_USERS,       /* char **, NULL terminated  */
	RESERVATION_DEL_USERS,       /* char **, NULL terminated  */
	RESERVATION_GROUPLIST,       /* char **, NULL terminated  */
	RESERVATION_ADD_GROUPS,      /* char **, NULL terminated  */
	RESERVATION_DEL_GROUPS,      /* char **, NULL terminated  */
	RESERVATION_MODE_SHARED,     /* int *; *data = 0 : Not Shared; *data = 1: Share  */
	RESERVATION_MODE_REMOVE_ON_IDLE,  /* int *; *data = 0 : Don't Remove; *data = 1 : Remove */
	RESERVATION_OWNER,           /* char *      */
	RESERVATION_GROUP,           /* char *      */
	RESERVATION_BY_BG_CNODE,     /* int *       */
	RESERVATION_BY_HOSTFILE,     /* char *      */
	RESERVATION_BINDING_METHOD,  /* int* 	*/
	RESERVATION_EXPIRATION,      /* char* 	*/
	RESERVATION_RECURRENCE,	 /* LL_contab_time* 	*/
	RESERVATION_OCCURRENCE,	 /* int* 	*/
	RESERVATION_NOTIFICATION_PROGRAM,      /* char *                    */
	RESERVATION_NOTIFICATION_PROGRAM_ARGS, /* char *                    */
	RESERVATION_FLOATING_RESOURCE          /* char **, NULL terminated  */
};

/* structure used by ll_change_reservation */
typedef struct {
	enum LL_reservation_data type;
	void *data;
}
LL_reservation_change_param;

/* structure for crontab time */
typedef struct {
	int *minutes; 			/* 0-59 */
	int *hours; 			/* 0-23 */
	int *dom; 			/* days of the month: 1-31 */
	int *months; 			/* 1 (January) through 12 (December) */
	int *dow; 			/* days of the week: 0 (Sunday) through 6 (Saturday) */
}
LL_crontab_time;

/* LL_reservation_type identifies the type of reservation */
enum LL_reservation_type {
	RESERVATION_TYPE_DEFAULT,
	RESERVATION_TYPE_FLEXIBLE
};
/* structure used by ll_make_reservation */
typedef struct {
	char **ID;                      /* -> to output string reservation id */
	char *start_time;               /* [mm/dd[/[yy]yy] ]HH:MM format start time */
	int duration;                   /* length of reservation in minutes   */
	enum LL_reservation_data data_type; /* how nodes should be reserved   */
	void *data;                     /* -> to data specifying the nodes    */
	int mode;                       /* shared/remove_on_idle one/both/neither*/
	char **users;                   /* NULL terminated array of user ids  */
	char **groups;                  /* NULL terminated array of LL groups */
	char *group;                    /* group which owns the reservation   */
	int unused1;										/* placeholder 												*/
	char *expiration;		/* recurring reservation expiration date */
	LL_crontab_time *recurrence;	/* crontab time structor */
	enum LL_reservation_type reservation_type_requested; /* type of reservation requested */
	char *notification_program; 	/* complete path for notification program */
	char *notification_program_args; /* argument to be passed to the notification program */
	char *floating_resources; 	/* floating consumable resources */
	int down_nodes_flag;		/* flag of reserving down nodes */
}
LL_reservation_param;

/* structure for binding jobsteps to a reservation */
typedef struct {
	char **jobsteplist;		/* host.jobid.stepid, null terminated */
	char *ID;			/* -> reservation id, NULL for unbind */
	int unbind;			/* TRUE = unbind, FALSE to bind       */
	int binding_method;		/* Default binding method for this reservation */
}
LL_bind_param;

typedef struct {
	char **IDs;			/* NULL-terminated array of reservation identifiers */
	char **user_list;		/* NULL-terminated array of user IDs that own reservation */
	char **group_list;	/* NULL-terminated array of LoadLeveler groups that own reservation */
	char **host_list;		/* NULL-terminated array of machine names. One member of all indicates all nodes in the cluster */
	char **base_partition_list;	/* NULL-terminated array of machine names. One member of all indicates all Blue Gene base partition. */
	char *begin;			/* The beginning of an interval for partial cancellation of occurrences of a recurring reservation */
	char *end;			/* The end of an interval for partial cancellation of occurrences of a recurring reservation */
} LL_remove_reservation_param;

/***********************************************************************
 * Status codes to support Advance Reservation
 **********************************************************************/
#define RESERVATION_OK                      0 /* Success */
#define RESERVATION_LIMIT_EXCEEDED         -1 /* Exceeds max # of reservations allowed in the LoadLeveler cluster */
#define RESERVATION_TOO_CLOSE              -2 /* Reservation is being made within the minimum advance time */
#define RESERVATION_NO_STORAGE             -3 /* The system cannot allocate memory  */
#define RESERVATION_CONFIG_ERR             -4 /* Errors encountered processing the LL admin/config files */
#define RESERVATION_CANT_TRANSMIT          -5 /* A data transmission failure occurred */
#define RESERVATION_GROUP_LIMIT_EXCEEDED   -6 /* Exceeds max # of reservations for the group */
#define RESERVATION_USER_LIMIT_EXCEEDED    -7 /* Exceeds max # of reservations for the user */
#define RESERVATION_SCHEDD_CANT_CONNECT    -8 /* Schedd cannot connect to CM */
#define RESERVATION_API_CANT_CONNECT       -9 /* API cannot connect to the Schedd or CM */
#define RESERVATION_JOB_SUBMIT_FAILED      -10 /* Submit of job command file failed */
#define RESERVATION_NO_MACHINE             -11 /* One or more machines in the hostlist are not in the LoadLeveler cluster */
#define RESERVATION_WRONG_MACHINE          -12 /* Reservations not permitted on one or more machines in the hostlist  */
#define RESERVATION_NO_RESOURCE            -13 /* Insufficient resources in the LoadLeveler cluster */
#define RESERVATION_NOT_SUPPORTED          -14 /* The scheduler in use does not support Advance Reservation  */
#define RESERVATION_NO_JOBSTEP             -15 /* The jobstep used for node selection doesn't exist */
#define RESERVATION_WRONG_JOBSTEP          -16 /* The jobstep used for node selection isn't in the right state */
#define RESERVATION_NOT_EXIST              -17 /* The reservation does not exist */
#define RESERVATION_REQUEST_DATA_NOT_VALID -18 /* Invalid input */
#define RESERVATION_NO_PERMISSION          -19 /* Permission cannot be granted */
#define RESERVATION_TOO_LONG               -20 /* Duration exceeds max allowed */
#define RESERVATION_WRONG_STATE            -21 /* Reservation state is not right for the requested operation */
#define RESERVATION_NO_DCE_CRED            -30 /* DCE is enabled, the user has no credentials */
#define RESERVATION_INSUFFICIENT_DCE_CRED  -31 /* DCE is enabled, credential lifetime < 5 minutes */
#define RESERVATION_COSCHEDULE_NOT_ALLOWED -32 /* Coschedule jobstep is not allowed for node selection */
#define RESERVATION_HOSTFILE_ERR           -33 /* Error encountered when processing the host file */
#define RESERVATION_WRONG_BINDINGMETHOD           -34 /*When customer submit make reservation requeset with wrong binding method. Such as not firm or soft method. */
#define RESERVATION_OCCURRENCE_ID_NOT_ALLOWED           -35 /*The occurrence ID can not be specified by llchres */
#define RESERVATION_IMMEDIATELY_EXPIRE_NOT_ALLOWED           -36 /* The expiration of a recurring reservation cannot be modified in a way that would cause the reservation to immediately expire.*/
#define RESERVATION_EXPIRE_TOO_LONG        -37 /* The requested expiration exceeds the max_reservation_expiration */
#define RESERVATION_VS_ERR                 -38 /* Errors encountered processing VS */
#define RESERVATION_OCCURRENCE_OVERLAP     -39 /* Recurring reservation's occurrence overlap with itself */
#define RESERVATION_RECURRING_SOFT_NOT_ALLOWED  -40 /* Recurring step can not use to make reservation */
#define RESERVATION_SCALE_ACROSS_NOT_ALLOWED -41 /* Scale-across jobstep is not allowed for making or changing a reservation */
#define RESERVATION_INCOMPATIBLE_STARTTIME_RECURRENCE -42 /* Incompatible config for start_time and recurrence*/
#define RESERVATION_FLEXIBLE_DSTG_NOT_ALLOWED -43 /* Data staging flexible job is not allowed for making or changing a reservation */
/***********************************************************************
 * Support for Blue Gene
 **********************************************************************/
typedef enum bg_bp_state_t {
	BG_BP_UP,
	BG_BP_DOWN,
	BG_BP_MISSING,
	BG_BP_ERROR,
	BG_BP_SOME_DOWN=BG_BP_ERROR,
	BG_BP_NAV
} BgBPState_t;

typedef enum bg_partition_state_t {
	BG_PARTITION_FREE,
	BG_PARTITION_CONFIGURING,
	BG_PARTITION_READY,
	BG_PARTITION_BUSY,
	BG_PARTITION_DEALLOCATING,
	BG_PARTITION_ERROR,
	BG_PARTITION_NAV,
	BG_PARTITION_UNINITIALIZED
} BgPartitionState_t;

typedef enum bg_connection_t {
	MESH  = 0,
	TORUS = 1,
	BG_NAV,
	PREFER_TORUS
} BgConnection_t;

typedef enum bg_node_mode_t {
	COPROCESSOR,
	VIRTUAL_NODE
} BgNodeMode_t;

typedef enum bg_partition_type_t {
	HPC,
	HTC_SMP,
	HTC_DUAL,
	HTC_VN,
	HTC_LINUX_SMP,
	PTYPE_NAV
} BgPartitionType_t;

typedef enum bg_port_t {
	BG_PORT_PLUS_X,
	BG_PORT_MINUS_X,
	BG_PORT_PLUS_Y,
	BG_PORT_MINUS_Y,
	BG_PORT_PLUS_Z,
	BG_PORT_MINUS_Z,
	BG_PORT_S0,
	BG_PORT_S1,
	BG_PORT_S2,
	BG_PORT_S3,
	BG_PORT_S4,
	BG_PORT_S5,
	BG_PORT_NAV
} BgPort_t;

typedef enum bg_switch_state_t {
	BG_SWITCH_UP,
	BG_SWITCH_DOWN,
	BG_SWITCH_MISSING,
	BG_SWITCH_ERROR,
	BG_SWITCH_NAV
} BgSwitchState_t;

typedef enum bg_switch_dimension_t {
	BG_DIM_X = 0,
	BG_DIM_Y = 1,
	BG_DIM_Z = 2,
	BG_DIM_NAV
} BgSwitchDimension_t;

typedef enum bg_wire_state_t {
	BG_WIRE_UP,
	BG_WIRE_DOWN,
	BG_WIRE_MISSING,
	BG_WIRE_ERROR,
	BG_WIRE_NAV
} BgWireState_t;

typedef enum bg_node_card_state_t {
	BG_NODE_CARD_UP,
	BG_NODE_CARD_DOWN,
	BG_NODE_CARD_MISSING,
	BG_NODE_CARD_ERROR,
	BG_NODE_CARD_NAV
} BgNodeCardState_t;

typedef enum bg_quarter_t {
	BG_QUARTER_Q1 = 0,
	BG_QUARTER_Q2 = 1,
	BG_QUARTER_Q3 = 2,
	BG_QUARTER_Q4 = 3,
	BG_QUARTER_Q_NAV
} BgQuarter_t;

typedef enum bg_job_state_t {
	BG_JOB_IDLE,
	BG_JOB_STARTING,
	BG_JOB_RUNNING,
	BG_JOB_TERMINATED,
	BG_JOB_KILLED,
	BG_JOB_ERROR,
	BG_JOB_DYING,
	BG_JOB_DEBUG,
	BG_JOB_LOAD,
	BG_JOB_LOADED,
	BG_JOB_BEGIN,
	BG_JOB_ATTACH,
	BG_JOB_NAV
} BgJobState_t;

typedef enum bg_bp_computenode_memory_t {
	BG_BP_COMPUTENODE_MEMORY_256M,
	BG_BP_COMPUTENODE_MEMORY_512M,
	BG_BP_COMPUTENODE_MEMORY_1G,
	BG_BP_COMPUTENODE_MEMORY_2G,
	BG_BP_COMPUTENODE_MEMORY_4G,
	BG_BP_COMPUTENODE_MEMORY_NAV
} BgBPComputenodeMemory_t;

/***********************************************************************
 * Support for ll_fair_share API.
 **********************************************************************/
enum {FAIR_SHARE_RESET, FAIR_SHARE_SAVE};

typedef struct {
	int operation;   /* either FAIR_SHARE_RESET or FAIR_SHARE_SAVE             */
	char *savedir;   /* directory to save data for FAIR_SHARE_SAVE option     */
	char *savedfile; /* optional saved file to restore to for FAIR_SHARE_RESET */
}
LL_fair_share_param;

/***********************************************************************
 * Support for Affinity
 **********************************************************************/
enum AffinitySupport { MCM_AFFINITY, USER_DEFINED_RSET, NO_AFFINITY };

/***********************************************************************
 * Support for SMT
 **********************************************************************/
typedef enum SmtStateType { SMT_DISABLED, SMT_ENABLED, SMT_NOT_SUPPORT, SMT_SMT2 } SmtStateType_t;
typedef enum SMTRequiredState { SMT_OFF, SMT_ON, SMT_AS_IS } SMTRequiredState_t;

/***********************************************************************
 * Support for ll_cluster_auth
 **********************************************************************/
typedef enum LL_cluster_auth_op {
	CLUSTER_AUTH_GENKEY
} ClusterAuthOp_t;

typedef struct {
	ClusterAuthOp_t type;
}
LL_cluster_auth_param;


/***********************************************************************
 * Support for ll_move_spool
 **********************************************************************/
/* Values for the members of the spool directory to transfer */
typedef enum move_spool_data {
	LL_MOVE_SPOOL_JOBS
} SpoolData_t;

typedef struct {
	char **schedd_host;
	char *spool_directory;
	SpoolData_t data;
}
LL_move_spool_param;

/***********************************************************************
 * Support for llconfig API
 **********************************************************************/
enum LL_config_op { CONFIG_I, CONFIG_C, CONFIG_D, CONFIG_STANZA_C, CONFIG_STANZA_D, CONFIG_STANZA_A,CONFIG_STANZA_R,CONFIG_STANZA_F};

#define LLCONFIG_INITIALIZE_OK 0
#define LLCONFIG_CHANGE_OK 0
#define LLCONFIG_DISPLAY_OK 0
#define LLCONFIG_DB_CANNOT_CONNECT -1
#define LLCONFIG_INITIALIZE_RUN_SQL_ERROR -2
#define LLCONFIG_CHECK_FAILED -3		/* check utility found error configures */
#define LLCONFIG_INITIALIZE_PUT_TO_DB_ERROR -4
#define LLCONFIG_INITIALIZE_CM_NEEDED -5 // tbd
#define LLCONFIG_INITIALIZE_SCHEDD_NEEDED -6 // tbd
#define LLCONFIG_CHANGE_AUTH_FAILED -7
#define LLCONFIG_CHANGE_RECONFIG_ERROR -8
#define LLCONFIG_NOT_VALID_INPUT -9
#define LLCONFIG_LOAD_ODBC_LIBRARY_FAILED -10
#define LLCONFIG_LOADLDB_NOT_CONFIG -11
#define LLCONFIG_CONFIG_OPT_UNSUPPORTED -12
#define LLCONFIG_INITIALIZE_PARSE_FILE_FAILED -13 // global config or admin file parse failed
#define LLCONFIG_NO_UPDATE_NEEDED -14

typedef enum Stanza_op {
	Stanza_All,
	Stanza_One,
	Stanza_Sub,
} Stanza_op_t;

/***********************************************************************
 * Function Declaration statements
 **********************************************************************/
#ifdef __cplusplus
extern "C" {
#endif
	void llfree_job_info(LL_job *, int);
	int llinit(int);
	int llinitiate(LL_job *, int);
	int llwait(LL_job **, LL_job_step **, int);
	int ll_start_job_ext(LL_start_job_info_ext *ptr);
	int ll_terminate_job(LL_terminate_job_info *);
	int llsubmit(char *, char *, char *, LL_job *, int);
	int GetHistory(char*, int (*)(LL_job *), int);

	LL_element *llpd_allocate(void);
	int ll_update(LL_element*,enum LL_Daemon);
	int ll_fetch(LL_element*,enum LLAPI_Specification, void *);

	LL_element *ll_query(enum QueryType);
	int ll_deallocate(LL_element *);
	int ll_set_request(LL_element *,enum QueryFlags,char **,enum DataFilter);
	int ll_reset_request(LL_element *);
	LL_element *ll_get_objs(LL_element *,enum LL_Daemon,char *,int *,int *);
	LL_element *ll_next_obj(LL_element *);
	int ll_free_objs(LL_element *);

	int ll_init_job(LL_element **);
	void ll_deallocate_job(LL_element *);
	int ll_parse_string(LL_element *,char *,LL_element **,int,char *,LL_element **);
	int ll_parse_file(LL_element *,char *,LL_element **,int,char *,LL_element **);
	int ll_parse_verify(LL_element *,LL_element *,LL_element **);
	int ll_request(LL_element *,LL_element *);
	int ll_spawn(LL_element *,LL_element *,LL_element *,char *);
	enum EventType ll_event(LL_element *,int,LL_element **,LL_element *);
	int ll_get_job(LL_element *,LL_element **);
	int ll_close(LL_element *);
	int ll_get_data(LL_element *,enum LLAPI_Specification, void *);
	int ll_set_data(LL_element *,enum LLAPI_Specification, void *);
	const char* ll_version(void);
	int ll_control(int, enum LL_control_op, char **, char **, char **, char **, int);
	int ll_spawn_task(LL_element *,LL_element *,char *, LL_element *, int flags);
	int ll_spawn_mpich_task(char *, char *,char *, int);
	int ll_spawn_mpich_error(char *);
	void free_crontab(LL_crontab_time * crontab);
	void ll_pe_rm_close_nullFP(void);

#if !defined(__linux__) || defined(LL_LINUX_CKPT) || defined(METACLUSTER_CKPT)
	int ll_ckpt(LL_ckpt_info *);
#endif /* !__linux__ || LL_LINUX_CKPT || METACLUSTER_CKPT */

	int ll_modify(int, LL_element **, LL_modify_param **, char **);
	int ll_run_scheduler(int, LL_element **);
	int ll_preempt(int, LL_element **, char *, enum LL_preempt_op);
	int ll_preempt_jobs(int, LL_element **, LL_preempt_param **);
	int ll_move_job(int, LL_element **, LL_move_job_param **);
	int ll_cluster(int, LL_element **, LL_cluster_param *);

#if !defined(__linux__) || defined(LL_LINUX_CKPT)
	enum CkptStart ll_metacluster_ckpt_start(time_t *, LL_element *, char **);
	time_t ll_metacluster_ckpt_complete(int, time_t, int, char *, LL_element *, char **);
	enum CkptStart ll_local_ckpt_start(time_t *);
	time_t ll_local_ckpt_complete(int, time_t, int, char *);
	void ckpt();			/* in support of old checkpoint code */
	int ll_init_ckpt(LL_ckpt_info *);
	time_t ll_ckpt_complete(LL_element *, int, cr_error_t *, time_t, int);
	int ll_set_ckpt_callbacks(callbacks_t *);
	int ll_unset_ckpt_callbacks(int);
#endif /* !__linux__ || LL_LINUX_CKPT */
	int ll_task_inst_pid_update(int*, int);

	char *ll_error(LL_element **, int);
	int ll_make_reservation(int, LL_element **, LL_reservation_param **);
	int ll_change_reservation(int, LL_element **, char **, LL_reservation_change_param **);
	int ll_bind(int, LL_element **, LL_bind_param **);
	int ll_remove_reservation(int, LL_element **, char **, char **, char **, char**,char**);
	int ll_init_reservation_param(int, LL_element **, LL_reservation_param **);
	int ll_spawn_connect(int, LL_element *,LL_element *,LL_element *,char *, LL_element **);
	int ll_spawn_connect_ext(int, LL_element **, LL_spawn_connect_param *, LL_element **);
	int ll_spawn_write(int, int, LL_element *,LL_element **);
	int ll_spawn_read(int, int, LL_element *,LL_element **);
	int ll_fair_share(int, LL_element **, LL_fair_share_param *);
	int ll_config_changed(void);
	int ll_read_config(LL_element **);

	int ll_remove_reservation_xtnd(int version, LL_element** errorObj, LL_remove_reservation_param* param);
	int ll_init_remove_reservation_param(int version, LL_remove_reservation_param *param);

	int ll_cluster_auth(int, LL_element **, LL_cluster_auth_param **);
	int ll_move_spool(int, LL_move_spool_param **, int (*func)(char *, int, LL_element **),LL_element **);
	int ll_config(int version, enum LL_config_op config_op, char **input, char ***output, char **hostlist, LL_element **err_obj, int no_reconfig);
#ifdef __cplusplus
}
#endif /* __cpluplus */

#endif /* _llapi_h__ */
