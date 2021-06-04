#ifdef AIX


/*============================================================================
  prints the wallclock time used by steps
   
   Usage "wall_clock_time stepid1 stepid2 stepid3"
============================================================================ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "llapi.h"

void get_wall_clock_used(char *stepid,char **hostname,int *wu,int *wh,int *ws);

#ifdef DEBUG
int  main(int argc, char **argv)
{	
	char *hostname;
	char *stepid;
	int wused, whard, wsoft;
	int rc, i, step_count;

	step_count = argc - 1;
	printf("%25s %15s \n"," ","WallClockTime");
	printf("%-15s %-20s %-5s %-5s %-5s \n","Hostname","StepId", "Used", "Hard", "Soft");

	/* loop through the steplist  to find wallclock time 
	   used for each step */

	for ( i = 0; i < step_count; i++ ) {
		stepid = argv[i+1];
		hostname = NULL;
		get_wall_clock_used(stepid, &hostname, &wused, &whard, &wsoft);
		printf ("%-15s %-20s %-5d %-5d %-5d \n", strtok(hostname,"."), stepid, wused, whard, wsoft);
		if (hostname != NULL) {
			free(hostname);
		}
	}	
}
#endif	

/* 
  This function finds the actual wall clock time used by a running step
  excluding the time spent while checkpointing or preempted by 
	suspend.

		params          type         descr
    stepid         char*     step id of the step whose master node to be found
    hostname       char**    on return on hostname will have the hostname for
                              machine where master task for this step runs
	  wu             int*      on return wallclock used
	  wh             int*      on return wallclock  hard limit
	  ws             int*      on return wallclock  soft limit

*/
	
 void get_wall_clock_used(char *stepid, char **hostname, int *wu, int *wh, int *ws )
{
	LL_element *queryObject; 
	LL_element *job = NULL; 
	LL_element *step = NULL;
	LL_element *machine = NULL;
	int rc, err_code ;
	int  obj_count ;
	char *steplist[2];
	steplist[0] = stepid;
	steplist[1] = NULL;

	*wu = 0;
	*wh = 0; 
	*ws = 0;
	
	/* Initialize the query: Job query  */
	queryObject = ll_query(JOBS);
	if (queryObject == NULL) { 
		printf("Query JOBS: ll_query() returns NULL.\n");
		exit(1);
	}

	/* Request information of job steps. */
	rc = ll_set_request(queryObject, QUERY_STEPID, steplist, ALL_DATA);
	if (rc != 0) {
		printf("ll_set_request() - QUERY_STEPID -  RC = %d\n", rc);
		exit(1);
	}

	/* Get the requested job objects from the central manager.      */
	/* NOTE: if running the API scheduler then the job query object */
	/*       must be requested from the LoadL_startd (LL_STARTD) on */
	/*       the machine running the master process                 */
	job = ll_get_objs(queryObject, LL_CM, NULL, &obj_count, &err_code);
	if (job == NULL) {
		*hostname = (char *)strdup("Not Found");
		return;
	}  

	ll_get_data(job, LL_JobGetFirstStep, &step);
	if (step != NULL) {
		ll_get_data(step, LL_StepGetFirstMachine, &machine);
		if ( machine ) {
			ll_get_data(machine, LL_MachineName, hostname);
			ll_get_data(step, LL_StepWallClockUsed, wu);
			ll_get_data(step, LL_StepWallClockLimitHard, wh);
			ll_get_data(step, LL_StepWallClockLimitSoft, ws);
		} else { 
			*hostname = (char *)strdup("Idle");
		}
	} else {
		*hostname = (char *)strdup("Not Found");
		return;
	}
	ll_free_objs(queryObject);
	ll_deallocate(queryObject);
}
#else
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
void get_wall_clock_used(char *stepid,char **hostname,int *wu,int *wh,int *ws)
{
  time_t curtime=time(NULL);
  char *JobStartTime=getenv("JobStartTime");
  char *JobTimeLimit=getenv("JobTimeLimit");
  int StartTime, TimeLimit;
  StartTime=curtime;
  TimeLimit=1800;
  if( JobStartTime != NULL )
    {
    sscanf(JobStartTime,"%d",&StartTime);
    }

  if( JobTimeLimit != NULL )
    {
    sscanf(JobTimeLimit,"%d",&TimeLimit);
    }
  *wh = *ws = TimeLimit;
  *wu = curtime - StartTime;
  return;
}
#ifdef DEBUG
int  main(int argc, char **argv)
{
        char *hostname;
        char *stepid;
        int wused, whard, wsoft;

        get_wall_clock_used(stepid, &hostname, &wused, &whard, &wsoft);
        fprintf(stderr,"used=%d, hard limit=%d, soft limit=%d\n",wused,whard,wsoft);
}
#endif
#endif

void c_get_my_resident_time(int *wu, int *wh, int *ws )
{
 char *hostname;
 *wu=0 ; *wh=1800 ; *ws=1800 ;
#ifdef AIX
 get_wall_clock_used(getenv("LOADL_STEP_ID"), &hostname, wu, wh, ws);
#endif
}

f_get_my_resident_time(int *wu, int *wh, int *ws )
{
 char *hostname;
 *wu=0 ; *wh=1800 ; *ws=1800 ;
#ifdef AIX
 get_wall_clock_used(getenv("LOADL_STEP_ID"), &hostname, wu, wh, ws);
#endif
}
f_get_my_resident_time_(int *wu, int *wh, int *ws )
{
 char *hostname;
 *wu=0 ; *wh=1800 ; *ws=1800 ;
#ifdef AIX
 get_wall_clock_used(getenv("LOADL_STEP_ID"), &hostname, wu, wh, ws);
#endif
}
f_get_my_resident_time__(int *wu, int *wh, int *ws )
{
 char *hostname;
 *wu=0 ; *wh=1800 ; *ws=1800 ;
#ifdef AIX
 get_wall_clock_used(getenv("LOADL_STEP_ID"), &hostname, wu, wh, ws);
#endif
}
