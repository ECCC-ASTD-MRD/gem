/* RMNLIB - Library of useful routines for C and FORTRAN programming
 * Copyright (C) 1975-2001  Division de Recherche en Prevision Numerique
 *                          Environnement Canada
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <rpnmacros.h>

typedef struct {
   pthread_cond_t condition;
   pthread_mutex_t mutex;
   int value;
} event;

/*
   INTEGER FUNCTION CREATE_LOCK (LOCK)
   INTEGER *8 LOCK
*/
/*! Create a user lock */
wordint f77name(create_lock)(pthread_mutex_t** lockno)
{
   if( ( *lockno = (pthread_mutex_t*) malloc(sizeof(pthread_mutex_t)) ) == NULL ) return(1);
   pthread_mutex_init(*lockno,NULL);      /* initialize  mutex */
   return(0);
}


pthread_mutex_t* c_create_lock()
{
   pthread_mutex_t *lockno;
   if( ( lockno = (pthread_mutex_t *) malloc(sizeof(pthread_mutex_t)) ) == NULL ) return(NULL) ;
   return(lockno);
}

/*
   INTEGER FUNCTION DESTROY_LOCK (LOCK)
   INTEGER *8 LOCK
*/
/*! Destroy a user lock */
wordint f77name(destroy_lock)(pthread_mutex_t** lockno)
{
   pthread_mutex_destroy(*lockno) ;
   free(*lockno);
   *lockno = NULL;
   return(0);
}

int c_destroy_lock(pthread_mutex_t* lockno)
{
   pthread_mutex_destroy(lockno) ;
   free(lockno);
   return( 0 );
}

/*
   INTEGER FUNCTION ACQUIRE_LOCK (LOCK)
   INTEGER *8 LOCK
*/
/*! Acquire a user lock */
wordint f77name(acquire_lock)(pthread_mutex_t** lockno)
{
   return( pthread_mutex_lock(*lockno) );
}

int c_acquire_lock(pthread_mutex_t* lockno)
{
   return( pthread_mutex_lock(lockno) );
}

/*
   INTEGER FUNCTION RELEASE_LOCK (LOCK)
   INTEGER *8 LOCK
*/
/*! Release a user lock */
wordint f77name(release_lock)(pthread_mutex_t** lockno)
{
   return( pthread_mutex_unlock(*lockno) );
}

int c_release_lock(pthread_mutex_t *lockno)
{
   return( pthread_mutex_unlock(lockno) );
}

/*
   INTEGER FUNCTION CREATE_EVENT (EVENT)
   INTEGER *8 EVENT
*/
/*! Create a user event */
wordint f77name(create_event)(event** user_event)
{
   event *the_event;

   if( ( the_event = (event *) malloc(sizeof(event)) ) == NULL )
      return(1);
   *user_event = the_event ;

   pthread_mutex_init(&(the_event->mutex),NULL);
   the_event->value = -1;
   return( pthread_cond_init(&(the_event->condition),NULL) );
}

event* c_create_event()
{
   event* user_event;
   if( f77name(create_event)(&user_event) ) return(NULL);
   /*
   printf("c_create_event: event address=%x, lock=%x,cond=%x\n",user_event,&(user_event->mutex),&(user_event->condition));
   */
   return(user_event);
}

/*
   INTEGER FUNCTION POST_EVENT (EVENT,VALUE) ! sets the value associated with event
   INTEGER *8 EVENT                          ! returns the value associated with event
   INTEGER VALUE
*/
wordint f77name(post_event)(event** user_event, wordint* value)
{
   event* the_event = *user_event;

   /*
   printf("post_event: value=%d,lock=%x,cond=%x\n",*value,&(the_event->mutex),&(the_event->condition));
   */
   pthread_mutex_lock( &(the_event->mutex) );
   the_event->value = *value;
   pthread_cond_broadcast( &(the_event->condition) );
   pthread_mutex_unlock( &(the_event->mutex) );
   return(*value);
}

int c_post_event(event* user_event, int value)
{
   event* the_event = user_event;
   wordint f_value = value;

   /*
   printf("c_post_event: event address=%x, lock=%x,cond=%x\n",the_event,&(the_event->mutex),&(the_event->condition));
   */
   return( f77name(post_event)(&the_event, &f_value) );
}

/*
   INTEGER FUNCTION CHECK_EVENT (EVENT) ! returns the value associated with event
   INTEGER *8 EVENT
*/
wordint f77name(check_event)(event **user_event)
{
   event *the_event = *user_event;

   return( the_event->value );
}

int c_check_event(event *user_event)
{
   event *the_event=user_event;
   return( f77name(check_event)(&the_event) );
}

/*
   INTEGER FUNCTION WAIT_EVENT (EVENT,VALUE) ! waits for the value associated with event
   INTEGER *8 EVENT                          ! to take a specific value
   INTEGER VALUE                             ! returns the value associated with event
*/
wordint f77name(wait_event)(event **user_event, wordint *value)
{
   event *the_event = *user_event;

   /*
   printf("wait_event: value=%d\n",*value);
   */
   pthread_mutex_lock( &(the_event->mutex) );
   while (the_event->value != *value) pthread_cond_wait( &(the_event->condition) , &(the_event->mutex) );
   pthread_mutex_unlock( &(the_event->mutex) );
   return(*value);
}

int c_wait_event(event *user_event, int value)
{
   event *the_event=user_event;
   wordint f_value=value;
   return( f77name(wait_event)(&the_event,&f_value) );
}

/*
   INTEGER FUNCTION DESTROY_EVENT (EVENT)    ! destroys a user event, returns 0
   INTEGER *8 EVENT
*/
wordint f77name(destroy_event)(event **user_event)
{
   event *the_event = *user_event;

   pthread_mutex_destroy( &(the_event->mutex) );
   pthread_cond_destroy( &(the_event->condition) );
   free(*user_event);
   *user_event=NULL;
   return(0);
}

int c_destroy_event(event *the_event)
{
   pthread_mutex_destroy( &(the_event->mutex) );
   pthread_cond_destroy( &(the_event->condition) );
   free(the_event);
   return(0);
}

/*! Ccreate a new thread (task) */
wordint f77name(create_thread)(void *(*function)(void *), void *arg)
{
   pthread_t thread_id;
   if ( pthread_create(&thread_id, NULL, function, arg) ) {
      return(-1);
   } else {
      return(thread_id);
   }
}

int c_create_thread(void *(*function)(void *), void *arg)
{
   return( f77name(create_thread)(function,arg) );
}

/*! Wait for a thread to terminate */
wordint f77name(join_thread)(wordint *id)
{
   pthread_t thread_id = *id;
   return( pthread_join(thread_id, NULL) );
}

int c_join_thread(int id)
{
   wordint f_id = id;
   return( f77name(join_thread)( &f_id ) );
}

/*! Get current thread id */
wordint f77name(id_thread)()
{
   return(pthread_self());
}

/*! Get current thread id */
int c_id_thread()
{
   return(pthread_self());
}

/*! Store address of an argument into a list */
void f77name(store_argument)(void **dest, void *arg)
{
   *dest = arg;
}

/*! Retrieve address of an argument from a list */
void f77name(fetch_argument)(void **dest, void **arg)
{
   *dest = *arg;
}
