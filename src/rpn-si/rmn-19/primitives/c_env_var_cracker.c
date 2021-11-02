#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pwd.h>
#include <rpnmacros.h>

#define null 0
#define MAXSTRINGS 32
#define MAXLINELENGTH 255

#define NONMATCHING -1

#define QUOTES     50
#define NOQUOTES  -2


#define _cptofcd(a,b) (a)

int check_start_end_char(char *var, int length);
static void convert_toupper(char *string);
static void trimleft(char *string);
static void trimright(char *string);
static void remove_quotes(char * string);
static void call_user_function(char *name, int index, char *value, char *lang, void (*user_function)());

void c_env_var_cracker(char *fstoption, void (*user_function)(), char *lang)
{
  char *token, *option;
  char *delimiter1 = ";", *delimiter2 = "=", *delimiter3 = ", ";
    
  char strings[MAXSTRINGS][MAXLINELENGTH];
  
  char *name, *values, *value;
  int length = 0, index = 0, i = 0, result = 0;

   if((option = getenv(fstoption)) == NULL)
    {
      return;
    }

   if(strlen(option = getenv(fstoption)) == 0)
    {
      return;
    }
   
#ifdef DEBUG
  fprintf(stderr, "option: %s\n", option);
#endif  
  
  /* parse options (extract variable name and values) using ";" separator */
  token = strtok(option, delimiter1);

  strncpy(strings[length], token, strlen(token));
  strings[length][strlen(token)] = '\0';
  length++;
 

  name = (char *)malloc(MAXLINELENGTH);
  values = (char *)malloc(MAXLINELENGTH);
  value = (char *)malloc(MAXLINELENGTH);

  while((token = strtok(null, delimiter1)) != null)
    {
      strncpy(strings[length], token, strlen(token));
      strings[length][strlen(token)] = '\0';
      length++;
    }
  
  for(i = 0; i < length; i++)
    {
      /* parse variable name /value using "=" separator */
      token = strtok(strings[i], delimiter2);
      strncpy(name, token, strlen(token));
            
      name[strlen(token)] = '\0'; 
      
      /* convert key to upper case */
      convert_toupper(name);
      /* remove leading whitespace */
      trimleft(name);
      /* remove trailing whitespace */
      trimright(name);

#ifdef DEBUG
      fprintf(stderr, "name: %s\n", name);
#endif
     
      /* parse variable name and values using "=" separator */
      while((token = strtok(null, delimiter2)) != NULL)
	{  

	  trimleft(token);
	  if(strlen(token) == 0)
	    {
	      continue;
	    }

	  strncpy(values, token, strlen(token));
	  values[strlen(token)] = '\0';
	  

#ifdef DEBUG
	  fprintf(stderr, "Values: %s\n", values);
#endif
	  
	  /* remove leading whitespace */
	  trimleft(values);
	  /* remove trailing whitespace */
	  trimright(values);

	  /* check opening and closing parenthesis/quotes */ 
	  result = check_start_end_char(values, strlen(values));

#ifdef DEBUG
	  fprintf(stderr, "white_space_count = %d\n", white_space_count);
#endif

	  if(result == NONMATCHING)
	    {
	      return;
	    }  
	  
	  if(result >= QUOTES)
	    {
	      /* remove quotes */ 
	      remove_quotes(values);
	    }

#ifdef DEBUG
	  fprintf(stderr, "Value(s) without white space and parentheses: %s, result = %d \n", values, result);
#endif
	  index++;

	  if(result >= QUOTES && strchr(values, ',') == NULL)
	    {
	      call_user_function(name, index, values, lang, user_function);
   	    }
	  
	  else /* case: multiple values per variable */
	    {
	      index = 1;

	      /* parse values using "," separator */
	      token = strtok(values, delimiter3);
	      strncpy(value, token, strlen(token));
	      value[strlen(token)] = '\0';
 
#ifdef DEBUG	      
	      fprintf(stderr, "Value, without white space and parentheses: %s, token = %s\n", value, token);
#endif

	      /* remove leading whitespace */
	      trimleft(value);
	      /* remove trailing whitespace */
	      trimright(value);


	      call_user_function(name, index, value, lang, user_function);

	      while((token = strtok(null, delimiter3)) != null)
		{
		  index++;
		  strncpy(value, token, strlen(token));
		  value[strlen(token)] = '\0';

		  /* remove leading whitespace */
		  trimleft(value);
		  /* remove trailing whitespace */
		  trimright(value);

		  call_user_function(name, index, value, lang, user_function);
			  
		}
	    }
	  
	}

    }

  if(name)
    free(name);
  if(values)
    free(values);
  if(value)
    free(value);
  
}


static void call_user_function(char *name, int index, char *value, char *lang, void (*user_function)())
{
  char *fp1, *fp2;

  if(strcmp(lang, "C") == 0) /* option: user function in C */
    {
#ifdef DEBUG
      fprintf(stderr, " Call user C function\n");
#endif
      user_function(name, index, value);
      
    }
  else if (strcmp(lang, "F") == 0) /* option: user function in Fortran */
    {
#ifdef DEBUG
      fprintf(stderr, "Call user Fortran function\n");
#endif
      fp1 = _cptofcd(name, strlen(name));
      fp2 = _cptofcd(value, strlen(value));
      user_function(fp1, &index, fp2, strlen(fp1), strlen(fp2));
            
    }
}


/* check values for matching/non matching quotes parenthesis */
int check_start_end_char(char *string, int length)
{
  int space_count = 0;

  /* count leading whitespaces in string */  
  while(isspace(string[space_count]))
    {
      space_count++;
    }

  if(string[space_count] == '[')
    {
      if(string[length - 1] == ']')
	return QUOTES;
      else
	return NONMATCHING;
    }
  
  if(string[length - 1] == ']')
    {
      if(string[space_count] == '[')

	return QUOTES;
      else
	return NONMATCHING;
    }
  
  if(string[length - 1] == '(')
    {
      if(string[space_count] == ')')

	return QUOTES;
      else
	return NONMATCHING;
    }
  
  if(string[space_count] == '(' )
    {
      if(string[length - 1] == ')')
	return QUOTES;
      else
	return NONMATCHING;
    }
  
  if(string[length - 1] == '}')
    {
      if(string[space_count] == '{')
	return QUOTES;
      else
	return NONMATCHING;
    }
  
  if(string[space_count] == '{')
    {
      if(string[length - 1] == '}')
	return QUOTES;
      else
	return NONMATCHING;
    }
  
  if(string[length - 1] == '\"')
    {
      if(string[space_count] == '\"')
	return QUOTES;
      else
	return NONMATCHING;
    }
  
  if(string[space_count] == '\"')
    {
      if(string[length - 1] == '\"')
	return QUOTES;
      else
	return NONMATCHING;
      
    }
  
  if(string[length - 1] == '\'')
    {
      if(string[space_count] == '\'')
	return QUOTES;
      else
	return NONMATCHING;
    }
  
  if(string[space_count] == '\'')
    {
      if(string[length - 1] == '\'')
	return QUOTES;
      else
	return NONMATCHING;
    }
  
  if(space_count >= 0)
    return space_count;

  return NOQUOTES;
}

/* convert a string to upper case */
static void convert_toupper(char *string)
{
  char *key;
  for(key = string; *key != '\0'; key++)
    *key = toupper(*key);
}


/* strip all leading whitesapce from a string */
static void trimleft(char *string)
{
  char *string_dest;

  if (string != NULL)
    {
      string_dest = string;
      while (*string != '\0' && isspace(*string))
	string++;
      while (*string != '\0')
	*string_dest++ = *string++;
      *string_dest = '\0';
    }
}

/*  strip all trailing whitesapce from a string*/
static void trimright(char *string)
{
  char *string_dest;

  if (string != NULL && string[0] != '\0')
    {
      string_dest = string + strlen(string) - 1;
      while(isspace(*string_dest))
	*string_dest-- = '\0';

    }

}

static void remove_quotes(char * string)
{
  int counter = 0;
  for(counter = 0; counter < strlen(string) - 2; counter++)
    string[counter] = string[counter + 1 ]; 

  string[strlen(string) - 2] = '\0';
  
}
