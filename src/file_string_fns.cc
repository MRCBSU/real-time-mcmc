#include "string_fns.h"

#define READ_ERROR_MESSAGE   printf("No instances of string:\"%s\" in file %s. Instance expected\n", varname, filename); \
      exit(2);								\


using namespace std;

// LOADS THE ENTIRE CONTENTS OF THE FILE filename INTO THE
// CHARACTER ARRAY POINTED TO BY OUT_STRING
// Throws std::ios_base::failure on error
void fn_load_file(string *out_string, const char *filename)
{

  vector<string> text;
  string line;
  ifstream textstream(filename);
  textstream.exceptions(ifstream::badbit);
  while(getline(textstream, line)){
    text.push_back(line + "\n");
  }
  textstream.close();

  for(int i = 0; i < text.size(); i++)
    *out_string += text[i];


}


// RETURNS THE NUMBER OF INSTANCES OF STRING sep_string IN FILE filename
int count_instances_in_file(const char *filename, const char *sep_string)
{

  string input;

  fn_load_file(&input, filename);

  return count_instances_in_string(input, sep_string);

}

// RETURNS A POINTER, str_ptr_out, TO THE n_instance^TH OCCURENCE OF THE CHARACTER ARRAY
// POINTED TO BY sep_string IN THE FILE NAMED filename
// MEMORY USED BY str_ptr_out IS ALLOCATED AND NEEDS TO BE FREED BY THE CALLING ROUTINE
void nth_instance(char *&str_ptr_out, const char *filename, const char *sep_string, const int n_instance)
{

  string input;
  int tot_appearances = 0;

  // LOAD FILE INTO A STRING
  fn_load_file(&input, filename);

  // COUNT THE NUMBER OF INSTANCES
  tot_appearances = count_instances_in_file(filename, sep_string);

  // CHECK THERE ARE SUFFICIENT INSTANCES IN THE FILE
  if(n_instance > tot_appearances)
    { // ERROR
      printf("Insufficient occurences of %s in %s\n", sep_string, filename);
      exit(2);
    }

  int substring_position = nth_instance_in_string(input, sep_string, tot_appearances);

  if(substring_position > -1)
    {
      // REALLOC THE MEMORY OF str_ptr_out (OF SIMPLY ALLOC IF NULL)
      if(str_ptr_out == NULL)
	str_ptr_out = (char *) calloc(strlen(input.c_str() + substring_position) + 1, sizeof(char)); // +1 TO ALLOW FOR THE NULL TERMINATING CHARACTER
      else
	str_ptr_out = (char *) realloc(str_ptr_out, (strlen(input.c_str() + substring_position) + 1) * sizeof(char)); // +1 TO ALLOW FOR THE NULL TERMINATING CHARACTER

      strcpy(str_ptr_out, input.c_str() + substring_position);

    } else str_ptr_out = NULL;
  
}


// FINDS THE LAST (DOUBLE) ASSIGNMENT OF VARIABLE varname IN FILE filename
void read_double_input(double& outval, const char* filename, const char* varname)
{

  char *tempchar = NULL;
  int num = count_instances_in_file(filename, varname);

  if(num <= 0){
    READ_ERROR_MESSAGE;}

  nth_instance(tempchar, filename, varname, num);

  if(tempchar != NULL){
    outval = atof(tempchar + strlen(varname));
    free(tempchar);
  } // else : do nothing, leave the value to which outval points unchanged

}

// FINDS THE LAST (INTEGER) ASSIGNMENT OF VARIABLE varname IN FILE filename
void read_int_input(int& outval, const char* filename, const char* varname)
{

  char *tempchar = NULL;
  int num = count_instances_in_file(filename, varname);

  if(num <= 0){
    READ_ERROR_MESSAGE;}

  nth_instance(tempchar, filename, varname, num);

  if(tempchar != NULL){
    outval = atoi(tempchar + strlen(varname));
    free(tempchar);
  } // else : do nothing, leave the value to which outval points unchanged

}


// FINDS THE LAST (CHARACTER STRING OF MAX LENGTH max_string_length)
// ASSIGNMENT OF VARIABLE varname IN FILE filename
// LEAVES out_string UNCHANGED IF varname NOT PRESENT IN filename
void read_char_input(char* out_string, const char* filename, const char* varname, const int max_string_length)
{

  // out_string is assumed to be pre-allocated
  char *tempchar = NULL;
  int num = count_instances_in_file(filename, varname);

  nth_instance(tempchar, filename, varname, num);

  if(tempchar != NULL){
    sscanf(tempchar + strlen(varname), "%s", out_string);

    if(strlen(out_string) >= max_string_length)
      { // ERROR
	printf("Input string exceeds max number of characters (length = %d)\n", (int) strlen(out_string));
	exit(2);
      }

    free(tempchar);

  } // else: do nothing, leave the string unchanged

}



