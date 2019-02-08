#include "../include/log.hh"

int return_with_mesg(const char *mesg, int return_code) {
  fprintf(stderr,"[ERROR] %s\n", mesg);
  return return_code;
}

