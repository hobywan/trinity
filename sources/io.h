
#ifndef IO_H
#define IO_H
/* ------------------------------------*/
#include "tools.h"
#include "mesh.h"
/* ------------------------------------*/

namespace trinity {
  namespace io {
    int find(const std::string key, std::ifstream& file);
  };
}
#endif
