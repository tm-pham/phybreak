#include <Rcpp.h>
#include <vector>
#include <iterator>
#include <cmath>
using namespace Rcpp;

// Calculation of the path to the root in a tree,
// given a parent vector and node/host ID

std::vector<int> ptroot(std::vector<int> pars, int ID) {
  std::vector<int> ans;
  ans.push_back(ID);
  int nextID = pars[ID - 1];
  while(nextID > 0) {
    ans.push_back(nextID);
    nextID = pars[nextID - 1];
  }
  return(ans);
}




// [[Rcpp::export(name=".CCtranstree")]]
std::vector<int> CCtranstree(std::vector<int> pars, std::vector<int> dims) {
  std::vector<bool> cladearray(dims[0] * dims[0] * dims[1]); //dims[0] = obs, dims[1] = samplesize
  std::vector<int> parset(dims[0]);
  std::vector<int> pos;

  for(int i = 0; i < dims[1]; ++i) {
    for(int ii = 0; ii< dims[0]; ++ii) {
      parset[ii] = pars[i*dims[0] + ii];
    }
    for(int j = 0; j < dims[0]; ++j) {
      pos = ptroot(parset, j + 1);
      for(unsigned int k = 0; k < pos.size(); ++k) {
        cladearray[i*dims[0]*dims[0] + (pos[k]-1)*dims[0] + j] = true;
      }
    }
  }

  std::vector<int> uniquecladepos;
  std::vector<int> uniquecladecount;
  std::vector<int> whichclade(dims[0] * dims[1]);
  bool cladeequal = false;
  bool nodeequal = true;
  unsigned int nextclade = 0;
  int nextnode = 0;
  uniquecladepos.push_back(0);
  uniquecladecount.push_back(1);
  for(int i = 1; i < dims[0] * dims[1]; ++i) {
    while(!cladeequal && nextclade < uniquecladepos.size()) {
      while(nodeequal && nextnode < dims[0]) {
        if(cladearray[i*dims[0] + nextnode] !=
           cladearray[uniquecladepos[nextclade]*dims[0] + nextnode]) {
          nodeequal = false;
        }
        ++nextnode;
      }
      cladeequal = nodeequal;
      nodeequal = true;
      ++nextclade;
      nextnode = 0;
    }
    if(!cladeequal) {
      uniquecladepos.push_back(i);
      uniquecladecount.push_back(1);
      whichclade[i] = uniquecladepos.size() - 1;
    } else {
      ++uniquecladecount[nextclade - 1];
      whichclade[i] = nextclade - 1;
    }
    cladeequal = false;
    nodeequal = true;
    nextclade = 0;
    nextnode = 0;
  }

  std::vector<double> cladescores(dims[0] * dims[1]);
  for(int i = 0; i < dims[0]*dims[1]; ++i) {
    cladescores[i] = std::log(uniquecladecount[whichclade[i]]);
  }

  int postsample = 0;

  std::vector<double> treescores(dims[1]);
  for(int i = 0; i < dims[1]; ++i) {
    for(int j = 0; j < dims[0]; ++j) {
      treescores[i] += cladescores[i * dims[0] + j];
    }
    if(treescores[i] > treescores[postsample]) {
      postsample = i;
    }
  }



 std::vector<int> result(2 * dims[0]);
  for(int i = 0; i < dims[0]; ++i) {
    result[i] = pars[dims[0] * postsample + i];
    result[i + dims[0]] = uniquecladecount[whichclade[dims[0] * postsample + i]];
  }




  return result;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

