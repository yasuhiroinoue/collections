/*
Yasuhiro Inoue
inoue.yasuhiro.4n@kyoto-u.ac.jp
*/
#ifndef RCLASS_H
#define RCLASS_H

#include "vec.h"
#include <map>
using namespace std;


class _cell_base{
  public:
	int count;
	_vec<int> loc;
	double mas;
	_vec<double> vel;
	double ene;
	double rot[2];
	_cell_base(void);
	
};


class _ptcl_base{
  public:
	int col;
	double mas;
	_vec<double> loc;
	_vec<double> vel;
	_ptcl_base(void);
};

class _list{
  public:
	multimap<int, _ptcl_base* > pmap;
	_list* lptr[8];
	double number_color(int i) const;
};

class _node{
	public:
	int count;
	double col;
	double mas;
	_vec<double> mom;
	double s[3];
	_node(void);
};



class _total{
  public:
  double   mas;
  _vec<double> mom;
  double   ene;
	_total(void);
};

inline _cell_base :: _cell_base(void){}
inline _ptcl_base :: _ptcl_base(void){}
inline _total :: _total(void){}
inline _node :: _node(void){}
inline double _list::number_color(int i) const{
	double color;
	color = (double)pmap.count(i);
	return color;
}
#endif
