/* 
Yasuhiro Inoue
inoue.yasuhiro.4n@kyoto-u.ac.jp
*/
#ifndef VEC_H

#define VEC_H

#include <stdio.h>

#include <math.h>



template <class TEMPLATE>

class _vec{

  public:

	TEMPLATE x;

	TEMPLATE y;

	//member functions
	//constractor

	_vec(void);

	_vec(const TEMPLATE&, const TEMPLATE& );

	_vec(const _vec<TEMPLATE>&);

	

	//operator

	void IN(const TEMPLATE& , const TEMPLATE& );

	void IN(const _vec<TEMPLATE>&);

	_vec<int> icast(void) const;

	_vec<TEMPLATE>& operator = (const _vec<TEMPLATE>&);

	_vec<TEMPLATE>& operator += (const _vec<TEMPLATE>&);

	_vec<TEMPLATE>& operator -= (const _vec<TEMPLATE>&);

	_vec<TEMPLATE>& operator *= (const _vec<TEMPLATE>&);

	_vec<TEMPLATE>& operator /= (const _vec<TEMPLATE>&);

	_vec<TEMPLATE>& operator *= (const TEMPLATE&);

	_vec<TEMPLATE>& operator /= (const TEMPLATE&);

	_vec<TEMPLATE>& rot( const TEMPLATE&);

	TEMPLATE norm(void) const;

};









//member functions

//constractor

template <class TEMPLATE> inline _vec<TEMPLATE>::_vec(void){

	x = 0;

	y = 0;

}



template <class TEMPLATE> inline _vec<TEMPLATE>::_vec(const TEMPLATE& a, const TEMPLATE& b){

	x = a;

	y = b;

}



template <class TEMPLATE> inline _vec<TEMPLATE>::_vec(const _vec<TEMPLATE>& v){

	x = v.x;

	y = v.y;

}



//operator

template <class TEMPLATE> inline void _vec<TEMPLATE>::IN(const TEMPLATE& a, const TEMPLATE& b){

	x = a; y = b;

}



template <class TEMPLATE> inline void _vec<TEMPLATE>::IN(const _vec<TEMPLATE>& v){

	x = v.x; y = v.y;

}



template<class TEMPLATE> inline _vec<int> _vec<TEMPLATE>::icast(void) const

{

  return _vec<int>((int)x, (int)y);

}



template <class TEMPLATE> inline _vec<TEMPLATE>& _vec<TEMPLATE>::operator =  ( const _vec<TEMPLATE>& v ){

  x = v.x;    y = v.y;

  return *this;

}



template <class TEMPLATE> inline _vec<TEMPLATE>& _vec<TEMPLATE>::operator +=  ( const _vec<TEMPLATE>& v ){

  x += v.x;    y += v.y;

  return *this;

}



template <class TEMPLATE> inline _vec<TEMPLATE>& _vec<TEMPLATE>::operator -=  ( const _vec<TEMPLATE>& v ){

  x -= v.x;    y -= v.y;

  return *this;

}



template <class TEMPLATE> inline _vec<TEMPLATE>& _vec<TEMPLATE>::operator *= ( const TEMPLATE& a){

	x *= a; y *= a;

	return *this;

}



template <class TEMPLATE> inline _vec<TEMPLATE>& _vec<TEMPLATE>::operator /= ( const TEMPLATE& a){

	x /= a; y /= a;

	return *this;

}



template <class TEMPLATE> inline _vec<TEMPLATE>& _vec<TEMPLATE>::operator *= ( const _vec<TEMPLATE>& v){

	x *= v.x; y *= v.y;

	return *this;

}



template <class TEMPLATE> inline _vec<TEMPLATE>& _vec<TEMPLATE>::operator /= ( const _vec<TEMPLATE>& v){

	x /= v.x; y /= v.y;

	return *this;

}



//scalor vs vector
template <class TEMPLATE> inline _vec<TEMPLATE> operator*( const TEMPLATE& a, const _vec<TEMPLATE> v  ){

	_vec<TEMPLATE> vec(a*v.x,a*v.y);

	return vec;

}



template <class TEMPLATE> inline _vec<TEMPLATE> operator*( const _vec<TEMPLATE> v, const TEMPLATE& a  ){

	_vec<TEMPLATE> vec(v.x*a,v.y*a);

	return vec;

}





template <class TEMPLATE> inline _vec<TEMPLATE> operator/( const TEMPLATE& a, const _vec<TEMPLATE> v  ){

	_vec<TEMPLATE> vec(a/v.x, a/v.y);

	

	return vec;

}



template <class TEMPLATE> inline _vec<TEMPLATE> operator/(const _vec<TEMPLATE> v, const TEMPLATE& a  ){

	_vec<TEMPLATE> vec(v.x/a,v.y/a);

	return vec;

}



template <class TEMPLATE> inline _vec<TEMPLATE> operator+(const _vec<TEMPLATE> v, const TEMPLATE& a  ){

	_vec<TEMPLATE> vec(v.x+a,v.y+a);

	return vec;

}



template <class TEMPLATE> inline _vec<TEMPLATE> operator+(const TEMPLATE& a ,const _vec<TEMPLATE> v ){

	_vec<TEMPLATE> vec(a+v.x,a+v.y);

	return vec;

}



template <class TEMPLATE> inline _vec<TEMPLATE> operator-(const _vec<TEMPLATE> v, const TEMPLATE& a  ){

	_vec<TEMPLATE> vec(v.x-a,v.y-a);

	return vec;

}



template <class TEMPLATE> inline _vec<TEMPLATE> operator-(const TEMPLATE& a ,const _vec<TEMPLATE> v ){

	_vec<TEMPLATE> vec(a-v.x,a-v.y);

	return vec;

}



//vector vs vector
template <class TEMPLATE> inline TEMPLATE operator*(const _vec<TEMPLATE>& v1, const _vec<TEMPLATE>& v2){

	return (v1.x * v2.x + v1.y * v2.y);

}



template <class TEMPLATE> inline TEMPLATE operator%(const _vec<TEMPLATE>& v1, const _vec<TEMPLATE>& v2){

	return (v1.x * v2.y - v1.y * v2.x);

}



template <class TEMPLATE> inline _vec<TEMPLATE> operator+(const _vec<TEMPLATE>& v1, const _vec<TEMPLATE>& v2){

	

	return _vec<TEMPLATE>(v1.x+v2.x,v1.y+v2.y);

}



template <class TEMPLATE> inline _vec<TEMPLATE> operator-(const _vec<TEMPLATE>& v1, const _vec<TEMPLATE>& v2){

	return _vec<TEMPLATE>(v1.x-v2.x,v1.y-v2.y);

}







//rotation vector

template <class TEMPLATE> inline _vec<TEMPLATE>& _vec<TEMPLATE>::rot( const TEMPLATE& angle )

{

  TEMPLATE c=cos(angle), s=sin(angle);

  TEMPLATE t;

  t = c*x - s*y;

  y = s*x + c*y;

  x = t;



  return *this;  

}



//norm of vector

template <class TEMPLATE> inline TEMPLATE _vec<TEMPLATE>::norm(void) const{

	TEMPLATE s = x*x + y*y;

	return sqrt(s);

}



#endif

