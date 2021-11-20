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
	//メンバ関数
	//コンストラクタ
	_vec(void);
	_vec(const TEMPLATE&, const TEMPLATE& );
	_vec(const _vec<TEMPLATE>&);
	
	//オペレータ
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




//メンバ関数実体
//コンストラクタ
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

//オペレータ
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

//スカラ対ベクトル
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

//ベクトル対ベクトル
template <class TEMPLATE> inline TEMPLATE operator*(const _vec<TEMPLATE>& v1, const _vec<TEMPLATE>& v2){
	return v1.x * v2.x + v1.y * v2.y;
}

template <class TEMPLATE> inline TEMPLATE operator%(const _vec<TEMPLATE>& v1, const _vec<TEMPLATE>& v2){
	return v1.x * v2.y - v1.y * v2.x;
}

template <class TEMPLATE> inline _vec<TEMPLATE> operator+(const _vec<TEMPLATE>& v1, const _vec<TEMPLATE>& v2){
	
	return _vec<TEMPLATE>(v1.x+v2.x,v1.y+v2.y);
}

template <class TEMPLATE> inline _vec<TEMPLATE> operator-(const _vec<TEMPLATE>& v1, const _vec<TEMPLATE>& v2){
	return _vec<TEMPLATE>(v1.x-v2.x,v1.y-v2.y);
}



//ベクトル回転
template <class TEMPLATE> inline _vec<TEMPLATE>& _vec<TEMPLATE>::rot( const TEMPLATE& angle )
{
  TEMPLATE c=cos(angle), s=sin(angle);
  TEMPLATE t;
  t = c*x - s*y;
  y = s*x + c*y;
  x = t;

  return *this;  
}

//ベクトルノルム
template <class TEMPLATE> inline TEMPLATE _vec<TEMPLATE>::norm(void) const{
	TEMPLATE s = x*x + y*y;
	return sqrt(s);
}

#endif
