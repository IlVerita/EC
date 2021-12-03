#include "Vettore.h"
#include <cmath>
#include <iostream>


double Vettore::X() const{
  return m_v[0];
}

double Vettore::Y() const{
  return m_v[1];
}

double Vettore::Z() const{
  return m_v[2];
}

void Vettore::X(double x){
  m_v[0] = x;
}

void Vettore::Y(double y){
  m_v[1] = y;
}

void Vettore::Z(double z){
  m_v[2] = z;
}

Vettore Vettore::operator-(){
  Vettore res;
  for (int i=0;i<3;i++){
    res.m_v[i] = -m_v[i];
  }
  return res;
}


Vettore Vettore::operator-( Vettore b){
  Vettore res;
  for (int i=0;i<3;i++){
    res.m_v[i] = m_v[i] - b.m_v[i];
  }
  return res;
}

Vettore Vettore::operator+( Vettore b){
  Vettore res;
  for (int i=0;i<3;i++){
    res.m_v[i] = m_v[i] + b.m_v[i];
  }
  return res;
}

double Vettore::operator*( Vettore b){
  double prod=0;
  for (int i=0;i<3;i++){
    prod += m_v[i]*b.m_v[i];
  }
  return prod;
}

Vettore Vettore::Cross( Vettore b){
  Vettore res;
  res.m_v[0] = m_v[1]*b.m_v[2]-b.m_v[1]*m_v[2];
  res.m_v[1] = m_v[2]*b.m_v[0]-b.m_v[2]*m_v[0];
  res.m_v[2] = m_v[0]*b.m_v[1]-b.m_v[0]*m_v[1];
  return res;
}

Vettore Vettore::operator*(double f){
  Vettore res;
  for (int i=0;i<3;i++){
    res.m_v[i] = m_v[i]*f;
  }
  return res;
}

double Vettore::Mod() {
  return sqrt((*this)*(*this));
}

Vettore Vettore::Vers() {
  return (*this)*(1/Mod());
}

Vettore operator*(double f,  Vettore b){
  return b*f;
}

std::ostream& operator<<(std::ostream& os,  Vettore b){
  os << "(";
  os << b.X() << ",";
  os << b.Y() << ",";
  os << b.Z() << ")";
  return os;
}

std::istream& operator>>(std:: istream& is, Vettore& b){
  double x,y,z;
  is >> x >> y >> z;
  b  = Vettore(x,y,z);
  return is;
}
