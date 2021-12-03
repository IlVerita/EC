#include "Vettore.h"
#include "OdeSolver.h"

void OdeSolver::Punto(PuntoMateriale tmp){
  m_p.push_back(tmp);
}

vector<PuntoMateriale> OdeSolver::Punti(){
  return m_p;
}

PuntoMateriale OdeSolver::Punto(unsigned int i){
  return m_p[i];
}

unsigned int OdeSolver::N(){
  return m_p.size();
}

void OdeSolver::T(double t0){
  m_t=t0;
}

double OdeSolver::T(){
  return      m_t;
}

void OdeSolver::Passo(double h){
  m_h = h;
}

double OdeSolver::Passo(){
  return   m_h;
}


Vettore OdeSolver::m_eqDiff(unsigned int i, double t, vector<PuntoMateriale> p){
  //Calcolo dell'accelerazione dovuta a forze interne e forze esterne
  double mi=p[i].Massa();
  Vettore A_fake;
  Vettore A;
  for(unsigned int j=0;j<p.size();j++){
    if (i!=j)
      A_fake=A_fake+fInterna(i,j,t,p)*(1/mi);
      }
  A=A_fake+fEsterna(i,t,p)*(1/mi);
  return A;
}


void OdeSolver::Cinematica(){

  if (m_method=="Eulero"){
    cout <<"Eulero!"<< endl;
    vector<Vettore>  k1(m_p.size());
    vector<Vettore>  w1(m_p.size());
    for (unsigned int i=0;i<m_p.size();i++){
      k1[i] = m_h*m_p[i].V();
      w1[i] = m_h*m_eqDiff(i,m_t,m_p);
    }

    for (unsigned int i=0;i<m_p.size();i++){
      m_p[i].R(m_p[i].R() + k1[i]);
      m_p[i].V(m_p[i].V() + w1[i]);
    }

  } else if (m_method=="Rk2"){
    vector<Vettore> k1(m_p.size());
    vector<Vettore> k2(m_p.size());
    vector<Vettore> w1(m_p.size());
    vector<Vettore> w2(m_p.size());
    
    for (unsigned int i=0;i<m_p.size();i++){
      k1[i]=m_h*m_p[i].V();
      w1[i]=m_h*m_eqDiff(i,m_t,m_p);
    }
    auto m_temp(m_p);
    for(unsigned int i=0;i<m_p.size();i++){
         m_temp[i].R(m_p[i].R()+k1[i]*(0.5));
      m_temp[i].V(m_p[i].V()+w1[i]*(0.5));
    }
    for (unsigned int i=0;i<m_p.size();i++){
      k2[i]=m_h*m_temp[i].V();
      w2[i]=m_h*m_eqDiff(i,m_t+m_h*0.5,m_temp);
    }
    for (unsigned int i=0;i<m_p.size();i++){
      m_p[i].R(m_p[i].R()+k2[i]);
      m_p[i].V(m_p[i].V()+w2[i]);
    }
    
  }else if(m_method=="Rk4"){
    vector<Vettore> k1(m_p.size());
    vector<Vettore> k2(m_p.size());
    vector<Vettore> w1(m_p.size());
    vector<Vettore> w2(m_p.size());
    vector<Vettore> k3(m_p.size());
    vector<Vettore> k4(m_p.size());
    vector<Vettore> w3(m_p.size());
    vector<Vettore> w4(m_p.size());

    for (unsigned int i=0;i<m_p.size();i++){
      k1[i]=m_h*m_p[i].V();
      w1[i]=m_h*m_eqDiff(i,m_t,m_p);
    }

    auto m_temp0(m_p);
    for(unsigned int i=0;i<m_p.size();i++){
    m_temp0[i].R(m_p[i].R()+k1[i]*(0.5));
    m_temp0[i].V(m_p[i].V()+w1[i]*(0.5));
    }

    for(unsigned int i=0;i<m_p.size();i++){
      k2[i]=m_h*m_temp0[i].V();
      w2[i]=m_h*m_eqDiff(i,m_t+m_h*0.5,m_temp0);
    }

    auto m_temp1(m_temp0);
    for(unsigned int i=0;i<m_p.size();i++){
      m_temp1[i].R(m_p[i].R()+k2[i]*(0.5));
      m_temp1[i].V(m_p[i].V()+w2[i]*(0.5));      
    }

    for(unsigned int i=0;i<m_p.size();i++){
      k3[i]=m_h*m_temp1[i].V();
      w3[i]=m_h*m_eqDiff(i,m_t+m_h*(0.5),m_temp1);
    }

    auto m_temp2(m_temp1);
    for(unsigned int i=0;i<m_p.size();i++){
      m_temp2[i].R(m_p[i].R()+k3[i]);
      m_temp2[i].V(m_p[i].V()+w3[i]);
     }
    

    for(unsigned int i=0;i<m_p.size();i++){
      k4[i]=m_h*m_temp2[i].V();
      w4[i]=m_h*m_eqDiff(i,m_t+m_h*(0.5),m_temp2);
	
    }

    for (unsigned int i=0;i<m_p.size();i++){
      m_p[i].R(m_p[i].R()+k1[i]*(1/6.)+k2[i]*(1/3.)+k3[i]*(1/3.)+k4[i]*(1/6.));
      m_p[i].V(m_p[i].V()+w1[i]*(1/6.)+w2[i]*(1/3.)+w3[i]*(1/3.)+w4[i]*(1/6.));
    }

  }
  else if(m_method=="Verlet-Velocity"){
    auto present(m_p); //posizioni e velocità al tempo iniziale
    auto future(m_p);
    for (unsigned int i=0;i<m_p.size();i++){
      m_p[i].R(m_p[i].R()+m_h*m_p[i].V()+(m_h*m_h)*m_eqDiff(i,m_t,present)*(1/2.));
      auto future(m_p); //L'assegnazione poteva essere fatta fuori dal for 
    }
    for (unsigned int i=0;i<m_p.size();i++){
      m_p[i].V(m_p[i].V()+m_h*(1/2.)*(m_eqDiff(i,m_t+m_h,future)+m_eqDiff(i,m_t,present)));
    }
      
  }

  else if(m_method=="Eulero-Cromer"){
    for (unsigned int i=0;i<m_p.size();i++){
      m_p[i].V(m_p[i].V()+m_h*m_eqDiff(i,m_t,m_p));
      m_p[i].R(m_p[i].R()+m_h*m_p[i].V());
      
    }
  }
  else{
    cout <<"Metodo scelto non valido"<< endl;
  }
  m_t += m_h;

}

