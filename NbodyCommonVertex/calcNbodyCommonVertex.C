////////////////////////////////////////////////////////////
// Author: Xianguo Lu, Xianguo.Lu@physics.ox.ac.uk        //
////////////////////////////////////////////////////////////

#include "stdlib.h"
#include <iostream>

#include "TVector3.h"
#include "TMatrixD.h"

using namespace std;

class LRay
{
public:
  LRay(const TVector3 aa, const TVector3 vv):fStartingPoint(aa), fUnitDirection(vv.Unit()) {}
  void Print(){ fStartingPoint.Print(); fUnitDirection.Print(); }

  TMatrixD GetMatrixA();
  TMatrixD GetMatrixB();
  const double GetDCA2(const TVector3 xx);

private:
  const TVector3 fStartingPoint;
  const TVector3 fUnitDirection;
};

TMatrixD LRay::GetMatrixA()
{
  TMatrixD aa(3,3);

  for(int ii=0; ii<3; ii++){
    for(int jj=0; jj<3; jj++){
      aa[ii][jj]=fUnitDirection[ii]*fUnitDirection[jj];
    }
  }

  return aa;
}

TMatrixD LRay::GetMatrixB()
{
  const TVector3 bv = fStartingPoint - ( fStartingPoint.Dot(fUnitDirection) )*fUnitDirection;

  TMatrixD bb(3,1);

  for(int ii=0; ii<3; ii++){
    bb[ii][0] = bv[ii];
  }

  return bb;
}

const double LRay::GetDCA2(const TVector3 xx)
{
  return (fUnitDirection.Cross(xx-fStartingPoint)).Mag2();
}

double GetDCA2(const int nray, LRay rs[], TVector3 *vertex=0x0)
{
  const int kprint = 1;

  if(kprint){
    cout<<"nray: "<<nray<<endl;
    for(int ii=0; ii<nray; ii++){
      cout<<"ray: "<<ii<<endl;
      rs[ii].Print();
      rs[ii].GetMatrixA().Print();
      rs[ii].GetMatrixB().Print();
      cout<<endl;
    }
  }

  TMatrixD lhs(3,3);
  for(int kk=0; kk<3; kk++){
    lhs[kk][kk]=nray;
  }
  for(int ii=0; ii<nray; ii++){
    lhs -= rs[ii].GetMatrixA();
  }
  if(kprint){
    cout<<"LHS: "<<endl;
    lhs.Print();
    cout<<endl;
  }

  TMatrixD rhs(3,1);
  for(int ii=0; ii<nray; ii++){
    rhs += rs[ii].GetMatrixB();
  }
  if(kprint){
   cout<<"RHS: "<<endl;
   rhs.Print();
   cout<<endl;
  }

  TMatrixD xx(3,1);
  xx=lhs.Invert()*rhs;
  const TVector3 Xvec(xx[0][0],xx[1][0],xx[2][0]);
  if(kprint){
    cout<<"X: "<<endl;
    Xvec.Print();
    cout<<endl;
  }

  if(vertex){
    (*vertex)=Xvec;
  }

  double dca2= 0;
  for(int ii=0; ii<nray; ii++){
    dca2 += rs[ii].GetDCA2(Xvec);
  }
  return dca2;
}

void calcNbodyCommonVertex()
{
  //-->
  //this will give vertex (x,y,z)=(0.500000,0.500000,0.500000), DCA2: 1.5
  const TVector3 a0(0,1,0);
  const TVector3 v0(10,0,0);
  const LRay r0(a0, v0);

  const TVector3 a1(0,0,1);
  const TVector3 v1(0,7,0);
  const LRay r1(a1, v1);

  const TVector3 a2(1,0,0);
  const TVector3 v2(0,0,90);
  const LRay r2(a2, v2);
  //<---

  LRay rs[]={r0, r1, r2};
  const int nray = sizeof(rs)/sizeof(LRay);

  const double dca2 = GetDCA2(nray, rs);
  cout<<"DCA2: "<<dca2<<endl<<endl;
}

int main(int argc, char* argv[])
{
  calcNbodyCommonVertex();
  return 0;
}
