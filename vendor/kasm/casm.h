//---------------------------------------------------------------------------
typedef double Real;
extern const Real UNDEF;   // testowanie czy wartosc niezdefiniowana
//---------------------------------------------------------------------------
class rs
{public:
 static void rs1 (int n, Real* x[2], Real* y[2], const Real& scl=1.);
 private:
 static const Real del;
 static const Real dela;    // do liczenia wag
};
//---------------------------------------------------------------------------
typedef Real* pReal;
//---------------------------------------------------------------------------
//
//  [  c   s ]  [ a ]  [ sqrt(a*a+b*b) ]
//  [ -s   c ]  [ b ]  [     0         ]
//
class giv
{public:
 giv () : c(1.), s(0.) {}
 giv& operator () (Real& a, Real& b);
 Real rot (const Real& a, const Real& b);
 private:
 Real c,s;
};
//---------------------------------------------------------------------------
class slv
{public:
 slv (int n);
 void solve1 (const Real& del, const Real& dela, pReal* x, pReal* y);
 ~slv ();
 private:
 void sol (pReal f, pReal x);
 void sol (const Real& rl, Real* f[2], Real* x[2]);
 void rx (pReal f, pReal x, const Real& xl, const Real& sa, const Real& sb);
 void wg (pReal y, const Real& dela);
 void wg (pReal* y, const Real& dela);
// Real ext (int n, pReal x);
 int n;
 pReal d,p,q,rl;
};
//---------------------------------------------------------------------------
void fit (int n, pReal x);
//---------------------------------------------------------------------------
