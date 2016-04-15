#ifndef CLASSES_H
#define CLASSES_H

class GS
{
 private:
  double u;
  double v;
 public:
  GS() {u=0.0; v=0.0;}

  ~GS() {}

  GS(GS &gs)
  {
   u = gs.u;
   v = gs.v;
  } 

  void edit_u(double uu) {u = uu;}
  void edit_v(double vv) {v = vv;}
  double get_u() {return u;}
  double get_v() {return v;}
};

#endif

