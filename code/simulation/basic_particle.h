#include <string>

class basic_particle
{
  double x;
  double y;
  double y_first; //y_first is the initial y position of the first parent, tracked for plotting reasons
  int first_parent;

 public:
  //Creator, destructor
  basic_particle();
  basic_particle (double pos_x, double pos_y, double pos_yfirst, int first);
  ~basic_particle();

  //Use private variables
  void set_x (double pos_x);
  void set_y (double pos_y);
  void set_yfirst (double pos_yfirst);
  void set_firstparent (int first);
  double get_x();
  double get_y();
  double get_yfirst();
  int get_firstparent();

  //Movement of each particle
  void diffusion(double Delta,double Lmax); //new values of x and y after particle random walk
  void check_boundaries(double Lmax); //checking that particles remain in a Lmax x Lmax grid
  void pierrehumbert_flow(double U_tot, double k, double phi, double theta, double Lmax); //new values of x and y after advection
};
